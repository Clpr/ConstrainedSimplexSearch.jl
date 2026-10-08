# Solver


# ------------------------------------------------------------------------------
"""
    RawEvaluation

Store constraint-first evaluation data for one candidate point. This internal
record exists so `optimize` can compare admissible and non-admissible points
without ever inventing an objective value for a point outside the admissible
set. Only `_raw_evaluate` constructs these records.
"""
struct RawEvaluation
    point::Vector{Float64}
    constraint_values::Vector{Float64}
    constraint_violation::Float64
    admissible::Bool
    objective_value::Float64
    signed_objective_value::Float64
end


# ------------------------------------------------------------------------------
"""
    ScoredEvaluation

Pair a raw guarded evaluation with its context-dependent search score. This
internal wrapper is used throughout move selection because non-admissible
points are ranked by violation while admissible points are ranked by the signed
objective.
"""
struct ScoredEvaluation
    raw::RawEvaluation
    score::Float64
end


# ------------------------------------------------------------------------------
"""
    ConstrainedSimplex()

Select the ConstrainedSimplexSearch algorithm in the optional `Optim.jl`-style
frontend. The marker carries no settings; the extension maps `Optim.Options`
and the remaining keywords into the native guarded `optimize` pipeline.
"""
struct ConstrainedSimplex end


# ------------------------------------------------------------------------------
"""
    ConstrainedSimplexOptimResult

Wrap the native result dictionary for the optional `Optim.jl` accessors. The
extension uses this internal-compatible container to implement
`Optim.minimizer`, `Optim.minimum`, `Optim.converged`, and `Optim.iterations`
without changing the native return format.
"""
struct ConstrainedSimplexOptimResult
    raw_result::Dict{String,Any}
end


# ------------------------------------------------------------------------------
"""
    fine_tune(prob, simplex0; true_control, ...)

Tune simplex hyperparameters against a benchmark with a known solution. The
method is supplied by the optional `Optim.jl` extension because tuning requires
an unconstrained optimizer, while ordinary constrained simplex searches do not
need that dependency. Load `Optim` before calling this function.
"""
function fine_tune end


# ------------------------------------------------------------------------------
"""
    _evaluate_constraints(prob, point)

Evaluate and validate the inequality callback at one in-box point. This helper
centralizes the callback contract for public admissibility checks and the
solver's guarded evaluation path. It never calls the objective.
"""
function _evaluate_constraints(
    prob::ConstrainedSimplexSearch,
    point::Vector{Float64},
)
    raw_values = prob.function_constraints(copy(point))
    raw_values isa AbstractVector ||
        throw(ArgumentError("function_constraints must return an AbstractVector."))
    length(raw_values) == prob.nconstraints ||
        throw(DimensionMismatch("function_constraints must return nconstraints values."))
    values = Vector{Float64}(raw_values)
    all(isfinite, values) ||
        throw(ArgumentError("function_constraints returned a non-finite value."))
    return values
end


# ------------------------------------------------------------------------------
"""
    _evaluate_objective(prob, point)

Evaluate and validate the scalar objective at a known-admissible point. This
function exists to enforce scalar, real, finite returns and is called only from
`_raw_evaluate` after all inequality constraints have passed.
"""
function _evaluate_objective(
    prob::ConstrainedSimplexSearch,
    point::Vector{Float64},
)
    raw_value = prob.function_objective(copy(point))
    raw_value isa Real || throw(ArgumentError("function_objective must return a real scalar."))
    value = Float64(raw_value)
    isfinite(value) ||
        throw(ArgumentError("function_objective returned a non-finite value at an admissible point."))
    return value
end


# ------------------------------------------------------------------------------
"""
    _raw_evaluate(prob, point, use_maximize)

Clamp a candidate into the box, evaluate constraints, and call the objective
only when the candidate is admissible. Every objective call made by `optimize`,
including final result selection, passes through this safety gate.

For maximization, `signed_objective_value` is the negated objective so the move
logic can remain a minimization algorithm. Non-admissible candidates retain
`NaN` objective fields because no objective call occurs.
"""
function _raw_evaluate(
    prob::ConstrainedSimplexSearch,
    point::AbstractVector,
    use_maximize::Bool,
)
    length(point) == prob.ncontrols ||
        throw(DimensionMismatch("candidate point must have length ncontrols."))
    checked = clamp.(Vector{Float64}(point), prob.lower_bounds, prob.upper_bounds)
    constraint_values = _evaluate_constraints(prob, checked)
    violation = constraintviolation(
        checked,
        constraint_values,
        prob.lower_bounds,
        prob.upper_bounds,
    )
    admissible = isadmissible(
        checked,
        constraint_values,
        prob.lower_bounds,
        prob.upper_bounds,
    )

    objective_value = NaN
    signed_objective_value = NaN
    if admissible
        objective_value = _evaluate_objective(prob, checked)
        signed_objective_value = use_maximize ? -objective_value : objective_value
    end

    return RawEvaluation(
        checked,
        constraint_values,
        violation,
        admissible,
        objective_value,
        signed_objective_value,
    )
end


# ------------------------------------------------------------------------------
"""
    _score_raw_items(raw_items, context=nothing)

Assign comparable scores to guarded evaluations. With no admissible context,
violation alone guides the simplex, so the first admissible candidate receives
a zero score regardless of its objective scale. Once the context contains an
admissible point, invalid points rank after the worst admissible objective by
adding their violation and valid points use signed objective values.

`optimize` uses this one policy for current vertices and every candidate move,
which keeps the objective guard and ranking behavior coherent across branches.
"""
function _score_raw_items(
    raw_items::Vector{RawEvaluation},
    context::Union{Nothing,Vector{ScoredEvaluation}} = nothing,
)
    context_raw = isnothing(context) ? raw_items : [item.raw for item in context]
    admissible_scores = [
        item.signed_objective_value for item in context_raw if item.admissible
    ]
    worst_admissible = isempty(admissible_scores) ? nothing : maximum(admissible_scores)

    scored = Vector{ScoredEvaluation}(undef, length(raw_items))
    for i in eachindex(raw_items)
        item = raw_items[i]
        score = if isnothing(worst_admissible)
            item.constraint_violation
        elseif item.admissible
            item.signed_objective_value
        else
            worst_admissible + item.constraint_violation
        end
        scored[i] = ScoredEvaluation(item, Float64(score))
    end
    return scored
end


# ------------------------------------------------------------------------------
"""
    _score_simplex(prob, vertices, use_maximize, parallel_mode, max_workers)

Guardedly evaluate and score every current simplex vertex. `optimize` calls this
at each iteration and for final selection. `:serial` is deterministic and safe
for ordinary callbacks; `:thread` uses Julia threads and therefore requires
thread-safe callbacks. `max_workers` optionally limits the number of spawned
worker tasks. Process callback transport is intentionally unsupported.
"""
function _score_simplex(
    prob::ConstrainedSimplexSearch,
    vertices::Matrix{Float64},
    use_maximize::Bool,
    parallel_mode::Symbol,
    max_workers::Union{Nothing,Int},
)
    raw_items = Vector{RawEvaluation}(undef, size(vertices, 1))
    if parallel_mode == :serial
        for i in axes(vertices, 1)
            raw_items[i] = _raw_evaluate(prob, view(vertices, i, :), use_maximize)
        end
    elseif parallel_mode == :thread
        if isnothing(max_workers)
            Threads.@threads for i in axes(vertices, 1)
                raw_items[i] = _raw_evaluate(prob, view(vertices, i, :), use_maximize)
            end
        else
            worker_count = min(max_workers, size(vertices, 1))
            @sync for worker in 1:worker_count
                Threads.@spawn for i in worker:worker_count:size(vertices, 1)
                    raw_items[i] = _raw_evaluate(prob, view(vertices, i, :), use_maximize)
                end
            end
        end
    else
        throw(ArgumentError("parallel_mode=:process is not implemented because callback transport is not available."))
    end
    return _score_raw_items(raw_items)
end


# ------------------------------------------------------------------------------
"""
    _score_candidate(prob, point, context, use_maximize)

Guardedly evaluate one move candidate and score it against the current simplex.
Reflection, expansion, and both contraction branches all use this helper, so no
candidate can bypass the constraint-first objective rule.
"""
function _score_candidate(
    prob::ConstrainedSimplexSearch,
    point::AbstractVector,
    context::Vector{ScoredEvaluation},
    use_maximize::Bool,
)
    raw = _raw_evaluate(prob, point, use_maximize)
    return only(_score_raw_items([raw], context))
end


# ------------------------------------------------------------------------------
"""
    _validated_hyperparam(hyperparam)

Return a checked floating-point copy of the five simplex coefficients. The
solver validates this mutable user-facing dictionary at the start of every run
so manual tuning cannot create undefined move geometry.
"""
function _validated_hyperparam(hyperparam::Dict{String,Float64})
    Set(keys(hyperparam)) == Set(keys(DEFAULT_HYPERPARAM)) ||
        throw(ArgumentError("hyperparam must contain exactly the default keys."))
    values = copy(hyperparam)
    0.0 < values["reflection_factor"] < Inf ||
        throw(ArgumentError("reflection_factor must be in (0, Inf)."))
    1.0 < values["expansion_factor"] < Inf ||
        throw(ArgumentError("expansion_factor must be in (1, Inf)."))
    0.0 < values["contraction_factor_outside"] < 0.5 ||
        throw(ArgumentError("contraction_factor_outside must be in (0, 0.5)."))
    0.0 < values["contraction_factor_inside"] < 0.5 ||
        throw(ArgumentError("contraction_factor_inside must be in (0, 0.5)."))
    0.0 < values["shrink_factor"] < 1.0 ||
        throw(ArgumentError("shrink_factor must be in (0, 1)."))
    return values
end


# ------------------------------------------------------------------------------
"""
    optimize(prob, simplex0; ...)

Solve a box-bounded problem with arbitrary nonlinear inequality constraints
using constrained simplex search. The algorithm ranks non-admissible points by
positive-part constraint violation and admissible points by objective value.
Every candidate is clipped into the box, and the objective is called only after
the inequality callback confirms admissibility.

Set `use_maximize=true` to maximize without changing the callback. Convergence
can be triggered by objective-score spread, simplex edge length, or both. The
returned `Dict{String,Any}` always reports convergence and admissibility
separately; rich returns additionally include traces, timing, settings, and move
counters.

`parallel_mode` accepts `:serial` and `:thread`. Thread mode applies only to the
current simplex because candidate branches are sequential, and callbacks must
be thread-safe. `:process` raises a clear `ArgumentError`.
"""
function optimize(
    prob::ConstrainedSimplexSearch,
    simplex0::AbstractInitialSimplex;
    use_maximize::Bool = false,
    max_iteration::Int = 1000,
    objective_tolerance::Float64 = 1E-4,
    control_tolerance::Float64 = 1E-4,
    verbose::Bool = false,
    showevery::Int = 2,
    rich_returns::Bool = true,
    parallel_mode::Symbol = :serial,
    max_workers::Union{Nothing,Int} = nothing,
)
    max_iteration > 0 || throw(ArgumentError("max_iteration must be positive."))
    objective_tolerance > 0.0 ||
        throw(ArgumentError("objective_tolerance must be positive."))
    control_tolerance > 0.0 ||
        throw(ArgumentError("control_tolerance must be positive."))
    showevery > 0 || throw(ArgumentError("showevery must be positive."))
    parallel_mode in (:serial, :thread, :process) ||
        throw(ArgumentError("parallel_mode must be :serial, :thread, or :process."))
    parallel_mode == :process &&
        throw(ArgumentError("parallel_mode=:process is not implemented because callback transport is not available."))
    isnothing(max_workers) || max_workers > 0 ||
        throw(ArgumentError("max_workers must be positive when provided."))

    hyperparam = _validated_hyperparam(prob.hyperparam)
    vertices = _validated_simplex(simplex0, prob)
    start_time = time_ns()

    lower_trace = Float64[]
    upper_trace = Float64[]
    objective_error_trace = Float64[]
    control_error_trace = Float64[]
    counters = Dict{String,Int}(
        "counter_reflection" => 0,
        "counter_expansion" => 0,
        "counter_outside_contraction" => 0,
        "counter_inside_contraction" => 0,
        "counter_shrink" => 0,
    )

    exit_iteration = 0
    exit_objective_error = Inf
    exit_control_error = Inf
    convergence_status = "not_converged"

    for iteration in 1:max_iteration
        exit_iteration = iteration
        scored_vertices = _score_simplex(
            prob,
            vertices,
            use_maximize,
            parallel_mode,
            max_workers,
        )
        scores = [item.score for item in scored_vertices]
        order = sortperm(scores)

        best_score = scores[order[1]]
        worst_score = scores[order[end]]
        second_worst_score = scores[order[end - 1]]
        exit_objective_error = abs(worst_score - best_score)
        exit_control_error = maxedgelen(vertices)

        push!(lower_trace, best_score)
        push!(upper_trace, worst_score)
        push!(objective_error_trace, exit_objective_error)
        push!(control_error_trace, exit_control_error)

        if verbose && iteration % showevery == 0
            @printf(
                "iter = %d, score in (%.3e, %.3e), score_gap = %.3e, max_edge = %.3e\n",
                iteration,
                best_score,
                worst_score,
                exit_objective_error,
                exit_control_error,
            )
        end

        objective_converged = exit_objective_error <= objective_tolerance
        control_converged = exit_control_error <= control_tolerance
        if objective_converged && control_converged
            convergence_status = "converged_both"
            break
        elseif objective_converged
            convergence_status = "converged_objective"
            break
        elseif control_converged
            convergence_status = "converged_control"
            break
        elseif iteration == max_iteration
            convergence_status = "reached_max_iteration"
            break
        end

        best_index = order[1]
        worst_index = order[end]
        best_point = copy(vertices[best_index, :])
        worst_point = copy(vertices[worst_index, :])
        remaining = [i for i in axes(vertices, 1) if i != worst_index]
        center = centroid(vertices[remaining, :])

        reflected = clamp.(
            reflect(center, worst_point, hyperparam["reflection_factor"]),
            prob.lower_bounds,
            prob.upper_bounds,
        )
        reflected_eval = _score_candidate(prob, reflected, scored_vertices, use_maximize)
        chosen_move = "reflection"

        if best_score <= reflected_eval.score < second_worst_score
            vertices[worst_index, :] = reflected
            counters["counter_reflection"] += 1
        elseif reflected_eval.score < best_score
            expanded = clamp.(
                expand(center, reflected, hyperparam["expansion_factor"]),
                prob.lower_bounds,
                prob.upper_bounds,
            )
            expanded_eval = _score_candidate(prob, expanded, scored_vertices, use_maximize)
            if expanded_eval.score < reflected_eval.score
                vertices[worst_index, :] = expanded
                counters["counter_expansion"] += 1
                chosen_move = "expansion"
            else
                vertices[worst_index, :] = reflected
                counters["counter_reflection"] += 1
            end
        elseif second_worst_score <= reflected_eval.score <= worst_score
            contracted = clamp.(
                contract_out(
                    center,
                    reflected,
                    hyperparam["contraction_factor_outside"],
                ),
                prob.lower_bounds,
                prob.upper_bounds,
            )
            contracted_eval = _score_candidate(prob, contracted, scored_vertices, use_maximize)
            if contracted_eval.score < worst_score
                vertices[worst_index, :] = contracted
                counters["counter_outside_contraction"] += 1
                chosen_move = "outside_contraction"
            else
                for i in axes(vertices, 1)
                    vertices[i, :] = clamp.(
                        shrink(
                            view(vertices, i, :),
                            best_point,
                            hyperparam["shrink_factor"],
                        ),
                        prob.lower_bounds,
                        prob.upper_bounds,
                    )
                end
                counters["counter_shrink"] += 1
                chosen_move = "shrink"
            end
        elseif reflected_eval.score > worst_score
            contracted = clamp.(
                contract_in(
                    center,
                    worst_point,
                    hyperparam["contraction_factor_inside"],
                ),
                prob.lower_bounds,
                prob.upper_bounds,
            )
            contracted_eval = _score_candidate(prob, contracted, scored_vertices, use_maximize)
            if contracted_eval.score < worst_score
                vertices[worst_index, :] = contracted
                counters["counter_inside_contraction"] += 1
                chosen_move = "inside_contraction"
            else
                for i in axes(vertices, 1)
                    vertices[i, :] = clamp.(
                        shrink(
                            view(vertices, i, :),
                            best_point,
                            hyperparam["shrink_factor"],
                        ),
                        prob.lower_bounds,
                        prob.upper_bounds,
                    )
                end
                counters["counter_shrink"] += 1
                chosen_move = "shrink"
            end
        else
            for i in axes(vertices, 1)
                vertices[i, :] = clamp.(
                    shrink(
                        view(vertices, i, :),
                        best_point,
                        hyperparam["shrink_factor"],
                    ),
                    prob.lower_bounds,
                    prob.upper_bounds,
                )
            end
            counters["counter_shrink"] += 1
            chosen_move = "shrink"
        end

        if verbose && iteration % showevery == 0
            println("  move: $chosen_move")
        end
    end

    final_scored = _score_simplex(
        prob,
        vertices,
        use_maximize,
        parallel_mode,
        max_workers,
    )
    admissible_indices = [i for i in eachindex(final_scored) if final_scored[i].raw.admissible]
    if isempty(admissible_indices)
        best_index = argmin(item.score for item in final_scored)
        best_item = final_scored[best_index]
        admissibility_status = "non_admissible"
        optimal_objective = NaN
    else
        local_index = argmin(final_scored[i].score for i in admissible_indices)
        best_item = final_scored[admissible_indices[local_index]]
        admissibility_status = "admissible"
        optimal_objective = best_item.raw.objective_value
    end

    result = Dict{String,Any}(
        "optimal_control" => copy(best_item.raw.point),
        "optimal_objective" => optimal_objective,
        "convergence_status" => convergence_status,
        "admissibility_status" => admissibility_status,
        "exit_objective_error" => exit_objective_error,
        "exit_control_error" => exit_control_error,
        "exit_iteration" => exit_iteration,
    )

    if rich_returns
        merge!(
            result,
            Dict{String,Any}(
                "problem_type" => use_maximize ? "maximization" : "minimization",
                "initial_simplex_type" => string(nameof(typeof(simplex0))),
                "objective_lower_bound_trace" => lower_trace,
                "objective_upper_bound_trace" => upper_trace,
                "error_trace_objective" => objective_error_trace,
                "error_trace_control" => control_error_trace,
                "elapsed_walltime_seconds" => (time_ns() - start_time) / 1.0E9,
                "max_iteration" => max_iteration,
                "objective_tolerance" => objective_tolerance,
                "control_tolerance" => control_tolerance,
                "parallel_mode" => parallel_mode,
                "max_workers" => max_workers,
            ),
        )
        merge!(result, counters)
    end
    return result
end
