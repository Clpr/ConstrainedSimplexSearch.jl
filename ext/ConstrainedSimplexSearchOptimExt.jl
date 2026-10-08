module ConstrainedSimplexSearchOptimExt
# ==============================================================================
import Optim
import ConstrainedSimplexSearch
import ConstrainedSimplexSearch: fine_tune, optimize
import LinearAlgebra: norm

const CSS = ConstrainedSimplexSearch


# ------------------------------------------------------------------------------
"""
    _option_tolerance(absolute, relative, fallback)

Map Optim's absolute and relative tolerance fields to the single positive
tolerance used by constrained simplex search. This adapter helper exists because
current Optim releases retain separate tolerances and use zero to mean disabled;
the native solver instead requires one positive value.
"""
function _option_tolerance(absolute::Real, relative::Real, fallback::Float64)
    selected = max(Float64(absolute), Float64(relative))
    return selected > 0.0 ? selected : fallback
end


# ------------------------------------------------------------------------------
"""
    optimize(objective, lower_bounds, upper_bounds, ConstrainedSimplex(); ...)

Run ConstrainedSimplexSearch through an Optim-style call. This extension builds
the native problem, infers the inequality count by evaluating constraints at
the box center, maps `Optim.Options`, and delegates to the same guarded solver
used by the package API. The objective is therefore still evaluated only at
admissible points.

The returned wrapper supports `Optim.minimizer`, `Optim.minimum`,
`Optim.converged`, and `Optim.iterations`. `initial_simplex` defaults to a
`MaxVolumeSimplex` when omitted.
"""
function optimize(
    objective_function::Function,
    lower_bounds::AbstractVector,
    upper_bounds::AbstractVector,
    ::CSS.ConstrainedSimplex;
    inequality_constraints::Function = c -> Float64[],
    initial_simplex::Union{Nothing,CSS.AbstractInitialSimplex} = nothing,
    options::Optim.Options = Optim.Options(),
    maximize::Bool = false,
    rich_returns::Bool = true,
)
    length(lower_bounds) == length(upper_bounds) ||
        throw(DimensionMismatch("lower_bounds and upper_bounds must have equal length."))
    center = CSS.boxcenter(lower_bounds, upper_bounds)
    sample_constraints = inequality_constraints(copy(center))
    sample_constraints isa AbstractVector ||
        throw(ArgumentError("inequality_constraints must return an AbstractVector."))

    prob = CSS.ConstrainedSimplexSearch(
        ncontrols = length(lower_bounds),
        nconstraints = length(sample_constraints),
        function_objective = objective_function,
        function_constraints = inequality_constraints,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
    )
    simplex0 = isnothing(initial_simplex) ? CSS.MaxVolumeSimplex() : initial_simplex
    objective_tolerance = _option_tolerance(
        options.f_abstol,
        options.f_reltol,
        1E-4,
    )
    control_tolerance = _option_tolerance(
        options.x_abstol,
        options.x_reltol,
        1E-4,
    )

    raw_result = optimize(
        prob,
        simplex0;
        use_maximize = maximize,
        max_iteration = options.iterations,
        objective_tolerance = objective_tolerance,
        control_tolerance = control_tolerance,
        verbose = options.show_trace,
        showevery = options.show_every,
        rich_returns = rich_returns,
    )
    return CSS.ConstrainedSimplexOptimResult(raw_result)
end


# ------------------------------------------------------------------------------
"""
    Optim.minimizer(result::ConstrainedSimplexOptimResult)

Return the optimal control vector from an Optim-style constrained simplex
result. This accessor lets downstream code consume the adapter result through a
familiar Optim interface.
"""
function Optim.minimizer(result::CSS.ConstrainedSimplexOptimResult)
    return result.raw_result["optimal_control"]
end


# ------------------------------------------------------------------------------
"""
    Optim.minimum(result::ConstrainedSimplexOptimResult)

Return the original objective value from an Optim-style result. For a run with
`maximize=true`, this is the maximum value rather than the internally negated
search score.
"""
function Optim.minimum(result::CSS.ConstrainedSimplexOptimResult)
    return result.raw_result["optimal_objective"]
end


# ------------------------------------------------------------------------------
"""
    Optim.converged(result::ConstrainedSimplexOptimResult)

Return whether the native convergence status starts with `"converged"`. This
reports numerical termination only; callers that need the separate
admissibility status can inspect `result.raw_result`.
"""
function Optim.converged(result::CSS.ConstrainedSimplexOptimResult)
    return startswith(result.raw_result["convergence_status"], "converged")
end


# ------------------------------------------------------------------------------
"""
    Optim.iterations(result::ConstrainedSimplexOptimResult)

Return the native solver's exit iteration. This accessor completes the small
Optim-compatible result surface needed by common frontend workflows.
"""
function Optim.iterations(result::CSS.ConstrainedSimplexOptimResult)
    return result.raw_result["exit_iteration"]
end


# ------------------------------------------------------------------------------
"""
    _sigmoid(value)

Map an unconstrained real value into `(0, 1)`. Hyperparameter tuning uses this
transformation to keep contraction and shrink factors inside their valid ranges
while Optim works in an unconstrained parameter space.
"""
function _sigmoid(value::Real)
    checked = clamp(Float64(value), -700.0, 700.0)
    return 1.0 / (1.0 + exp(-checked))
end


# ------------------------------------------------------------------------------
"""
    _pack_hyperparam(hyperparam)

Transform valid simplex hyperparameters into an unconstrained five-vector.
`fine_tune` uses this inverse parameterization as the initial point for
`Optim.optimize`, preserving the current user settings.
"""
function _pack_hyperparam(hyperparam::Dict{String,Float64})
    values = CSS._validated_hyperparam(hyperparam)
    return [
        log(values["reflection_factor"]),
        log(values["expansion_factor"] - 1.0),
        log(values["contraction_factor_outside"] /
            (0.5 - values["contraction_factor_outside"])),
        log(values["contraction_factor_inside"] /
            (0.5 - values["contraction_factor_inside"])),
        log(values["shrink_factor"] / (1.0 - values["shrink_factor"])),
    ]
end


# ------------------------------------------------------------------------------
"""
    _unpack_hyperparam(parameters)

Transform an unconstrained five-vector into valid simplex hyperparameters. This
function implements the release's documented exponential and logistic maps and
is used for every benchmark solve performed by `fine_tune`.
"""
function _unpack_hyperparam(parameters::AbstractVector)
    length(parameters) == 5 || throw(DimensionMismatch("tuning parameters must have length 5."))
    return Dict{String,Float64}(
        "reflection_factor" => exp(parameters[1]),
        "expansion_factor" => 1.0 + exp(parameters[2]),
        "contraction_factor_outside" => 0.5 * _sigmoid(parameters[3]),
        "contraction_factor_inside" => 0.5 * _sigmoid(parameters[4]),
        "shrink_factor" => _sigmoid(parameters[5]),
    )
end


# ------------------------------------------------------------------------------
"""
    fine_tune(prob, simplex0; true_control, ...)

Tune the five simplex coefficients on a representative problem with a known
solution. Optim minimizes the distance between each guarded solver result and
`true_control` in an unconstrained transformed parameter space; non-admissible
solver exits receive a large added loss.

The tuned dictionary is assigned to `prob.hyperparam` and returned with
diagnostics. If tuning throws, the original dictionary is restored. This
routine is intended for offline calibration, not for every disposable solve in
a dynamic-programming loop.
"""
function fine_tune(
    prob::CSS.ConstrainedSimplexSearch,
    simplex0::CSS.AbstractInitialSimplex;
    true_control::AbstractVector,
    use_maximize::Bool = false,
    max_iteration::Int = 1000,
    tolerance::Float64 = 1E-4,
    verbose::Bool = false,
    showevery::Int = 2,
    method = Optim.BFGS(),
)
    length(true_control) == prob.ncontrols ||
        throw(DimensionMismatch("true_control must have length ncontrols."))
    all(isfinite, true_control) || throw(ArgumentError("true_control must be finite."))
    max_iteration > 0 || throw(ArgumentError("max_iteration must be positive."))
    tolerance > 0.0 || throw(ArgumentError("tolerance must be positive."))
    showevery > 0 || throw(ArgumentError("showevery must be positive."))

    target = Vector{Float64}(true_control)
    original_hyperparam = copy(prob.hyperparam)
    start_parameters = _pack_hyperparam(original_hyperparam)
    call_counter = Ref(0)

    calibration_loss = function(parameters)
        call_counter[] += 1
        prob.hyperparam = _unpack_hyperparam(parameters)
        solved = optimize(
            prob,
            simplex0;
            use_maximize = use_maximize,
            max_iteration = max_iteration,
            objective_tolerance = tolerance,
            control_tolerance = tolerance,
            rich_returns = false,
        )
        loss = norm(solved["optimal_control"] .- target)
        solved["admissibility_status"] == "admissible" || (loss += 1.0E6)
        if verbose && call_counter[] % showevery == 0
            println("fine_tune call = $(call_counter[]), loss = $loss")
        end
        return loss
    end

    try
        tuning_result = Optim.optimize(
            calibration_loss,
            start_parameters,
            method,
            Optim.Options(
                iterations = max_iteration,
                f_abstol = tolerance,
                x_abstol = tolerance,
                show_trace = false,
            ),
        )
        tuned_hyperparam = _unpack_hyperparam(Optim.minimizer(tuning_result))
        prob.hyperparam = copy(tuned_hyperparam)
        final_solve = optimize(
            prob,
            simplex0;
            use_maximize = use_maximize,
            max_iteration = max_iteration,
            objective_tolerance = tolerance,
            control_tolerance = tolerance,
            rich_returns = false,
        )
        minimized_error = norm(final_solve["optimal_control"] .- target)
        diagnostics = Dict{String,Any}(
            "optimal_control" => final_solve["optimal_control"],
            "true_control" => target,
            "minimized_error" => minimized_error,
            "converged" => Optim.converged(tuning_result) || minimized_error <= tolerance,
            "admissible" => final_solve["admissibility_status"] == "admissible",
            "exit_iteration" => final_solve["exit_iteration"],
            "max_iteration" => max_iteration,
            "method" => string(method),
        )
        return tuned_hyperparam, diagnostics
    catch
        prob.hyperparam = original_hyperparam
        rethrow()
    end
end

# ==============================================================================
end # module ConstrainedSimplexSearchOptimExt
