# Problem definition


# ------------------------------------------------------------------------------
const DEFAULT_HYPERPARAM = Dict{String,Float64}(
    "reflection_factor"          => 1.0,
    "expansion_factor"           => 1.5,
    "contraction_factor_outside" => 0.4,
    "contraction_factor_inside"  => 0.4,
    "shrink_factor"              => 0.5,
)


# ------------------------------------------------------------------------------
"""
    ConstrainedSimplexSearch

Describe a box-bounded optimization problem with nonlinear inequality
constraints. This type is the central problem object used by `optimize`, the
initial-simplex builders, and the optional `Optim.jl` extension.

The constraint callback must return `nconstraints` finite values written as
`g(c) <= 0` and must be defined throughout the box. The objective callback only
needs to be defined at admissible points because the solver always checks the
constraints before calling it.

Use the keyword constructor rather than constructing fields directly so bounds,
dimensions, and callbacks are validated and the default hyperparameters are
copied for the new problem.
"""
mutable struct ConstrainedSimplexSearch{F,G}
    ncontrols::Int
    nconstraints::Int
    function_objective::F
    function_constraints::G
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    hyperparam::Dict{String,Float64}
end


# ------------------------------------------------------------------------------
"""
    ConstrainedSimplexSearch(;
        ncontrols,
        nconstraints,
        function_objective,
        function_constraints,
        lower_bounds,
        upper_bounds,
    )

Construct a validated inequality-constrained optimization problem. This
constructor exists to establish the dimensional and finite-bound assumptions
used by every simplex builder and solver evaluation.

`function_objective(c)` must return one finite real value whenever `c` is
admissible. `function_constraints(c)` must return exactly `nconstraints` finite
values at every point in the box. Each lower bound must be strictly below its
matching upper bound so all controls have a searchable interval.
"""
function ConstrainedSimplexSearch(;
    ncontrols::Int,
    nconstraints::Int,
    function_objective::F,
    function_constraints::G,
    lower_bounds::AbstractVector,
    upper_bounds::AbstractVector,
) where {F<:Function,G<:Function}
    ncontrols > 0 || throw(ArgumentError("ncontrols must be positive."))
    nconstraints >= 0 || throw(ArgumentError("nconstraints must be non-negative."))
    length(lower_bounds) == ncontrols ||
        throw(DimensionMismatch("lower_bounds must have length ncontrols."))
    length(upper_bounds) == ncontrols ||
        throw(DimensionMismatch("upper_bounds must have length ncontrols."))

    lower = Vector{Float64}(lower_bounds)
    upper = Vector{Float64}(upper_bounds)
    all(isfinite, lower) || throw(ArgumentError("lower_bounds must be finite."))
    all(isfinite, upper) || throw(ArgumentError("upper_bounds must be finite."))
    all(lower .< upper) ||
        throw(ArgumentError("each lower bound must be strictly smaller than its upper bound."))

    return ConstrainedSimplexSearch{F,G}(
        ncontrols,
        nconstraints,
        function_objective,
        function_constraints,
        lower,
        upper,
        copy(DEFAULT_HYPERPARAM),
    )
end


# ------------------------------------------------------------------------------
"""
    Base.show(io, prob::ConstrainedSimplexSearch)

Print a compact mathematical summary of a problem. This display is used in the
REPL to make dimensions and the distinction between box feasibility and
inequality admissibility visible without evaluating either callback.
"""
function Base.show(io::IO, prob::ConstrainedSimplexSearch)
    println(io, "ConstrainedSimplexSearch")
    println(io, " min/max f(c), c in R^$(prob.ncontrols)")
    println(io, " subject to $(prob.nconstraints) inequalities g(c) <= 0")
    print(io, " and lower_bounds <= c <= upper_bounds")
end
