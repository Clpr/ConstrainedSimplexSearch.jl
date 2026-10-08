# Geometry and constraint helpers


# ------------------------------------------------------------------------------
"""
    boxcenter(lower_bounds, upper_bounds)
    boxcenter(prob)

Return the midpoint of every box interval. This helper is used when a neutral
starting point is needed for a `CenteredSimplex` or for inspecting a problem's
search domain. It performs no callback evaluation.
"""
function boxcenter(lower_bounds::AbstractVector, upper_bounds::AbstractVector)
    length(lower_bounds) == length(upper_bounds) ||
        throw(DimensionMismatch("bound vectors must have equal length."))
    return (Vector{Float64}(lower_bounds) .+ Vector{Float64}(upper_bounds)) ./ 2.0
end

function boxcenter(prob::ConstrainedSimplexSearch)
    return boxcenter(prob.lower_bounds, prob.upper_bounds)
end


# ------------------------------------------------------------------------------
"""
    inbox(point, lower_bounds, upper_bounds)
    inbox(point, prob)

Return whether every coordinate lies inside or on the box boundary. The simplex
builders and public admissibility check use this helper before any objective
evaluation because box feasibility is a prerequisite for admissibility.
"""
function inbox(
    point::AbstractVector,
    lower_bounds::AbstractVector,
    upper_bounds::AbstractVector,
)
    length(point) == length(lower_bounds) == length(upper_bounds) || return false
    return all(lower_bounds .<= point) && all(point .<= upper_bounds)
end

function inbox(point::AbstractVector, prob::ConstrainedSimplexSearch)
    return inbox(point, prob.lower_bounds, prob.upper_bounds)
end


# ------------------------------------------------------------------------------
"""
    centroid(vertices)

Return the coordinate-wise mean of row-wise simplex vertices. The solver uses
this operation on all vertices except the current worst point to construct
reflection and contraction candidates.
"""
function centroid(vertices::AbstractMatrix)
    size(vertices, 1) > 0 || throw(ArgumentError("vertices must not be empty."))
    return vec(sum(vertices; dims = 1)) ./ size(vertices, 1)
end


# ------------------------------------------------------------------------------
"""
    reflect(center, worst, factor)

Reflect `worst` through `center` by a positive factor. This is the first
candidate move considered in each `optimize` iteration; the solver clips and
guardedly evaluates the returned point before deciding whether to accept it.
"""
function reflect(center::AbstractVector, worst::AbstractVector, factor::Real)
    return Vector{Float64}(center .+ factor .* (center .- worst))
end


# ------------------------------------------------------------------------------
"""
    expand(center, reflected, factor)

Move farther from `center` through an improving reflected point. `optimize`
uses this geometry only after reflection improves on the current best point and
then applies the same constraint-first objective guard to the candidate.
"""
function expand(center::AbstractVector, reflected::AbstractVector, factor::Real)
    return Vector{Float64}(center .+ factor .* (reflected .- center))
end


# ------------------------------------------------------------------------------
"""
    contract_out(center, reflected, factor)

Move from `center` partway toward a mediocre reflected point. This outside
contraction is used by `optimize` when reflection is no better than the
second-worst simplex score.
"""
function contract_out(center::AbstractVector, reflected::AbstractVector, factor::Real)
    return Vector{Float64}(center .+ factor .* (reflected .- center))
end


# ------------------------------------------------------------------------------
"""
    contract_in(center, worst, factor)

Move from `center` partway toward the current worst point. This inside
contraction is used by `optimize` after a reflection scores worse than the
current worst vertex.
"""
function contract_in(center::AbstractVector, worst::AbstractVector, factor::Real)
    return Vector{Float64}(center .+ factor .* (worst .- center))
end


# ------------------------------------------------------------------------------
"""
    shrink(point, best, factor)

Move one point toward the best simplex point. `optimize` applies this operation
to every vertex when both reflection and contraction fail, then clips the new
simplex into the box before its next guarded evaluation.
"""
function shrink(point::AbstractVector, best::AbstractVector, factor::Real)
    return Vector{Float64}(best .+ factor .* (point .- best))
end


# ------------------------------------------------------------------------------
"""
    maxedgelen(vertices)

Return the largest infinity-norm edge length in a row-wise simplex. `optimize`
records this quantity as its control-space error and compares it with
`control_tolerance` when determining convergence.
"""
function maxedgelen(vertices::AbstractMatrix)
    edge_max = 0.0
    for i in axes(vertices, 1), j in (i + 1):size(vertices, 1)
        edge_max = max(edge_max, norm(view(vertices, i, :) .- view(vertices, j, :), Inf))
    end
    return edge_max
end


# ------------------------------------------------------------------------------
"""
    simplexvolume(vertices)

Return the geometric volume of a row-wise N-dimensional simplex. Initial
simplex validation uses this determinant formula to reject degenerate geometry
before constraints or objectives are evaluated.
"""
function simplexvolume(vertices::AbstractMatrix)
    ncontrols = size(vertices, 2)
    size(vertices, 1) == ncontrols + 1 ||
        throw(DimensionMismatch("a simplex must have ncontrols + 1 rows."))
    edges = vertices[2:end, :] .- vertices[1, :]'
    return abs(det(edges)) / factorial(ncontrols)
end


# ------------------------------------------------------------------------------
"""
    constraintviolation(point, constraint_values, lower_bounds, upper_bounds)
    constraintviolation(point, prob)

Return the L1 positive-part violation of inequality and box constraints. The
solver uses this score to guide a wholly non-admissible simplex toward the
admissible set without calling the objective at any invalid point.

The problem overload evaluates only the constraint callback. Constraints must
therefore be finite throughout the feasible box, while the objective remains
guarded and is not used by this function.
"""
function constraintviolation(
    point::AbstractVector,
    constraint_values::AbstractVector,
    lower_bounds::AbstractVector,
    upper_bounds::AbstractVector,
)
    inequality = sum((max(0.0, value) for value in constraint_values); init = 0.0)
    below = sum((max(0.0, lower_bounds[i] - point[i]) for i in eachindex(point)); init = 0.0)
    above = sum((max(0.0, point[i] - upper_bounds[i]) for i in eachindex(point)); init = 0.0)
    return Float64(inequality + below + above)
end

function constraintviolation(point::AbstractVector, prob::ConstrainedSimplexSearch)
    length(point) == prob.ncontrols ||
        throw(DimensionMismatch("point must have length ncontrols."))
    checked = clamp.(Vector{Float64}(point), prob.lower_bounds, prob.upper_bounds)
    values = _evaluate_constraints(prob, checked)
    return constraintviolation(point, values, prob.lower_bounds, prob.upper_bounds)
end


# ------------------------------------------------------------------------------
"""
    isadmissible(point, constraint_values, lower_bounds, upper_bounds)
    isadmissible(point, prob)

Return whether a point lies in the feasible box and satisfies every inequality
`g(c) <= 0`. The solver's guarded evaluation path uses the value-based method
after checking callback shape and finiteness; the problem overload is provided
for users and never calls the objective.
"""
function isadmissible(
    point::AbstractVector,
    constraint_values::AbstractVector,
    lower_bounds::AbstractVector,
    upper_bounds::AbstractVector,
)
    return inbox(point, lower_bounds, upper_bounds) && all(constraint_values .<= 0.0)
end

function isadmissible(point::AbstractVector, prob::ConstrainedSimplexSearch)
    length(point) == prob.ncontrols ||
        throw(DimensionMismatch("point must have length ncontrols."))
    inbox(point, prob) || return false
    values = _evaluate_constraints(prob, Vector{Float64}(point))
    return all(values .<= 0.0)
end
