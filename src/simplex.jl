# Initial simplex strategies


# ------------------------------------------------------------------------------
"""
    AbstractInitialSimplex

Supertype for initial-simplex construction strategies. `optimize` accepts this
interface so callers can choose a local simplex, a box-spanning simplex, or an
explicit matrix without changing the solver pipeline.
"""
abstract type AbstractInitialSimplex end


# ------------------------------------------------------------------------------
"""
    CenteredSimplex

Describe an axis-aligned simplex from the base point `c0`. This strategy is used
when a caller has a meaningful initial guess and wants the first search moves to
remain local relative to the box width on each axis.
"""
struct CenteredSimplex <: AbstractInitialSimplex
    c0::Vector{Float64}
    radius::Float64
    towards::Symbol
end


# ------------------------------------------------------------------------------
"""
    CenteredSimplex(; c0, radius=0.5, towards=:upper)

Construct a centered-simplex strategy after validating its scalar options.
`radius` is the fraction of the distance from `c0` to the selected box side,
and `towards` must be `:upper` or `:lower`. `build` uses this object to create
the row-wise vertex matrix passed to `optimize`.
"""
function CenteredSimplex(;
    c0::AbstractVector,
    radius::Real = 0.5,
    towards::Symbol = :upper,
)
    0.0 < radius <= 1.0 || throw(ArgumentError("radius must be in (0, 1]."))
    towards in (:upper, :lower) ||
        throw(ArgumentError("towards must be :upper or :lower."))
    return CenteredSimplex(Vector{Float64}(c0), Float64(radius), towards)
end


# ------------------------------------------------------------------------------
"""
    MaxVolumeSimplex

Describe an axis-aligned simplex rooted at the lower-bound corner. This
strategy is used when no strong starting point is available and a broad initial
search is preferable to a local one.
"""
struct MaxVolumeSimplex <: AbstractInitialSimplex
    radius::Float64
end


# ------------------------------------------------------------------------------
"""
    MaxVolumeSimplex(; radius=0.7)

Construct a box-spanning simplex strategy with a validated relative radius.
`build` uses the radius to place one vertex along each coordinate axis, spanning
that fraction of the corresponding box width.
"""
function MaxVolumeSimplex(; radius::Real = 0.7)
    0.0 < radius <= 1.0 || throw(ArgumentError("radius must be in (0, 1]."))
    return MaxVolumeSimplex(Float64(radius))
end


# ------------------------------------------------------------------------------
"""
    ExplicitSimplex

Store a user-provided row-wise simplex matrix. This strategy exists for callers
who need complete control over all initial vertices. `optimize` validates shape,
finiteness, box bounds, and nonzero volume before evaluating callbacks.
"""
struct ExplicitSimplex <: AbstractInitialSimplex
    vertices::Matrix{Float64}

    """
        ExplicitSimplex(vertices::Matrix{Float64})

    Store a private copy of an already converted vertex matrix. This inner
    constructor is used by the public abstract-matrix constructor so later
    mutations of caller-owned geometry cannot change a configured strategy.
    """
    function ExplicitSimplex(vertices::Matrix{Float64})
        return new(copy(vertices))
    end
end


# ------------------------------------------------------------------------------
"""
    ExplicitSimplex(vertices)

Construct an explicit-simplex strategy from any matrix. A floating-point copy is
retained so later changes to the caller's matrix cannot silently alter an
optimization setup.
"""
function ExplicitSimplex(vertices::AbstractMatrix)
    return ExplicitSimplex(Matrix{Float64}(vertices))
end


# ------------------------------------------------------------------------------
"""
    build(simplex0::CenteredSimplex, prob::ConstrainedSimplexSearch)

Create the `(ncontrols + 1) x ncontrols` vertex matrix for a local centered
simplex. `optimize` calls this method before common simplex validation. It does
not evaluate constraints or the objective.
"""
function build(simplex0::CenteredSimplex, prob::ConstrainedSimplexSearch)
    length(simplex0.c0) == prob.ncontrols ||
        throw(DimensionMismatch("c0 must have length ncontrols."))
    all(isfinite, simplex0.c0) || throw(ArgumentError("c0 must be finite."))
    inbox(simplex0.c0, prob) || throw(ArgumentError("c0 must lie inside the box."))

    vertices = repeat(simplex0.c0', prob.ncontrols + 1, 1)
    for axis in 1:prob.ncontrols
        if simplex0.towards == :upper
            offset = simplex0.radius * (prob.upper_bounds[axis] - simplex0.c0[axis])
        else
            offset = -simplex0.radius * (simplex0.c0[axis] - prob.lower_bounds[axis])
        end
        vertices[axis + 1, axis] += offset
    end
    return vertices
end


# ------------------------------------------------------------------------------
"""
    build(simplex0::MaxVolumeSimplex, prob::ConstrainedSimplexSearch)

Create a box-spanning row-wise simplex rooted at the lower-bound corner. This is
the broad default used by examples and the Optim frontend when the caller does
not provide an initial simplex. No callback is evaluated during construction.
"""
function build(simplex0::MaxVolumeSimplex, prob::ConstrainedSimplexSearch)
    vertices = repeat(prob.lower_bounds', prob.ncontrols + 1, 1)
    boxwidth = prob.upper_bounds .- prob.lower_bounds
    for axis in 1:prob.ncontrols
        vertices[axis + 1, axis] += simplex0.radius * boxwidth[axis]
    end
    return vertices
end


# ------------------------------------------------------------------------------
"""
    build(simplex0::ExplicitSimplex, prob::ConstrainedSimplexSearch)

Return a validated-size copy of explicit vertices. This method catches shape
errors close to the user-facing strategy while `optimize` performs the remaining
common checks for finite, in-box, nondegenerate geometry.
"""
function build(simplex0::ExplicitSimplex, prob::ConstrainedSimplexSearch)
    expected = (prob.ncontrols + 1, prob.ncontrols)
    size(simplex0.vertices) == expected ||
        throw(DimensionMismatch("explicit simplex must have shape $expected."))
    return copy(simplex0.vertices)
end


# ------------------------------------------------------------------------------
"""
    _validated_simplex(simplex0, prob)

Build and validate an initial simplex before any user callback is evaluated.
This shared gate exists so every strategy obeys the same shape, finiteness, box,
and nonzero-volume requirements used by `optimize`.
"""
function _validated_simplex(
    simplex0::AbstractInitialSimplex,
    prob::ConstrainedSimplexSearch,
)
    vertices = build(simplex0, prob)
    expected = (prob.ncontrols + 1, prob.ncontrols)
    size(vertices) == expected || throw(DimensionMismatch("simplex must have shape $expected."))
    all(isfinite, vertices) || throw(ArgumentError("simplex vertices must be finite."))
    all(prob.lower_bounds' .<= vertices .<= prob.upper_bounds') ||
        throw(ArgumentError("simplex vertices must lie inside the box."))
    simplexvolume(vertices) > 0.0 ||
        throw(ArgumentError("simplex volume must be strictly nonzero."))
    return Matrix{Float64}(vertices)
end
