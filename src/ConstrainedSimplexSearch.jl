module ConstrainedSimplexSearch
# ==============================================================================
import LinearAlgebra: det, norm
import Printf: @printf

export ConstrainedSimplexSearch
export CenteredSimplex
export MaxVolumeSimplex
export ExplicitSimplex
export ConstrainedSimplex
export optimize
export fine_tune
export isadmissible
export constraintviolation
export boxcenter
export inbox
export centroid
export reflect
export expand
export contract_out
export contract_in
export shrink
export maxedgelen
export simplexvolume
export DEFAULT_HYPERPARAM


# ------------------------------------------------------------------------------
include("problem.jl")
include("simplex.jl")
include("math.jl")
include("solve.jl")

# ==============================================================================
end # module ConstrainedSimplexSearch
