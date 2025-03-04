module ConstrainedSimplexSearch
# ==============================================================================
import LinearAlgebra: norm
import Printf: @printf, @sprintf

using StaticArrays


# ------------------------------------------------------------------------------
# alias: std elementary type
const F64 = Float64
const I64 = Int64
const Str = String
const Sym = Symbol

# alias: std collections
const V64 = Vector{Float64}
const M64 = Matrix{Float64}
const Dict64 = Dict{Symbol,Float64}

const Vec{T} = Vector{T}
const Mat{T} = Matrix{T}

# alias: std abstract types
const AbsV = AbstractVector
const AbsM = AbstractMatrix
const AbsVM = AbstractVecOrMat

# alias: StaticArrays.jl
const SV64{D}   = SVector{D,Float64}
const SM64{D,K} = SMatrix{D,K,Float64}


# ------------------------------------------------------------------------------
const Point{D}    = SV64{D}




# ------------------------------------------------------------------------------
include("problem.jl")

include("math.jl")

include("solve.jl")







































# ==============================================================================
end # module ConstrainedSimplexSearch
