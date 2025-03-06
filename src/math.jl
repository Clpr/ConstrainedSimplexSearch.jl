# math helpers
# ------------------------------------------------------------------------------
export boxcenter
export inbox
export modify

export centroid
export reflect
export expand
export contract_out
export contract_in
export shrink

export cv
export isadmissible



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pure geometric & algorithmatic functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function boxcenter(lb::Point{N}, ub::Point{N})::Point{N} where {N}
    return (lb .+ ub) ./ 2.0
end
boxcenter(mp::MinimizeProblem) = boxcenter(mp.lb, mp.ub)
# ------------------------------------------------------------------------------
function inbox(x::Point{N}, lb::Point{N}, ub::Point{N})::Bool where {N}
    return all(lb .<= x) && all(x .<= ub)
end
# ------------------------------------------------------------------------------
function inbox(x::Point{N}, mp::MinimizeProblem{N,P,Q})::Bool where {N,P,Q}
    return inbox(x, mp.lb, mp.ub)
end
# ------------------------------------------------------------------------------
function modify(x::Point{N}, i::Int, v::Real)::Point{N} where {N}
    return SV64{N}([(j == i) ? v : x[j] for j in 1:N])
end
# ------------------------------------------------------------------------------
"""
    centroid(xs::AbsV{Point{N}})::Point{N}

Compute the centroid of a set of points.
"""
function centroid(xs::AbsV{Point{N}})::Point{N} where {N}
    return sum(xs) / length(xs)
end
# ------------------------------------------------------------------------------
"""
    centroid(xs::AbsV{Point{N}}, i::Int)::Point{N}

Compute the centroid of a set of points, exoluding the i-th point.
"""
function centroid(xs::AbsV{Point{N}}, i::Int)::Point{N} where {N}
    n = length(xs)
    return sum(xs[1:i-1;i+1:n]) / (n - 1)
end
# ------------------------------------------------------------------------------
"""
    reflect(xo::Point{N}, xw::Point{N}, α::Real)::Point{N}

Compute the reflection of worst point `xw` w.r.t. the centroid `xo` of a simplex
, where `α` > 0 is the reflection coefficient.

Formula: `xr = xo + α * (xo - xw)`
"""
function reflect(xo::Point{N}, xw::Point{N}, α::Real)::Point{N} where {N}
    return xo .+ α .* (xo .- xw)
end
# ------------------------------------------------------------------------------
"""
    expand(xo::Point{N}, xr::Point{N}, γ::Real)::Point{N}

Compute the expansion of point `xr` w.r.t. the centroid `xo` of a simplex, where
`γ` > 1 is the expansion coefficient.

Formula: `xe = xo + γ * (xr - xo)`
"""
function expand(xo::Point{N}, xr::Point{N}, γ::Real)::Point{N} where {N}
    return xo .+ γ .* (xr .- xo)
end
# ------------------------------------------------------------------------------
"""
    contract_out(xo::Point{N}, xr::Point{N}, β::Real)::Point{N}

Compute the outside contraction of the reflection point `xr` w.r.t. the centroid
`xo` of a simplex, where `0 < ρ <= 0.5` is the contraction coefficient.

Formula: `xc_out = xo + β * (xr - xo)`
"""
function contract_out(xo::Point{N}, xr::Point{N}, ρ::Real)::Point{N} where {N}
    return xo .+ ρ .* (xr .- xo)
end
# ------------------------------------------------------------------------------
"""
    contract_in(xo::Point{N}, xw::Point{N}, ρ::Real)::Point{N} where {N}

Compute the inside contraction of the worst point `xw` w.r.t. the centroid `xo`
of a simplex, where `0 < ρ <= 0.5` is the contraction coefficient.
"""
function contract_in(xo::Point{N}, xw::Point{N}, ρ::Real)::Point{N} where {N}
    return xo .+ ρ .* (xw .- xo)
end
# ------------------------------------------------------------------------------
"""
    shrink(x::Point{N}, xbest::Point{N}, σ::Real)::Point{N} where {N}

Shrink the point `x` towards the best point `xbest` by a factor `σ` ∈ (0,1).
"""
function shrink(x::Point{N}, xbest::Point{N}, σ::Real)::Point{N} where {N}
    return xbest .+ σ .* (x .- xbest)
end
# ------------------------------------------------------------------------------
"""
    shrink!(spl::Vec{Point{N}}, ibest::Int, σ::Real)::Nothing where {N}

Shrink all points in the simplex `spl` towards the best point `spl[ibest]` by a
factor `σ` ∈ (0,1).
"""
function shrink!(
    spl::Vec{Point{N}},
    ibest::Int,
    σ::Real,
) where {N}
    for j in 1:N+1
        if j != ibest
            spl[j] = shrink(spl[j], spl[ibest], σ)
        end
    end
    return nothing
end
# ------------------------------------------------------------------------------
"""
    maxedgelen(spl::Vec{Point{N}})::F64

Compute the maximum edge length of a simplex `spl`.
"""
function maxedgelen(spl::Vec{Point{N}})::F64 where {N}
    lenMax = -Inf
    for i in 1:N+1, j in 1:N+1
        if i != j
            Δ = norm(spl[i] - spl[j], Inf)
            lenMax = max(lenMax, Δ)
        end
    end
    return lenMax
end
# ------------------------------------------------------------------------------
"""
    rebound!(spl::Vec{Point{N}}, lb::Point{N}, ub::Point{N})::Nothing where {N}

Rebound all points in the simplex `spl` within the box constraints `lb` and `ub`
"""
function rebound!(
    spl::Vec{Point{N}},
    lb::Point{N},
    ub::Point{N},
)::Nothing where {N}
    for i in 1:N+1
        spl[i] = clamp.(spl[i], lb, ub)
    end
    return nothing
end









# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constraint violation functions & feasibility check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function _cv_g(gvals::AbsV)::F64
    return sum(Float64, max.(0.0, gvals))
end
function _cv_h(hvals::AbsV, δ::Real, R::Real)::F64
    term1 = sum( Float64, max.(0.0, abs.(hvals) .- δ) )
    term2 = R * Float64(sum(abs2, hvals))
    return term1 + term2
end
function _cv_lub(x::Point{N}, lb::Point{N}, ub::Point{N})::F64 where {N}
    term1 = max.(0.0, lb .- x) |> sum
    term2 = max.(0.0, x .- ub) |> sum
    return term1 + term2
end
# ------------------------------------------------------------------------------
function cv(
    x    ::Point{N},
    gvals::AbsV,
    hvals::AbsV,
    lb   ::Point{N},
    ub   ::Point{N},
    δ    ::Real,
    R    ::Real,
)::F64 where {N}
    res = 0.0    
    if length(gvals) > 0
        _val = _cv_g(gvals)
        @assert !isnan(_val) "NaN found in g(x) <=0 constraint function"
        res += _val
    end
    if length(hvals) > 0
        _val = _cv_h(hvals, δ, R)
        @assert !isnan(_val) "NaN found in h(x) == 0 constraint function"
        res += _val
    end
    res += _cv_lub(x, lb, ub)
    return res
end
# ------------------------------------------------------------------------------
"""
    cv(x::Point{N}, mp::MinimizeProblem{N,P,Q})::F64 where {N,P,Q}

Computes the value of constraint violation (CV) function at point `x` for the
non-linear constraint function `g(x) <= 0`; equality constraint function
`h(x) = 0`; and box constraints `lb <= x <= ub`.
"""
function cv(x::Point{N}, mp::MinimizeProblem{N,P,Q})::F64 where {N,P,Q}
    return cv(x, mp.g(x), mp.h(x), mp.lb, mp.ub, mp.δ, mp.R)
end
# ------------------------------------------------------------------------------
function isadmissible(
    x    ::Point{N},
    gvals::AbsV,
    hvals::AbsV,
    lb   ::Point{N},
    ub   ::Point{N},
    δ    ::Real,
)::Bool where {N}
    if any(gvals .> 0.0)
        return false
    end
    if any(abs.(hvals) .> δ)
        return false
    end
    if any(x .< lb) || any(x .> ub)
        return false
    end
    return true
end
# ------------------------------------------------------------------------------
"""
    isadmissible(x::Point{N}, mp::MinimizeProblem{N,P,Q},δ::Real)::Bool

Check if point `x` is feasible w.r.t. the constraints of the minimization
problem `mp`.
"""
function isadmissible(
    x::Point{N}, 
    mp::MinimizeProblem{N,P,Q}, 
    δ::Real
)::Bool where {N,P,Q}
    return isadmissible(x, mp.g(x), mp.h(x), mp.lb, mp.ub, δ)
end
# ------------------------------------------------------------------------------










 