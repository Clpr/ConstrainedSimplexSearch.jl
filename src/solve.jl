# Solver
export initsimplex
export solve



# ------------------------------------------------------------------------------
"""
    initsimplex(x0, mp; radius=0.5, towards=:upper)

Initialize a simplex around the initial point `x0` for the minimization problem 
`mp`. The simplex is initialized by moving `radius` towards the box boundaries,
where `radius` is a relative radius regarding the distance between `x0` and the
box boundaries. The `towards` parameter specifies towards which boundary to
build the simplex (:upper or :lower)

Returns a simlex which is a vector of N+1 Points{N}.

## Example
```julia
x0 = css.boxcenter(mp)
initsimplex(x0, mp)
```
"""
function initsimplex(
    x0::Point{N},
    mp::MinimizeProblem{N,P,Q} ;

    radius::F64 = 0.5, # simplex (relative radius) towards the box boundaries
    towards::Symbol = :upper, # towards which boundary to build the simplex
)::Vec{Point{N}} where {N,P,Q}
    
    @assert inbox(x0,mp) "x0 must be interior wrt the box constraints"
    @assert (0.0 < radius <= 1.0) "radius must be in (0,1]"

    spl = Point{N}[x0, ]

    for i in 1:N
        Δxi ::F64 = if towards == :upper
            radius * (mp.ub[i] - x0[i])
        elseif towards == :lower
            - radius * (x0[i] - mp.lb[i])
        else
            error("invalid `towards` parameter: $toward")
        end

        push!(spl, modify(x0, i, x0[i] + Δxi))
    end
    return spl
end







# ------------------------------------------------------------------------------
"""
    solve(
        mp::MinimizeProblem{N,P,Q} ;

        x0    ::Point{N} = boxcenter(mp),
        radius::F64      = 0.5,

        verbose   ::Bool = false, # print the trace
        showevery ::Int  = 1,     # print the trace every `showevery` iterations

        δ   ::F64 = 1E-4, # tolerance for the equality constraint violation
        R   ::F64 = 1.0 , # penalty factor for the equality constraint violation
        α   ::F64 = 1.0,  # reflection factor, (0,∞)
        γ   ::F64 = 2.0,  # expansion factor, (1,∞)
        ρout::F64 = 0.5,  # outside contraction factor, (0,0.5]
        ρin ::F64 = 0.5,  # inside contraction factor, (0,0.5]
        σ   ::F64 = 0.5,  # shrink factor, (0,1)

        ftol   ::F64 = 1E-5, # tol for the function value change at centroids
        xtol   ::F64 = 1E-5, # tol for the max simplex edge length/size
        maxiter::Int = 1000, # maximum number of iterations
    ) where {N,P,Q}

Solves the minimization problem `mp` using the Constrained Simplex Search
algorithm. The algorithm is a variant of the Nelder-Mead simplex search
algorithm that can handle box constraints and non-linear constraints.

## Reference
Vivek Kumar Mehta & Bhaskar Dasgupta (2012) A constrained optimization algorithm
based on the simplex search method, Engineering Optimization, 44:5, 537-550

DOI: 10.1080/0305215X.2011.598520

## Example
```julia
css = include("src/ConstrainedSimplexSearch.jl")

mp = css.MinimizeProblem(
    x -> sum(x .^ 2), 
    x -> [], 
    x -> [], 
    3, 0, 0
)

res = css.solve(
    mp,

    x0 = css.boxcenter(mp),
    radius = 0.5,

    verbose = true,
    showevery = 1,

    maxiter = 10,
)


```
"""
function solve(
    mp::MinimizeProblem{N,P,Q} ;

    x0    ::AbstractVector = boxcenter(mp),
    radius::F64      = 0.5,

    verbose   ::Bool = false, # print the progress
    showevery ::Int  = 1,     # print the progress every `showevery` iterations

    δ   ::F64 = 1E-4, # tolerance for the equality constraint violation
    R   ::F64 = 2.0 , # penalty factor for the equality constraint violation
    α   ::F64 = 1.0,  # reflection factor, (0,∞)
    γ   ::F64 = 1.5,  # expansion factor, (1,∞)
    ρout::F64 = 0.4,  # outside contraction factor, (0,0.5]
    ρin ::F64 = 0.4,  # inside contraction factor, (0,0.5]
    σ   ::F64 = 0.5,  # shrink factor, (0,1)

    ftol   ::F64 = 1E-5, # tol for the function value gap at vertexes
    xtol   ::F64 = 1E-5, # tol for the max simplex edge length/size
    maxiter::Int = 1000, # maximum number of iterations
) where {N,P,Q}

    @assert showevery > 0 "showevery must be positive but got: $showevery"

    @assert 0.0 < δ < Inf "δ must be positive & finite"
    @assert 0.0 < R < Inf "R must be positive & finite"
    @assert 0.0 < α < Inf "must be finite positive: α, reflection coef"
    @assert 1.0 < γ < Inf "must be finite positive: γ, expansion coef"
    @assert 0.0 < ρout <= 0.5 "must be in (0,0.5]: ρout, out contraction"
    @assert 0.0 < ρin <= 0.5 "must be in (0,0.5]: ρin, in contraction"
    @assert 0.0 < σ < 1.0 "must be in (0,1): σ, shrink coefficient"

    @assert 0.0 < ftol < Inf "ftol must be positive & finite"
    @assert 0.0 < xtol < Inf "xtol must be positive & finite"
    @assert maxiter > 0 "maxiter must be positive"

    
    # step: initialize a simplex
    spl = initsimplex(Point{N}(x0), mp, radius=radius) ::Vec{Point{N}}
    @assert length(spl) == N+1 "invalid simplex initialization: $(length(spl))"

    # malloc: iteration trace indicators
    lastRes = Dict(
        :fo       => Inf,
        :edgeLen  => maxedgelen(spl),
        :ftrace   => F64[],
        :centroid_feasibility => false,
    )

    # iteration
    for t in 1:maxiter

        # ----------------------------------------------------------------------
        # Step: eval the objective according to the constraint violations of the
        #       current simplex
        # ----------------------------------------------------------------------
        # eval the constraints at the simplex vertexes
        gvalsAll = V64[mp.g(xi) for xi in spl]
        hvalsAll = V64[mp.h(xi) for xi in spl]

        # check feasibility of all vertexes, determine the status of the simplex
        feasAll = Bool[
            isfeasible(xi, gvals, hvals, mp.lb, mp.ub, δ)
            for (xi, gvals, hvals) in zip(spl, gvalsAll, hvalsAll)
        ]
        splLoc ::Symbol = if all(feasAll)
            :interior
        elseif any(feasAll)
            :boundary
        elseif !any(feasAll)
            :exterior
        else
            error("unreachable code: simplex location. what happened?")
        end

        # eval the functions at all feasible simplex vertexes (in case of the
        # boundary simplex
        fvalsAll = fill(NaN, N+1)
        if splLoc == :interior
            # the same as the standard Nelder-Mead simplex search
            for i in 1:N+1
                fvalsAll[i] = mp.f(spl[i])
            end
        elseif splLoc == :exterior
            # fill all vertexes with the CV value
            for i in 1:N+1
                fvalsAll[i] = cv(
                    spl[i],
                    gvalsAll[i],
                    hvalsAll[i],
                    mp.lb,
                    mp.ub,
                    δ,
                    R,
                )
            end
        elseif splLoc == :boundary
            # first, fill the feasible vertexes with the function value
            fmax = -Inf
            for i in 1:N+1
                if feasAll[i]
                    fvalsAll[i] = mp.f(spl[i])
                    fmax = max(fmax, fvalsAll[i])
                end
            end
            # then, fill the infeasible vertexes with fmax + CV values
            for i in 1:N+1
                if !feasAll[i]
                    fvalsAll[i] = fmax + cv(
                        spl[i],
                        gvalsAll[i],
                        hvalsAll[i],
                        mp.lb,
                        mp.ub,
                        δ,
                        R,
                    )
                end
            end # i
        end


        # ----------------------------------------------------------------------        
        # sort the simplex vertexes, and determine the best, worst, 2nd worst
        isort = sortperm(fvalsAll)
        xw  = spl[isort[end]]   # worst

        fb  = fvalsAll[isort[1]]
        fw  = fvalsAll[isort[end]]
        f2w = fvalsAll[isort[end-1]]

        # compute the centroid of the simplex except the worst vertex
        xo = centroid(spl[isort[1:end-1]])
        fo = mp.f(xo)

        # ----------------------------------------------------------------------
        # Choose one move from: reflection, expansion, contraction, shrinkage

        # must-do: find the reflection point of the worst
        xr = reflect(xo, xw, α)
        fr = mp.f(xr)

        # branching
        if fb <= fr < f2w
            # accept the reflection point (to replace the worst)
            spl[isort[end]] = xr
        elseif fr < fb
            # get: expansion point (try more radical move)
            xe = expand(xo, xr, γ)
            fe = mp.f(xe)
            if fe < fr
                # accept the expansion point
                spl[isort[end]] = xe
            else
                # accept the reflection point (no radical move but still good)
                spl[isort[end]] = xr
            end
        elseif f2w <= fr <= fw
            # get: OUTSIDE contraction point
            xc_out = contract_out(xo, xr, ρout)
            fc_out = mp.f(xc_out)
            if fc_out < fw
                # accept the outside contraction point
                spl[isort[end]] = xc_out
            else
                # shrink the simplex towards the best
                shrink!(spl, isort[1], σ)
            end
        elseif fr > fw
            # get: INSIDE contraction point
            xc_in = contract_in(xo, xw, ρin)
            fc_in = mp.f(xc_in)
            if fc_in < fw
                # accept the inside contraction point
                spl[isort[end]] = xc_in
            else
                # shrink the simplex towards the best
                shrink!(spl, isort[1], σ)
            end
        else
            error("unreachable code. what happened?")
        end # if branching

        # rebound the simplex vertexes to the box boundaries
        rebound!(spl, mp.lb, mp.ub)
        
        # error statistics
        maxEdgeLen = maxedgelen(spl)
        maxFvalGap = abs(fw - fb)
        lastRes[:fo]      = fo
        lastRes[:edgeLen] = maxEdgeLen
        push!(lastRes[:ftrace], fo)

        if verbose & (t % showevery == 0)
            @printf(
                "iter = %d, status = %s, f ∈ (%.2e, %.2e), f(centroid) = %.4e",
                t, splLoc, fb, fw, fo
            )
            @printf(
                ", |f(worst)-f(best)| = %.2e, max(edge) = %.2e\n",
                maxFvalGap, maxEdgeLen
            )
        end

        status_centroid = isfeasible(xo, mp, δ)
        lastRes[:centroid_feasibility] = status_centroid

        # convergence check
        if (maxFvalGap < ftol) | (maxEdgeLen < xtol)
            if verbose
                @printf(
                    "Converged: iter = %d, status = %s, centroid feasible = %s, ",
                    t, splLoc, status_centroid,
                )
                @printf(
                    "f ∈ (%.2e, %.2e), f(centroid) = %.4e, max(edge) = %.2e\n",
                    fb, fw, fo, maxEdgeLen
                )
                status_centroid || @warn "centroid is infeasible"
            end
            break
        elseif t == maxiter
            if verbose
                @printf(
                    "reached maxiter: iter = %d, status = %s, centroid feasible = %s, ",
                    t, splLoc, status_centroid,
                )
                @printf(
                    "f ∈ (%.2e, %.2e), f(centroid) = %.4e, max(edge) = %.2e\n",
                    fb, fw, fo, maxEdgeLen
                )
                status_centroid || @warn "centroid is infeasible"
            end
            break
        else
            continue
        end # if


    end # t

    return (
        x = centroid(spl),
        f = lastRes[:fo],
        ftrace = lastRes[:ftrace],
        centroid_feasibility = lastRes[:centroid_feasibility],
    )
end # solve








#===============================================================================
if x_best <= xr < x_2nd_best
    # accept the reflection point (to replace the worst)
    x_worst = xr
elseif xr < x_best
    # get: expansion point (try more radical move)
    xe = expand(xc, xr, γ)
    if xe < xr
        # accept the expansion point
        x_worst = xe
    else
        # accept the reflection point (ok, no radical move and still good)
        x_worst = xr
    end
elseif x_2nd_best <= xr <= x_worst
    # get: OUTSIDE contraction point
    xc_out = contract_out(xc, xr, β)
    if xc_out < x_worst
        # accept the outside contraction point
        x_worst = xc_out
    else
        # shrink the simplex towards the best
        shrink(spl, x_best)
    end
elseif xr > x_worst
    # get: INSIDE contraction point
    xc_in = contract_in(xc, x_worst, β)
    if xc_in < x_worst
        # accept the inside contraction point
        x_worst = xc_in
    else
        # shrink the simplex towards the best
        shrink(spl, x_best)
    end
else
    error("unreachable code. what happened?")
end


Reference:
https://alexdowad.github.io/visualizing-nelder-mead/
===============================================================================#