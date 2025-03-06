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
    towards::Symbol  = :upper, # or :lower, the direction of simplex init

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
    spl ::Vec{Point{N}} = initsimplex(
        Point{N}(x0), 
        mp, 
        radius  = radius, 
        towards = towards,
    )

    # malloc: intermediate results
    fvalsAll = fill(Inf, N+1)
    gvalsAll = [zeros(P) for _ in 1:N+1]
    hvalsAll = [zeros(Q) for _ in 1:N+1]
    feasAll  = fill(false, N+1)
    splLoc   = :interior
    isort    = Vector{Int}(undef, N+1)

    flagConverged = false

    # iteration
    for t in 1:maxiter

        # ----------------------------------------------------------------------
        # Step: eval the objective according to the constraint violations of the
        #       current simplex
        # ----------------------------------------------------------------------
        # eval the constraints at the simplex vertexes
        for i in 1:N+1
            gvalsAll[i] .= mp.g(spl[i]) ::AbstractVector
            hvalsAll[i] .= mp.h(spl[i]) ::AbstractVector

            # check admissibility of the vertex
            feasAll[i] = isadmissible(
                spl[i], 
                gvalsAll[i], 
                hvalsAll[i], 
                mp.lb, mp.ub, 
                δ
            )
        end
        splLoc = if all(feasAll)
            :interior
        elseif any(feasAll)
            :boundary
        elseif !any(feasAll)
            :exterior
        else
            error("unreachable code: simplex location. what happened?")
        end

        # eval the functions at all feasible simplex vertexes (in case of the
        # boundary simplex)
        fmaxAdmissible = NaN

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
            # finally, save this value for the reflection point
            fmaxAdmissible = fmax
        end


        # ----------------------------------------------------------------------        
        # sort the simplex vertexes, and determine the best, worst, 2nd worst
        isort .= sortperm(fvalsAll)

        xw  = spl[isort[end]]   # worst

        fb  = fvalsAll[isort[1]]     # best
        fw  = fvalsAll[isort[end]]   # worst
        f2w = fvalsAll[isort[end-1]] # 2nd worst

        @assert all(!isnan, [fb,fw,f2w]) "NaN found in the function values"


        # compute the centroid of the simplex except the worst vertex
        xo = centroid(spl[isort[1:end-1]])

        # ----------------------------------------------------------------------
        # Choose one move from: reflection, expansion, contraction, shrinkage

        # must-do: find the reflection point of the worst
        # must-do: restrict the reflection point to the box boundaries, to avoid
        #          the potential undefined behvaior of the objective function.
        # caution: the reflection point should respect the location status of
        #          the current simplex (interior, exterior, or boundary)
        xr = clamp.(reflect(xo, xw, α), mp.lb, mp.ub)
        fr = if splLoc == :interior
            mp.f(xr)
        elseif splLoc == :exterior
            cv(
                xr,
                mp.g(xr),
                mp.h(xr),
                mp.lb,
                mp.ub,
                δ,
                R,
            )
        elseif splLoc == :boundary
            _gval = mp.g(xr)
            _hval = mp.h(xr)
            _flag = isadmissible(xr, _gval, _hval, mp.lb, mp.ub, δ)
            if _flag
                mp.f(xr)
            else
                fmaxAdmissible + cv(xr, _gval, _hval, mp.lb, mp.ub, δ, R)
            end
        else
            error("unreachable code: simplex location. what happened?")
        end

        # logging the chosen move this iteration
        chosenMove = :reflection

        # branching
        if fb <= fr < f2w
            # accept the reflection point (to replace the worst)
            spl[isort[end]] = xr

            chosenMove = :reflection
        elseif fr < fb
            # get: expansion point (try more radical move)
            xe = clamp.(expand(xo, xr, γ), mp.lb, mp.ub)
            fe = mp.f(xe)
            if fe < fr
                # accept the expansion point
                spl[isort[end]] = xe

                chosenMove = :expansion
            else
                # accept the reflection point (no radical move but still good)
                spl[isort[end]] = xr

                chosenMove = :reflection
            end
        elseif f2w <= fr <= fw
            # get: OUTSIDE contraction point
            xc_out = clamp.(contract_out(xo, xr, ρout), mp.lb, mp.ub)
            fc_out = mp.f(xc_out)
            if fc_out < fw
                # accept the outside contraction point
                spl[isort[end]] = xc_out

                chosenMove = :outside_contraction
            else
                # shrink the simplex towards the best
                shrink!(spl, isort[1], σ)

                # rebound the simplex vertexes to the box boundaries (if shrink)
                rebound!(spl, mp.lb, mp.ub)

                chosenMove = :shrink
            end
        elseif fr > fw
            # get: INSIDE contraction point
            xc_in = clamp.(contract_in(xo, xw, ρin), mp.lb, mp.ub)
            fc_in = mp.f(xc_in)
            if fc_in < fw
                # accept the inside contraction point
                spl[isort[end]] = xc_in

                chosenMove = :inside_contraction
            else
                # shrink the simplex towards the best
                shrink!(spl, isort[1], σ)

                rebound!(spl, mp.lb, mp.ub)

                chosenMove = :shrink
            end
        else
            @info "Diagnostic:" f_reflection=fr f_best=fb f_worst=fw
            error("unreachable code: move selection. what happened?")
        end # if branching
        
        # error statistics
        maxEdgeLen = maxedgelen(spl) # wrt: xtol
        maxFvalGap = abs(fw - fb)    # wrt: ftol

        if verbose & (t % showevery == 0)
            @printf(
                "iter = %d, simplex status = %s, f ∈ (%.2e, %.2e)",
                t, splLoc, fb, fw
            )
            @printf(
                ", gap(f) = %.2e, max|edge| = %.2e, move: %s\n",
                maxFvalGap, maxEdgeLen, chosenMove
            )
        end

        # convergence check
        if (maxFvalGap < ftol) | (maxEdgeLen < xtol)
            flagConverged = true
            if verbose
                @printf(
                    "Converged: iter = %d, simplex status = %s, ",
                    t, splLoc,
                )
                @printf(
                    "f ∈ (%.2e, %.2e), max|edge| = %.2e\n",
                    fb, fw, maxEdgeLen
                )
            end
            break
        elseif t == maxiter
            if verbose
                @printf(
                    "reached maxiter: iter = %d, simplex status = %s, ",
                    t, splLoc, 
                )
                @printf(
                    "f ∈ (%.2e, %.2e), max(edge) = %.2e\n",
                    fb, fw, maxEdgeLen
                )
            end
            break
        else
            continue
        end # if


    end # t


    # post-iter: decide which vertex to return
    #=
    rule: check all N+1 vertexes. return the best admissible one among them.
    =#
    xOpt, fOpt = if any(feasAll)
        # if any feasible vertex or centroid, return the best admissible one
        _f, _i = findmin(fvalsAll[feasAll])
        spl[feasAll][_i], _f
    else
        # if no admissible vertex or centroid, return the centroid with a NaN
        # (because in this exterior status, sthe objective function is not
        # defined and its value was not stored, so nothing to return)
        centroid(spl), NaN
    end

    return (
        x = xOpt,
        f = fOpt,
        converged  = flagConverged,
        admissible = any(feasAll),
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