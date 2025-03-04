using Test
import ConstrainedSimplexSearch as css


# Test: defining the problem
@test try
    
    # Objective function: receives an N-vector and returns a scalar
    function f(x::AbstractVector)
        return sum(x .^ 2)
    end

    # Inequality constraint g(x): receives an N-vector, returns a P-scalar
    function g(x::AbstractVector)
        return [
            - x[1]       , # x >= 0
            x[1]^2 - x[2], # y >= x^2
        ]
    end

    # Equality constraint h(x): receives an N-vector, returns a Q-vector
    function h(x::AbstractVector)
        return []
    end

    # Define a box/rectangular region for searching
    lb = [-2.0, -3.0]
    ub = [1.0, 5.0]

    # Define the minimization problem
    mp = css.MinimizeProblem(
        f, g, h,
        2, 2, 0, # (N,P,Q), report the dimensionalities of the problem
        lb = lb,
        ub = ub,
    )

    true
catch
    false
end



# Test: solving the problem
@test try

    # Objective function: receives an N-vector and returns a scalar
    function f(x::AbstractVector)
        return sum(x .^ 2)
    end

    # Inequality constraint g(x): receives an N-vector, returns a P-scalar
    function g(x::AbstractVector)
        return [
            - x[1]       , # x >= 0
            x[1]^2 - x[2], # y >= x^2
        ]
    end

    # Equality constraint h(x): receives an N-vector, returns a Q-vector
    function h(x::AbstractVector)
        return []
    end

    # Define a box/rectangular region for searching
    lb = [-2.0, -3.0]
    ub = [1.0, 5.0]

    # Define the minimization problem
    mp = css.MinimizeProblem(
        f, g, h,
        2, 2, 0, # (N,P,Q), report the dimensionalities of the problem
        lb = lb,
        ub = ub,
    )

    # Run the algorithm
    @time res = css.solve(
        mp,
        x0 = [-1,-2], # initial interior point to establish the simplex, not necessary to be feasible; `boxcenter()` computes the very mid point of the box constraints (default)
        radius = 0.5, # how large the initial simplex should be, it is the relative distance from x0 to the boundaries of the box. raidus in (0,1)

        verbose = false, # if to print the iteration details
        showevery = 1, # every how many iterations to print if verbose == true

        maxiter = 1000, # max iteration

        δ    = 1E-6, # tolerance for the equality constraint violation
        R    = 10.0 , # penalty factor for the equality constraint violation
        
        α    = 1.0,  # reflection factor, (0,∞)
        γ    = 1.5,  # expansion factor, (1,∞)
        ρout = 0.4,  # outside contraction factor, (0,0.5]
        ρin  = 0.4,  # inside contraction factor, (0,0.5]
        σ    = 0.5,  # shrink factor, (0,1)

        ftol = 1E-3, # tol for the function value gap at vertexes
        xtol = 1E-5, # tol for the max simplex edge length/size
    )

    true
catch
    false
end

