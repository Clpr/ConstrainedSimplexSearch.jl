# ConstrainedSimplexSearch.jl

This repository contains the implementation of the Constrained Simplex Search algorithm
which is generalized from the unconstrained Nelder-Mead search method.
The algorithm is designed to solve optimization problems with constraints.

$$
\begin{align}
& \min_{x \in \mathbb{R}^N} f(x) \\
\text{s.t. }& g_i(x) \leq 0, i = 1,\dots,P \\
& h_j(x) = 0, j = 1,\dots,Q \\
& \text{lb}_k \leq x_k \leq \text{ub}_k, k = 1,\dots, N
\end{align}
$$

## Features

- Implements the Constrained Simplex Search algorithm
- Handles generic non-linear constraints using penalty which guides the simplex to move to the feasible region.


## Installation

To install the required dependencies, run:

```bash
pkg> add ConstrainedSimplexSearch

# or directly install from GitHub
```

## Usage

To use the Constrained Simplex Search algorithm, import the necessary modules and call the functions as demonstrated below:

$$
\begin{aligned}
& \min_{x,y} x^2 + y^2 \\
\text{s.t. }& x \geq 0 \\
& y \geq x^2 \\
& -2 \leq x \leq 1 \\
& -3 \leq y \leq 5
\end{aligned}
$$

which has a known corner solution $(0,0)$.

![](asset/feasible_region_with_x0.svg)


```julia
import ConstrainedSimplexSearch as css

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
    x0 = [-1, -2], # within-box point to initialize the simplex

    verbose = true, # if to print the iteration details
    showevery = 1, # every how many iterations to print if verbose == true

    maxiter = 1000, # max iteration

    δ    = 1E-4, # tolerance for the equality constraint violation
    R    = 2.0 , # penalty factor for the equality constraint violation
    
    α    = 1.0,  # reflection factor, (0,∞)
    γ    = 1.5,  # expansion factor, (1,∞)
    ρout = 0.4,  # outside contraction factor, (0,0.5]
    ρin  = 0.4,  # inside contraction factor, (0,0.5]
    σ    = 0.5,  # shrink factor, (0,1)

    ftol = 1E-5, # tol for the function value gap at vertexes
    xtol = 1E-5, # tol for the max simplex edge length/size
)

# what are included in the result (a named tuple)
res.x ::SVector{N,Float64}     # the final simplex's centroid point
res.f ::Float64                # the objective function value at the centroid
res.ftrace::Vector{Float64}    # the trace of f(centroid) by iteration
res.centroid_feasibility::Bool # if the final point is feasible wrt (g, h, lb, ub)

```
The algorithm should converged at iteration 31, at approximately $(0.024,0.028)$, with the objective function value about $0.0014$. One can further improve the solution precision by reducing `ftol` and `xtol`.

**Notes**
1. The keyword argument `x0::AbstractVector=` specifies a point where the simplex is initialized. Such as point must be interior regarding the box constraints `[lb,ub]`, but not necessary to satisfy the inequality constraint $g(x)$ and equality constraints $h(x)$.
    - By default, the center point of the box space (say `(lb + ub)/2`) is applied. One may call `css.boxcenter(mp)` to do this job
    - In the above example, the illustrative `x0` $(-1,-2)$ is within the box $[-2,1]\times [-3,5]$, but does not satisfy the inequality constraints
2. Detailed docstrings are availble by using `?function_name`.
3. Fine tuning is necessary for quick convergence. Try different strategies regarding your problems.



## Pros & Cons

**Pros**
- As robust as Nelder-Mead
- Can handle box constraints and inequality constraints effectively
- Works with high-dimensional approximations (e.g. sparse grid interpolations)
in which the derivatives are oscilliating or unstable across the space
- Very useful in quantitative macroeconomic modeling, esp. working with intra-period portfolio optimization problems.

**Cons**
- Can handle equality constraints but very **_not_** recommended. In general, the simplex-based methods do not work well with equality constraints which in fact "reduce" the dimensionality of the searching space. Please consider substituting out the equality constraints by re-defining the control variables.
- Less efficient if the function evaluations are expensive




## Reference

- Mehta, Vivek Kumar, and Bhaskar Dasgupta. "A constrained optimization algorithm based on the simplex search method." _Engineering Optimization_ 44, no. 5 (2012): 537-550.

- Nealder-Mead implementation: `https://alexdowad.github.io/visualizing-nelder-mead/`


## License

MIT license