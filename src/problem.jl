# Defining the problem instance




# ------------------------------------------------------------------------------
"""
    MinimizeProblem{N,P,Q}

Minimization problem definition:

min_x f(x):R^N -> R

s.t. 

g_i(x) <= 0, i = 1,...,P

h_j(x) = 0, j = 1,...,Q

lb_k <= x_k <= ub_k, k = 1,...,N

where:
- `f` is the objective function
- `g` is the vector of inequality constraint functions
- `h` is the vector of equality constraint functions
- `lb` is the lower bound vector
- `ub` is the upper bound vector

## Constraint penalty parameters
- `δ` is the tolerance for the equality constraint violation
- `R` is the penalty factor for the equality constraint violation

## Simplex search parameters
- `α` is the reflection factor, (0,∞)
- `γ` is the expansion factor, (1,∞)
- `ρout` is the outside contraction factor, (0,0.5]
- `ρin` is the inside contraction factor, (0,0.5]
- `σ` is the shrink factor, (0,1)


## Constructor

    MinimizeProblem(
        f, g, h, n, p, q; 
        lb=, ub=, δ=, R=, α=, γ=, ρout=, ρin=, σ=
    )

Create a new MinimizeProblem instance.

### Arguments
- `f::Function`: the objective function, receives a N-dim abstract vector &
returns a scalar
- `g::Function`: the inequality constraint function, receives a N-dim 
abstract vector and returns a P-dim abstract vector
- `h::Function`: the equality constraint function, receives a N-dim abstract
vector and returns a Q-dim abstract vector
- `n::Int`: the dimensionality of the problem, N
- `p::Int`: the number of inequality constraints, P
- `q::Int`: the number of equality constraints, Q

### Optional Arguments
- `lb::AbsV`: the lower bound vector, default is zeros(n)
- `ub::AbsV`: the upper bound vector, default is one(n)
- `δ::Real`: the tolerance for the equality constraint violation, default is
1E-5
- `R::Real`: the penalty factor for the equality constraint violation, 
default is 1.0
- `α::Real`: the reflection factor, (0,∞), default is 1.0
- `γ::Real`: the expansion factor, (1,∞), default is 2.0
- `ρout::Real`: the outside contraction factor, (0,0.5], default is 0.5
- `ρin::Real`: the inside contraction factor, (0,0.5], default is 0.5
- `σ::Real`: the shrink factor, (0,1), default is 0.5



## Example
```julia
css = include("src/ConstrainedSimplexSearch.jl")

# create an unconstrained minimization problem
mp = css.MinimizeProblem(
    x -> sum(x .^ 2), 
    x -> [], 
    x -> [], 
    3, 0, 0
)




```
"""
mutable struct MinimizeProblem{N,P,Q}
    f ::Function # R^N -> R
    g ::Function # R^N -> R^P
    h ::Function # R^N -> R^Q
    lb::Point{N} # control lower bound
    ub::Point{N} # control upper bound

    # --------------------------------------------------------------------------
    function MinimizeProblem(
        f::Function, # min f(x), x ∈ R^N
        g::Function, # g(x) <= 0, i = 1,...,P
        h::Function, # h(x) == 0, j = 1,...,Q
        n::Int, # user-report dimensionality of the problem
        p::Int, # number of inequality constraints
        q::Int; # number of equality constraints

        lb ::AbsV = zeros(n),
        ub ::AbsV = ones(n),
    )
        @assert n > 0 "dimensionality of the problem must be positive: $n"
        @assert n == length(lb) == length(ub) "dimensionality mismatch: lb,ub"
        @assert p >= 0 "# of inequality constraints must be non-negative: $p"
        @assert q >= 0 "# of equality constraints must be non-negative: $q"
        @assert all(isfinite, lb) "lb must be finite"
        @assert all(isfinite, ub) "ub must be finite"
        @assert all(lb .<= ub) "lb must be less than or equal to ub"

        new{n,p,q}( f, g, h, Point{n}(lb), Point{n}(ub) )
    end
end




# ------------------------------------------------------------------------------
function Base.show(io::IO, mp::MinimizeProblem{N,P,Q}) where {N,P,Q}
    println(io, "MinimizeProblem{N=$N,P=$P,Q=$Q}")
    println(io, " min_{x ∈ R^$N} f(x)")
    println(io, " s.t.")
    println(io, " g(x) <= 0, i = 1,...,$P")
    println(io, " h(x) == 0, j = 1,...,$Q")
    println(io, " lb <= x <= ub")
end