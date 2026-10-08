using Test
using Optim
import ConstrainedSimplexSearch as css


# This testset checks constructor validation, field storage, and independent
# copies of the documented default hyperparameters.
@testset "Problem construction and defaults" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 1,
        function_objective = c -> sum(abs2, c),
        function_constraints = c -> [c[1] + c[2] - 1.0],
        lower_bounds = [-1.0, -2.0],
        upper_bounds = [2.0, 3.0],
    )

    @test prob.ncontrols == 2
    @test prob.nconstraints == 1
    @test prob.lower_bounds == [-1.0, -2.0]
    @test prob.upper_bounds == [2.0, 3.0]
    @test prob.hyperparam == css.DEFAULT_HYPERPARAM
    @test prob.hyperparam !== css.DEFAULT_HYPERPARAM
    @test_throws ArgumentError css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 0,
        function_objective = c -> sum(c),
        function_constraints = c -> Float64[],
        lower_bounds = [0.0, 0.0],
        upper_bounds = [1.0, 0.0],
    )
end


# This testset verifies that a centered strategy builds row-wise vertices at the
# requested local scale and produces nondegenerate simplex geometry.
@testset "Centered simplex" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 0,
        function_objective = c -> sum(abs2, c),
        function_constraints = c -> Float64[],
        lower_bounds = [-2.0, -4.0],
        upper_bounds = [2.0, 4.0],
    )
    strategy = css.CenteredSimplex(c0 = [0.0, 0.0], radius = 0.5, towards = :upper)
    vertices = css.build(strategy, prob)

    @test size(vertices) == (3, 2)
    @test vertices == [0.0 0.0; 1.0 0.0; 0.0 2.0]
    @test css.simplexvolume(vertices) > 0.0
end


# This testset checks that the broad strategy remains in the box and spans every
# coordinate direction by a positive amount.
@testset "Maximum-volume simplex" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 3,
        nconstraints = 0,
        function_objective = c -> sum(abs2, c),
        function_constraints = c -> Float64[],
        lower_bounds = [-2.0, -1.0, 0.0],
        upper_bounds = [1.0, 3.0, 5.0],
    )
    vertices = css.build(css.MaxVolumeSimplex(radius = 0.7), prob)

    @test size(vertices) == (4, 3)
    @test all(prob.lower_bounds' .<= vertices .<= prob.upper_bounds')
    @test css.simplexvolume(vertices) > 0.0
end


# This testset ensures explicit vertices reject wrong dimensions and degenerate
# geometry before any objective or constraint callback can run.
@testset "Explicit simplex validation" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 0,
        function_objective = c -> sum(abs2, c),
        function_constraints = c -> Float64[],
        lower_bounds = [-1.0, -1.0],
        upper_bounds = [1.0, 1.0],
    )

    vertices = [0.0 0.0; 0.5 0.0; 0.0 0.5]
    strategy = css.ExplicitSimplex(vertices)
    vertices[1, 1] = 0.75
    @test strategy.vertices[1, 1] == 0.0

    @test_throws DimensionMismatch css.optimize(
        prob,
        css.ExplicitSimplex([0.0 0.0; 1.0 0.0]),
    )
    @test_throws ArgumentError css.optimize(
        prob,
        css.ExplicitSimplex([0.0 0.0; 0.5 0.0; 1.0 0.0]),
    )
end


# This testset solves the package's illustrative nonlinear corner problem and
# checks both numerical accuracy and the separate admissibility report.
@testset "Known corner solution" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 2,
        function_objective = c -> c[1]^2 + c[2]^2,
        function_constraints = c -> [-c[1], c[1]^2 - c[2]],
        lower_bounds = [-2.0, -3.0],
        upper_bounds = [1.0, 5.0],
    )
    result = css.optimize(
        prob,
        css.MaxVolumeSimplex(radius = 0.7);
        max_iteration = 2000,
        objective_tolerance = 1E-10,
        control_tolerance = 1E-8,
    )

    @test result["admissibility_status"] == "admissible"
    @test result["optimal_control"] ≈ [0.0, 0.0] atol = 2E-3
    @test result["optimal_objective"] ≈ 0.0 atol = 1E-5
end


# This testset confirms ordinary objective ordering after the simplex has entered
# an admissible region with a known interior minimizer.
@testset "Known interior solution" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 1,
        function_objective = c -> (c[1] - 0.25)^2 + (c[2] + 0.5)^2,
        function_constraints = c -> [c[1] + c[2] - 1.0],
        lower_bounds = [-2.0, -2.0],
        upper_bounds = [2.0, 2.0],
    )
    result = css.optimize(
        prob,
        css.CenteredSimplex(c0 = [0.0, 0.0], radius = 0.5);
        objective_tolerance = 1E-10,
        control_tolerance = 1E-8,
    )

    @test result["optimal_control"] ≈ [0.25, -0.5] atol = 2E-3
    @test result["admissibility_status"] == "admissible"
end


# This testset checks that the search can stop on an active linear inequality
# instead of requiring an artificial interior offset.
@testset "Known boundary solution" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 2,
        nconstraints = 1,
        function_objective = c -> (c[1] - 1.0)^2 + (c[2] - 1.0)^2,
        function_constraints = c -> [c[1] + c[2] - 1.0],
        lower_bounds = [-1.0, -1.0],
        upper_bounds = [2.0, 2.0],
    )
    result = css.optimize(
        prob,
        css.CenteredSimplex(c0 = [0.0, 0.0], radius = 0.5);
        max_iteration = 2000,
        objective_tolerance = 1E-10,
        control_tolerance = 1E-8,
    )

    @test sum(result["optimal_control"]) <= 1.0 + 1E-8
    @test result["optimal_control"] ≈ [0.5, 0.5] atol = 3E-3
end


# This testset ensures one-dimensional callbacks always receive an ordinary
# Vector{Float64} of length one throughout constraints and objective evaluation.
@testset "One-dimensional callback types" begin
    objective_calls = Ref(0)
    constraint_calls = Ref(0)
    objective = function(c)
        @test c isa Vector{Float64}
        @test length(c) == 1
        objective_calls[] += 1
        return (c[1] - 0.2)^2
    end
    constraints = function(c)
        @test c isa Vector{Float64}
        @test length(c) == 1
        constraint_calls[] += 1
        return Float64[]
    end
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 1,
        nconstraints = 0,
        function_objective = objective,
        function_constraints = constraints,
        lower_bounds = [-1.0],
        upper_bounds = [1.0],
    )
    result = css.optimize(prob, css.MaxVolumeSimplex(); control_tolerance = 1E-8)

    @test result["optimal_control"] ≈ [0.2] atol = 2E-3
    @test objective_calls[] > 0
    @test constraint_calls[] >= objective_calls[]
end


# This is the central safety regression test: the objective throws outside the
# admissible half-line, while the initial simplex deliberately contains an
# invalid point. A successful solve proves all candidate branches are guarded.
@testset "Objective evaluation guard" begin
    objective = function(c)
        c[1] >= 0.0 || error("objective evaluated outside the admissible set")
        return 1.0E9 + (c[1] - 0.25)^2
    end
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 1,
        nconstraints = 1,
        function_objective = objective,
        function_constraints = c -> [-c[1]],
        lower_bounds = [-1.0],
        upper_bounds = [1.0],
    )
    result = css.optimize(
        prob,
        css.MaxVolumeSimplex(radius = 0.7);
        max_iteration = 1000,
        objective_tolerance = 1E-10,
        control_tolerance = 1E-8,
    )

    @test result["admissibility_status"] == "admissible"
    @test result["optimal_control"] ≈ [0.25] atol = 2E-3
    @test result["optimal_objective"] ≈ 1.0E9 atol = 1E-5
end


# This testset verifies sign handling for maximization while confirming the
# returned objective remains in the user's original, non-negated convention.
@testset "Maximization" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 1,
        nconstraints = 0,
        function_objective = c -> -(c[1] - 0.3)^2 + 2.0,
        function_constraints = c -> Float64[],
        lower_bounds = [-1.0],
        upper_bounds = [1.0],
    )
    result = css.optimize(
        prob,
        css.MaxVolumeSimplex();
        use_maximize = true,
        objective_tolerance = 1E-10,
        control_tolerance = 1E-8,
    )

    @test result["optimal_control"] ≈ [0.3] atol = 2E-3
    @test result["optimal_objective"] ≈ 2.0 atol = 1E-5
end


# This testset exercises threaded vertex evaluation and the explicit rejection
# of process mode, whose callback transport is not implemented in v0.2.0.
@testset "Parallel mode contract" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 1,
        nconstraints = 0,
        function_objective = c -> c[1]^2,
        function_constraints = c -> Float64[],
        lower_bounds = [-1.0],
        upper_bounds = [1.0],
    )
    threaded = css.optimize(
        prob,
        css.MaxVolumeSimplex();
        parallel_mode = :thread,
        max_workers = 1,
    )

    @test threaded["parallel_mode"] == :thread
    @test threaded["max_workers"] == 1
    @test_throws ArgumentError css.optimize(
        prob,
        css.MaxVolumeSimplex();
        parallel_mode = :process,
    )
end


# This testset verifies the optional Optim frontend delegates to the native
# solver and exposes the familiar result accessors expected by Optim users.
@testset "Optim integration" begin
    objective = c -> (c[1] - 0.4)^2
    constraints = c -> [-c[1]]
    options = Optim.Options(
        iterations = 1000,
        f_abstol = 1E-10,
        x_abstol = 1E-8,
        show_trace = false,
    )
    result = css.optimize(
        objective,
        [-1.0],
        [1.0],
        css.ConstrainedSimplex();
        inequality_constraints = constraints,
        initial_simplex = css.MaxVolumeSimplex(),
        options = options,
    )

    @test Optim.minimizer(result) ≈ [0.4] atol = 2E-3
    @test Optim.minimum(result) ≈ 0.0 atol = 1E-5
    @test Optim.converged(result)
    @test Optim.iterations(result) > 0
end


# This lightweight calibration test checks the transformation ranges, returned
# diagnostics, and assignment of the tuned dictionary without demanding a
# particular optimizer path through the nonsmooth benchmark loss.
@testset "Hyperparameter fine tuning" begin
    prob = css.ConstrainedSimplexSearch(
        ncontrols = 1,
        nconstraints = 0,
        function_objective = c -> (c[1] - 0.2)^2,
        function_constraints = c -> Float64[],
        lower_bounds = [-1.0],
        upper_bounds = [1.0],
    )
    tuned, diagnostics = css.fine_tune(
        prob,
        css.MaxVolumeSimplex();
        true_control = [0.2],
        max_iteration = 8,
        tolerance = 0.1,
        method = Optim.NelderMead(),
    )

    @test prob.hyperparam == tuned
    @test Set(keys(tuned)) == Set(keys(css.DEFAULT_HYPERPARAM))
    @test tuned["reflection_factor"] > 0.0
    @test tuned["expansion_factor"] > 1.0
    @test 0.0 < tuned["contraction_factor_outside"] < 0.5
    @test 0.0 < tuned["contraction_factor_inside"] < 0.5
    @test 0.0 < tuned["shrink_factor"] < 1.0
    @test haskey(diagnostics, "minimized_error")
    @test haskey(diagnostics, "admissible")
end
