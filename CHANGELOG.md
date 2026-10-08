# Changelog

## [0.2.0]

### Added

- Added the new `ConstrainedSimplexSearch(...)` problem constructor with explicit
  objective, inequality-constraint, and box-bound callbacks.
- Added `CenteredSimplex`, `MaxVolumeSimplex`, and `ExplicitSimplex` initial
  simplex strategies with common geometry validation.
- Added `optimize(prob, simplex0; ...)`, returning a `Dict{String,Any}` with
  convergence, admissibility, trace, timing, settings, and operation-counter
  diagnostics.
- Added optional threaded evaluation of current simplex vertices.
- Added an optional `Optim.jl` frontend with `ConstrainedSimplex()` and the
  familiar `Optim.minimizer`, `Optim.minimum`, `Optim.converged`, and
  `Optim.iterations` accessors.
- Added hyperparameter fine tuning for representative benchmark problems with
  known solutions.

### Changed

- Changed the package API to a modern pipeline.
- Changed solver evaluation so the objective function is called only after a
  candidate point has been confirmed admissible.
- Changed result reporting from the old compact named tuple to a richer
  dictionary-based result.
- Updated package compatibility for Julia 1.10+ and made `Optim.jl` an optional
  extension dependency.
- Rewrote the README motivation, usage guide, example, and API documentation.

### Removed

- Removed equality-constraint support. This is a breaking change in v0.2.0.
- Removed the old `MinimizeProblem(f, g, h, n, p, q; ...)` constructor and
  `solve` workflow.
- Removed the `StaticArrays.jl` dependency from the public data pipeline.

### Fixed

- Fixed candidate evaluation paths that could call the objective at feasible
  but non-admissible expansion or contraction points.
- Fixed one-dimensional callback handling so callbacks consistently receive a
  `Vector{Float64}`.
