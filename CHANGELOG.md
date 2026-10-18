# Changelog

All notable changes to CorePeriphery.jl are documented here. The project follows
semantic versioning within the Julia package ecosystem.

## [0.2.0] - 2026-07-13

### Scientific correctness

- Implemented Lip's published exact degree-prefix loss minimizer.
- Implemented Cucuringu LowRank-Core rank-two reconstruction and Find-Cut.
- Replaced approximate Surprise scoring with the exact multivariate-hypergeometric
  joint-tail statistic in log space.
- Corrected shared objective, degeneracy, dyad, self-loop, symmetry, and binary input
  contracts.
- Reworked Kojaku--Masuda objectives and move gains with consistent ER and
  configuration-null conventions.
- Added exhaustive and direct mathematical validation oracles.

### Della Rossa profiling

- Implemented exact equation-(5) persistence and equation-(8) cp-centralization.
- Added the paper's stochastic initial and final tie choices through
  `tie_break=:paper` and an explicit RNG.
- Added strongly connected directed profiling using stationary random-walk flow and
  total in-plus-out strength.
- Added lazy, direct linear, and automatic stationary solvers, validated supplied
  stationary vectors, and solver residual/iteration diagnostics.
- Added reproducible ensemble uncertainty summaries through `rossa_profile_ensemble`.

### New features

- Added `CPDirectedResult`, `minres_svd_directed`, and the distinct
  `minres_symmetric` estimator.
- Added penalized pair-count model selection and candidate diagnostics to
  `CPMultiResult`.
- Added Monte Carlo significance with ER, undirected degree-configuration, directed
  in/out-degree configuration, and weight-permutation null models.
- Added explicit configuration-switch completion diagnostics and strict shortfall
  handling.
- Added optional Graphs.jl dispatch through a package extension.

### Performance

- Added sufficient-statistic Pearson kernels, cached Rombach and Surprise updates,
  incremental Rossa boundary flows, adjacency-list multi-pair moves, analytic MINRES
  residuals, and sparse partial eigenpairs.
- Added deterministic threaded significance scheduling with reusable buffers.
- Added allocation and scaling gates for Julia 1.10 and 1.12 and a versioned benchmark
  report generator.

### Documentation and validation

- Added algorithm/result contracts, development record, migration guide, scientific
  validation matrix, cross-package comparison, causal divergence analysis, performance
  ablation, and optimization reports.
- Compared nine common configurations with cpnet commit
  `6aad458a6d434a3617d33e74f7163d514a27fecb` on identical planted networks.

### Breaking and compatibility changes

- Raised the Julia compatibility floor from 1.6 to 1.10.
- Algorithm-specific quality values and several optimizer results differ from 0.1
  because published statistics and objectives replaced earlier approximations.
- Configuration-null significance now errors when requested switching is incomplete
  unless the caller explicitly selects `swap_shortfall=:warn` or `:accept`.
- `CPMultiResult` and `CPSignificanceResult` contain additional diagnostics fields.
- Rossa `n_walks` and `walk_length` keywords remain accepted but are deprecated and
  ignored; profiling is analytic rather than simulated.

[0.2.0]: https://github.com/simoneSantoni/CorePeriphery.jl/compare/a2226832508c57efdeab71cde1f055549ca028e7...v0.2.0
