# Changelog

All notable changes to CorePeriphery.jl are documented here. The project follows
semantic versioning within the Julia package ecosystem.

## [0.2.0] - Unreleased

### Scientific correctness

- Implemented Lip's published exact degree-prefix loss minimizer.
- Implemented Cucuringu LowRank-Core rank-two reconstruction and Find-Cut.
- Replaced approximate Surprise scoring with the exact multivariate-hypergeometric
  joint-tail statistic in log space.
- Corrected shared objective, degeneracy, dyad, self-loop, symmetry, and binary input
  contracts.
- Reworked Kojaku--Masuda objectives and move gains with consistent ER and
  configuration-null conventions.
- **Corrected the Kojaku--Masuda configuration objective.** The `:configuration`
  objective of `multiple_cp_pairs` omitted the `i = j` null term of Kojaku & Masuda
  (2018) eq. (13), which lowered it by `Σ_{i ∈ core} d_i² / (4m²)` (about 17% on a
  30-node planted test). `_joint_expected` and every move delta now carry the term,
  the objective reduces to modularity when every node is core as the paper states, and
  literal eq. (3)/(13) oracles check both nulls. `benchmark/DIVERGENCE.md` no longer
  attributes this difference to cpnet.
- **Corrected the Rombach template.** Rombach et al. (2017) eq. (4.9) is implemented
  exactly with the `⌊βN⌋` boundary (the former template used ranks `i/n`), `beta` is
  the periphery fraction and accepts the closed interval `[0, 1]`, the core vector is
  normalised to sum to one as in the paper, `quality` is the full double sum `R` of eq.
  (4.10), and the default `partition=:beta` reports the `N − ⌊βN⌋` core ranks.
- Added `rombach_core_score`, the paper's aggregate core score (eq. 4.11) averaged over
  an `(alpha, beta)` grid.
- Added the Borgatti--Everett missing-block discrete model (`cp_block=:missing`, BE 2000
  eq. 4) to `borgatti_everett_discrete`, `label_switching_cp`, `core_quality`, and
  `ideal_cp_matrix`, and directed input (`directed=true`, correlation over all ordered
  dyads as in UCINET) to the Borgatti--Everett family.
- Replaced the `>=` median split of the continuous detectors with a strict `>` rule
  (ties at the median are periphery) and made the rule explicit through
  `partition=:median|:none|:threshold|:core_size` with `threshold`/`core_size`.
- Isolates always receive coreness `0` and are never core (pair label `0` in
  `multiple_cp_pairs`); an edgeless network returns an empty core; detectors warn once
  per session on input with more than one non-trivial connected component.
- Added exhaustive and direct mathematical validation oracles, including a `BigInt`
  Surprise oracle, a brute-force optimum for every graph with `n ≤ 6` under both
  Borgatti--Everett patterns, and a Rombach template oracle.

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
- Added `n_runs`/`rng` restarts and `method=:exhaustive` (Gray-code enumeration of
  every core set, `n ≤ 24`, global optimum) to `borgatti_everett_discrete`.
- Added `block_densities(A, result_or_core_nodes)`, the UCINET-style density table of
  the core-core, core-periphery, periphery-core, and periphery-periphery blocks.
- Added a labelled input path: `adjacency_to_matrix(edges; directed, combine,
  sort_labels)` accepts arbitrary node identifiers and returns `(A, labels)`;
  `adjacency_to_matrix(edges, n)` gained `directed=true` and `combine=:sum|:max|:min|:last`.
- Added `null_mean`, `null_std`, and `zscore` effect sizes to `CPSignificanceResult`,
  made `n_samples=999` and `pass_rng=true` the defaults, and gave non-binary input with
  a binary null an error that names `:weight_permutation`.
- Added the exported `ConvergenceError`, raised by ARPACK in sparse LowRank-Core (with
  the solver settings in the message), by Rossa stationary solvers that exhaust their
  budget, and by incomplete configuration-null switching under `swap_shortfall=:error`
  (the message explains stars and threshold graphs and suggests `:er`,
  `max_swap_attempts`, or `swap_shortfall=:accept`).
- Added structural `==`, `isequal`, `hash`, and compact `show` for every result type
  and the documented `CorePeriphery.NullDiagnostics`.
- Relaxed numeric keyword types to `Real` and `ideal_cp_matrix` to any real vector;
  `Symmetric`/`Hermitian`/`Transpose`/`Adjoint` wrappers of sparse matrices are
  unwrapped at entry; symmetry is tested with tolerance `1e-12 · max(1, max|A|)`.
- `Base.depwarn` for the ignored Rossa `n_walks`/`walk_length` and Lip `max_iter`
  keywords; `rossa_profile_ensemble` errors on forwarded `tie_break`/`rng`;
  `multiple_cp_pairs_config` errors on a forwarded `null_model`.

### Performance

- Added sufficient-statistic Pearson kernels, cached Rombach and Surprise updates,
  incremental Rossa boundary flows, adjacency-list multi-pair moves, analytic MINRES
  residuals, and sparse partial eigenpairs.
- Surprise evaluates the hypergeometric tail from its mode outwards with rigorous
  geometric truncation instead of the `O(m²)` double sum: a 400-node ER graph at
  density 0.05 went from seconds to tens of milliseconds and a planted `n = 512` graph
  from 11 s to 6 ms in development measurements.
- Sparse-aware significance: null samples are drawn into a `_NullWorkspace` holding
  either a dense matrix or a `SparseMatrixCSC`, chosen from the input storage, so
  sparse input stays sparse through the detectors (`n = 2048` ER, 1.0 s to 0.32 s in
  development measurements; identical results to the dense path).
- Sparse LowRank-Core uses ARPACK with a finite tolerance (`1e-10`), `ncv = 48`, and
  `maxiter = max(300, 20n)` instead of `tol = 0` (a 2000-node ring went from 13 s to
  about 3 s in development measurements); the thresholded reconstruction is pruned by
  a per-row bound; ARPACK calls are serialised through a lock because Arpack.jl is not
  reentrant.
- Threaded significance uses a `Channel` work queue (`_run_null_samples!`) with one
  task per thread; `thread_schedule=:auto` maps to `:greedy` on every Julia version and
  the greedy, static, and serial schedules are bit-for-bit identical (a 144,000-run
  stress showed no mismatches).
- Added allocation and scaling gates for Julia 1.10, 1.12, and 1.13 and a versioned
  benchmark report generator; the gate now uses the minimum over 5 (quick) or 7 (full)
  repetitions and reports minimum/median/maximum columns.

### Bug fixes

- Fixed intermittent non-determinism of threaded `:greedy` significance sampling:
  buffer names assigned inside the task closures were shared across tasks; the
  per-task workspace is now created inside each task and declared `local`.
- `cp_significance` with sparse input no longer falls back to dense `n × n` buffers
  per worker.
- `cp_significance` no longer runs the disconnected-input check on null samples. Random
  null graphs are often disconnected, so the check emitted spurious warnings and added
  an `O(n^2)` component search per sample on dense input. The observed network is still
  checked, and isolate handling still applies to every sample.
- `lip_discrete` derives isolates from its degree vector instead of a second full
  matrix scan, restoring dense significance-testing throughput.

### Documentation and validation

- Added algorithm/result contracts, development record, migration guide, scientific
  validation matrix, cross-package comparison, causal divergence analysis, performance
  ablation, and optimization reports.
- Compared nine common configurations with cpnet commit
  `6aad458a6d434a3617d33e74f7163d514a27fecb` on identical planted networks.
- Corrected documentation claims found by an expert review: `multiple_cp_pairs`
  accepts only the `:er` and `:configuration` nulls; Rombach `beta` is the fraction of
  nodes in the periphery (paper semantics); Lip `quality` is the discrete Pearson
  correlation of the returned partition rather than the Lip loss; LowRank-Core `beta`
  is the Find-Cut lower bound rather than a core fraction; the Surprise optimiser and
  the multi-pair pair-count selection are package variants of the published
  procedures.
- Replaced the tutorial's Rombach grid search by raw quality with a sensitivity sweep,
  fixed the crashing Jaccard example (now a doctest), and added "Which null model
  answers which question" guidance to the tutorial and algorithms pages.
- Added a shared "Choosing an algorithm" table to the README, full references with
  DOIs, a `CITATION.cff`, and a documented statement that the pre-refactor benchmark
  snapshot is historical and not reproducible from this repository's history.

### Tooling

- CI tests Julia 1.10, 1.12, and 1.13 with one and two threads, plus an
  allowed-to-fail nightly leg; the documentation workflow builds on Julia 1.13.
- Added CompatHelper and TagBot workflows, Aqua.jl quality checks in the test target,
  and a split of `src/CorePeriphery.jl` into one file per algorithm family
  (`errors.jl`, `results.jl`, `validation.jl`, `quality.jl`, `borgatti_everett.jl`,
  `lip.jl`, `rombach.jl`, `rossa.jl`, `surprise.jl`, `label_switching.jl`,
  `spectral.jl`, `directed.jl`, `multipair.jl`, `significance.jl`).
- Test hygiene: no global `Random.seed!`; explicit `rng=` on stochastic calls; the
  one-time disconnected-input warning is captured where planted fixtures are
  disconnected.

### Breaking and compatibility changes

- Raised the Julia compatibility floor from 1.6 to 1.10.
- Algorithm-specific quality values and several optimizer results differ from 0.1
  because published statistics and objectives replaced earlier approximations.
- Continuous detectors split at the median with a strict `>` (ties are periphery);
  `rombach_continuous` defaults to `partition=:beta`, its `coreness` sums to one, and
  its `quality` is the full double sum of eq. (4.10).
- `multiple_cp_pairs` defaults to `pair_selection=:penalized`; its configuration-null
  objective includes the eq. (13) `i = j` term.
- `cp_significance` defaults to `n_samples=999` and `pass_rng=true`; binary nulls
  reject weighted input.
- `adjacency_to_matrix` raises on conflicting duplicate dyads (use `combine=`) instead
  of letting the last weight win.
- Isolates are never core and an edgeless network has an empty core.
- Configuration-null significance now raises `ConvergenceError` when requested
  switching is incomplete unless the caller explicitly selects `swap_shortfall=:warn`
  or `:accept`.
- `CPMultiResult` and `CPSignificanceResult` contain additional diagnostics and
  effect-size fields; `ideal_cp_matrix(...; cp_block=:missing)` returns a
  `Matrix{Union{Missing,Float64}}`.
- Rossa `n_walks` and `walk_length` keywords and Lip's `max_iter` remain accepted but
  are deprecated, ignored, and emit `Base.depwarn`; profiling is analytic rather than
  simulated. `rossa_profile_ensemble` rejects forwarded `tie_break`/`rng`.

<!-- Replace `main` with `v0.2.0` once the release is tagged. -->
[0.2.0]: https://github.com/simoneSantoni/CorePeriphery.jl/compare/d5f2c33...main
