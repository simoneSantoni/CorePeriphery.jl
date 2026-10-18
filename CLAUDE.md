# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

CorePeriphery.jl is a Julia package (v0.2.0, Julia ≥ 1.10) implementing ten published
core-periphery detection algorithms plus Monte Carlo significance testing. The package's
defining concern is **scientific fidelity to the source papers**: each algorithm has a
documented contract (input requirements, estimand, quality definition) and an oracle test
suite. The papers are cited by DOI in `docs/src/algorithms.md`; they are not stored in
the repository.

Julia 1.10 is the compatibility floor; Julia 1.13 is the primary development,
regression-gate, and documentation runtime; CI also tests 1.12 and an informational
nightly leg. The versions are installed via juliaup (`julia +1.10`, `julia +1.12`,
`julia +1.13.2`).

## Commands

The project has no committed Manifest (gitignored). On a fresh checkout run
`julia --project=. -e 'using Pkg; Pkg.instantiate()'` first, otherwise loading fails on
the missing Arpack dependency.

```sh
# Full test suite (~546,000 checks, 541,852 of them from the brute-force Borgatti-Everett
# optimum over every graph with n ≤ 6; also pulls in Graphs for the extension tests,
# Aqua for package-quality checks, and Logging). Run with two threads so the threaded
# significance schedules are exercised.
julia +1.13.2 --threads=2 --project=. -e 'using Pkg; Pkg.test()'
# Repeat on +1.12 (threads=2) and +1.10 before a release.

# Single test file: every file under test/ carries its own `using` lines and can be
# run directly (keep it that way when adding tests; runtests.jl just includes them)
julia +1.13.2 --project=. test/significance.jl
# test/graphs_extension.jl additionally needs Graphs.jl, test/aqua.jl needs Aqua.jl

# Strict docs build (Documenter 1: checkdocs=:exports, doctests on, warnings fatal).
# Needs Julia >= 1.11: docs/Project.toml uses a [sources] entry that 1.10's Pkg
# cannot resolve. Instantiate the docs env once per checkout.
julia +1.13.2 --project=docs -e 'using Pkg; Pkg.instantiate()'
julia +1.13.2 --project=docs docs/make.jl
# .github/workflows/Documentation.yml runs this on Julia 1.13. Pull requests only
# build; pushes to main and tags deploy to the gh-pages branch (dev/ and stable/).
# .github/workflows/CI.yml runs Pkg.test() on Julia 1.10, 1.12, and 1.13 with 1 and 2
# threads (nightly allowed to fail), uploads coverage from the 1.13/2-thread leg, and
# runs the quick scaling/allocation gate on all three versions.

# Dependency-free scaling + allocation regression gate (thresholds in
# benchmark/regression_targets.jl; statistic is the minimum over 5 quick / 7 full
# repetitions; TSV columns are minimum/median/maximum)
julia +1.13.2 --threads=2 --project=. benchmark/scaling.jl --quick --check

# Julia 1.10/1.12/1.13.2 versioned performance snapshot (writes TSV + Markdown; pass
# scratch paths unless you intend to replace the committed benchmark/ snapshot)
benchmark/run_versioned_performance.sh [out.tsv out.md] [--quick]

# Example script
julia +1.13.2 --project=. examples/basic_usage.jl
```

Pre-release checks recorded in `docs/src/development.md` also include `git diff --check`
and a YAML parse of `.github/workflows/*.yml`.

The cross-package comparison against Python `cpnet` (`benchmark/fit_cpnet.py`,
`analyze.py`, `diagnose_divergence.py`) needs a manually assembled Python env pinned to
NumPy 1.x and a specific cpnet commit; see `benchmark/README.md`. It is not part of
`Pkg.test()`.

## Architecture

### Module layout

`src/CorePeriphery.jl` holds only the module docstring, `using` lines, the export list,
and the `include` order. Everything else lives in one file per concern:

| File | Contents | Public entry points |
|---|---|---|
| `src/errors.jl` | `ConvergenceError` | `ConvergenceError` |
| `src/results.jl` | All result structs, `NullDiagnostics`, `==`/`isequal`/`hash`/`show` for each, `coreness_scores` | `CPResult`, `CPDirectedResult`, `CPMultiResult`, `CPSignificanceResult`, `StationaryDistributionResult`, `RossaEnsembleResult`, `coreness_scores` |
| `src/validation.jl` | `_unwrap_adjacency` (sparse wrappers → `SparseMatrixCSC`), `_issymmetric` with scaled tolerance, `_validate_adjacency`, `_normalize_scores`, component helpers `_component_labels`/`_isolates`/`_warn_disconnected`/`_enforce_isolates!` | none |
| `src/quality.jl` | `adjacency_to_matrix` (indexed and labelled), `block_densities`, `ideal_cp_matrix`, `core_quality`, the Pearson sufficient-statistic kernels, the discrete `_DiscreteState` flip machinery shared by BE-discrete and label switching, and the `_partition` helper | `adjacency_to_matrix`, `block_densities`, `ideal_cp_matrix`, `core_quality` |
| `src/borgatti_everett.jl` | Continuous coordinate ascent; discrete greedy flips, restarts, Gray-code `_exhaustive_discrete` | `borgatti_everett_continuous`, `borgatti_everett_discrete` |
| `src/lip.jl` | Exact degree-prefix minimiser | `lip_discrete` |
| `src/rombach.jl` | Eq. (4.9) template (`_rombach_template`, `⌊βN⌋` boundary), pair-swap optimiser, eq. (4.11) aggregate | `rombach_continuous`, `rombach_core_score` |
| `src/rossa.jl` | Persistence profile, stationary solvers (`:lazy`, `:linear`, `:auto`), ensemble | `random_walker_profiling`, `rossa_stationary_distribution`, `rossa_profile_ensemble` |
| `src/surprise.jl` | Log-factorial table, mode-outward hypergeometric tail (`_log_bimodular_tail`), greedy optimiser | `surprise_cp` |
| `src/label_switching.jl` | Random-order greedy flips on the discrete BE objective | `label_switching_cp` |
| `src/spectral.jl` | LowRank-Core (`_lowrank_core`), ARPACK settings and `_ARPACK_LOCK`, pruned thresholded scores, Find-Cut | `spectral_method` |
| `src/directed.jl` | Symmetric `ww'` and directed `uv'` MINRES | `minres_symmetric`, `minres_svd_directed`, `minres_svd` (compatibility projection) |
| `src/multipair.jl` | Kojaku–Masuda joint pair/role label switching, ER and configuration null contexts (`_joint_diagonal_term` is the eq. 13 `i = j` term) | `multiple_cp_pairs`, `multiple_cp_pairs_config` |
| `src/significance.jl` | Null samplers, `_NullWorkspace` (dense or sparse), `_NullSpec`, `_run_null_samples!` scheduler, `cp_significance` | `cp_significance` |

`ext/CorePeripheryGraphsExt.jl` is a weak-dependency package extension: it defines
`adjacency_to_matrix(::AbstractGraph)` (sparse Float64) and metaprograms one method per
detector that converts a Graphs.jl graph and forwards. Graphs.jl is never a hard
dependency; the detector list in the extension must be kept in sync when adding a
public detector.

### Significance workspace and threading model

`cp_significance` draws `n_samples + 1` seeds from `rng` up front (one per sample plus
the observed fit), so every schedule consumes the same per-sample `MersenneTwister`.
`_NullSpec` holds the read-only null parameters; each worker owns a `_NullWorkspace`
whose `network` is a dense `Matrix{Float64}` or a `SparseMatrixCSC{Float64,Int}`
chosen from the input storage (`_assemble_csc!` builds valid CSC in place with a
counting sort). `_draw_null!` has one method per storage type and both consume the RNG
identically, which is what the dense/sparse `isequal` tests assert. `_run_null_samples!`
implements the schedules: `:serial`, `:greedy` (a `Channel` work queue with one
`Threads.@spawn` task per thread, each creating its own workspace inside the task), and
`:static` (round-robin `@threads :static`). `thread_schedule=:auto` is `:greedy` on
every Julia version. Names assigned inside the task closures must be declared `local`:
a bare assignment to a name that also exists in the enclosing function becomes a shared
boxed variable, which was the root cause of the 0.2.0 greedy race.

Arpack.jl keeps workspace in global Fortran state and is not reentrant, so sparse
LowRank-Core eigensolves go through `_ARPACK_LOCK` (`src/spectral.jl`); any new
ARPACK call must take the same lock.

### Cross-cutting conventions

- **Dense/sparse dual dispatch.** Many internal kernels have a generic
  `AbstractMatrix{<:Real}` method and a specialised `SparseMatrixCSC` method
  (`_validate_adjacency`, `_issymmetric`, `_adjacency_pattern_product`,
  `_node_block_links`, `_block_sum`, `_neighbor_lists`, `_positive_reachable`,
  `_edge_list`, `_joint_build_graph`, `_lowrank_eigenpairs`, `_draw_null!`, ...).
  Every public entry point starts with `A = _unwrap_adjacency(A)` so `Symmetric`,
  `Hermitian`, `Transpose`, and `Adjoint` wrappers around sparse matrices hit the sparse
  methods. Tests assert dense and sparse paths give identical results; preserve that.
- **Centralised input validation.** `_validate_adjacency` enforces square, finite,
  nonnegative, loop-free, and optionally symmetric/binary input; algorithms pass the
  `symmetric` / `binary` flags they need (Lip, LowRank-Core, Surprise, and the binary
  significance nulls require simple binary graphs). Symmetry uses the absolute
  tolerance `1e-12 · max(1, max|A|)`. Connectivity is not checked there:
  `random_walker_profiling` and `rossa_stationary_distribution` require connected /
  strongly connected input via `_positive_reachable`; every other detector except
  `surprise_cp` and `spectral_method` calls `_warn_disconnected(A, :detector)`.
- **Null samples skip the component check.** `cp_significance` evaluates each null
  sample inside `_with_component_checks_suppressed` (a task-local flag in
  `src/validation.jl`), so random null graphs neither warn nor pay the `O(n^2)`
  component search. Isolate handling still runs, keeping null and observed
  statistics on the same estimator.
- **`_warn_disconnected` fires at most once per detector per process.** It uses
  `@warn ... maxlog=1 _id=Symbol(:disconnected_input_, detector)`, so the global
  logger shows one warning per detector per Julia session (confirmed by a
  `TestLogger` test in `test/components_isolates.jl`). A fresh `TestLogger` (as
  installed by `@test_logs`) counts separately, which is how the per-detector warning
  tests work. Tests whose planted fixtures are disconnected capture the warning with
  `@test_logs (:warn, r"connected components") match_mode=:any` or run under
  `Logging.with_logger(Logging.NullLogger())`; never weaken the warning itself. Null
  samples drawn by `cp_significance` are random graphs and can trigger it too.
- **Isolate guarantee.** Nodes with zero in- and out-strength get coreness `0` and are
  never in `core_nodes` (`_enforce_isolates!`; pair label `0` in `multiple_cp_pairs`);
  an edgeless network returns an empty core. Rombach and the multi-pair engine drop
  isolates before fitting and fit the non-isolated subnetwork.
- **Partition rule.** Continuous scores are normalised to [0,1] by `_normalize_scores`
  (Rombach's instead sum to one) and split by `_partition` in `src/quality.jl`:
  `:median` is **strict** (`>` the median, ties are periphery), plus `:none`,
  `:threshold`, `:core_size`; Rombach adds `:beta` (its default: the `N − ⌊βN⌋` core
  ranks). LowRank-Core uses Find-Cut and the discrete detectors return their binary
  state. `core_nodes` of a continuous detector is never an estimate of the model.
- **Quality is algorithm-specific**, evaluated on the returned state, and documented per
  algorithm in `docs/src/contracts.md`. Do not assume it is a Pearson correlation. That
  page is also the source of truth for which detectors accept directed, weighted, or
  sparse input ("Directionality guide") and for every result container's field semantics.
- **Results compare with `isequal`.** All result types define `==`, `isequal`, and
  `hash`; `==` is `NaN`-sensitive (a degenerate null gives `zscore = NaN`, Rossa on two
  nodes gives `quality = NaN`), so tests and user code comparing such results must use
  `isequal`.
- **Errors.** Invalid input raises `ArgumentError`/`DimensionMismatch`; failure to
  converge (ARPACK, Rossa stationary solvers, incomplete configuration switching under
  `swap_shortfall=:error`) raises `ConvergenceError`. Deprecated-but-accepted keywords
  (`n_walks`, `walk_length`, Lip `max_iter`) go through `Base.depwarn`.
- **Reproducibility.** Stochastic algorithms take an explicit `rng` keyword; tests pass
  `rng=MersenneTwister(seed)` on every stochastic call and never call `Random.seed!`.
  Rossa's `tie_break=:deterministic` (lower node index) is the default; `tie_break=:paper`
  reproduces the paper's stochastic rule and needs an explicit RNG. `cp_significance`
  forwards a child RNG to detectors with an `rng` keyword (`pass_rng=true`).
- **Keyword types.** Numeric keywords are `::Real` (converted to `Float64` after
  validation); `ideal_cp_matrix` and `core_quality` accept any real vector.
- **Legacy API.** `minres_svd` exists only for 0.1 compatibility.

### Tests

`test/runtests.jl` holds the basic per-algorithm testset and then includes the
specialised files, each standalone:

- Scientific oracles: `scientific_validation.jl` (exhaustive Lip, directed Rossa),
  `scientific_contracts.jl`, `spectral_lowrank.jl`, `multipair_joint.jl` (literal KM eq.
  3/13 transcriptions), `directed_minres.jl`, `borgatti_everett.jl` (brute-force optimum
  for every graph with n ≤ 6 under both `cp_block` patterns, directed and exhaustive
  variants), `rombach.jl` (eq. 4.9 template oracle, `beta` semantics, core score),
  `surprise_kernel.jl` (`BigInt` tail oracle), `network_algorithms.jl`,
  `julia_algorithms.jl`.
- Behavioural contracts: `significance.jl` (defaults, effect sizes, sparse-vs-dense
  `isequal`, schedule stress test), `hpc_algorithms.jl`, `performance_kernels.jl`
  (sparse/dense equivalence), `result_types.jl` (`==`/`hash`/`show`),
  `generic_inputs.jl` (wrappers, `Real` keywords, tolerance), `partitions.jl`,
  `components_isolates.jl` (component labels, isolate guarantee, warning and its
  `maxlog`), `input_path.jl` (labelled edge lists, `combine`, `block_densities`),
  `graphs_extension.jl`, `aqua.jl`.

A change to an algorithm's objective or move rule must be checked against the paper and
the corresponding oracle test, and the resulting contract change recorded in
`docs/src/contracts.md`, `CHANGELOG.md`, and (if user-visible) `docs/src/migration.md`.
New exports must be added to `docs/src/api.md` (`checkdocs=:exports`) and, if they are
detectors, to the Graphs extension.

### Documentation and benchmark records

`docs/src/` is a Documenter site (index, tutorial, algorithms, contracts, performance,
development, migration, api). The Quick Start in `index.md` and the labelled-input and
block-density examples in `tutorial.md` are doctests; `README.md` mirrors the index
page's choosing table and Quick Start by hand. `docs/src/development.md` is the
audit/refactor record, including the October 2026 panel review. `benchmark/*.md`
(COMPARISON, DIVERGENCE, ABLATION, OPTIMIZATION, VALIDATION, VERSIONED_PERFORMANCE)
are generated or hand-maintained evidence reports; the TSV files beside them are their
machine-readable sources and are regenerated by the scripts named in
`benchmark/README.md`. The committed versioned snapshot predates the Surprise and
sparse-significance rewrites and is historical.
