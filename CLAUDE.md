# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

CorePeriphery.jl is a Julia package (v0.2.0, Julia ≥ 1.10) implementing ten published
core-periphery detection algorithms plus Monte Carlo significance testing. The package's
defining concern is **scientific fidelity to the source papers**: each algorithm has a
documented contract (input requirements, estimand, quality definition) and an oracle test
suite. The reference PDFs/HTML live in `algo/`.

Julia 1.10 is the compatibility floor; Julia 1.12 is the primary optimization and
regression-gate runtime. Both are installed via juliaup (`julia +1.10`, `julia +1.12`).

## Commands

The project has no committed Manifest (gitignored). On a fresh checkout run
`julia --project=. -e 'using Pkg; Pkg.instantiate()'` first, otherwise loading fails on
the missing Arpack dependency.

```sh
# Full test suite (~1,500 checks; also pulls in Graphs for the extension tests)
julia --project=. -e 'using Pkg; Pkg.test()'

# Single test file: every file under test/ except performance_kernels.jl is
# self-contained (own `using` lines) and can be run directly
julia --project=. test/significance.jl
# test/performance_kernels.jl relies on names loaded by runtests.jl (Statistics, etc.)
# test/graphs_extension.jl needs Graphs.jl available in the environment

# Strict docs build (Documenter 1: checkdocs=:exports, doctests on, warnings fatal).
# .github/workflows/Documentation.yml runs this on Julia 1.12 and deploys to the
# gh-pages branch on pushes to main and on tags (dev/ and stable/ URLs).
julia --project=docs docs/make.jl

# Dependency-free scaling + allocation regression gate (thresholds in benchmark/regression_targets.jl)
julia --threads=2 --project=. benchmark/scaling.jl --quick --check

# Julia 1.10/1.12 versioned performance snapshot (writes TSV + Markdown in benchmark/)
benchmark/run_versioned_performance.sh [out.tsv out.md] [--quick]

# Example script
julia --project=. examples/basic_usage.jl
```

Pre-release checks recorded in `docs/src/development.md` also include `git diff --check`.

The cross-package comparison against Python `cpnet` (`benchmark/fit_cpnet.py`,
`analyze.py`, `diagnose_divergence.py`) needs a manually assembled Python env pinned to
NumPy 1.x and a specific cpnet commit; see `benchmark/README.md`. It is not part of
`Pkg.test()`.

## Architecture

### Module layout

`src/CorePeriphery.jl` is the main module (~1,900 lines). It declares all exports, the
`CPResult` / `CPMultiResult` / `StationaryDistributionResult` / `RossaEnsembleResult`
structs, shared validation and kernels, and the Borgatti–Everett, Lip, Rombach, Della
Rossa, Surprise, and label-switching algorithms. Near the end it `include`s four files
and then defines thin public wrappers around them:

| File | Contents | Public entry points |
|---|---|---|
| `src/spectral.jl` | Cucuringu LowRank-Core rank-two reconstruction + Find-Cut (`_lowrank_core`); sparse path uses two ARPACK eigenpairs | `spectral_method` |
| `src/multipair.jl` | Kojaku–Masuda joint pair/role label switching with ER and configuration null contexts (`_multiple_cp_pairs_joint`) | `multiple_cp_pairs`, `multiple_cp_pairs_config` |
| `src/directed.jl` | `CPDirectedResult`; symmetric `ww'` and directed `uv'` MINRES rank-one estimators | `minres_symmetric`, `minres_svd_directed` (`minres_svd` is a `CPResult` compatibility projection) |
| `src/significance.jl` | `CPSignificanceResult`; ER, degree-configuration, directed in/out configuration, and weight-permutation null samplers; deterministic threaded sampling | `cp_significance` |

`ext/CorePeripheryGraphsExt.jl` is a weak-dependency package extension: it defines
`adjacency_to_matrix(::AbstractGraph)` (sparse Float64) and metaprograms one method per
detector that converts a Graphs.jl graph and forwards. Graphs.jl is never a hard
dependency; the detector list in the extension must be kept in sync when adding a
public detector.

### Cross-cutting conventions

- **Dense/sparse dual dispatch.** Many internal kernels have a generic
  `AbstractMatrix{<:Real}` method and a specialised `SparseMatrixCSC` method
  (`_validate_adjacency`, `_adjacency_pattern_product`, `_links_to_periphery`,
  `_positive_reachable`, `_edge_list`, `_joint_build_graph`, `_lowrank_eigenpairs`, ...).
  Tests assert dense and sparse paths give identical results; preserve that when touching
  either.
- **Centralised input validation.** `_validate_adjacency` enforces square, finite,
  nonnegative, loop-free, and optionally symmetric/binary input. Method-specific
  requirements (Lip, LowRank-Core, Surprise, and configuration-null significance require
  simple binary graphs; Rossa requires connected / strongly connected) are enforced there
  rather than silently reinterpreted.
- **Quality is algorithm-specific**, evaluated on the returned state, and documented per
  algorithm in `docs/src/contracts.md`. Do not assume it is a Pearson correlation. That
  page is also the source of truth for which detectors accept directed, weighted, or
  sparse input ("Directionality guide") and for every result container's field semantics.
- **Coreness normalisation and thresholds.** Continuous scores are passed through
  `_normalize_scores` to [0,1]; core/periphery split defaults to the median unless the
  algorithm defines its own partition.
- **Reproducibility.** Stochastic algorithms take an explicit `rng` keyword. Rossa's
  `tie_break=:deterministic` (lower node index) is the default; `tie_break=:paper`
  reproduces the paper's stochastic rule and needs an explicit RNG. Threaded significance
  uses deterministic child RNG streams so serial and threaded runs agree.
- **Legacy keywords.** Rossa `n_walks` / `walk_length` are accepted but ignored
  (profiling is analytic). `minres_svd` exists only for 0.1 compatibility.
- **Configuration-null shortfall policy.** `cp_significance` errors when requested
  edge switching does not complete unless `swap_shortfall=:warn` or `:accept` is passed.

### Tests

`test/runtests.jl` holds the basic per-algorithm testset and then includes the
specialised files. The scientifically important ones: `scientific_validation.jl`
(exhaustive/direct oracles, ~1,000 checks), `scientific_contracts.jl`,
`spectral_lowrank.jl`, `multipair_joint.jl`, `directed_minres.jl`, `significance.jl`,
`performance_kernels.jl` (sparse/dense equivalence), `graphs_extension.jl`.

A change to an algorithm's objective or move rule must be checked against the paper in
`algo/` and the corresponding oracle test, and the resulting contract change recorded in
`docs/src/contracts.md`, `CHANGELOG.md`, and (if user-visible) `docs/src/migration.md`.

### Documentation and benchmark records

`docs/src/` is a Documenter site (index, tutorial, algorithms, contracts, performance,
development, migration, api). `docs/src/development.md` is the audit/refactor record and
lists the release-roadmap issues. `benchmark/*.md` (COMPARISON, DIVERGENCE, ABLATION,
OPTIMIZATION, VALIDATION, VERSIONED_PERFORMANCE) are generated or hand-maintained
evidence reports; the TSV files beside them are their machine-readable sources and are
regenerated by the scripts named in `benchmark/README.md`.
