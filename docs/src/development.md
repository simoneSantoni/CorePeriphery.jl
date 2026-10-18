# Development Record and Roadmap

This page records the scientific audit, implementation refactor, performance work,
and cross-package validation carried out in July 2026. It also identifies the work
that should precede the next release. The implementation described here is development
work on top of baseline commit `a2226832508c57efdeab71cde1f055549ca028e7`; it should
not be treated as released until it has been reviewed, committed, and tagged.

## Scope and evidence

The work combined four perspectives: network-analysis definitions, Julia language and
API design, software-engineering correctness, and high-performance computing. Published
algorithm descriptions were checked against implementations, direct mathematical
oracles were added where practical, and CorePeriphery.jl was compared empirically with
Python `cpnet` commit `6aad458a6d434a3617d33e74f7163d514a27fecb`.

The evidence is retained in four reports:

- [`benchmark/COMPARISON.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/COMPARISON.md)
  records identical-data fits, recovery metrics, concordance, and warmed runtimes.
- [`benchmark/ABLATION.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/ABLATION.md)
  separates the pre-refactor implementation, current dense algorithms, and current
  sparse kernels.
- [`benchmark/DIVERGENCE.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/DIVERGENCE.md)
  traces lower-concordance algorithms one mechanism at a time.
- [`benchmark/OPTIMIZATION.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/OPTIMIZATION.md)
  records before/after hotspot measurements and the current regression gate.

These links target `main` and will resolve after the development work is merged.

## Completed scientific and API work

| Area | Completed work | Resulting contract |
|---|---|---|
| Shared input validation | Centralized square, finite, nonnegative, loop-free, symmetry, and binary checks, including sparse matrices | Invalid inputs fail consistently instead of being silently reinterpreted |
| Borgatti--Everett | Reworked objective evaluation around sufficient statistics and protected degenerate inputs | Reported quality is evaluated on the returned state; empty/constant cases are finite |
| Lip | Replaced the previous heuristic with the published degree-prefix minimizer | Simple undirected binary input; deterministic smallest minimizing prefix |
| LowRank-Core | Implemented rank-two reconstruction, thresholded row scores, and Find-Cut | Binary undirected input; sparse matrices use two ARPACK eigenpairs without a dense reconstruction |
| Rombach | Replaced coordinate-grid work with seeded template-permutation restarts and cached swap deltas | Explicit RNG and reproducible restarts; quality is the optimized template objective |
| MINRES | Added distinct symmetric `ww'` and directed `uv'` rank-one estimators | Symmetric comparisons use `minres_symmetric`; `CPDirectedResult` retains sender, receiver, combined scores, residual, and explained-weight quality |
| Della Rossa | Implemented equation-(5) persistence and equation-(8) centralization, then added the published stochastic tie mode and directed stationary-flow formulation | Deterministic compatibility default; `tie_break=:paper` with an RNG; directed input must be strongly connected |
| Surprise | Replaced the approximation with the multivariate-hypergeometric joint-tail statistic in log space | Exact statistic for simple undirected binary graphs |
| Multiple CP pairs | Added joint pair/role label switching plus compatibility and penalized pair-count selection | Seeded restarts, pair-size constraints, and candidate objective/selection diagnostics in `CPMultiResult` |
| Significance | Added ER, degree-switching, directed in/out-degree switching, and fixed-topology weight-permutation nulls | Reproducible serial/threaded sampling, strict switch completion, and per-sample diagnostics through `CPSignificanceResult` |
| Rossa robustness | Added lazy, constrained-linear, and automatic stationary solvers plus paper-mode ensembles | Solver residual/iteration metadata and serial/thread-reproducible uncertainty summaries |
| Ecosystem integration | Added an optional Graphs.jl extension and generic dense/sparse matrix support | Graph objects are converted lazily through the weak dependency; Graphs.jl is not a hard dependency |

Legacy `n_walks` and `walk_length` Rossa keywords remain accepted but are ignored,
because the implemented method computes persistence analytically rather than simulating
walks. Julia 1.10 is the compatibility floor; Julia 1.12 is the primary optimization,
documentation, and regression-gate runtime.

## Rossa paper audit

The Della Rossa audit deserves separate treatment because it changed the interpretation
of the earlier cross-package results.

For an undirected graph, CorePeriphery.jl incrementally evaluates the exact candidate
persistence ratio. For a directed graph it computes the stationary distribution of the
row-normalized transition matrix, then grows the peripheral set using stationary edge
flow. Directed tie strength is total in-plus-out strength, as specified by the paper.
A lazy stationary iteration preserves the target distribution and also converges on
periodic strongly connected graphs.

The API now distinguishes two selection policies:

- `tie_break=:deterministic` selects the lower node index after the paper's persistence
  and minimum-strength criteria. This preserves the package's reproducible default.
- `tie_break=:paper` randomly chooses the initial minimum-strength node and randomly
  resolves the paper's final ties. An explicit RNG makes a run repeatable.

The cpnet divergence is not just random tie variation. On noisy planted graphs, cpnet
minimizes a non-equivalent surrogate and omits the secondary minimum-strength filter.
It also reports `1 - mean(alpha)` rather than equation-(8) centralization and does not
implement the paper's stationary directed equation. The exact trace audit found a mean
of 11.385 and 16.950 paper-rule violations per cpnet run on the noisy single-pair and
two-pair graphs, respectively.

## Performance refactor

The main performance gains came from changing asymptotic work and data movement:

- sufficient statistics replaced repeated full Pearson-correlation scans;
- cached block counts replaced repeated Surprise rescans;
- incremental boundary weights replaced cubic Rossa set reconstruction;
- sparse LowRank-Core uses partial eigenpairs and streams rank-two scores;
- symmetric and directed MINRES reuse matrix-vector workspaces and evaluate residuals analytically;
- multi-pair moves use adjacency lists and cached pair aggregates;
- threaded significance allocates deterministic child RNG streams and reuses per-worker
  or per-task buffers.

Development measurements reported approximately 339x and 336x improvements for one
dense Borgatti--Everett discrete and continuous sweep, 32x for dense Rossa profiling,
more than 1,000x for the exact Surprise workload, 33--72x for multi-pair workloads,
74x for sparse directed MINRES, and 17x for sparse LowRank-Core. These are same-machine,
warmed development measurements, not portable performance guarantees.

The dense/sparse ablation is scientifically important: current dense and sparse paths
produced identical top-k sets, explicit core sets, and pair partitions in all 27 fitted
cases. Therefore storage-specific optimization does not explain the CorePeriphery.jl/
cpnet result divergence on this corpus. Sparse storage is a scalability path, not an
automatic speedup for tiny graphs; it was faster in 11 of the 27 small-network fits.

## Cross-package findings

Across three planted networks and nine comparable algorithm configurations,
CorePeriphery.jl was faster in 25 of 27 warmed rows, with a median cpnet/Julia runtime
ratio of 33.6x in the checked comparison. The comparison reflects each package's
configured estimator and restart policy, so it must not be interpreted as equal
low-level work.

Lower concordance was traced to algorithmic or output-definition differences:

- The headline MINRES row now fits the same symmetric `ww'` estimand in both packages.
  cpnet still selects the largest residual across restarts, uses a scorer inconsistent
  with its gradient, and stops loosely. Correct selection and tighter convergence raise
  mean rank concordance with Julia from 0.320 to 0.918.
- Multi-pair methods enter different local basins and use different state constraints;
  cpnet additionally has inconsistent ER density conventions and a missing node-degree
  factor in the configuration move gain.
- cpnet Rombach exits after its first improving sweep. Completing the swap search
  recovers Julia's ordering on the noisy and two-pair cases.
- LowRank-Core exposes different node quantities: Julia retains reconstruction scores,
  while cpnet returns the binary Find-Cut state as coreness. The latent reconstruction
  rankings agree on all three graphs.
- Rossa differs for the substantive reasons documented above.

Concordance is diagnostic, not a claim that cpnet is ground truth. Several identically
named methods expose different mathematical estimands.

## Verification completed

The 0.2 release candidate passes:

- all 1,525 package checks through `Pkg.test()` on Julia 1.12 (and 1,524 on Julia
  1.10), including 1,056
  exhaustive/direct scientific-validation checks, sparse/dense equivalence, Graphs.jl
  dispatch, directed Rossa, statistical contracts, and threading;
- a direct stationary linear-system oracle for weighted directed Rossa;
- seeded stochastic Rossa reproducibility and periodic directed-cycle tests;
- the strict Documenter build and doctests;
- `git diff --check`;
- the existing cpnet divergence replay after the Rossa extension;
- the dependency-free two-thread scaling and allocation gate on Julia 1.10 and 1.12.

## Remaining limitations and risks

1. The empirical cross-package corpus contains only three small planted networks. The
   independent oracle suite is much broader at small `n`, but larger licensed real,
   temporal, weighted, directed, and adversarial corpora remain future research work.
2. Edge switching has no universal finite mixing-time guarantee. CorePeriphery.jl now
   guarantees that the requested accepted-swap count completes or is explicitly
   reported; that is a completion guarantee, not proof that a chain has mixed.
3. The direct stationary solver intentionally materializes dense matrices and is for
   smaller problems. Sparse `:lazy` remains subject to slow mixing, but failures now
   report method, iterations, tolerance, and residual and users can supply a validated
   stationary vector.
4. Pair-count penalization controls the tested planted single-pair cases but is a
   model-selection heuristic, not a universal significance guarantee. Candidate scores
   are exposed so sensitivity to `pair_penalty` can be audited.
5. The weight-permutation null preserves topology and the global weight multiset, not
   node strengths. A strength-preserving weighted configuration null would answer a
   different question and remains a possible future feature.
6. The cross-package Python environment is manually assembled around the 2022 cpnet
   revision and NumPy 1.x; exact revisions are recorded, but this leg is not in CI.

## Release-roadmap resolution

The audit roadmap is implemented in 0.2.0. Each issue is closed only after the clean
supported-version CI and release checks complete.

| Priority | Delivered change | Tracking issue | Release status |
|---:|---|---|---|
| 1 | Explicit configuration-switch completion and diagnostics | [#3](https://github.com/simoneSantoni/CorePeriphery.jl/issues/3) | Implemented and tested |
| 2 | Alternative Rossa stationary solvers and convergence metadata | [#4](https://github.com/simoneSantoni/CorePeriphery.jl/issues/4) | Implemented and tested |
| 3 | Exhaustive/direct oracle and provenance validation matrix | [#5](https://github.com/simoneSantoni/CorePeriphery.jl/issues/5) | Implemented and tested |
| 4 | Penalized pair-count selection with candidate diagnostics | [#6](https://github.com/simoneSantoni/CorePeriphery.jl/issues/6) | Implemented and tested |
| 5 | Symmetric rank-one MINRES alongside directed MINRES | [#7](https://github.com/simoneSantoni/CorePeriphery.jl/issues/7) | Implemented and tested |
| 6 | Reproducible Rossa ensemble uncertainty summaries | [#8](https://github.com/simoneSantoni/CorePeriphery.jl/issues/8) | Implemented and tested |
| 7 | Directed configuration and weighted permutation significance | [#9](https://github.com/simoneSantoni/CorePeriphery.jl/issues/9) | Implemented and tested |
| 8 | Julia 1.10/1.12 versioned performance and allocation tracking | [#10](https://github.com/simoneSantoni/CorePeriphery.jl/issues/10) | Implemented and tested |
| 9 | Algorithm/result/directionality contract documentation | [#11](https://github.com/simoneSantoni/CorePeriphery.jl/issues/11) | Implemented; strict docs pass |
| 10 | Changelog, migration guide, version 0.2.0, and release review | [#12](https://github.com/simoneSantoni/CorePeriphery.jl/issues/12) | Awaiting clean CI and tag |

## Reproduction commands

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=docs docs/make.jl
julia --threads=2 --project=. benchmark/scaling.jl --quick --check
```

The full cross-package and ablation workflows, including pinned environment details,
are documented in `benchmark/README.md`.
