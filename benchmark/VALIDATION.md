# Scientific validation matrix

This report separates mathematical correctness checks, literature fixtures, recovery
experiments, and cross-package concordance. Concordance with another implementation is
never treated as a correctness oracle. The current deterministic CI subset is run by
`test/scientific_validation.jl`; broader planted-network results are retained in
`COMPARISON.md`, `ABLATION.md`, and `DIVERGENCE.md`.

## Provenance

| Evidence | Provenance | Version or seed | License/data status |
|---|---|---|---|
| Published formulas | Papers cited in package documentation | Bibliographic references in `docs/src/index.md` | Equations independently implemented in tests |
| Synthetic planted graphs | `benchmark/generate_networks.py` | Recorded generator seed and TSV inputs | Generated in-repository, no external data license |
| cpnet comparison | `skojaku/core-periphery-detection` | Commit `6aad458a6d434a3617d33e74f7163d514a27fecb` | External package used only at benchmark runtime |
| Pre-refactor ablation | CorePeriphery.jl | Commit `a2226832508c57efdeab71cde1f055549ca028e7` | Same repository |

No third-party network dataset is vendored. Canonical star, complete, ideal
core-periphery, directed-cycle, and planted-pair fixtures are generated directly from
their definitions.

## Independent checks by algorithm

| Algorithm | Independent evidence | Location |
|---|---|---|
| Borgatti--Everett continuous/discrete | Direct Pearson correlation from materialized dyad vectors; returned-state objective and degeneracy tests | `test/performance_kernels.jl`, `test/scientific_contracts.jl` |
| Lip | Brute-force loss over every nontrivial partition for all 1,024 labeled five-node simple graphs | `test/scientific_validation.jl` |
| Rombach | Direct template objective and monotone accepted-swap checks | `test/hpc_algorithms.jl` |
| LowRank-Core | Dense eigendecomposition, explicit rank-two reconstruction, thresholding, and exhaustive feasible Find-Cut scores | `test/spectral_lowrank.jl` |
| Della Rossa | Analytic star/complete profiles; direct stationary linear solve; naive equation-(5) set sums on weighted directed graphs | `test/network_algorithms.jl`, `test/scientific_validation.jl` |
| Symmetric MINRES | Explicit off-diagonal residual and permutation/dense/sparse invariance | `test/directed_minres.jl` |
| Directed MINRES | Explicit asymmetric residual, transpose duality, and permutation/dense/sparse invariance | `test/directed_minres.jl` |
| Surprise | Published exact statistic contract, log-space finite result, and planted star fixture | `test/scientific_contracts.jl` |
| Multiple CP pairs | Direct dyad objective, brute move deltas, aggregate rebuilds, and planted one-/two-pair model selection | `test/multipair_joint.jl` |
| Label switching | Direct discrete Pearson quality and seeded monotonic search checks | `test/julia_algorithms.jl` |
| Significance | Exact degree constraints, topology/weight preservation, plus-one p-value bounds, and serial/threaded replay | `test/significance.jl` |

## Recovery and sensitivity coverage

The versioned planted corpus contains an ideal single-pair graph, a noisy single-pair
graph, and a two-pair graph. Reports include Spearman rank correlation, planted-core
AUC, top-k Jaccard, pair adjusted Rand index, runtime, storage ablation, and seeded
mechanism-level diagnostics. Stochastic Rossa additionally reports profile counts,
centralization variation, and rank stability through `rossa_profile_ensemble`.

The corpus remains intentionally small enough for reproducible CI and cross-language
diagnostics. Larger real-network validation is a continuing research activity; any
future dataset must record source, license, checksum, preprocessing, and expected
statistic before inclusion.

## Reproduction

```sh
julia --project=. test/scientific_validation.jl
julia --project=. -e 'using Pkg; Pkg.test()'
python3 benchmark/generate_networks.py
```

See `benchmark/README.md` for the pinned Python environment and full comparison
commands.
