# Performance-refactor ablation

This experiment separates the pre-refactor implementation at commit `a2226832508c57efdeab71cde1f055549ca028e7`, the current algorithms on dense matrices, and the current optimized sparse paths. Each runtime is the median of three warmed fits with one Julia and one BLAS thread. The same seed is restored before every fit. Rombach is limited to ten iterations in all three legs because the old coordinate-grid implementation otherwise requires billions of allocations even on these small graphs.

## Causal conclusion

The current dense and sparse execution paths select the same cores and pair partitions in every fitted case. Storage-specific performance work therefore does **not** explain the observed CorePeriphery.jl/cpnet result divergence on this corpus. The broader refactor can explain part of it only where it deliberately changed the mathematical statistic, optimizer, tie convention, or returned score.

## Old, current, and cpnet concordance

Values are means over the three networks. Rank columns are Spearman correlations; top-k columns compare the highest-ranked nodes using the planted core size. AUC is measured independently against planted truth. The old package has no KM-configuration implementation.

| Algorithm | Rank O↔N | Rank O↔P | Rank N↔P | Top-k O↔N | Top-k O↔P | Top-k N↔P | AUC old | AUC new | AUC cpnet |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BE | 0.667 | 0.621 | 0.915 | 1.000 | 1.000 | 1.000 | 0.772 | 0.939 | 0.931 |
| KM_ER | 0.000 | 0.000 | 0.658 | 0.773 | 0.773 | 0.717 | 0.500 | 0.894 | 0.873 |
| KM_config | NA | NA | 0.652 | NA | NA | 0.794 | NA | 0.780 | 0.819 |
| Lip | 0.491 | 0.444 | 0.961 | 0.804 | 0.804 | 1.000 | 0.797 | 0.958 | 0.947 |
| LowRankCore | 0.808 | 0.634 | 0.743 | 0.833 | 0.717 | 0.806 | 0.972 | 0.977 | 0.817 |
| MINRES | 0.992 | 0.318 | 0.320 | 1.000 | 0.333 | 0.333 | 0.903 | 0.915 | 0.724 |
| Rombach | 0.668 | 0.727 | 0.729 | 0.667 | 0.717 | 0.751 | 0.989 | 0.885 | 0.914 |
| Rossa | 0.811 | 0.662 | 0.778 | 0.905 | 0.867 | 0.867 | 0.996 | 0.998 | 0.957 |
| Surprise | 0.925 | 0.978 | 0.947 | 0.949 | 1.000 | 0.949 | 0.956 | 0.995 | 0.972 |

Multi-pair partition concordance on the planted two-pair graph:

| Algorithm | Pair ARI O↔N | Pair ARI O↔P | Pair ARI N↔P |
|---|---:|---:|---:|
| KM_ER | 0.634 | 0.753 | 0.773 |
| KM_config | NA | NA | 1.000 |

O = pre-refactor CorePeriphery.jl, N = current CorePeriphery.jl on dense matrices, and P = cpnet. Concordance is diagnostic only: the old spectral implementation, current LowRank-Core, and several cpnet methods do not expose identical mathematical scores.

## Pre-refactor to current-dense

This comparison includes both scientific changes and optimization. It must not be interpreted as the effect of sparse kernels alone. Values are means over the three planted networks; speedup is the median old/current runtime ratio.

| Algorithm | Rank concordance | Top-k concordance | Old AUC | Current AUC | Pair ARI | Median speedup |
|---|---:|---:|---:|---:|---:|---:|
| BE | 0.667 | 1.000 | 0.772 | 0.939 | NA | 115.8× |
| KM_ER | 0.000 | 0.773 | 0.500 | 0.894 | 0.262 | 112.0× |
| Lip | 0.491 | 0.804 | 0.797 | 0.958 | NA | 4.9× |
| LowRankCore | 0.808 | 0.833 | 0.972 | 0.977 | NA | 1.1× |
| MINRES | 0.992 | 1.000 | 0.903 | 0.915 | NA | 260.4× |
| Rombach | 0.668 | 0.667 | 0.989 | 0.885 | NA | 124614.7× |
| Rossa | 0.811 | 0.905 | 0.996 | 0.998 | NA | 113.2× |
| Surprise | 0.925 | 0.949 | 0.956 | 0.995 | NA | 0.5× |

## Current dense to current sparse

This is the kernel/storage ablation: algorithm definitions, seeds, and parameters are held fixed. Speedup is dense/sparse, so values above one favor sparse storage. Maximum score and quality differences are maxima over all three networks.

| Algorithm | Rank concordance | Top-k concordance | Core-set Jaccard | Pair ARI | Max score Δ | Max quality Δ | Median sparse speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| BE | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 1.53× |
| KM_ER | 1.000 | 1.000 | 1.000 | 1.000 | 0.000e+00 | 0.000e+00 | 1.02× |
| KM_config | 1.000 | 1.000 | 1.000 | 1.000 | 0.000e+00 | 0.000e+00 | 0.78× |
| Lip | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 1.25× |
| LowRankCore | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 0.45× |
| MINRES | 0.997 | 1.000 | 1.000 | NA | 3.331e-16 | 4.441e-16 | 0.87× |
| Rombach | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 0.29× |
| Rossa | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 1.12× |
| Surprise | 1.000 | 1.000 | 1.000 | NA | 0.000e+00 | 0.000e+00 | 0.81× |

## Interpretation

- Current dense and sparse paths produced identical top-k sets, explicit core sets, and pair partitions in 27 of 27 fits.
- Sparse storage was faster in 11 of 27 tiny-network fits, with a median dense/sparse runtime ratio of 0.87×. Sparse storage is therefore an output-preserving scalability path, not a universal latency win at 20–48 vertices.
- The pre-refactor spectral method was leading-eigenvector centrality, not LowRank-Core, and is included only to measure public-API drift.
- The old package had only an ER-style multi-pair detector, so there is no pre-refactor KM-configuration row.
- Rossa, Surprise, Rombach, Lip, LowRank-Core, and MINRES received scientific or optimizer changes in addition to faster kernels. Their old/current divergence cannot be attributed to performance mechanics.
- MINRES is especially diagnostic: old/current rank concordance is about 0.99 while current Julia/cpnet concordance is much lower. Its cross-package divergence is therefore principally an implementation/model difference, not a consequence of the Julia performance refactor.
- The pre-refactor BE routine returns a degenerate all-core result on the ideal graph; the current objective correction, rather than sparse storage, accounts for that particular change.
- Exact Surprise is slower than the old approximation on two of the three tiny graphs. That cost buys a different, published statistic and should not be presented as a pure optimization speedup.
- Rombach's very large old/current runtime ratio is descriptive only: the objective and search were changed, and both legs use the stated ten-iteration budget.
- Dense/sparse discrepancies, if any, isolate numerical or tie sensitivity in the current implementation and are listed in `ablation_kernels.tsv`.
- This dense/sparse leg isolates storage-specific execution paths. Algebraic sufficient-statistic, LowRank reconstruction, and MINRES residual kernels are additionally checked against direct oracles in the Julia test suite; the ablation does not manufacture a second slow implementation of every current optimizer.
