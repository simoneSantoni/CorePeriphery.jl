# CorePeriphery.jl vs cpnet empirical comparison

## Environment

```text
CorePeriphery.jl 0.2.0
Julia 1.12.6
Julia threads 1
BLAS threads 1
cpnet 6aad458
upstream commit 6aad458a6d434a3617d33e74f7163d514a27fecb
Python 3.12.3
NumPy 1.26.4
SciPy 1.14.1
Numba 0.66.0
```

Quality values are reported by each package but are not directly comparable when objective definitions or scaling differ. Runtime is the median of five warmed fits; cpnet's Numba compilation is excluded. NumPy, Python, and Numba RNG state is reset before every cpnet fit. Both packages receive sparse matrices and BLAS is pinned to one thread by the benchmark launcher.

| Dataset | Algorithm | Spearman | Top-k Jaccard | Julia AUC | cpnet AUC | Pair ARI | Julia truth ARI | cpnet truth ARI | Julia ms | cpnet ms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ideal_single | BE | 0.882 | 1.000 | 1.000 | 0.967 | 1.000 | 1.000 | 1.000 | 0.004 | 0.139 |
| ideal_single | KM_ER | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.069 | 1.039 |
| ideal_single | KM_config | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.066 | 1.313 |
| ideal_single | Lip | 0.882 | 1.000 | 1.000 | 0.967 | 1.000 | 1.000 | 1.000 | 0.004 | 0.066 |
| ideal_single | LowRankCore | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.080 | 2.690 |
| ideal_single | MINRES | 0.230 | 0.429 | 1.000 | 0.653 | 1.000 | 1.000 | 1.000 | 0.007 | 4.160 |
| ideal_single | Rombach | 0.591 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.064 | 0.646 |
| ideal_single | Rossa | 0.977 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.004 | 1.600 |
| ideal_single | Surprise | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.007 | 12.910 |
| noisy_single | BE | 0.928 | 1.000 | 0.900 | 0.950 | 1.000 | 1.000 | 1.000 | 0.008 | 1.172 |
| noisy_single | KM_ER | 0.888 | 0.818 | 0.917 | 0.883 | 0.551 | 0.000 | 0.000 | 0.180 | 4.236 |
| noisy_single | KM_config | -0.142 | 0.538 | 0.717 | 0.567 | 0.000 | 0.000 | 1.000 | 0.234 | 5.127 |
| noisy_single | Lip | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.006 | 0.058 |
| noisy_single | LowRankCore | 0.556 | 0.818 | 0.993 | 0.700 | 1.000 | 1.000 | 1.000 | 0.343 | 3.458 |
| noisy_single | MINRES | 0.581 | 0.429 | 0.993 | 0.823 | 1.000 | 1.000 | 1.000 | 0.022 | 4.273 |
| noisy_single | Rombach | 0.935 | 0.818 | 0.970 | 0.987 | 1.000 | 1.000 | 1.000 | 2.378 | 1.341 |
| noisy_single | Rossa | 0.768 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.015 | 2.495 |
| noisy_single | Surprise | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.932 | 48.760 |
| two_pairs | BE | 0.936 | 1.000 | 0.917 | 0.875 | 1.000 | 0.000 | 0.000 | 0.008 | 1.539 |
| two_pairs | KM_ER | 0.510 | 0.333 | 0.931 | 0.736 | 0.774 | 1.000 | 0.774 | 0.178 | 4.850 |
| two_pairs | KM_config | 0.374 | 0.500 | 0.736 | 0.889 | 1.000 | 1.000 | 1.000 | 0.155 | 5.342 |
| two_pairs | Lip | 1.000 | 1.000 | 0.875 | 0.875 | 1.000 | 0.000 | 0.000 | 0.005 | 0.049 |
| two_pairs | LowRankCore | 0.673 | 0.600 | 0.938 | 0.750 | 1.000 | 0.000 | 0.000 | 0.216 | 4.311 |
| two_pairs | MINRES | 0.134 | 0.143 | 0.752 | 0.697 | 1.000 | 0.000 | 0.000 | 0.113 | 5.699 |
| two_pairs | Rombach | 0.987 | 0.846 | 0.748 | 0.755 | 1.000 | 0.000 | 0.000 | 3.295 | 2.247 |
| two_pairs | Rossa | 0.589 | 0.600 | 0.995 | 0.870 | 1.000 | 0.000 | 0.000 | 0.017 | 3.783 |
| two_pairs | Surprise | 0.842 | 0.846 | 0.986 | 0.917 | 1.000 | 0.000 | 0.000 | 0.285 | 43.590 |

## Aggregate description

CorePeriphery.jl was faster in 25 of 27 warmed rows; the median cpnet/Julia runtime ratio was 33.6×. This ratio includes each implementation's configured search procedure: cpnet uses five random starts for BE, MINRES, Surprise, Rombach, and both KM methods, whereas CorePeriphery.jl uses five starts only for Rombach and KM and deterministic searches for BE and Surprise. The runtime ratio therefore describes these fitted estimators, not equal low-level work per restart.

## Interpretation

- Lip agrees exactly on the noisy and two-pair core rankings; its ideal-case difference is a tied boundary convention in cpnet.
- LowRank-Core agrees exactly on the ideal graph. On noisy graphs, CorePeriphery.jl retains rank-two reconstruction scores whereas cpnet exposes the binary Find-Cut result as coreness.
- Rombach, Rossa, and discrete BE show strong agreement on planted single-pair structure. The headline MINRES comparison now fits the same symmetric off-diagonal `w*w'` model in both packages. Its remaining divergence is optimizer behavior: cpnet selects the largest residual across restarts, uses inconsistent update/scoring equations, and stops loosely.
- On the planted two-pair graph, configuration-null pair recovery agrees exactly in this seeded fit; ER-null recovery differs more. Both packages can over-partition single-pair graphs without a significance-filtering stage.
- Surprise agrees on single-pair graphs but uses different statistics and search procedures, so raw quality values—including their signs—must not be compared.
- CorePeriphery.jl is faster in 25 of 27 rows. Only Rombach is runtime-competitive here: cpnet is faster on the noisy and two-pair Rombach fits. The largest ratios occur for deterministic Julia fits compared with cpnet's five-start stochastic fits. These tiny-network timings are descriptive microbenchmarks, not hardware-independent claims.
