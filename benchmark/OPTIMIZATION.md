# Optimization results

These development measurements were collected on the same workstation with
warmed Julia specializations and one BLAS thread. They are regression evidence,
not portable performance guarantees; network size, density, convergence, CPU,
and Julia version all affect absolute times.

## Before/after hotspot measurements

| Workload | Before | After | Change |
|---|---:|---:|---:|
| BE discrete, dense n=256, one sweep | 57.0 ms | 0.168 ms | about 339× faster |
| BE continuous, dense n=96, one sweep | 18.8 ms | 0.056 ms | about 336× faster |
| Random-walker profile, dense n=256 | 8.79 ms | 0.273 ms | about 32× faster |
| Exact Surprise, dense n=40 | 469 ms | 0.354 ms | over 1,000× faster |
| Multi-pair configuration, dense n=256 | 26.44 ms | 0.805 ms | about 33× faster |
| Multi-pair configuration, sparse n=256 | 73.95 ms | 1.021 ms | about 72× faster |
| Directed MINRES, sparse n=300, 2 iterations | 4.31 ms | 0.058 ms | about 74× faster |
| LowRank-Core, sparse n=300 | 11.75 ms | 0.709 ms | about 17× faster |

Multi-pair allocation fell about 16×; sparse LowRank-Core allocation fell 92%.
Configuration-null significance sampling fell from roughly 694 KiB to 49 KiB
per sample, about 14× lower.

The large greedy-method improvements come from changing complexity, not loop
annotations: sufficient statistics replace full Pearson recomputation; cached
block counts replace Surprise rescans; and incremental links-to-set replace the
cubic persistence implementation. Sparse spectral and directed improvements come
from partial eigenpairs, streamed reconstruction scores, `mul!`, and analytic
residuals.

## Current regression-gate sample

The two-thread quick gate on Julia 1.12.4 reported zero steady-state allocations
for dense and sparse `core_quality`; at n=256 it measured 0.045 ms dense and
0.010 ms sparse. The exact command was:

```sh
julia --threads=2 --project=. benchmark/scaling.jl --quick --check
```

The gate checks dense/sparse objective equality, deterministic threaded Monte
Carlo output, allocation ceilings, and size-doubling ratios. See
`regression_targets.jl` for its deliberately hardware-independent limits.
