# Performance

CorePeriphery.jl has separate dense and sparse execution paths. Use a
`SparseMatrixCSC` when the network is sparse; validation, Pearson quality,
LowRank-Core, directed MINRES, random-walker profiling, and multi-pair fitting
then operate on stored edges or sparse matrix-vector products instead of scanning
all possible dyads.

For a directed Rossa profile, the stationary distribution is computed by a lazy
matrix-vector iteration that remains sparse and converges even for periodic strongly
connected graphs. `stationary_tol` trades precision for runtime; increase
`stationary_max_iter` only when a slowly mixing network does not reach the requested
tolerance, or use `stationary_method=:linear` for a dense constrained solve.

```julia
using CorePeriphery, SparseArrays

A = sparse(adjacency_to_matrix([(1, 2), (1, 3), (2, 3), (1, 4)], 4))
result = random_walker_profiling(A)
```

## Parallel significance tests

Monte Carlo samples can run concurrently. Start Julia with multiple threads and
set `threaded=true`:

```julia
result = cp_significance(
    A,
    lip_discrete;
    null_model=:er,
    n_samples=1_000,
    threaded=true,
)
```

The package assigns a deterministic child RNG to each sample before scheduling,
so a fixed input RNG produces the same sample sequence in serial and threaded
runs. Julia 1.12 uses a benchmark-informed mix of greedy task-local scheduling
and static workers by default; pass `thread_schedule=:static` for uniform
workloads or explicitly use `:greedy` for irregular null-model/detector runtimes. Julia 1.10 retains the
static generic path. When nesting threaded significance tests around algorithms that use BLAS,
set BLAS to one thread to avoid oversubscription:

```julia
using LinearAlgebra
BLAS.set_num_threads(1)
```

## Cross-package benchmark

CorePeriphery.jl was compared with Python
[`cpnet`](https://github.com/skojaku/core-periphery-detection) commit
`6aad458a6d434a3617d33e74f7163d514a27fecb` on identical sparse matrices. The
comparison covers every algorithm configuration with a counterpart in both packages,
including discrete Borgatti--Everett (BE). The three fixtures are small planted
correctness cases:

| Dataset | Nodes | Undirected edges | Density |
|:--------|------:|-----------------:|--------:|
| `ideal_single` | 20 | 85 | 44.7% |
| `noisy_single` | 40 | 162 | 20.8% |
| `two_pairs` | 48 | 124 | 11.0% |

Each package was warmed once, after which the reported runtime was the median of five
fits. Numba compilation was excluded, RNG state was reset before each cpnet fit, both
packages received sparse matrices, and BLAS was restricted to one thread. In the table,
each fixture cell is **CorePeriphery.jl ms / cpnet ms (cpnet/Julia ratio)**; a ratio
above 1 means that the Julia call completed sooner.

| Algorithm | Fit budget (Julia / cpnet) | `ideal_single` | `noisy_single` | `two_pairs` | Median ratio | Julia faster |
|:----------|:----------------------------|---------------:|---------------:|------------:|-------------:|-------------:|
| BE | 1 deterministic fit / 5 starts | 0.004 / 0.139 (36.2×) | 0.008 / 1.172 (141.8×) | 0.008 / 1.539 (204.0×) | 141.8× | 3/3 |
| Lip | 1 deterministic fit / 1 deterministic fit | 0.004 / 0.066 (17.7×) | 0.006 / 0.058 (10.3×) | 0.005 / 0.049 (9.5×) | 10.3× | 3/3 |
| LowRank-Core | 1 deterministic fit / 1 deterministic fit | 0.080 / 2.690 (33.6×) | 0.343 / 3.458 (10.1×) | 0.216 / 4.311 (20.0×) | 20.0× | 3/3 |
| Rombach | 5 starts / 5 starts | 0.064 / 0.646 (10.0×) | 2.378 / 1.341 (0.6×) | 3.295 / 2.247 (0.7×) | 0.7× | 1/3 |
| Rossa | 1 fit / 1 fit | 0.004 / 1.600 (416.1×) | 0.015 / 2.495 (171.7×) | 0.017 / 3.783 (222.6×) | 222.6× | 3/3 |
| MINRES | 1 deterministic fit / 5 starts | 0.007 / 4.160 (623.8×) | 0.022 / 4.273 (196.8×) | 0.113 / 5.699 (50.5×) | 196.8× | 3/3 |
| Surprise | 1 deterministic fit / 5 starts | 0.007 / 12.910 (1834.7×) | 0.932 / 48.760 (52.3×) | 0.285 / 43.590 (152.9×) | 152.9× | 3/3 |
| KM-ER | 5 starts / 5 starts | 0.069 / 1.039 (15.1×) | 0.180 / 4.236 (23.6×) | 0.178 / 4.850 (27.2×) | 23.6× | 3/3 |
| KM-configuration | 5 starts / 5 starts | 0.066 / 1.313 (19.9×) | 0.234 / 5.127 (21.9×) | 0.155 / 5.342 (34.4×) | 21.9× | 3/3 |

Across these 27 fitted calls, CorePeriphery.jl was faster in 25 and the median of the
27 per-fixture cpnet/Julia ratios was 33.6×. Rombach was the exception: cpnet was faster
on the two larger fixtures. This is a comparison of each package's configured estimator,
not equal optimizer work. In particular, cpnet uses five starts for BE, MINRES, and
Surprise where CorePeriphery.jl uses a deterministic fit. The largest ratios are
therefore not per-restart speedups.

The comparison also has three important limits:

- these `n = 20--48` measurements are descriptive microbenchmarks, and the many
  sub-millisecond Julia entries are sensitive to timer noise and host load;
- the median ratio in each algorithm row is the median of three paired ratios, not a
  ratio of runtimes pooled across differently sized fixtures;
- identically named methods can expose different objectives or node scores. The MINRES
  row deliberately uses [`minres_symmetric`](@ref) in Julia to match cpnet's symmetric
  off-diagonal `ww'` model; it does not benchmark [`minres_svd_directed`](@ref).

See the
[`COMPARISON.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/COMPARISON.md)
report for recovery and concordance metrics, environment versions, and the complete
27-row result table. The raw fitted outputs and reproduction scripts are retained in
the repository's
[`benchmark/`](https://github.com/simoneSantoni/CorePeriphery.jl/tree/main/benchmark)
directory.

## Julia scaling and regression benchmark

The independent Julia harness separates dense and sparse scaling, warms compilation
before timing, records allocations, and pins BLAS to one thread. Run it from the
repository root as described in
[`benchmark/README.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/README.md).

For reliable measurements, use a fresh Julia process, report allocations as well as
elapsed time, and compare scaling over several network sizes and densities. GPU
execution is not currently exposed: most greedy algorithms are sequential, while the
matrix-free sparse CPU paths avoid transfer and launch overhead that would dominate
typical network sizes.

CorePeriphery.jl supports Julia 1.10 and newer. Julia 1.12 is the primary benchmark and
documentation runtime; Julia 1.10 is the compatibility floor. Run
`benchmark/run_versioned_performance.sh` to regenerate the Julia 1.10/1.12 TSV and
Markdown snapshot. The current
[`VERSIONED_PERFORMANCE.md`](https://github.com/simoneSantoni/CorePeriphery.jl/blob/main/benchmark/VERSIONED_PERFORMANCE.md)
report includes seven-repetition ranges through `n = 512`, source/CPU/thread/BLAS
provenance, allocations, dense/sparse LowRank-Core, Surprise, directed Rossa, symmetric
and directed MINRES, penalized multi-pair fitting, and significance scheduling.
