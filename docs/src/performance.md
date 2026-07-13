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

## Benchmarking

The `benchmark/` directory contains the cross-package comparison and the
performance regression harness. Benchmarks separate dense and sparse scaling,
warm up compilation before timing, and pin BLAS to one thread. Run the harness
from the repository root as described in `benchmark/README.md`.

For reliable measurements, use a fresh Julia process, report allocations as well
as elapsed time, and compare scaling over several network sizes and densities.
GPU execution is not currently exposed: most greedy algorithms are sequential,
while the matrix-free sparse CPU paths avoid the transfer and launch overhead that
would dominate typical network sizes.

CorePeriphery.jl supports Julia 1.10 and newer. Julia 1.12 is the primary
benchmark and documentation runtime; Julia 1.10 is the compatibility floor.

Run `benchmark/run_versioned_performance.sh` to regenerate the Julia 1.10/1.12 TSV
and Markdown snapshot. The current report is `benchmark/VERSIONED_PERFORMANCE.md`;
it includes seven-repetition ranges through `n=512`, source/CPU/thread/BLAS provenance,
dense/sparse LowRank-Core, Surprise, directed Rossa, symmetric and directed MINRES,
penalized multi-pair fitting, and significance scheduling.
