# Cross-package empirical comparison

The independent correctness-oracle and provenance matrix is documented in
`VALIDATION.md`. Cross-package concordance below is complementary evidence, not a
substitute for mathematical validation.

This benchmark fits CorePeriphery.jl and
[skojaku/core-periphery-detection](https://github.com/skojaku/core-periphery-detection)
to identical ideal, noisy single-pair, and planted two-pair networks.

Run:

```sh
python3 benchmark/generate_networks.py
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=. benchmark/fit_coreperiphery.jl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python3 benchmark/fit_cpnet.py /path/to/core-periphery-detection
python3 benchmark/analyze.py
```

The Python environment must contain the upstream package's dependencies. The
tested 2022 cpnet revision requires NumPy 1.x (the benchmark used NumPy 1.26.4
and SciPy 1.14.1).
`analyze.py` writes `comparison.tsv` and `COMPARISON.md`. Raw quality values
remain package-specific; agreement is evaluated using rank correlation,
truth AUC, top-k Jaccard similarity, and adjusted Rand indices.
Both packages receive sparse matrices. The Python harness resets Python,
NumPy, and Numba RNG state before each fit, and the generated comparison records
the exact language, dependency, and upstream-commit versions.
The row named `MINRES` fits the symmetric off-diagonal `w*w'` model in both packages;
the directed Julia `u*v'` estimator is deliberately excluded from that row.

## Scaling and regression gate

Independently of Python, `scaling.jl` is a dependency-free, warmed benchmark for dense and sparse
quality/Lip kernels, LowRank-Core, symmetric and directed MINRES, Surprise,
directed Rossa, penalized multi-pair fitting, plus ER and configuration significance
sampling. It prints tab-separated warmed medians, ranges, and allocations and checks objective equivalence,
allocation ceilings, and generous size-doubling limits:

```sh
julia --project=. benchmark/scaling.jl
julia --threads=4 --project=. benchmark/scaling.jl --quick --check
```

The thresholds in `regression_targets.jl` deliberately gate complexity and
allocation regressions rather than machine-specific absolute timings. For
publication-quality performance work, run the full size grid on an idle host,
record Julia/CPU/BLAS/thread settings, and retain the TSV output with the commit
being measured.

Generate the checked Julia 1.10/1.12 full snapshot and Markdown summary with:

```sh
benchmark/run_versioned_performance.sh
```

Pass output paths followed by `--quick` for the smaller CI-sized snapshot.

The implementation-phase before/after measurements and their interpretation are
recorded in [`OPTIMIZATION.md`](OPTIMIZATION.md).

The algorithm-specific investigation of lower-concordance methods is in
[`DIVERGENCE.md`](DIVERGENCE.md). Reproduce its targeted numerical probes with:

```sh
/path/to/python-with-numpy-scipy-numba benchmark/diagnose_divergence.py \
    /path/to/core-periphery-detection
```

The script writes long-form measurements to `divergence_diagnostics.tsv`.

## Performance-refactor ablation

The ablation separates the committed pre-refactor implementation, the current
algorithms on dense matrices, and the current sparse paths. From the repository root:

```sh
git worktree add --detach /tmp/CorePeriphery-pre-refactor \
    a2226832508c57efdeab71cde1f055549ca028e7
julia --project=/tmp/CorePeriphery-pre-refactor benchmark/fit_ablation.jl \
    pre_refactor dense benchmark/ablation_pre_refactor.tsv
julia --project=. benchmark/fit_ablation.jl \
    current dense benchmark/ablation_current_dense.tsv
julia --project=. benchmark/fit_ablation.jl \
    current sparse benchmark/ablation_current_sparse.tsv
/path/to/python-with-numpy-and-scipy benchmark/analyze_ablation.py
```

The generated [`ABLATION.md`](ABLATION.md) reports three-way concordance, semantic
drift, and dense/sparse kernel equivalence. Machine-readable three-way measurements
are written to `ablation_three_way.tsv`. The pre-refactor spectral method and current
LowRank-Core implementation are not mathematically equivalent, and the old package
did not contain a configuration-null multi-pair detector; the report preserves these
limitations explicitly. All legs use three warmed timed fits and the same ten-iteration
Rombach budget; the bounded budget is necessary because the old coordinate-grid search
otherwise performs billions of allocations on the benchmark graphs.
