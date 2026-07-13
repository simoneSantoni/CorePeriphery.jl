# Hardware-independent guardrails for `benchmark/scaling.jl --check`. These are
# intentionally generous: CI gates complexity/allocation regressions, not small
# wall-clock fluctuations between machines.
const CORE_QUALITY_MAX_BYTES = 256
const CORE_QUALITY_MAX_DOUBLING_RATIO = 12.0
const ROSSA_MAX_DOUBLING_RATIO = 20.0
const SIGNIFICANCE_ER_MAX_BYTES_PER_SAMPLE = 100_000
const SIGNIFICANCE_CONFIGURATION_MAX_BYTES_PER_SAMPLE = 120_000
const SPARSE_DENSE_QUALITY_ATOL = 1e-12
