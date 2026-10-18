using CorePeriphery
using LinearAlgebra
using Random
using SparseArrays
using Statistics

BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "regression_targets.jl"))

const QUICK = "--quick" in ARGS
const CHECK = "--check" in ARGS
# Five repetitions in quick mode and the minimum over repetitions as the ratio
# statistic keep the doubling-ratio gates insensitive to shared-runner noise.
const REPETITIONS = QUICK ? 5 : 7

function planted_network(n::Int; seed::Int=20260713)
    rng = MersenneTwister(seed + n)
    core_size = max(2, round(Int, 0.15 * n))
    A = zeros(Float64, n, n)
    for j in 2:n, i in 1:(j - 1)
        probability = i <= core_size && j <= core_size ? 0.45 :
                      (i <= core_size || j <= core_size ? 0.08 : 0.01)
        if rand(rng) < probability
            A[i, j] = A[j, i] = 1.0
        end
    end
    # Connectivity makes the same fixture valid for persistence profiling.
    for i in 1:(n - 1)
        A[i, i + 1] = A[i + 1, i] = 1.0
    end
    return A
end

function directed_network(n::Int; seed::Int=20260713)
    rng = MersenneTwister(seed + 3 * n)
    A = zeros(Float64, n, n)
    for node in 1:n
        A[node, mod1(node + 1, n)] = 1.0
        A[node, mod1(node + 2, n)] = 0.5 + rand(rng)
    end
    for _ in 1:(2 * n)
        source, target = rand(rng, 1:n, 2)
        source == target && continue
        A[source, target] = 0.25 + rand(rng)
    end
    return A
end

function measure(factory; repetitions::Int=REPETITIONS)
    factory()() # compile and initialize the specialization outside measurements
    times = Vector{Float64}(undef, repetitions)
    bytes = Vector{Int}(undef, repetitions)
    for repetition in 1:repetitions
        GC.gc()
        measurement = @timed factory()()
        times[repetition] = measurement.time
        bytes[repetition] = measurement.bytes
    end
    # `seconds` (the gate statistic) is the minimum over repetitions: it is the
    # least noisy estimate of the intrinsic cost on a shared runner.
    return (
        seconds=minimum(times),
        median_seconds=median(times),
        maximum_seconds=maximum(times),
        bytes=round(Int, median(bytes)),
    )
end

function record!(rows, name, n, storage, measurement)
    push!(rows, (name=name, n=n, storage=storage,
                 minimum_milliseconds=1_000 * measurement.seconds,
                 median_milliseconds=1_000 * measurement.median_seconds,
                 maximum_milliseconds=1_000 * measurement.maximum_seconds,
                 bytes=measurement.bytes))
end

sizes = QUICK ? (128, 256) : (128, 256, 512)
rows = NamedTuple[]
quality_times = Dict{Tuple{Int,String},Float64}()
rossa_times = Dict{Tuple{Int,String},Float64}()
surprise_times = Dict{Tuple{Int,String},Float64}()
failures = String[]

for n in sizes
    dense = planted_network(n)
    sparse_matrix = sparse(dense)
    directed_dense = directed_network(n)
    directed_sparse = sparse(directed_dense)
    c = collect(range(0.0, stop=1.0, length=n))
    dense_quality = core_quality(dense, c)
    sparse_quality = core_quality(sparse_matrix, c)
    if !isapprox(dense_quality, sparse_quality; atol=SPARSE_DENSE_QUALITY_ATOL, rtol=0.0)
        push!(failures, "dense/sparse quality differs at n=$n")
    end

    for (storage, A) in (("dense", dense), ("sparse", sparse_matrix))
        quality = measure(() -> () -> core_quality(A, c))
        record!(rows, "core_quality", n, storage, quality)
        quality_times[(n, storage)] = quality.seconds
        quality.bytes <= CORE_QUALITY_MAX_BYTES ||
            push!(failures, "core_quality allocated $(quality.bytes) bytes at n=$n/$storage")

        lip = measure(() -> () -> lip_discrete(A))
        record!(rows, "lip_discrete", n, storage, lip)

        symmetric_minres = measure(
            () -> () -> minres_symmetric(A; max_iter=50, tol=1e-6))
        record!(rows, "minres_symmetric_50", n, storage, symmetric_minres)

        lowrank = measure(() -> () -> spectral_method(A; beta=0.1))
        record!(rows, "lowrank_core", n, storage, lowrank)

        surprise = measure(() -> () -> surprise_cp(A; max_iter=10))
        record!(rows, "surprise_10", n, storage, surprise)

        # Full greedy run with the mode-centred tail kernel: the per-flip cost
        # is governed by the hypergeometric peak widths, not O(m²).
        surprise_full = measure(() -> () -> surprise_cp(A))
        record!(rows, "surprise_full", n, storage, surprise_full)
        surprise_times[(n, storage)] = surprise_full.seconds
    end

    directed_dense_result = random_walker_profiling(
        directed_dense; stationary_tol=1e-9)
    directed_sparse_result = random_walker_profiling(
        directed_sparse; stationary_tol=1e-9)
    isapprox(directed_dense_result.coreness, directed_sparse_result.coreness;
             atol=1e-8, rtol=0.0) ||
        push!(failures, "dense/sparse directed Rossa differs at n=$n")
    for (storage, A) in (("dense", directed_dense), ("sparse", directed_sparse))
        rossa = measure(
            () -> () -> random_walker_profiling(A; stationary_tol=1e-9))
        record!(rows, "rossa_directed", n, storage, rossa)
        rossa_times[(n, storage)] = rossa.seconds

        directed_minres = measure(
            () -> () -> minres_svd_directed(A; max_iter=50, tol=1e-6))
        record!(rows, "minres_directed_50", n, storage, directed_minres)
    end
end

pair_n = QUICK ? 48 : 96
pair_A = planted_network(pair_n; seed=191)
pair_measurement = measure(() -> begin
    () -> multiple_cp_pairs(
        pair_A;
        max_pairs=3,
        max_iter=8,
        n_runs=2,
        pair_selection=:penalized,
        rng=MersenneTwister(717),
    )
end)
record!(rows, "multiple_pairs_penalized", pair_n, "dense", pair_measurement)

for storage in ("dense", "sparse")
    for index in 2:length(sizes)
        smaller, larger = sizes[index - 1], sizes[index]
        ratio = quality_times[(larger, storage)] / quality_times[(smaller, storage)]
        ratio <= CORE_QUALITY_MAX_DOUBLING_RATIO || push!(
            failures,
            "core_quality time ratio $ratio exceeds target for $storage $smaller->$larger",
        )
    end
end


for storage in ("dense", "sparse")
    for index in 2:length(sizes)
        smaller, larger = sizes[index - 1], sizes[index]
        ratio = rossa_times[(larger, storage)] / rossa_times[(smaller, storage)]
        ratio <= ROSSA_MAX_DOUBLING_RATIO || push!(
            failures,
            "directed Rossa time ratio $ratio exceeds target for $storage $smaller->$larger",
        )
    end
end

for storage in ("dense", "sparse")
    for index in 2:length(sizes)
        smaller, larger = sizes[index - 1], sizes[index]
        ratio = surprise_times[(larger, storage)] / surprise_times[(smaller, storage)]
        ratio <= SURPRISE_MAX_DOUBLING_RATIO || push!(
            failures,
            "surprise_cp time ratio $ratio exceeds target for $storage $smaller->$larger",
        )
    end
end

significance_n = QUICK ? 64 : 128
significance_samples = QUICK ? 12 : 30
significance_A = planted_network(significance_n; seed=991)
significance_sparse = sparse(significance_A)
for null_model in (:er, :configuration)
    # Null samples are drawn in the storage of the input: the sparse path never
    # materializes an n×n buffer and detectors run their sparse kernels.
    for (storage, A) in (("dense", significance_A), ("sparse", significance_sparse))
        factory = () -> begin
            rng = MersenneTwister(818)
            () -> cp_significance(
                A,
                lip_discrete;
                null_model=null_model,
                n_samples=significance_samples,
                n_swaps=100,
                rng=rng,
                threaded=Threads.nthreads() > 1,
            )
        end
        measurement = measure(factory)
        record!(rows, "significance_$(null_model)", significance_n, storage, measurement)
        bytes_per_sample = measurement.bytes / significance_samples
        target = storage == "sparse" ? SIGNIFICANCE_SPARSE_MAX_BYTES_PER_SAMPLE :
                 null_model === :er ? SIGNIFICANCE_ER_MAX_BYTES_PER_SAMPLE :
                 SIGNIFICANCE_CONFIGURATION_MAX_BYTES_PER_SAMPLE
        bytes_per_sample <= target || push!(
            failures,
            "significance $null_model/$storage allocated $bytes_per_sample bytes/sample",
        )
    end

    serial = cp_significance(
        significance_A, lip_discrete; null_model=null_model,
        n_samples=significance_samples, n_swaps=100,
        rng=MersenneTwister(919), threaded=false,
    )
    parallel = cp_significance(
        significance_A, lip_discrete; null_model=null_model,
        n_samples=significance_samples, n_swaps=100,
        rng=MersenneTwister(919), threaded=true,
    )
    # `isequal`, not `==`: a degenerate null has a NaN z-score.
    isequal(serial, parallel) ||
        push!(failures, "serial/threaded $null_model results differ")
    sparse_serial = cp_significance(
        significance_sparse, lip_discrete; null_model=null_model,
        n_samples=significance_samples, n_swaps=100,
        rng=MersenneTwister(919), threaded=false,
    )
    isequal(serial, sparse_serial) ||
        push!(failures, "dense/sparse $null_model significance results differ")

    if Threads.nthreads() > 1
        for schedule in (:static, :greedy)
            schedule_factory = () -> begin
                rng = MersenneTwister(818)
                () -> cp_significance(
                    significance_A,
                    lip_discrete;
                    null_model=null_model,
                    n_samples=significance_samples,
                    n_swaps=100,
                    rng=rng,
                    threaded=true,
                    thread_schedule=schedule,
                )
            end
            schedule_measurement = measure(schedule_factory)
            record!(rows, "significance_$(null_model)_$(schedule)",
                significance_n, "dense", schedule_measurement)
            scheduled = cp_significance(
                significance_A, lip_discrete; null_model=null_model,
                n_samples=significance_samples, n_swaps=100,
                rng=MersenneTwister(919), threaded=true, thread_schedule=schedule,
            )
            isequal(serial, scheduled) ||
                push!(failures, "serial/$schedule $null_model results differ")
        end
    end
end

println("benchmark\tn\tstorage\tminimum_ms\tmedian_ms\tmaximum_ms\tbytes")
for row in rows
    println(join((row.name, row.n, row.storage, row.minimum_milliseconds,
                  row.median_milliseconds, row.maximum_milliseconds, row.bytes), '\t'))
end
println("context\tjulia=$(VERSION)\tthreads=$(Threads.nthreads())" *
        "\tblas_threads=$(BLAS.get_num_threads())\tcpu=$(Sys.CPU_NAME)" *
        "\tword_size=$(Sys.WORD_SIZE)\tquick=$QUICK\trepetitions=$REPETITIONS")

if CHECK && !isempty(failures)
    error("benchmark regression gate failed:\n- " * join(failures, "\n- "))
elseif CHECK
    println("regression_gate\tpass")
end
