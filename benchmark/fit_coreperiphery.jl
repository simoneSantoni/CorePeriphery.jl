using CorePeriphery
using LinearAlgebra
using Random
using SparseArrays
using Statistics

const ROOT = @__DIR__
const DATASETS = ("ideal_single", "noisy_single", "two_pairs")
BLAS.set_num_threads(1)

function read_matrix(path)
    rows = [parse.(Float64, split(line, ',')) for line in readlines(path)]
    return reduce(vcat, permutedims.(rows))
end

encode(values) = join(Float64.(values), ';')

function factories(seed)
    return [
        "BE" => () -> A -> borgatti_everett_discrete(A),
        "Lip" => () -> A -> lip_discrete(A),
        "LowRankCore" => () -> A -> spectral_method(A; beta=0.1),
        "Rombach" => () -> A -> rombach_continuous(
            A; alpha=0.5, beta=0.8, n_runs=5, rng=MersenneTwister(seed + 3)),
        "Rossa" => () -> A -> random_walker_profiling(A),
        # cpnet.MINRES fits the symmetric off-diagonal w*w' model. Use the
        # matching Julia estimator here; minres_svd_directed is intentionally
        # reserved for directed-network comparisons.
        "MINRES" => () -> A -> minres_symmetric(A),
        "Surprise" => () -> A -> surprise_cp(A),
        "KM_ER" => () -> A -> multiple_cp_pairs(
            A; null_model=:er, n_runs=5, rng=MersenneTwister(seed + 7)),
        "KM_config" => () -> A -> multiple_cp_pairs(
            A; null_model=:configuration, n_runs=5, rng=MersenneTwister(seed + 8)),
    ]
end

function fit(factory, A)
    factory()(A)
    result = nothing
    timings = Float64[]
    for _ in 1:5
        detector = factory()
        push!(timings, @elapsed result = detector(A))
    end
    return result, median(timings) * 1000
end

open(joinpath(ROOT, "coreperiphery_results.tsv"), "w") do io
    println(io, "package\tdataset\talgorithm\tquality\truntime_ms\tcoreness\tpairs")
    for (dataset_index, dataset) in enumerate(DATASETS)
        A = sparse(read_matrix(joinpath(ROOT, "data", "$dataset.csv")))
        for (name, factory) in factories(1000 + 100 * (dataset_index - 1))
            try
                result, elapsed = fit(factory, A)
                pairs = result isa CPMultiResult ? result.pair_labels : ones(Int, size(A, 1))
                println(
                    io,
                    join(
                        ("CorePeriphery.jl", dataset, name, result.quality, elapsed,
                         encode(result.coreness), encode(pairs)),
                        '\t',
                    ),
                )
            catch error
                println(
                    io,
                    join(
                        ("CorePeriphery.jl", dataset, name, "NaN", "NaN", "",
                         "ERROR:$(typeof(error)):$(sprint(showerror, error))"),
                        '\t',
                    ),
                )
            end
        end
    end
end

open(joinpath(ROOT, "coreperiphery_environment.txt"), "w") do io
    println(io, "CorePeriphery.jl ", Base.pkgversion(CorePeriphery))
    println(io, "Julia ", VERSION)
    println(io, "Julia threads ", Threads.nthreads())
    println(io, "BLAS threads ", BLAS.get_num_threads())
end
