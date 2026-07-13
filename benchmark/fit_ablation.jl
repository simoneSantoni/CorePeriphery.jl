using CorePeriphery
using LinearAlgebra
using Random
using Statistics

const ROOT = @__DIR__
const DATASETS = ("ideal_single", "noisy_single", "two_pairs")
const ALGORITHMS = (
    "BE", "Lip", "LowRankCore", "Rombach", "Rossa", "MINRES", "Surprise",
    "KM_ER", "KM_config",
)

length(ARGS) == 3 || error(
    "usage: fit_ablation.jl PRE_REFACTOR|CURRENT DENSE|SPARSE OUTPUT.tsv",
)

const VARIANT = Symbol(lowercase(ARGS[1]))
const STORAGE = Symbol(lowercase(ARGS[2]))
const OUTPUT = abspath(ARGS[3])
VARIANT in (:pre_refactor, :current) || error("unknown implementation variant: $(ARGS[1])")
STORAGE in (:dense, :sparse) || error("unknown storage variant: $(ARGS[2])")
VARIANT === :pre_refactor && STORAGE === :sparse &&
    error("the pre-refactor API accepted only Matrix{Float64}")

BLAS.set_num_threads(1)

function read_matrix(path)
    rows = [parse.(Float64, split(line, ',')) for line in readlines(path)]
    return reduce(vcat, permutedims.(rows))
end

function matrix_for_storage(A)
    STORAGE === :dense && return A
    sparse_arrays = getfield(CorePeriphery, :SparseArrays)
    return sparse_arrays.sparse(A)
end

encode(values) = join(Float64.(values), ';')
encode_indices(values) = join(Int.(values), ';')

function factory(name, seed)
    if name == "BE"
        return () -> A -> borgatti_everett_discrete(A)
    elseif name == "Lip"
        return () -> A -> lip_discrete(A)
    elseif name == "LowRankCore"
        # Before the refactor this public name was a leading-eigenvector method;
        # the current implementation is Cucuringu et al.'s LowRank-Core.
        return VARIANT === :pre_refactor ?
            (() -> A -> spectral_method(A)) :
            (() -> A -> spectral_method(A; beta=0.1))
    elseif name == "Rombach"
        return VARIANT === :pre_refactor ?
            (() -> A -> rombach_continuous(
                A; alpha=0.5, beta=0.8, max_iter=10, n_runs=5,
            )) :
            (() -> A -> rombach_continuous(
                A; alpha=0.5, beta=0.8, max_iter=10, n_runs=5,
                rng=MersenneTwister(seed),
            ))
    elseif name == "Rossa"
        return () -> A -> random_walker_profiling(A)
    elseif name == "MINRES"
        return () -> A -> minres_svd(A)
    elseif name == "Surprise"
        return () -> A -> surprise_cp(A)
    elseif name == "KM_ER"
        return VARIANT === :pre_refactor ?
            (() -> A -> multiple_cp_pairs(A)) :
            (() -> A -> multiple_cp_pairs(
                A; null_model=:er, n_runs=5, rng=MersenneTwister(seed),
            ))
    elseif name == "KM_config"
        VARIANT === :pre_refactor && return nothing
        return () -> A -> multiple_cp_pairs(
            A; null_model=:configuration, n_runs=5, rng=MersenneTwister(seed),
        )
    end
    error("unknown algorithm: $name")
end

function fit(make_detector, A, seed; repetitions=3)
    Random.seed!(seed)
    make_detector()(A)
    result = nothing
    timings = Float64[]
    for _ in 1:repetitions
        Random.seed!(seed)
        detector = make_detector()
        push!(timings, @elapsed result = detector(A))
    end
    return result, median(timings) * 1000
end

mkpath(dirname(OUTPUT))
open(OUTPUT, "w") do io
    println(
        io,
        "variant\tstorage\tdataset\talgorithm\tquality\truntime_ms\t" *
        "coreness\tpairs\tcore_nodes",
    )
    for (dataset_index, dataset) in enumerate(DATASETS)
        A = matrix_for_storage(read_matrix(joinpath(ROOT, "data", "$dataset.csv")))
        for (algorithm_index, name) in enumerate(ALGORITHMS)
            seed = 1000 + 100 * (dataset_index - 1) + algorithm_index
            make_detector = factory(name, seed)
            make_detector === nothing && continue
            try
                result, elapsed = fit(make_detector, A, seed)
                pairs = result isa CPMultiResult ? result.pair_labels : ones(Int, size(A, 1))
                core_nodes = hasproperty(result, :core_nodes) ? result.core_nodes : Int[]
                println(
                    io,
                    join(
                        (VARIANT, STORAGE, dataset, name, result.quality, elapsed,
                         encode(result.coreness), encode_indices(pairs),
                         encode_indices(core_nodes)),
                        '\t',
                    ),
                )
                flush(io)
            catch error
                println(
                    io,
                    join(
                        (VARIANT, STORAGE, dataset, name, "NaN", "NaN", "", "",
                         "ERROR:$(typeof(error)):$(sprint(showerror, error))"),
                        '\t',
                    ),
                )
                flush(io)
            end
        end
    end
end

environment_path = replace(OUTPUT, r"\.tsv$" => "_environment.txt")
open(environment_path, "w") do io
    println(io, "variant ", VARIANT)
    println(io, "storage ", STORAGE)
    println(io, "CorePeriphery.jl ", Base.pkgversion(CorePeriphery))
    println(io, "Julia ", VERSION)
    println(io, "Julia threads ", Threads.nthreads())
    println(io, "BLAS threads ", BLAS.get_num_threads())
    println(io, "active project ", Base.active_project())
end
