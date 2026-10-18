using CorePeriphery
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Statistics
using Test

# Null samples are random graphs and are frequently disconnected, so the
# detectors' one-time disconnected-input warning fires on them while a test
# evaluates null samples. That warning is legitimate but not under test in this
# file. This logger drops exactly that record (by its `_id`) and forwards every
# other message to the enclosing logger, so unexpected warnings still surface.
struct DropDisconnectedLogger <: Logging.AbstractLogger
    parent::Logging.AbstractLogger
end
Logging.min_enabled_level(logger::DropDisconnectedLogger) =
    Logging.min_enabled_level(logger.parent)
Logging.catch_exceptions(logger::DropDisconnectedLogger) =
    Logging.catch_exceptions(logger.parent)
function Logging.shouldlog(logger::DropDisconnectedLogger, level, _module, group, id)
    startswith(string(id), "disconnected_input_") && return false
    return Logging.shouldlog(logger.parent, level, _module, group, id)
end
Logging.handle_message(logger::DropDisconnectedLogger, args...; kwargs...) =
    Logging.handle_message(logger.parent, args...; kwargs...)
without_disconnected_warnings(f) =
    Logging.with_logger(f, DropDisconnectedLogger(Logging.current_logger()))

function significance_planted_network(n::Int; seed::Int=20260713)
    rng = MersenneTwister(seed + n)
    core_size = max(2, round(Int, 0.15 * n))
    A = zeros(Float64, n, n)
    for j in 2:n, i in 1:(j - 1)
        probability = i <= core_size && j <= core_size ? 0.45 :
                      (i <= core_size || j <= core_size ? 0.08 : 0.01)
        rand(rng) < probability && (A[i, j] = A[j, i] = 1.0)
    end
    for i in 1:(n - 1)
        A[i, i + 1] = A[i + 1, i] = 1.0
    end
    return A
end

# Every testset in this file evaluates detectors on random null samples; see
# `DropDisconnectedLogger` above. The closing `end` is at the bottom of the file.
without_disconnected_warnings() do

@testset "Core-periphery significance" begin
    A = adjacency_to_matrix(
        [
            (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
            (6, 7), (7, 8), (8, 1), (1, 5), (2, 6),
        ],
        8,
    )

    first_result = cp_significance(
        A,
        lip_discrete;
        n_samples=12,
        rng=MersenneTwister(91),
        n_swaps=20,
    )
    second_result = cp_significance(
        A,
        lip_discrete;
        n_samples=12,
        rng=MersenneTwister(91),
        n_swaps=20,
    )
    @test first_result == second_result
    @test length(first_result.null_qualities) == 12
    @test 1 / 13 <= first_result.pvalue <= 1
    @test first_result.null_model == :configuration
    @test first_result.null_diagnostics isa CorePeriphery.NullDiagnostics
    @test first_result.null_diagnostics.complete
    @test first_result.null_diagnostics.accepted_swaps == fill(20, 12)
    @test first_result.null_diagnostics.preserved == :degree_sequence
    @test hash(first_result) == hash(second_result)
    @test cp_significance(
        A, lip_discrete; n_samples=12, rng=MersenneTwister(91), n_swaps=20, alpha=1//20,
    ) == first_result
    @test_throws ArgumentError cp_significance(A, lip_discrete; alpha=0)
    @test_throws ArgumentError cp_significance(A, lip_discrete; alpha=1)

    rewired = CorePeriphery._sample_configuration(
        A,
        MersenneTwister(12);
        n_swaps=50,
    )
    @test vec(sum(rewired, dims=2)) == vec(sum(A, dims=2))

    er_result = cp_significance(
        A,
        borgatti_everett_discrete;
        null_model=:er,
        n_samples=5,
        rng=MersenneTwister(3),
    )
    @test er_result.null_model == :er

    for null_model in (:er, :configuration)
        serial = cp_significance(
            A,
            label_switching_cp;
            null_model=null_model,
            n_samples=16,
            rng=MersenneTwister(404),
            n_swaps=20,
            detector_kwargs=(max_iter=3, n_runs=2),
            pass_rng=true,
            threaded=false,
        )
        parallel = cp_significance(
            A,
            label_switching_cp;
            null_model=null_model,
            n_samples=16,
            rng=MersenneTwister(404),
            n_swaps=20,
            detector_kwargs=(max_iter=3, n_runs=2),
            pass_rng=true,
            threaded=true,
        )
        @test isequal(serial, parallel)

        static = cp_significance(
            A,
            label_switching_cp;
            null_model=null_model,
            n_samples=16,
            rng=MersenneTwister(404),
            n_swaps=20,
            detector_kwargs=(max_iter=3, n_runs=2),
            pass_rng=true,
            threaded=true,
            thread_schedule=:static,
        )
        @test isequal(serial, static)
        # The work-queue schedule is available on every supported Julia version.
        greedy = cp_significance(
            A,
            label_switching_cp;
            null_model=null_model,
            n_samples=16,
            rng=MersenneTwister(404),
            n_swaps=20,
            detector_kwargs=(max_iter=3, n_runs=2),
            pass_rng=true,
            threaded=true,
            thread_schedule=:greedy,
        )
        @test isequal(serial, greedy)
    end

    # In-place samplers reset their buffers and therefore reproduce the public
    # allocating wrappers exactly for the same RNG stream.
    density = length(CorePeriphery._edge_list(A)) / (8 * 7 / 2)
    er_buffer = fill(9.0, 8, 8)
    er_inplace = CorePeriphery._sample_er!(er_buffer, density, MersenneTwister(8))
    er_allocating = CorePeriphery._sample_er(A, MersenneTwister(8))
    @test er_inplace == er_allocating

    base_edges = CorePeriphery._edge_list(A)
    config_buffer = fill(9.0, 8, 8)
    config_edges = copy(base_edges)
    config_inplace = CorePeriphery._sample_configuration!(
        config_buffer, A, config_edges, base_edges, MersenneTwister(13), 30,
    )
    config_allocating = CorePeriphery._sample_configuration(
        A, MersenneTwister(13); n_swaps=30,
    )
    @test config_inplace == config_allocating
    @test vec(sum(config_inplace, dims=2)) == vec(sum(A, dims=2))
    @test CorePeriphery._sample_er(sparse(A), MersenneTwister(8)) == er_allocating
    @test CorePeriphery._sample_configuration(
        sparse(A), MersenneTwister(13); n_swaps=30,
    ) == config_allocating

    # Some degree sequences admit no simple double-edge switch. The strict
    # default rejects a scientifically ineffective null, while explicit
    # acceptance returns complete diagnostics instead of failing silently.
    star = adjacency_to_matrix([(1, node) for node in 2:7], 7)
    @test_throws ConvergenceError CorePeriphery._sample_configuration(
        star, MersenneTwister(1); n_swaps=1,
    )
    @test_throws ConvergenceError cp_significance(
        star, lip_discrete; n_samples=2, n_swaps=1, max_swap_attempts=25,
    )
    short = CorePeriphery._sample_configuration(
        star,
        MersenneTwister(1);
        n_swaps=1,
        max_swap_attempts=25,
        swap_shortfall=:accept,
        return_stats=true,
    )
    @test short.network == star
    @test short.stats.requested == 1
    @test short.stats.accepted == 0
    @test short.stats.attempts == 25
    @test !CorePeriphery._swap_complete(short.stats)

    # Exercise both ends of the switchability spectrum. A four-cycle has only
    # four edges, while K6 minus a perfect matching is nearly complete. Both
    # must either finish the requested chain or fail explicitly; with the
    # declared budget these fixtures finish and preserve a simple graph.
    low_edge = adjacency_to_matrix([(1, 2), (2, 3), (3, 4), (4, 1)], 4)
    nearly_complete = ones(Float64, 6, 6) - Matrix{Float64}(I, 6, 6)
    for (left, right) in ((1, 2), (3, 4), (5, 6))
        nearly_complete[left, right] = nearly_complete[right, left] = 0.0
    end
    for constrained in (low_edge, nearly_complete)
        sampled = CorePeriphery._sample_configuration(
            constrained,
            MersenneTwister(2026);
            n_swaps=10,
            max_swap_attempts=10_000,
            return_stats=true,
        )
        @test sampled.stats.accepted == 10
        @test vec(sum(sampled.network, dims=2)) ==
              vec(sum(constrained, dims=2))
        @test all(iszero, diag(sampled.network))
        @test all(value -> value in (0.0, 1.0), sampled.network)
        @test issymmetric(sampled.network)
        @test CorePeriphery._sample_configuration(
            sparse(constrained),
            MersenneTwister(2026);
            n_swaps=10,
            max_swap_attempts=10_000,
        ) == sampled.network
    end

    directed = zeros(Int, 8, 8)
    for node in 1:8
        directed[node, mod1(node + 1, 8)] = 1
        directed[node, mod1(node + 3, 8)] = 1
    end
    directed_null = CorePeriphery._sample_directed_configuration(
        directed, MersenneTwister(4); n_swaps=30, return_stats=true,
    )
    @test directed_null.stats.accepted == 30
    @test vec(sum(directed_null.network, dims=2)) == vec(sum(directed, dims=2))
    @test vec(sum(directed_null.network, dims=1)) == vec(sum(directed, dims=1))
    @test all(iszero, diag(directed_null.network))
    @test all(value -> value in (0.0, 1.0), directed_null.network)
    @test CorePeriphery._sample_directed_configuration(
        sparse(directed), MersenneTwister(4); n_swaps=30,
    ) == directed_null.network

    directed_significance = cp_significance(
        directed,
        minres_svd_directed;
        null_model=:directed_configuration,
        n_samples=6,
        n_swaps=15,
        rng=MersenneTwister(81),
    )
    @test directed_significance.null_diagnostics.complete
    @test directed_significance.null_diagnostics.preserved ==
          :in_out_degree_sequences

    weighted = adjacency_to_matrix(
        [(1, 2, 0.5), (1, 3, 1.5), (2, 3, 2.5),
         (2, 4, 3.5), (3, 5, 4.5)],
        5,
    )
    weighted_null = CorePeriphery._sample_weight_permutation(
        weighted, MersenneTwister(12),
    )
    @test issymmetric(weighted_null)
    @test (weighted_null .!= 0.0) == (weighted .!= 0.0)
    @test sort(nonzeros(sparse(triu(weighted_null, 1)))) ==
          sort(nonzeros(sparse(triu(weighted, 1))))
    weighted_significance = cp_significance(
        weighted,
        borgatti_everett_continuous;
        null_model=:weight_permutation,
        n_samples=6,
        rng=MersenneTwister(90),
        detector_kwargs=(max_iter=3,),
    )
    @test weighted_significance.null_diagnostics.complete
    @test weighted_significance.null_diagnostics.preserved ==
          :topology_and_weight_multiset

    for (network, detector, model, detector_kwargs) in (
        (directed, minres_svd_directed, :directed_configuration, NamedTuple()),
        (weighted, borgatti_everett_continuous, :weight_permutation, (max_iter=3,)),
    )
        serial_extended = cp_significance(
            network,
            detector;
            null_model=model,
            n_samples=8,
            n_swaps=10,
            rng=MersenneTwister(707),
            detector_kwargs=detector_kwargs,
            threaded=false,
        )
        threaded_extended = cp_significance(
            network,
            detector;
            null_model=model,
            n_samples=8,
            n_swaps=10,
            rng=MersenneTwister(707),
            detector_kwargs=detector_kwargs,
            threaded=true,
        )
        # `isequal`: degenerate nulls give a NaN z-score, which is not `==` to itself.
        @test isequal(serial_extended, threaded_extended)
    end

    @test_throws ArgumentError cp_significance(A, lip_discrete; n_samples=0)
    @test_throws ArgumentError cp_significance(A, lip_discrete; null_model=:bad)
    @test_throws ArgumentError cp_significance(
        A, lip_discrete; thread_schedule=:invalid,
    )
    @test_throws ArgumentError cp_significance(
        A, lip_discrete; null_model=:configuration, n_swaps=-1,
    )
    @test_throws ArgumentError cp_significance(
        A, lip_discrete; max_swap_attempts=-1,
    )
    @test_throws ArgumentError cp_significance(
        A, lip_discrete; swap_shortfall=:invalid,
    )
    incompatible = try
        cp_significance(
            directed,
            lip_discrete;
            null_model=:directed_configuration,
            n_samples=2,
        )
        nothing
    catch error
        error
    end
    @test incompatible isa ArgumentError
    @test occursin("incompatible with null_model=:directed_configuration",
                   sprint(showerror, incompatible))
end

@testset "Significance defaults, effect sizes, and messages" begin
    A = adjacency_to_matrix(
        [
            (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
            (6, 7), (7, 8), (8, 1), (1, 5), (2, 6),
        ],
        8,
    )

    # Defaults: 999 samples and the child RNG forwarded to detectors that
    # accept one. `label_switching_cp` with pass_rng (the default) is therefore
    # reproducible from `rng` alone.
    default_result = cp_significance(A, lip_discrete; rng=MersenneTwister(5), n_swaps=10)
    @test length(default_result.null_qualities) == 999
    @test 1 / 1000 <= default_result.pvalue <= 1
    repeat_a = cp_significance(
        A, label_switching_cp; n_samples=20, n_swaps=10, rng=MersenneTwister(77),
        detector_kwargs=(max_iter=3, n_runs=3),
    )
    repeat_b = cp_significance(
        A, label_switching_cp; n_samples=20, n_swaps=10, rng=MersenneTwister(77),
        detector_kwargs=(max_iter=3, n_runs=3),
    )
    @test isequal(repeat_a, repeat_b)

    # Effect sizes follow the null sample moments.
    result = cp_significance(A, lip_discrete; n_samples=30, rng=MersenneTwister(9), n_swaps=10)
    @test result.null_mean ≈ mean(result.null_qualities)
    @test result.null_std ≈ std(result.null_qualities)
    @test result.zscore ≈ (result.observed_quality - result.null_mean) / result.null_std
    @test propertynames(result) == (
        :observed_quality, :null_qualities, :pvalue, :significant, :alpha,
        :null_model, :null_diagnostics, :null_mean, :null_std, :zscore,
    )
    shown = sprint(show, MIME("text/plain"), result)
    @test occursin("null mean", shown)
    @test occursin("z-score", shown)
    single = cp_significance(A, lip_discrete; n_samples=1, rng=MersenneTwister(9), n_swaps=10)
    @test isnan(single.null_std)
    @test isnan(single.zscore)
    @test single.null_mean == single.null_qualities[1]
    # A degenerate null (every sample identical) gives a NaN z-score, not ±Inf.
    star = adjacency_to_matrix([(1, node) for node in 2:7], 7)
    degenerate = cp_significance(
        star, lip_discrete; null_model=:configuration, n_samples=3, n_swaps=1,
        max_swap_attempts=5, swap_shortfall=:accept, rng=MersenneTwister(1),
    )
    @test degenerate.null_std == 0.0
    @test isnan(degenerate.zscore)

    # pass_rng: package detectors without an rng keyword are silently called
    # without it; foreign detectors warn once.
    @test_logs min_level = Base.CoreLogging.Warn cp_significance(
        A, lip_discrete; n_samples=2, n_swaps=5, rng=MersenneTwister(2),
    )
    foreign_detector = network -> lip_discrete(network)
    @test_logs (:warn, r"does not accept an `rng` keyword") match_mode = :any cp_significance(
        A, foreign_detector; n_samples=2, n_swaps=5, rng=MersenneTwister(2),
    )
    @test_logs min_level = Base.CoreLogging.Warn cp_significance(
        A, network -> lip_discrete(network); n_samples=2, n_swaps=5,
        rng=MersenneTwister(2), pass_rng=false,
    )

    # Actionable errors name the null model and the way out.
    weighted = adjacency_to_matrix([(1, 2, 0.5), (1, 3, 1.5), (2, 3, 2.5)], 3)
    binary_error = try
        cp_significance(weighted, lip_discrete; n_samples=2)
        nothing
    catch error
        error
    end
    @test binary_error isa ArgumentError
    @test occursin("null_model=:configuration", sprint(showerror, binary_error))
    @test occursin(":weight_permutation", sprint(showerror, binary_error))
    shortfall_error = try
        cp_significance(star, lip_discrete; n_samples=2, n_swaps=1, max_swap_attempts=25)
        nothing
    catch error
        error
    end
    @test shortfall_error isa ConvergenceError
    @test occursin("null_model=:er", sprint(showerror, shortfall_error))
    @test occursin("swap_shortfall=:accept", sprint(showerror, shortfall_error))
    @test occursin("no valid double-edge swap", sprint(showerror, shortfall_error))
end

@testset "Significance: sparse input uses sparse null samples" begin
    planted = significance_planted_network(48)
    sparse_planted = sparse(planted)
    for (null_model, detector, kwargs) in (
        (:er, lip_discrete, NamedTuple()),
        (:configuration, borgatti_everett_discrete, NamedTuple()),
        (:configuration, label_switching_cp, (max_iter=3, n_runs=2)),
    )
        dense_result = cp_significance(
            planted, detector; null_model=null_model, n_samples=10, n_swaps=50,
            rng=MersenneTwister(31), detector_kwargs=kwargs,
        )
        sparse_result = cp_significance(
            sparse_planted, detector; null_model=null_model, n_samples=10, n_swaps=50,
            rng=MersenneTwister(31), detector_kwargs=kwargs,
        )
        @test isequal(dense_result, sparse_result)
        @test dense_result.null_qualities == sparse_result.null_qualities
        @test dense_result.null_diagnostics == sparse_result.null_diagnostics
    end

    directed = zeros(Int, 10, 10)
    for node in 1:10
        directed[node, mod1(node + 1, 10)] = 1
        directed[node, mod1(node + 4, 10)] = 1
    end
    dense_directed = cp_significance(
        directed, minres_svd_directed; null_model=:directed_configuration,
        n_samples=6, n_swaps=20, rng=MersenneTwister(8),
    )
    sparse_directed = cp_significance(
        sparse(directed), minres_svd_directed; null_model=:directed_configuration,
        n_samples=6, n_swaps=20, rng=MersenneTwister(8),
    )
    @test isequal(dense_directed, sparse_directed)

    weighted = adjacency_to_matrix(
        [(1, 2, 0.5), (1, 3, 1.5), (2, 3, 2.5), (2, 4, 3.5), (3, 5, 4.5), (4, 5, 1.0)],
        5,
    )
    dense_weighted = cp_significance(
        weighted, borgatti_everett_continuous; null_model=:weight_permutation,
        n_samples=6, rng=MersenneTwister(90), detector_kwargs=(max_iter=3,),
    )
    sparse_weighted = cp_significance(
        sparse(weighted), borgatti_everett_continuous; null_model=:weight_permutation,
        n_samples=6, rng=MersenneTwister(90), detector_kwargs=(max_iter=3,),
    )
    @test isequal(dense_weighted, sparse_weighted)

    # The sparse samplers produce exactly the dense samples, as valid CSC
    # matrices with sorted row indices.
    n = size(planted, 1)
    base_edges = CorePeriphery._edge_list(sparse_planted)
    for null_model in (:er, :configuration)
        spec = CorePeriphery._NullSpec(
            null_model, sparse_planted, base_edges, Float64[],
            length(base_edges) / (n * (n - 1) ÷ 2), 50, 5_000, :error, true,
        )
        dense_spec = CorePeriphery._NullSpec(
            null_model, planted, base_edges, Float64[],
            length(base_edges) / (n * (n - 1) ÷ 2), 50, 5_000, :error, true,
        )
        sparse_ws = CorePeriphery._sparse_workspace(n, base_edges, Float64[])
        dense_ws = CorePeriphery._dense_workspace(n, base_edges, Float64[])
        for seed in 1:4
            CorePeriphery._draw_null!(sparse_ws, spec, MersenneTwister(seed))
            CorePeriphery._draw_null!(dense_ws, dense_spec, MersenneTwister(seed))
            @test sparse_ws.network isa SparseMatrixCSC{Float64,Int}
            @test Matrix(sparse_ws.network) == dense_ws.network
            @test issymmetric(sparse_ws.network)
            rows = rowvals(sparse_ws.network)
            @test all(issorted(rows[nzrange(sparse_ws.network, column)])
                      for column in 1:n)
            @test nnz(sparse_ws.network) == count(!iszero, dense_ws.network)
        end
    end
end

@testset "Significance: schedules are bit-for-bit identical" begin
    # Regression test for the 0.2.0 greedy race: repeated threaded runs must
    # reproduce the serial result exactly for every null model and schedule.
    planted = significance_planted_network(40; seed=991)
    inputs = (planted, sparse(planted))
    for (null_model, detector, kwargs) in (
        (:er, lip_discrete, NamedTuple()),
        (:configuration, lip_discrete, NamedTuple()),
        (:configuration, label_switching_cp, (max_iter=3, n_runs=2)),
    ), network in inputs
        serial = cp_significance(
            network, detector; null_model=null_model, n_samples=12, n_swaps=60,
            rng=MersenneTwister(919), detector_kwargs=kwargs, threaded=false,
        )
        mismatches = 0
        for repetition in 1:25, schedule in (:static, :greedy, :auto)
            threaded = cp_significance(
                network, detector; null_model=null_model, n_samples=12, n_swaps=60,
                rng=MersenneTwister(919), detector_kwargs=kwargs, threaded=true,
                thread_schedule=schedule,
            )
            isequal(threaded, serial) || (mismatches += 1)
        end
        @test mismatches == 0
    end
end

end # without_disconnected_warnings
