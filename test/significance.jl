using CorePeriphery
using LinearAlgebra
using Random
using SparseArrays
using Test

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
    @test first_result.null_diagnostics.complete
    @test first_result.null_diagnostics.accepted_swaps == fill(20, 12)

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
        @test serial == parallel

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
        @test serial == static
        if VERSION >= v"1.12"
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
            @test serial == greedy
        end
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
    @test_throws ErrorException CorePeriphery._sample_configuration(
        star, MersenneTwister(1); n_swaps=1,
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
        @test serial_extended == threaded_extended
    end

    @test_throws ArgumentError cp_significance(A, lip_discrete; n_samples=0)
    @test_throws ArgumentError cp_significance(A, lip_discrete; null_model=:bad)
    @test_throws ArgumentError cp_significance(
        A, lip_discrete; thread_schedule=:invalid,
    )
    if VERSION < v"1.12"
        @test_throws ArgumentError cp_significance(
            A, lip_discrete; threaded=true, thread_schedule=:greedy,
        )
    end
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
