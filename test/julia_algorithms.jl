using CorePeriphery
using Random
using Test

@testset "Julia algorithm corrections" begin
    @testset "Label switching optimizes discrete correlation" begin
        A = adjacency_to_matrix(
            [(1, 2), (1, 3), (1, 4), (1, 5),
             (2, 3), (2, 4), (2, 5), (3, 4),
             (6, 1), (7, 2)],
            7,
        )

        degrees = vec(sum(A, dims=2))
        order = sortperm(degrees; rev=true)
        initial = zeros(7)
        initial[order[1:3]] .= 1.0
        initial_quality = core_quality(A, initial; discrete=true)

        result = label_switching_cp(A; rng=MersenneTwister(42), n_runs=5)
        repeated = label_switching_cp(A; rng=MersenneTwister(42), n_runs=5)

        @test result.coreness == repeated.coreness
        @test result.quality == repeated.quality
        @test result.quality == core_quality(A, result.coreness; discrete=true)
        @test result.quality >= initial_quality
        @test 1 <= length(result.core_nodes) <= size(A, 1) - 2
        @test length(result.periphery_nodes) >= 2

        @test_throws ArgumentError label_switching_cp(zeros(2, 2))
        @test_throws ArgumentError label_switching_cp(A; max_iter=-1)
        @test_throws ArgumentError label_switching_cp(A; n_runs=0)
    end

    @testset "Multiple-pair constraints" begin
        A = zeros(8, 8)
        for offset in (0, 4)
            core = (offset + 1):(offset + 2)
            pair = (offset + 1):(offset + 4)
            for i in core, j in pair
                i == j && continue
                A[i, j] = A[j, i] = 1.0
            end
        end

        capped = multiple_cp_pairs(A; max_pairs=2, min_pair_size=1, max_iter=0)
        @test capped.n_pairs <= 2
        @test sort(unique(capped.pair_labels)) == collect(1:capped.n_pairs)

        sized = multiple_cp_pairs(A; max_pairs=4, min_pair_size=3, max_iter=0)
        @test sized.n_pairs <= 4
        @test all(count(==(pair), sized.pair_labels) >= 3 for pair in 1:sized.n_pairs)

        @test_throws ArgumentError multiple_cp_pairs(A; max_pairs=0)
        @test multiple_cp_pairs(A; max_pairs=9, min_pair_size=1).n_pairs <= 8
        @test_throws ArgumentError multiple_cp_pairs(A; min_pair_size=0)
        @test_throws ArgumentError multiple_cp_pairs(A; min_pair_size=9)
        @test_throws ArgumentError multiple_cp_pairs(A; max_iter=-1)
    end
end
