using CorePeriphery
using Test

@testset "Scientific and input contracts" begin
    @testset "Borgatti-Everett objective and degeneracy" begin
        A = adjacency_to_matrix(
            [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)],
            5,
        )
        init = [0.8, 0.7, 0.4, 0.2, 0.1]
        initial_quality = core_quality(A, init)
        result = borgatti_everett_continuous(A; init=init, step=0.1)
        @test result.quality >= initial_quality
        @test result.quality == core_quality(A, result.coreness)

        empty_result = borgatti_everett_continuous(zeros(Int, 4, 4))
        @test empty_result.coreness == fill(0.5, 4)
        @test empty_result.quality == 0.0

        discrete = borgatti_everett_discrete(zeros(Int, 4, 4))
        @test !isempty(discrete.core_nodes)
        @test !isempty(discrete.periphery_nodes)
    end

    @testset "Directed quality uses ordered dyads" begin
        A = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
        c = [1.0, 0.4, 0.0]
        q = core_quality(A, c; directed=true)
        @test q ≈ core_quality(transpose(A), c; directed=true)

        permutation = [3, 1, 2]
        @test q ≈ core_quality(
            A[permutation, permutation],
            c[permutation];
            directed=true,
        ) atol=1e-14
    end

    @testset "Shared validation and generic matrices" begin
        A_int = [0 1 1; 1 0 0; 1 0 0]
        @test length(spectral_method(A_int).coreness) == 3
        @test_throws DimensionMismatch borgatti_everett_continuous(zeros(2, 3))
        @test_throws ArgumentError spectral_method([0.0 1.0; 0.0 0.0])
        @test_throws ArgumentError borgatti_everett_discrete([0.0 -1.0; -1.0 0.0])
        @test_throws ArgumentError surprise_cp([0.0 0.5; 0.5 0.0])
        @test_throws ArgumentError minres_svd([1.0 0.0; 0.0 0.0])
    end

    @testset "Exact surprise statistic" begin
        star = adjacency_to_matrix([(1, i) for i in 2:6], 6)
        result = surprise_cp(star)
        @test result.coreness[1] == 1.0
        @test result.quality > 0.0
        @test isfinite(result.quality)
    end
end
