"""
Unit tests for CorePeriphery module.
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using CorePeriphery
using Test
using LinearAlgebra
using Random

Random.seed!(123)

@testset "CorePeriphery Tests" begin

    @testset "Utility Functions" begin
        # Test adjacency_to_matrix
        edges = [(1, 2), (2, 3), (1, 3)]
        A = adjacency_to_matrix(edges, 3)
        @test size(A) == (3, 3)
        @test A[1, 2] == 1.0
        @test A[2, 1] == 1.0
        @test A[1, 3] == 1.0
        @test A[2, 3] == 1.0
        @test A[1, 1] == 0.0  # No self-loops

        # Test weighted edges
        weighted_edges = [(1, 2, 2.5), (2, 3, 1.5)]
        A_w = adjacency_to_matrix(weighted_edges, 3)
        @test A_w[1, 2] == 2.5
        @test A_w[2, 3] == 1.5

        # Test ideal_cp_matrix continuous
        c = [1.0, 0.5, 0.0]
        Δ = ideal_cp_matrix(c)
        @test Δ[1, 2] == 0.5  # 1.0 * 0.5
        @test Δ[1, 3] == 0.0  # 1.0 * 0.0
        @test Δ[2, 3] == 0.0  # 0.5 * 0.0
        @test Δ[1, 1] == 0.0  # Diagonal is zero

        # Test ideal_cp_matrix discrete
        Δ_d = ideal_cp_matrix(c; discrete=true)
        @test Δ_d[1, 2] == 1.0  # max(1.0, 0.5)
        @test Δ_d[1, 3] == 1.0  # max(1.0, 0.0)
        @test Δ_d[2, 3] == 0.5  # max(0.5, 0.0)
    end

    @testset "CPResult Structure" begin
        result = CPResult(
            [0.8, 0.9, 0.2, 0.1],
            [1, 2],
            [3, 4],
            0.75,
            "Test Algorithm"
        )
        @test length(result.coreness) == 4
        @test result.core_nodes == [1, 2]
        @test result.periphery_nodes == [3, 4]
        @test result.quality == 0.75
        @test result.algorithm == "Test Algorithm"

        # Test coreness_scores function
        scores = coreness_scores(result)
        @test scores == [0.8, 0.9, 0.2, 0.1]
    end

    @testset "Borgatti-Everett Continuous" begin
        # Create simple CP network
        # Nodes 1, 2 are core; nodes 3, 4, 5 are periphery
        edges = [
            (1, 2),  # Core-core
            (1, 3), (1, 4), (1, 5),  # Core-periphery
            (2, 3), (2, 4), (2, 5),  # Core-periphery
        ]
        A = adjacency_to_matrix(edges, 5)

        result = borgatti_everett_continuous(A)

        @test length(result.coreness) == 5
        @test result.quality > 0  # Should find positive correlation
        @test result.algorithm == "Borgatti-Everett Continuous"

        # Core nodes should have higher coreness
        @test result.coreness[1] > result.coreness[3]
        @test result.coreness[2] > result.coreness[4]
    end

    @testset "Borgatti-Everett Discrete" begin
        edges = [
            (1, 2), (1, 3), (1, 4),
            (2, 3), (2, 4),
            (5, 1)  # Periphery connection
        ]
        A = adjacency_to_matrix(edges, 5)

        result = borgatti_everett_discrete(A)

        @test all(x -> x == 0.0 || x == 1.0, result.coreness)  # Binary
        @test !isempty(result.core_nodes)
        @test !isempty(result.periphery_nodes)
        @test length(result.core_nodes) + length(result.periphery_nodes) == 5
    end

    @testset "Lip Discrete" begin
        # Star graph: center is core, leaves are periphery
        n = 6
        edges = [(1, i) for i in 2:n]
        A = adjacency_to_matrix(edges, n)

        result = lip_discrete(A)

        @test length(result.coreness) == n
        @test result.algorithm == "Lip Discrete"

        # Node 1 (center) should be core
        @test 1 in result.core_nodes
    end

    @testset "Spectral Method" begin
        # Complete graph on core, sparse on periphery
        A = zeros(Float64, 6, 6)
        # Core: nodes 1, 2, 3 (complete)
        for i in 1:3
            for j in (i+1):3
                A[i, j] = A[j, i] = 1.0
            end
        end
        # Periphery connections
        A[4, 1] = A[1, 4] = 1.0
        A[5, 2] = A[2, 5] = 1.0
        A[6, 3] = A[3, 6] = 1.0

        result = spectral_method(A)

        @test length(result.coreness) == 6
        @test result.algorithm == "Spectral Method"
        @test result.quality >= -1.0 && result.quality <= 1.0  # Valid correlation
    end

    @testset "Random Walker Profiling" begin
        # Hub-and-spoke: hub should be core
        n = 5
        A = zeros(Float64, n, n)
        # Node 1 is hub
        for i in 2:n
            A[1, i] = A[i, 1] = 1.0
        end

        result = random_walker_profiling(A; n_walks=500, walk_length=5)

        @test length(result.coreness) == n
        @test result.algorithm == "Random Walker Profiling"

        # Hub (node 1) should have highest coreness
        @test argmax(result.coreness) == 1
    end

    @testset "Rombach Continuous" begin
        # Simple test with small network
        edges = [(1, 2), (1, 3), (2, 3), (4, 1)]
        A = adjacency_to_matrix(edges, 4)

        result = rombach_continuous(A; n_runs=2, max_iter=50)

        @test length(result.coreness) == 4
        @test result.algorithm == "Rombach Continuous"
        @test all(x -> 0.0 <= x <= 1.0, result.coreness)
    end

    @testset "MINRES/SVD" begin
        # Test with symmetric network
        edges = [
            (1, 2), (1, 3), (1, 4),
            (2, 3), (2, 4),
            (5, 1)
        ]
        A = adjacency_to_matrix(edges, 5)

        result = minres_svd(A)

        @test length(result.coreness) == 5
        @test result.algorithm == "MINRES/SVD"
        @test all(x -> 0.0 <= x <= 1.0, result.coreness)

        # Test with asymmetric network (directed)
        A_asym = zeros(Float64, 4, 4)
        A_asym[1, 2] = 1.0
        A_asym[1, 3] = 1.0
        A_asym[2, 3] = 1.0
        A_asym[4, 1] = 1.0

        result_asym = minres_svd(A_asym)
        @test length(result_asym.coreness) == 4
    end

    @testset "Multiple CP Pairs" begin
        # Create network with two distinct CP structures
        A = zeros(Float64, 8, 8)
        # First CP pair: nodes 1-4
        for i in 1:2
            for j in (i+1):4
                A[i, j] = A[j, i] = 1.0
            end
        end
        A[3, 1] = A[1, 3] = 1.0
        A[4, 2] = A[2, 4] = 1.0

        # Second CP pair: nodes 5-8
        for i in 5:6
            for j in (i+1):8
                A[i, j] = A[j, i] = 1.0
            end
        end
        A[7, 5] = A[5, 7] = 1.0
        A[8, 6] = A[6, 8] = 1.0

        result = multiple_cp_pairs(A)

        @test length(result.pair_labels) == 8
        @test length(result.coreness) == 8
        @test result.n_pairs >= 1
        @test result.algorithm == "Multiple CP Pairs"
    end

    @testset "Surprise CP" begin
        # Star graph: hub should be core
        n = 6
        A = zeros(Float64, n, n)
        for i in 2:n
            A[1, i] = A[i, 1] = 1.0
        end

        result = surprise_cp(A)

        @test length(result.coreness) == n
        @test result.algorithm == "Surprise CP"
        @test all(x -> x == 0.0 || x == 1.0, result.coreness)  # Binary

        # Hub should be core
        @test result.coreness[1] == 1.0
    end

    @testset "Label Switching CP" begin
        # Simple CP network
        edges = [
            (1, 2), (1, 3), (1, 4), (1, 5),
            (2, 3), (2, 4), (2, 5),
            (3, 4),
            (6, 1), (7, 2)
        ]
        A = adjacency_to_matrix(edges, 7)

        result = label_switching_cp(A)

        @test length(result.coreness) == 7
        @test result.algorithm == "Label Switching CP"
        @test all(x -> x == 0.0 || x == 1.0, result.coreness)  # Binary
        @test !isempty(result.core_nodes)
        @test !isempty(result.periphery_nodes)
    end

    @testset "Core Quality Function" begin
        # Perfect core-periphery structure
        c = [1.0, 1.0, 0.0, 0.0]
        A = zeros(Float64, 4, 4)
        # Core-core edges
        A[1, 2] = A[2, 1] = 1.0
        # Core-periphery edges
        A[1, 3] = A[3, 1] = 1.0
        A[2, 4] = A[4, 2] = 1.0
        # No periphery-periphery edges

        quality = core_quality(A, c; discrete=true)
        @test quality > 0.5  # Should be high quality

        # Poor assignment (swap labels)
        c_bad = [0.0, 0.0, 1.0, 1.0]
        quality_bad = core_quality(A, c_bad; discrete=true)
        @test quality > quality_bad  # Good assignment should be better
    end

    @testset "Edge Cases" begin
        # Empty network
        A_empty = zeros(Float64, 3, 3)
        result = spectral_method(A_empty)
        @test length(result.coreness) == 3

        # Single node
        A_single = zeros(Float64, 1, 1)
        result = borgatti_everett_continuous(A_single)
        @test length(result.coreness) == 1

        # Complete graph (no CP structure)
        A_complete = ones(Float64, 4, 4) - I(4)
        result = borgatti_everett_discrete(A_complete)
        @test length(result.coreness) == 4
    end

    @testset "Consistency Checks" begin
        # Same network should give consistent results
        Random.seed!(999)
        edges = [(1, 2), (1, 3), (2, 3), (1, 4), (2, 5)]
        A = adjacency_to_matrix(edges, 5)

        result1 = borgatti_everett_continuous(A)
        result2 = borgatti_everett_continuous(A)

        # With same initialization, should get same result
        @test result1.coreness ≈ result2.coreness atol=1e-6
    end

end

println("\nAll tests passed!")
