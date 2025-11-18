using Test
using Graphs
using Random

# Include the module
include("../src/CorePeriphery.jl")
using .CorePeriphery

@testset "CorePeriphery.jl Tests" begin

    @testset "Basic functionality" begin
        # Test with a simple graph that has clear core-periphery structure
        # Create a graph: core (nodes 1-3 densely connected) and periphery (nodes 4-6 connected to core)
        g = SimpleGraph(6)

        # Core edges (complete among 1,2,3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)

        # Core-periphery edges
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 6)

        config = KMConfig(n_runs=5, n_random_networks=50, seed=42)
        result = detect_core_periphery(g; config=config)

        @test result.total_quality isa Float64
        @test length(result.node_assignments) == 6
        @test length(result.node_roles) == 6
        @test all(x -> x ∈ [0, 1], result.node_roles)
    end

    @testset "Empty graph" begin
        g = SimpleGraph(5)
        config = KMConfig(n_runs=1, n_random_networks=10, seed=42)
        result = detect_core_periphery(g; config=config)

        @test result.total_quality == 0.0
        @test length(result.pairs) == 0
        @test length(result.residual_nodes) == 5
    end

    @testset "Complete graph" begin
        g = complete_graph(5)
        config = KMConfig(n_runs=3, n_random_networks=20, seed=42)
        result = detect_core_periphery(g; config=config)

        @test result.total_quality isa Float64
        @test length(result.node_assignments) == 5
    end

    @testset "Two-community structure" begin
        # Create graph with two clear communities that should be detected as
        # two core-periphery pairs
        g = SimpleGraph(10)

        # Community 1: nodes 1-5
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 5)

        # Community 2: nodes 6-10
        add_edge!(g, 6, 7)
        add_edge!(g, 6, 8)
        add_edge!(g, 7, 8)
        add_edge!(g, 6, 9)
        add_edge!(g, 7, 10)

        # Weak inter-community connection
        add_edge!(g, 3, 6)

        config = KMConfig(n_runs=5, n_random_networks=30, seed=42)
        result = detect_core_periphery(g; config=config)

        @test length(result.pairs) >= 1
        @test result.total_quality isa Float64
    end

    @testset "Star graph" begin
        # Star graph: one hub connected to all others
        g = star_graph(6)
        config = KMConfig(n_runs=3, n_random_networks=20, seed=42)
        result = detect_core_periphery(g; config=config)

        @test result.total_quality isa Float64
        @test length(result.node_assignments) == 6
    end

    @testset "Quality calculation" begin
        # Test that quality is properly calculated
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)

        config = KMConfig(n_runs=3, n_random_networks=20, seed=42)
        result = detect_core_periphery(g; config=config)

        @test quality(result) == result.total_quality

        if !isempty(result.pairs)
            @test pair_quality(result.pairs[1]) == result.pairs[1].quality
        end
    end

    @testset "Utility functions" begin
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 6)

        config = KMConfig(n_runs=3, n_random_networks=20, seed=42)
        result = detect_core_periphery(g; config=config)

        @test n_significant_pairs(result) >= 0
        @test n_significant_pairs(result) <= length(result.pairs)

        if !isempty(result.pairs)
            core = get_core_nodes(result, 1)
            periphery = get_periphery_nodes(result, 1)
            @test core isa Vector{Int}
            @test periphery isa Vector{Int}
        end

        # Test with invalid index
        @test isempty(get_core_nodes(result, 999))
        @test isempty(get_periphery_nodes(result, 0))
    end

    @testset "Karate club network" begin
        # Test on a classic benchmark: Zachary's karate club
        # Using the known structure
        g = SimpleGraph(34)

        # Karate club edges
        edges = [
            (1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,11), (1,12),
            (1,13), (1,14), (1,18), (1,20), (1,22), (1,32), (2,3), (2,4), (2,8),
            (2,14), (2,18), (2,20), (2,22), (2,31), (3,4), (3,8), (3,9), (3,10),
            (3,14), (3,28), (3,29), (3,33), (4,8), (4,13), (4,14), (5,7), (5,11),
            (6,7), (6,11), (6,17), (7,17), (9,31), (9,33), (9,34), (10,34),
            (14,34), (15,33), (15,34), (16,33), (16,34), (19,33), (19,34),
            (20,34), (21,33), (21,34), (23,33), (23,34), (24,26), (24,28),
            (24,30), (24,33), (24,34), (25,26), (25,28), (25,32), (26,32),
            (27,30), (27,34), (28,34), (29,32), (29,34), (30,33), (30,34),
            (31,33), (31,34), (32,33), (32,34), (33,34)
        ]

        for (i, j) in edges
            add_edge!(g, i, j)
        end

        config = KMConfig(n_runs=5, n_random_networks=50, seed=42)
        result = detect_core_periphery(g; config=config)

        @test nv(g) == 34
        @test ne(g) == length(edges)
        @test result.total_quality > 0  # Should find some structure
        @test length(result.pairs) >= 1

        # The karate club should have at least one core-periphery pair
        total_in_pairs = sum(length(p.core_nodes) + length(p.periphery_nodes)
                            for p in result.pairs)
        @test total_in_pairs > 0
    end

    @testset "Configuration model sampling" begin
        # Test internal configuration model sampling
        degrees = [4, 3, 3, 2, 2, 2]

        # Run multiple times to check consistency
        for _ in 1:5
            g_rand = CorePeriphery._configuration_model_sample(degrees)
            @test nv(g_rand) == 6
            # Due to self-loop removal, actual edges may differ
            @test ne(g_rand) >= 0
            @test ne(g_rand) <= sum(degrees) ÷ 2
        end
    end

    @testset "Reproducibility with seed" begin
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 5)
        add_edge!(g, 3, 6)

        config1 = KMConfig(n_runs=3, n_random_networks=20, seed=12345)
        config2 = KMConfig(n_runs=3, n_random_networks=20, seed=12345)

        result1 = detect_core_periphery(g; config=config1)
        result2 = detect_core_periphery(g; config=config2)

        @test result1.total_quality == result2.total_quality
        @test result1.node_assignments == result2.node_assignments
        @test result1.node_roles == result2.node_roles
    end

    @testset "Data structures" begin
        # Test CorePeripheryPair
        pair = CorePeripheryPair([1, 2], [3, 4, 5], 0.5, true, 0.01)
        @test length(pair.core_nodes) == 2
        @test length(pair.periphery_nodes) == 3
        @test pair.quality == 0.5
        @test is_significant(pair)
        @test pair.p_value == 0.01

        # Test KMConfig
        config = KMConfig()
        @test config.n_runs == 10
        @test config.n_random_networks == 500
        @test config.significance_level == 0.05
        @test config.seed === nothing

        config2 = KMConfig(n_runs=5, seed=42)
        @test config2.n_runs == 5
        @test config2.seed == 42
    end

end

println("\nAll tests passed!")
