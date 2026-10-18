using CorePeriphery
using LinearAlgebra
using Random
using SparseArrays
using Test

@testset "adjacency_to_matrix: duplicate dyads and direction" begin
    # Repeating a dyad with the same weight is harmless.
    @test adjacency_to_matrix([(1, 2), (2, 1)], 2) == [0.0 1.0; 1.0 0.0]
    @test adjacency_to_matrix([(1, 2, 2.0), (2, 1, 2.0)], 2) == [0.0 2.0; 2.0 0.0]
    # Conflicting weights on an undirected dyad are an error, not last-wins.
    @test_throws ArgumentError adjacency_to_matrix([(1, 2, 5.0), (2, 1, 1.0)], 3)
    @test_throws ArgumentError adjacency_to_matrix([(1, 2, 5.0), (1, 2, 1.0)], 3)
    @test_throws ArgumentError adjacency_to_matrix([(1, 2), (1, 2, 3.0)], 3)
    @test_throws ArgumentError adjacency_to_matrix(Any[(1, 2, 5.0), (2, 1)], 3)
    # Explicit resolution rules.
    edges = [(1, 2, 5.0), (2, 1, 1.0), (2, 3, 2.0)]
    @test adjacency_to_matrix(edges, 3; combine=:sum)[1, 2] == 6.0
    @test adjacency_to_matrix(edges, 3; combine=:max)[1, 2] == 5.0
    @test adjacency_to_matrix(edges, 3; combine=:min)[1, 2] == 1.0
    @test adjacency_to_matrix(edges, 3; combine=:last)[1, 2] == 1.0
    @test adjacency_to_matrix(edges, 3; combine=:last)[2, 1] == 1.0
    @test adjacency_to_matrix(edges, 3; combine=:sum)[2, 3] == 2.0
    @test adjacency_to_matrix([(1, 2, 2.0), (1, 2, 2.0)], 2; combine=:sum)[1, 2] == 4.0
    @test_throws ArgumentError adjacency_to_matrix(edges, 3; combine=:mean)
    # Directed: (i, j) and (j, i) are different dyads.
    D = adjacency_to_matrix([(1, 2, 5.0), (2, 1, 1.0), (3, 1)], 3; directed=true)
    @test D == [0.0 5.0 0.0; 1.0 0.0 0.0; 1.0 0.0 0.0]
    @test adjacency_to_matrix([(1, 2), (3, 1)], 3; directed=true) ==
          [0.0 1.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0]
    @test_throws ArgumentError adjacency_to_matrix(
        [(1, 2, 5.0), (1, 2, 1.0)], 3; directed=true)
    @test adjacency_to_matrix([(1, 2, 5.0), (1, 2, 1.0)], 3; directed=true, combine=:sum)[1, 2] == 6.0
    @test_throws ArgumentError adjacency_to_matrix([(1, 1)], 2; directed=true)
    @test_throws BoundsError adjacency_to_matrix([(1, 3)], 2; directed=true)
    @test_throws ArgumentError adjacency_to_matrix([(1, 2, -1.0)], 2; directed=true)
    @test minres_svd_directed(D) isa CPDirectedResult
end

@testset "adjacency_to_matrix: labelled edge lists" begin
    A, labels = adjacency_to_matrix([("b", "a"), ("a", "c", 2.0), ("c", "b")])
    @test labels == ["b", "a", "c"]
    @test labels isa Vector{String}
    @test A == [0.0 1.0 1.0; 1.0 0.0 2.0; 1.0 2.0 0.0]
    @test issymmetric(A)

    sorted_A, sorted_labels = adjacency_to_matrix(
        [("b", "a"), ("a", "c", 2.0), ("c", "b")]; sort_labels=true)
    @test sorted_labels == ["a", "b", "c"]
    @test sorted_A == [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0]

    S, symbols = adjacency_to_matrix([(:x, :y), (:y, :z)]; directed=true)
    @test symbols == [:x, :y, :z]
    @test S == [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]

    # Integer identifiers are labels too (first appearance, not value).
    I_A, ints = adjacency_to_matrix([(10, 20), (20, 5)])
    @test ints == [10, 20, 5]
    @test I_A == [0.0 1.0 0.0; 1.0 0.0 1.0; 0.0 1.0 0.0]
    @test adjacency_to_matrix([(10, 20), (20, 5)]; sort_labels=true)[2] == [5, 10, 20]

    @test_throws ArgumentError adjacency_to_matrix([("a", "a")])
    @test_throws ArgumentError adjacency_to_matrix([("a", "b", 2.0), ("b", "a", 3.0)])
    @test adjacency_to_matrix([("a", "b", 2.0), ("b", "a", 3.0)]; combine=:max)[1][1, 2] == 3.0
    @test_throws ArgumentError adjacency_to_matrix([("a", "b", -1.0)])
    @test_throws ArgumentError adjacency_to_matrix(Tuple{String,String}[])
    @test_throws ArgumentError adjacency_to_matrix([("a", "b", "c", "d")])

    # Labelled input flows straight into a detector (the matrix is weighted,
    # so use a weight-aware detector).
    result = borgatti_everett_discrete(A)
    @test labels[result.core_nodes] isa Vector{String}
    @test length(labels[result.core_nodes]) == length(result.core_nodes)
end

@testset "block_densities" begin
    # Core 1:2 complete; node 1 -- 3, node 2 -- 4; no periphery-periphery edge.
    A = adjacency_to_matrix([(1, 2), (1, 3), (2, 4)], 4)
    result = CPResult([1.0, 1.0, 0.0, 0.0], [1, 2], [3, 4], 0.0, "test")
    densities = block_densities(A, result)
    @test densities == (core_core=1.0, core_periphery=0.5, periphery_core=0.5,
                        periphery_periphery=0.0)
    @test block_densities(A, [1, 2]) == densities
    @test block_densities(sparse(A), result) == densities
    @test block_densities(A, [2, 1]) == densities
    @test keys(densities) == (:core_core, :core_periphery, :periphery_core, :periphery_periphery)

    # Detector results plug in directly, including the directed container.
    fit = lip_discrete(A)
    @test block_densities(A, fit) == densities
    @test block_densities(A, minres_svd_directed(A)) isa NamedTuple

    # Directed input distinguishes the two off-diagonal blocks.
    D = adjacency_to_matrix([(1, 2), (2, 1), (1, 3), (1, 4), (3, 2)], 4; directed=true)
    directed = block_densities(D, [1, 2])
    @test directed.core_core == 1.0
    @test directed.core_periphery == 0.5     # 1->3, 1->4 among 4 ordered dyads
    @test directed.periphery_core == 0.25    # 3->2 among 4 ordered dyads
    @test directed.periphery_periphery == 0.0

    # Weighted input reports mean weights; empty blocks are NaN.
    W = adjacency_to_matrix([(1, 2, 4.0), (1, 3, 2.0)], 3)
    weighted = block_densities(W, [1])
    @test isnan(weighted.core_core)
    @test weighted.core_periphery == 3.0
    @test weighted.periphery_core == 3.0
    @test weighted.periphery_periphery == 0.0
    everyone = block_densities(W, [1, 2, 3])
    @test isnan(everyone.core_periphery) && isnan(everyone.periphery_periphery)
    @test everyone.core_core == 2.0

    @test_throws ArgumentError block_densities(A, [0, 1])
    @test_throws ArgumentError block_densities(A, [5])
    @test_throws DimensionMismatch block_densities(zeros(2, 3), [1])
end
