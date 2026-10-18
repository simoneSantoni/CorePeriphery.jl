using CorePeriphery
using Graphs
using SparseArrays
using Test

@testset "Graphs.jl extension" begin
    graph = SimpleGraph(6)
    for edge in [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4)]
        add_edge!(graph, edge...)
    end

    A = adjacency_to_matrix(graph)
    @test A isa SparseMatrixCSC{Float64,Int}
    @test size(A) == (6, 6)
    @test lip_discrete(graph).coreness == lip_discrete(A).coreness
    @test spectral_method(graph; beta=0.2).coreness ==
          spectral_method(A; beta=0.2).coreness

    digraph = SimpleDiGraph(3)
    add_edge!(digraph, 1, 2)
    add_edge!(digraph, 1, 3)
    directed = minres_svd_directed(digraph)
    @test length(directed.out_coreness) == 3

    rossa_digraph = SimpleDiGraph(4)
    for edge in [(1, 2), (2, 3), (3, 4), (4, 1), (1, 3)]
        add_edge!(rossa_digraph, edge...)
    end
    rossa_matrix = adjacency_to_matrix(rossa_digraph)
    @test random_walker_profiling(rossa_digraph).coreness ≈
          random_walker_profiling(rossa_matrix).coreness
end
