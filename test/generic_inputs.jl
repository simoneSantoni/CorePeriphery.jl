using CorePeriphery
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Test

@testset "Generic inputs: element types, wrappers, and tolerances" begin
    edges = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 6), (4, 7)]
    A = adjacency_to_matrix(edges, 7)
    S = sparse(A)

    @testset "ideal_cp_matrix accepts any real vector" begin
        reference = ideal_cp_matrix([1.0, 0.5, 0.0])
        @test ideal_cp_matrix([1, 0, 1]) == ideal_cp_matrix([1.0, 0.0, 1.0])
        @test ideal_cp_matrix([1, 0, 1]) isa Matrix{Float64}
        @test ideal_cp_matrix(@view([1.0, 0.5, 0.0, 9.0][1:3])) == reference
        @test ideal_cp_matrix(Float32[1, 0.5, 0]) == reference
        @test ideal_cp_matrix([1//1, 1//2, 0//1]; discrete=true) ==
              ideal_cp_matrix([1.0, 0.5, 0.0]; discrete=true)
        @test ideal_cp_matrix(1:3) == ideal_cp_matrix([1.0, 2.0, 3.0])
    end

    @testset "Real keyword arguments" begin
        reference = borgatti_everett_continuous(A)
        @test borgatti_everett_continuous(A; tol=1//1_000_000, step=1//20) == reference
        @test borgatti_everett_continuous(A; tol=1.0f-6, step=0.05f0).coreness ≈
              reference.coreness atol=1e-6
        @test borgatti_everett_continuous(A; tol=0, step=1).coreness ==
              borgatti_everett_continuous(A; tol=0.0, step=1.0).coreness
        @test_throws ArgumentError borgatti_everett_continuous(A; tol=Inf)

        rombach_reference = rombach_continuous(
            A; alpha=0.5, beta=0.5, n_runs=2, rng=MersenneTwister(5))
        @test rombach_continuous(
            A; alpha=1//2, beta=1//2, tol=1//1_000_000, n_runs=2, rng=MersenneTwister(5),
        ) == rombach_reference
        @test rombach_continuous(
            A; alpha=0.5f0, beta=0.5f0, n_runs=2, rng=MersenneTwister(5),
        ).coreness ≈ rombach_reference.coreness
        @test rombach_continuous(
            A; alpha=0, beta=1//2, n_runs=1, rng=MersenneTwister(5),
        ).coreness isa Vector{Float64}

        # Null samples may be disconnected (one-time detector warning); not under test.
        significance = Logging.with_logger(Logging.NullLogger()) do
            cp_significance(
                A, lip_discrete; n_samples=3, rng=MersenneTwister(1), n_swaps=10, alpha=1//10)
        end
        @test significance.alpha == 0.1
        @test significance.alpha isa Float64
    end

    @testset "Element types" begin
        reference = lip_discrete(A)
        for converted in (round.(Int, A), Float32.(A), A .> 0, BitMatrix(A .> 0))
            @test lip_discrete(converted).core_nodes == reference.core_nodes
            @test borgatti_everett_discrete(converted).coreness ==
                  borgatti_everett_discrete(A).coreness
            @test random_walker_profiling(converted).coreness ≈
                  random_walker_profiling(A).coreness
        end
        @test core_quality(Float32.(A), Float32.(reference.coreness); discrete=true) ≈
              core_quality(A, reference.coreness; discrete=true) atol=1e-6
        @test core_quality(@view(A[:, :]), @view(reference.coreness[:]); discrete=true) ==
              core_quality(A, reference.coreness; discrete=true)
    end

    @testset "Sparse wrappers are unwrapped to SparseMatrixCSC" begin
        @test CorePeriphery._unwrap_adjacency(Symmetric(S)) isa SparseMatrixCSC
        @test CorePeriphery._unwrap_adjacency(Hermitian(S)) isa SparseMatrixCSC
        @test CorePeriphery._unwrap_adjacency(transpose(S)) isa SparseMatrixCSC
        @test CorePeriphery._unwrap_adjacency(S') isa SparseMatrixCSC
        @test CorePeriphery._unwrap_adjacency(S) === S
        @test CorePeriphery._unwrap_adjacency(A) === A
        @test CorePeriphery._unwrap_adjacency(Symmetric(A)) isa Symmetric
        @test CorePeriphery._unwrap_adjacency(Symmetric(S)) == A
        @test CorePeriphery._unwrap_adjacency(transpose(S)) == A

        # Only the upper triangle of an asymmetric parent defines Symmetric.
        asymmetric_parent = sparse([0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0])
        @test CorePeriphery._unwrap_adjacency(Symmetric(asymmetric_parent)) ==
              Matrix(Symmetric(asymmetric_parent))

        wrapped = (Symmetric(S), transpose(S), S')
        for W in wrapped
            @test borgatti_everett_continuous(W) == borgatti_everett_continuous(A)
            @test borgatti_everett_discrete(W) == borgatti_everett_discrete(A)
            @test lip_discrete(W) == lip_discrete(A)
            @test random_walker_profiling(W) == random_walker_profiling(A)
            @test surprise_cp(W) == surprise_cp(A)
            @test spectral_method(W; beta=0.2) == spectral_method(A; beta=0.2)
            @test minres_symmetric(W).coreness ≈ minres_symmetric(A).coreness atol=1e-10
            @test label_switching_cp(W; rng=MersenneTwister(1)) ==
                  label_switching_cp(A; rng=MersenneTwister(1))
            @test core_quality(W, lip_discrete(A).coreness; discrete=true) ==
                  core_quality(A, lip_discrete(A).coreness; discrete=true)
            @test multiple_cp_pairs(W; n_runs=2, rng=MersenneTwister(2)) ==
                  multiple_cp_pairs(A; n_runs=2, rng=MersenneTwister(2))
        end
        directed = sparse([0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0])
        @test random_walker_profiling(transpose(directed)).coreness ==
              random_walker_profiling(copy(transpose(directed))).coreness
        @test minres_svd_directed(transpose(directed)).out_coreness ≈
              minres_svd_directed(Matrix(transpose(directed))).out_coreness atol=1e-10
    end

    @testset "Symmetry tolerance" begin
        # Floating-point noise far below 1e-12 does not make a graph directed.
        noisy = copy(A)
        noisy[1, 2] += 1e-15
        @test !issymmetric(noisy)
        @test CorePeriphery._issymmetric(noisy)
        @test borgatti_everett_continuous(noisy).coreness ≈
              borgatti_everett_continuous(A).coreness atol=1e-12
        @test borgatti_everett_discrete(noisy).core_nodes ==
              borgatti_everett_discrete(A).core_nodes
        @test CorePeriphery._issymmetric(sparse(noisy))
        @test borgatti_everett_discrete(sparse(noisy)).core_nodes ==
              borgatti_everett_discrete(A).core_nodes
        @test minres_symmetric(noisy).coreness ≈ minres_symmetric(A).coreness atol=1e-10

        # The tolerance scales with the largest entry.
        scaled = 1e6 .* A
        scaled[1, 2] += 1e-7
        @test CorePeriphery._issymmetric(scaled)
        scaled[1, 2] += 1e-5
        @test !CorePeriphery._issymmetric(scaled)

        # Genuinely directed input still fails the undirected contract.
        directed = copy(A)
        directed[1, 2] = 0.0
        @test !CorePeriphery._issymmetric(directed)
        @test_throws ArgumentError borgatti_everett_continuous(directed)
        @test_throws ArgumentError lip_discrete(sparse(directed))
        @test_throws ArgumentError surprise_cp(directed)
        @test_throws ArgumentError spectral_method(directed)
        @test_throws ArgumentError multiple_cp_pairs(directed)
        one_edge = zeros(3, 3)
        one_edge[1, 2] = 1e-11
        @test !CorePeriphery._issymmetric(one_edge)
        @test !CorePeriphery._issymmetric(sparse(one_edge))
        @test CorePeriphery._issymmetric(zeros(2, 3)) == false
    end
end
