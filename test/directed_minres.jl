using CorePeriphery
using LinearAlgebra
using SparseArrays
using Test

# The production module will include this file from CorePeriphery.jl. Loading it
# explicitly keeps this focused test independent until that compatibility wiring
# is added.
if !isdefined(CorePeriphery, :CPDirectedResult)
    Base.include(CorePeriphery, joinpath(@__DIR__, "..", "src", "directed.jl"))
end

@testset "Symmetric rank-one MINRES" begin
    A = adjacency_to_matrix(
        [(1, 2, 3.0), (1, 3, 2.0), (1, 4, 1.0),
         (2, 3, 2.0), (2, 4, 1.0), (1, 5, 0.5)],
        5,
    )
    result = minres_symmetric(A)
    @test result isa CPResult
    @test result.algorithm == "MINRES Symmetric Rank One"
    @test all(0.0 .<= result.coreness .<= 1.0)
    @test 0.0 <= result.quality <= 1.0

    factor, work = CorePeriphery._minres_symmetric_factor(A, 5, 10_000, 1e-8, 0.5)
    analytic = CorePeriphery._directed_residual!(work, A, factor, factor)
    explicit = sum(abs2(A[i, j] - factor[i] * factor[j])
                   for i in axes(A, 1), j in axes(A, 2) if i != j)
    @test analytic ≈ explicit atol=1e-12
    @test result.quality ≈ 1 - explicit / sum(abs2, A) atol=1e-12

    permutation = [3, 5, 1, 4, 2]
    permuted = minres_symmetric(A[permutation, permutation])
    @test permuted.coreness ≈ result.coreness[permutation] atol=1e-10
    @test permuted.quality ≈ result.quality atol=1e-12

    sparse_result = minres_symmetric(sparse(A))
    @test sparse_result.coreness ≈ result.coreness atol=1e-10
    @test sparse_result.quality ≈ result.quality atol=1e-12
    @test minres_symmetric(transpose(A)).coreness ≈ result.coreness atol=1e-10

    @test_throws ArgumentError minres_symmetric([0.0 1.0; 0.0 0.0])
    @test_throws ArgumentError minres_symmetric(A; max_iter=-1)
    @test_throws ArgumentError minres_symmetric(A; tol=-1)
    @test_throws ArgumentError minres_symmetric(A; damping=0)
end

@testset "Directed MINRES result" begin
    A = [
        0.0 4.0 3.0 2.0 1.0;
        1.0 0.0 2.0 1.0 0.0;
        0.0 3.0 0.0 2.0 0.0;
        0.0 1.0 0.0 0.0 0.0;
        2.0 0.0 0.0 0.0 0.0
    ]

    result = CorePeriphery.minres_svd_directed(A)
    @test result isa CorePeriphery.CPDirectedResult
    @test length(result.out_coreness) == size(A, 1)
    @test length(result.in_coreness) == size(A, 1)
    @test result.coreness ≈ (result.out_coreness + result.in_coreness) / 2
    @test all(0.0 .<= result.out_coreness .<= 1.0)
    @test all(0.0 .<= result.in_coreness .<= 1.0)
    @test result.residual >= 0.0
    @test 0.0 <= result.quality <= 1.0
    @test sort(vcat(result.core_nodes, result.periphery_nodes)) == 1:size(A, 1)

    transposed = CorePeriphery.minres_svd_directed(transpose(A))
    @test transposed.out_coreness ≈ result.in_coreness atol=1e-10
    @test transposed.in_coreness ≈ result.out_coreness atol=1e-10
    @test transposed.coreness ≈ result.coreness atol=1e-10
    @test transposed.residual ≈ result.residual atol=1e-10
    @test transposed.quality ≈ result.quality atol=1e-10

    permutation = [3, 5, 1, 4, 2]
    permuted = CorePeriphery.minres_svd_directed(A[permutation, permutation])
    @test permuted.out_coreness ≈ result.out_coreness[permutation] atol=1e-10
    @test permuted.in_coreness ≈ result.in_coreness[permutation] atol=1e-10
    @test permuted.coreness ≈ result.coreness[permutation] atol=1e-10
    @test permuted.residual ≈ result.residual atol=1e-10

    integer_result = CorePeriphery.minres_svd_directed(round.(Int, A))
    @test integer_result.coreness ≈ result.coreness

    sparse_result = CorePeriphery.minres_svd_directed(sparse(A))
    @test sparse_result.out_coreness ≈ result.out_coreness atol=1e-10
    @test sparse_result.in_coreness ≈ result.in_coreness atol=1e-10
    @test sparse_result.residual ≈ result.residual atol=1e-10

    u = [0.2, 0.8, 0.5, 0.1, 0.3]
    v = [0.7, 0.4, 0.9, 0.2, 0.6]
    analytic = CorePeriphery._directed_residual!(similar(u), sparse(A), u, v)
    explicit = sum(abs2(A[i, j] - u[i] * v[j])
                   for i in axes(A, 1), j in axes(A, 2) if i != j)
    @test analytic ≈ explicit atol=1e-12

    @test_throws DimensionMismatch CorePeriphery.minres_svd_directed(zeros(2, 3))
    @test_throws ArgumentError CorePeriphery.minres_svd_directed([0.0 -1.0; 1.0 0.0])
    @test_throws ArgumentError CorePeriphery.minres_svd_directed([1.0 0.0; 0.0 0.0])
    @test_throws ArgumentError CorePeriphery.minres_svd_directed(A; max_iter=-1)
    @test_throws ArgumentError CorePeriphery.minres_svd_directed(A; tol=-1.0)
end
