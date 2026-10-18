using CorePeriphery
using Random
using Statistics
using Test

@testset "Streaming core quality" begin
    A = [0 1 0; 1 0 1; 0 1 0]
    c = [1, 0.5, 0]
    expected = cor([1.0, 0.0, 1.0], [0.5, 0.0, 0.0])

    @test core_quality(A, c) ≈ expected
    @test core_quality(@view(A[:, :]), @view(c[:])) ≈ expected
    @test core_quality(zeros(Int, 3, 3), ones(Int, 3)) == 0.0
    @test core_quality(zeros(2, 2), [1.0, 0.0]; discrete=true) == 0.0
    @test core_quality(zeros(1, 1), [1.0]) == 0.0
    @test_throws DimensionMismatch core_quality(zeros(2, 3), [1.0, 0.0])
    @test_throws DimensionMismatch core_quality(zeros(3, 3), [1.0, 0.0])
end

@testset "Rombach template permutation optimizer" begin
    A = zeros(5, 5)
    for leaf in 2:5
        A[1, leaf] = A[leaf, 1] = 1.0
    end

    first_result = rombach_continuous(
        A; alpha=0.8, beta=0.5, n_runs=4, max_iter=100,
        rng=MersenneTwister(2025),
    )
    second_result = rombach_continuous(
        A; alpha=0.8, beta=0.5, n_runs=4, max_iter=100,
        rng=MersenneTwister(2025),
    )

    @test first_result.coreness == second_result.coreness
    @test first_result.quality == second_result.quality
    @test first_result.coreness[1] == maximum(first_result.coreness)

    objective = sum(
        A[i, j] * first_result.coreness[i] * first_result.coreness[j]
        for i in 1:4 for j in (i + 1):5
    )
    @test first_result.quality ≈ objective

    sharp = rombach_continuous(
        A; alpha=1.0, beta=0.5, n_runs=2, max_iter=50,
        rng=MersenneTwister(7),
    )
    smooth = rombach_continuous(
        A; alpha=0.0, beta=0.5, n_runs=2, max_iter=50,
        rng=MersenneTwister(7),
    )
    @test sort(sharp.coreness) != sort(smooth.coreness)

    @test_throws ArgumentError rombach_continuous(A; alpha=-0.1)
    @test_throws ArgumentError rombach_continuous(A; alpha=1.1)
    @test_throws ArgumentError rombach_continuous(A; beta=0.0)
    @test_throws ArgumentError rombach_continuous(A; beta=1.0)
    @test_throws DimensionMismatch rombach_continuous(zeros(2, 3))
end
