using SparseArrays

function reference_quality(A, c; discrete=false, directed=false)
    adjacency = Float64[]
    pattern = Float64[]
    n = length(c)
    if directed
        for j in 1:n, i in 1:n
            i == j && continue
            push!(adjacency, A[i, j])
            push!(pattern, discrete ? max(c[i], c[j]) : c[i] * c[j])
        end
    else
        for j in 2:n, i in 1:(j - 1)
            push!(adjacency, A[i, j])
            push!(pattern, discrete ? max(c[i], c[j]) : c[i] * c[j])
        end
    end
    (var(adjacency) > 0 && var(pattern) > 0) || return 0.0
    return cor(adjacency, pattern)
end

@testset "Sufficient-statistic and sparse kernels" begin
    rng = MersenneTwister(9182)
    for directed in (false, true), discrete in (false, true)
        A = rand(rng, 9, 9)
        A[diagind(A)] .= 0.0
        directed || (A = Matrix(Symmetric(A, :U)))
        c = rand(rng, 9)
        expected = reference_quality(A, c; discrete=discrete, directed=directed)
        @test core_quality(A, c; discrete=discrete, directed=directed) ≈ expected
        @test core_quality(sparse(A), c; discrete=discrete, directed=directed) ≈ expected
    end

    A = adjacency_to_matrix([(1, 2), (1, 3), (2, 3), (1, 4), (2, 5)], 5)
    @test borgatti_everett_continuous(sparse(A); max_iter=4).quality ≈
          borgatti_everett_continuous(A; max_iter=4).quality
    @test borgatti_everett_discrete(sparse(A)).coreness ==
          borgatti_everett_discrete(A).coreness
    @test random_walker_profiling(sparse(A)).coreness ==
          random_walker_profiling(A).coreness
    @test surprise_cp(sparse(A); max_iter=5).coreness ==
          surprise_cp(A; max_iter=5).coreness

    @test_throws ArgumentError adjacency_to_matrix([(1, 1)], 2)
    @test_throws ArgumentError adjacency_to_matrix([(1, 2, -1.0)], 2)
    @test_throws BoundsError adjacency_to_matrix([(1, 3)], 2)
end
