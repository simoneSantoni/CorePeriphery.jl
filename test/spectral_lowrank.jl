using CorePeriphery
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Test

function ideal_cp_adjacency(n, core)
    A = zeros(Float64, n, n)
    for i in 1:(n - 1), j in (i + 1):n
        if i in core || j in core
            A[i, j] = A[j, i] = 1.0
        end
    end
    return A
end

function cp_density(A, core)
    n = size(A, 1)
    periphery = setdiff(1:n, core)
    cc = sum(A[i, j] for (position, i) in enumerate(core)
             for j in core[(position + 1):end]; init=0.0)
    cp = sum(A[i, j] for i in core for j in periphery; init=0.0)
    pp = sum(A[i, j] for (position, i) in enumerate(periphery)
             for j in periphery[(position + 1):end]; init=0.0)
    vcc = length(core) * (length(core) - 1) / 2
    vcp = length(core) * length(periphery)
    vpp = length(periphery) * (length(periphery) - 1) / 2
    return (iszero(vcc) ? 0.0 : cc / vcc) + cp / vcp -
           (iszero(vpp) ? 0.0 : pp / vpp)
end

@testset "Cucuringu LowRank-Core reconstruction and Find-Cut" begin
    A = ideal_cp_adjacency(8, 1:3)
    result = CorePeriphery._lowrank_core(A; beta=0.25)

    @test result.coreness == [8.0, 8.0, 8.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    @test result.core_nodes == collect(1:3)
    @test result.periphery_nodes == collect(4:8)
    @test result.quality == 2.0

    # Independent reconstruction oracle: top two eigenvalues by magnitude,
    # threshold strictly above 0.5, then take row sums.
    decomposition = eigen(Symmetric(A))
    selected = sortperm(abs.(decomposition.values); rev=true)[1:2]
    reconstructed = sum(
        decomposition.values[index] *
        decomposition.vectors[:, index] * decomposition.vectors[:, index]'
        for index in selected
    )
    oracle_scores = Float64.(vec(sum(reconstructed .> 0.5, dims=2)))
    @test result.coreness == oracle_scores

    order = sortperm(oracle_scores; rev=true, alg=Base.Sort.MergeSort)
    feasible = 2:6
    oracle_qualities = [cp_density(A, order[1:k]) for k in feasible]
    @test result.quality == maximum(oracle_qualities)
    @test length(result.core_nodes) == feasible[argmax(oracle_qualities)]
end

@testset "LowRank-Core invariance and matrix interfaces" begin
    A = ideal_cp_adjacency(8, 1:3)
    permutation = [5, 2, 8, 1, 7, 3, 6, 4]
    permuted = CorePeriphery._lowrank_core(A[permutation, permutation]; beta=0.25)
    original = CorePeriphery._lowrank_core(A; beta=0.25)

    @test permuted.coreness == original.coreness[permutation]
    @test Set(permutation[permuted.core_nodes]) == Set(original.core_nodes)
    @test permuted.quality == original.quality

    sparse_result = CorePeriphery._lowrank_core(sparse(A); beta=0.25)
    @test sparse_result == original

    dense_statistics = CorePeriphery._ordered_cut_statistics(A, permutation)
    sparse_statistics = CorePeriphery._ordered_cut_statistics(sparse(A), permutation)
    @test sparse_statistics == dense_statistics
end

@testset "LowRank-Core input contract" begin
    A = ideal_cp_adjacency(6, 1:2)
    @test_throws DimensionMismatch CorePeriphery._lowrank_core(A[1:5, :])
    @test_throws ArgumentError CorePeriphery._lowrank_core(A; beta=0.6)

    asymmetric = copy(A)
    asymmetric[1, 2] = 0.0
    @test_throws ArgumentError CorePeriphery._lowrank_core(asymmetric)

    weighted = copy(A)
    weighted[1, 2] = weighted[2, 1] = 0.5
    @test_throws ArgumentError CorePeriphery._lowrank_core(weighted)

    looped = copy(A)
    looped[1, 1] = 1.0
    @test_throws ArgumentError CorePeriphery._lowrank_core(looped)
    @test_throws ArgumentError CorePeriphery._lowrank_core(sparse(weighted))
    @test_throws ArgumentError CorePeriphery._lowrank_core(sparse(looped))
end

ring_graph(n) = sparse([1:n; 1:n], [mod1.(2:(n + 1), n); mod1.(0:(n - 1), n)], 1.0, n, n)
star_graph(n) = sparse([ones(Int, n - 1); 2:n], [2:n; ones(Int, n - 1)], 1.0, n, n)

@testset "LowRank-Core ARPACK convergence on clustered spectra" begin
    # Rings have eigenvalues 2cos(2πk/n) clustered near ±2; stars have a
    # degenerate pair ±sqrt(n - 1). Both must converge with the finite ARPACK
    # tolerance and agree with the dense LAPACK path.
    for (name, build, n) in (("ring", ring_graph, 600), ("star", star_graph, 600),
                             ("ring", ring_graph, 1001), ("star", star_graph, 1500))
        A = build(n)
        values, vectors = CorePeriphery._lowrank_eigenpairs(A)
        dense_values, _ = CorePeriphery._lowrank_eigenpairs(Matrix(A))
        @test abs.(values) ≈ abs.(dense_values) atol = 1e-8
        @test abs(values[1]) >= abs(values[2])
        sparse_result = spectral_method(A)
        dense_result = spectral_method(Matrix(A))
        @test sparse_result == dense_result
        if name == "star"
            # Hub scores 2(n - 1) - ish rank-two mass: it is the unique core.
            @test argmax(sparse_result.coreness) == 1
            @test sparse_result.core_nodes[1] == 1
        end
    end

    # Convergence failures surface as ConvergenceError carrying the ARPACK
    # settings, never as a bare ErrorException.
    @test CorePeriphery._LOWRANK_TOL > 0
    # Arpack.jl logs the Fortran XYAUPD failure at error level before it
    # throws; that record is expected here and is silenced.
    failure = Logging.with_logger(Logging.NullLogger()) do
        try
            CorePeriphery._lowrank_eigenpairs(ring_graph(400); ncv=3, maxiter=1)
            nothing
        catch caught
            caught
        end
    end
    @test failure isa ConvergenceError
    @test occursin("maxiter=1", sprint(showerror, failure))
    @test occursin("ncv=3", sprint(showerror, failure))
end

@testset "LowRank-Core pruned threshold scores" begin
    # The active-row pruning must reproduce the full O(n²) pair loop exactly.
    function reference_scores(values, vectors)
        n = size(vectors, 1)
        scores = zeros(Float64, n)
        for i in 1:n
            diagonal = values[1] * vectors[i, 1]^2 + values[2] * vectors[i, 2]^2
            diagonal > 0.5 && (scores[i] += 1.0)
            for j in (i + 1):n
                reconstructed = values[1] * vectors[i, 1] * vectors[j, 1] +
                                values[2] * vectors[i, 2] * vectors[j, 2]
                if reconstructed > 0.5
                    scores[i] += 1.0
                    scores[j] += 1.0
                end
            end
        end
        return scores
    end
    rng = MersenneTwister(77)
    for trial in 1:20
        n = rand(rng, 2:40)
        values = [10.0 * randn(rng), 10.0 * randn(rng)]
        vectors = randn(rng, n, 2)
        @test CorePeriphery._thresholded_rank_two_scores(values, vectors) ==
              reference_scores(values, vectors)
    end
    planted = ideal_cp_adjacency(30, 1:6)
    values, vectors = CorePeriphery._lowrank_eigenpairs(planted)
    @test CorePeriphery._thresholded_rank_two_scores(values, vectors) ==
          reference_scores(values, vectors)
end
