using LinearAlgebra
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
end
