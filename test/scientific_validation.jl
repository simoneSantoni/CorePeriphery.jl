using CorePeriphery
using LinearAlgebra
using Random
using Test

function reference_rossa_directed(A)
    n = size(A, 1)
    out_strength = vec(sum(A, dims=2))
    transition = A ./ out_strength
    system = Matrix{Float64}(I, n, n) - transpose(transition)
    system[end, :] .= 1.0
    right_hand_side = zeros(n)
    right_hand_side[end] = 1.0
    stationary = system \ right_hand_side
    total_strength = out_strength + vec(sum(A, dims=1))

    selected = falses(n)
    coreness = zeros(n)
    profile = zeros(n)
    selected[argmin(total_strength)] = true
    for step in 2:n
        best_node = 0
        best_alpha = Inf
        best_strength = Inf
        for candidate in 1:n
            selected[candidate] && continue
            candidate_set = copy(selected)
            candidate_set[candidate] = true
            nodes = findall(candidate_set)
            numerator = sum(
                stationary[source] * transition[source, target]
                for source in nodes for target in nodes
            )
            alpha = numerator / sum(stationary[nodes])
            if alpha < best_alpha ||
               (alpha == best_alpha &&
                (total_strength[candidate] < best_strength ||
                 (total_strength[candidate] == best_strength && candidate < best_node)))
                best_node = candidate
                best_alpha = alpha
                best_strength = total_strength[candidate]
            end
        end
        selected[best_node] = true
        profile[step] = best_alpha
        coreness[best_node] = best_alpha
    end
    centralization = 1 - 2 * sum(profile[1:(n - 1)]) / (n - 2)
    return coreness, centralization
end

@testset "Independent scientific validation" begin
    @testset "Exhaustive Lip loss on all five-node graphs" begin
        n = 5
        dyads = [(i, j) for i in 1:(n - 1) for j in (i + 1):n]
        for graph_mask in 0:(2^length(dyads) - 1)
            A = zeros(Int, n, n)
            for (bit, (source, target)) in enumerate(dyads)
                graph_mask & (1 << (bit - 1)) == 0 && continue
                A[source, target] = A[target, source] = 1
            end
            result = lip_discrete(A)
            function loss(core)
                return sum(
                    (source in core && target in core) ? 1 - A[source, target] :
                    (source ∉ core && target ∉ core) ? A[source, target] : 0
                    for (source, target) in dyads
                )
            end
            attained = loss(Set(result.core_nodes))
            optimum = minimum(
                loss(Set(node for node in 1:n if mask & (1 << (node - 1)) != 0))
                for mask in 1:(2^n - 2)
            )
            @test attained == optimum
        end
    end

    @testset "Directed Rossa equation-five oracle" begin
        rng = MersenneTwister(8128)
        for n in 4:7, _ in 1:4
            A = rand(rng, n, n)
            A[diagind(A)] .= 0.0
            expected_coreness, expected_centralization = reference_rossa_directed(A)
            result = random_walker_profiling(A; stationary_method=:linear)
            @test result.coreness ≈ expected_coreness atol=2e-12
            @test result.quality ≈ expected_centralization atol=2e-12
        end
    end
end
