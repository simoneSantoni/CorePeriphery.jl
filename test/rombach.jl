using CorePeriphery
using Random
using SparseArrays
using Test

const _CP = CorePeriphery

# Direct transcription of Rombach et al. (2017) eq. (4.9): ranks 1:⌊βN⌋ are
# peripheral, ranks ⌊βN⌋+1:N are core.
function rombach_template_oracle(N, α, β)
    b = floor(Int, β * N + 1e-9)
    return [i <= b ? i * (1 - α) / (2b) : (i - b) * (1 - α) / (2 * (N - b)) + (1 + α) / 2
            for i in 1:N]
end

function planted_rombach_fixture()
    # Core 1:5 is a clique; each peripheral node 6:20 attaches to three core
    # nodes; no periphery-periphery edges. Deterministic and connected.
    edges = Tuple{Int,Int}[]
    for i in 1:4, j in (i + 1):5
        push!(edges, (i, j))
    end
    for node in 6:20
        for offset in 0:2
            push!(edges, (mod1(node + offset, 5), node))
        end
    end
    return adjacency_to_matrix(edges, 20)
end

@testset "Rombach transition function, eq. (4.9)" begin
    # The reviewer's reference values for N = 5, α = β = 0.5.
    @test _CP._rombach_template(5, 0.5, 0.5) ≈ [0.125, 0.25, 5 / 6, 11 / 12, 1.0]
    # Hand-computed: N = 10, α = 0, β = 0.29 gives ⌊2.9⌋ = 2 peripheral ranks.
    @test _CP._rombach_template(10, 0.0, 0.29) ≈
          vcat([1 / 4, 2 / 4], [(i - 2) / 16 + 1 / 2 for i in 3:10])
    # α = 1 is the discrete classification: zeros then ones.
    @test _CP._rombach_template(6, 1.0, 0.5) == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    # β = 0: only the core case applies; β = 1: only the periphery case.
    @test _CP._rombach_template(4, 0.5, 0.0) ≈ [i * 0.5 / 8 + 0.75 for i in 1:4]
    @test _CP._rombach_template(4, 0.5, 1.0) ≈ [i * 0.5 / 8 for i in 1:4]
    @test _CP._rombach_boundary(100, 0.29) == 29
    @test _CP._rombach_boundary(5, 0.5) == 2
    @test _CP._rombach_boundary(7, 1.0) == 7

    for (N, α, β) in (
        (5, 0.5, 0.5), (7, 0.3, 0.4), (10, 0.0, 0.29), (10, 1.0, 0.5),
        (6, 0.8, 0.0), (6, 0.8, 1.0), (1, 0.5, 0.5), (13, 0.25, 0.37),
        (20, 0.9, 0.05), (3, 0.0, 0.0), (3, 0.0, 1.0),
    )
        template = _CP._rombach_template(N, α, β)
        @test template ≈ rombach_template_oracle(N, α, β)
        @test issorted(template)
        # The top rank scores 1 whenever a core case exists (β < 1); with
        # β = 1 every rank is peripheral and the top one scores (1 - α) / 2.
        @test template[end] ≈ (β < 1 ? 1.0 : (1 - α) / 2)
        @test all(x -> 0 <= x <= 1, template)
    end
    @test _CP._rombach_template(4, 1.0, 1.0) == zeros(4)
    @test sum(_CP._rombach_normalized_template(9, 0.3, 0.6)) ≈ 1.0
    @test _CP._rombach_normalized_template(4, 1.0, 1.0) == zeros(4)
end

@testset "rombach_continuous: normalisation, objective, beta" begin
    star = adjacency_to_matrix([(1, i) for i in 2:5], 5)
    result = rombach_continuous(star; alpha=0.5, beta=0.5, n_runs=4, rng=MersenneTwister(1))
    raw = _CP._rombach_template(5, 0.5, 0.5)
    @test sum(result.coreness) ≈ 1.0
    @test sort(result.coreness) ≈ raw ./ sum(raw)
    @test result.coreness[1] == maximum(result.coreness)
    # Quality is eq. (4.10): the full double sum over ordered pairs.
    @test result.quality ≈ sum(
        star[i, j] * result.coreness[i] * result.coreness[j] for i in 1:5, j in 1:5)
    # Default partition: N - ⌊βN⌋ = 3 core ranks.
    @test length(result.core_nodes) == 3
    @test 1 in result.core_nodes
    @test rombach_continuous(star; partition=:beta, n_runs=4, rng=MersenneTwister(1)) == result

    @test rombach_continuous(star; beta=0.0, n_runs=1).core_nodes == collect(1:5)
    @test isempty(rombach_continuous(star; beta=1.0, n_runs=1).core_nodes)
    @test length(rombach_continuous(star; beta=0.6, n_runs=1).core_nodes) == 2
    @test length(rombach_continuous(star; beta=0.84, n_runs=2,
                                    rng=MersenneTwister(1)).core_nodes) == 1
    @test_throws ArgumentError rombach_continuous(star; beta=-0.1)
    @test_throws ArgumentError rombach_continuous(star; beta=1.1)
    @test_throws ArgumentError rombach_continuous(star; partition=:beta, core_size=2)
    @test_throws ArgumentError rombach_continuous(star; partition=:beta, threshold=0.1)
    @test_throws ArgumentError rombach_continuous(star; partition=:unknown)
    @test isempty(rombach_continuous(star; partition=:none).core_nodes)
    @test rombach_continuous(star; partition=:core_size, core_size=1,
                             n_runs=2, rng=MersenneTwister(1)).core_nodes == [1]

    # α = β = 1: the template is identically zero.
    degenerate = rombach_continuous(star; alpha=1.0, beta=1.0, n_runs=1)
    @test degenerate.coreness == zeros(5)
    @test degenerate.quality == 0.0
    @test isempty(degenerate.core_nodes)

    @test rombach_continuous(sparse(star); n_runs=3, rng=MersenneTwister(4)) ==
          rombach_continuous(star; n_runs=3, rng=MersenneTwister(4))
end

@testset "rombach_core_score, eq. (4.11)" begin
    A = planted_rombach_fixture()
    alphas = 0.25:0.25:1.0
    betas = 0.25:0.25:1.0
    score = rombach_core_score(
        A; alphas=alphas, betas=betas, n_runs=2, rng=MersenneTwister(3))
    @test score.algorithm == "Rombach aggregate core score"
    @test maximum(score.coreness) == 1.0
    @test all(x -> 0 <= x <= 1, score.coreness)
    @test minimum(score.coreness[1:5]) > maximum(score.coreness[6:20])
    @test rombach_core_score(A; alphas=alphas, betas=betas, n_runs=2,
                             rng=MersenneTwister(3), partition=:core_size,
                             core_size=5).core_nodes == collect(1:5)

    # Oracle: CS(i) = Z Σ_γ C_i(γ) R_γ assembled from single-(α, β) fits that
    # consume the same random stream in the same grid order.
    rng = MersenneTwister(3)
    aggregate = zeros(20)
    total = 0.0
    for α in alphas, β in betas
        fit = rombach_continuous(A; alpha=α, beta=β, n_runs=2, rng=rng)
        aggregate .+= fit.coreness .* fit.quality
        total += fit.quality
    end
    @test score.coreness ≈ aggregate ./ maximum(aggregate)
    @test score.quality ≈ total / (length(alphas) * length(betas))

    # Isolates score zero and never join the core.
    with_isolate = zeros(21, 21)
    with_isolate[1:20, 1:20] .= A
    isolated = rombach_core_score(
        with_isolate; alphas=0.5:0.5:1.0, betas=0.5:0.5:1.0, n_runs=1,
        rng=MersenneTwister(8))
    @test isolated.coreness[21] == 0.0
    @test 21 ∉ isolated.core_nodes
    @test isolated.coreness[1:20] ≈ rombach_core_score(
        A; alphas=0.5:0.5:1.0, betas=0.5:0.5:1.0, n_runs=1,
        rng=MersenneTwister(8)).coreness

    @test_throws ArgumentError rombach_core_score(A; alphas=Float64[])
    @test_throws ArgumentError rombach_core_score(A; betas=[0.5, 1.5])
    @test_throws ArgumentError rombach_core_score(A; n_runs=0)
    @test_throws ArgumentError rombach_core_score(A; partition=:core_size)
end
