using CorePeriphery
using LinearAlgebra
using Random
using SparseArrays
using Test

@testset "Lip exact degree-prefix algorithm" begin
    star = adjacency_to_matrix([(1, i) for i in 2:6], 6)
    result = lip_discrete(star)

    # k=1 and k=2 tie on Lip's loss for a star; the documented smallest
    # minimizer consists solely of the center.
    @test result.core_nodes == [1]
    @test result.periphery_nodes == collect(2:6)

    # Exhaustively check the attained Lip loss on a nontrivial small graph.
    A = adjacency_to_matrix([(1, 2), (1, 3), (1, 4), (2, 3), (3, 5)], 5)
    lip_loss(core) = sum(
        (i in core && j in core) ? 1.0 - A[i, j] :
        (i ∉ core && j ∉ core) ? A[i, j] : 0.0
        for i in 1:4 for j in (i + 1):5
    )
    result = lip_discrete(A)
    attained = lip_loss(Set(result.core_nodes))
    optimum = minimum(
        lip_loss(Set(i for i in 1:5 if mask & (1 << (i - 1)) != 0))
        for mask in 1:(2^5 - 2)
    )
    @test attained == optimum

    weighted = copy(A)
    weighted[1, 2] = weighted[2, 1] = 0.5
    @test_throws ArgumentError lip_discrete(weighted)
    @test_throws DimensionMismatch lip_discrete(A[1:4, :])
end

@testset "Della Rossa persistence profile" begin
    star = adjacency_to_matrix([(1, i) for i in 2:6], 6)
    result = random_walker_profiling(star)
    @test result.coreness[1] == 1.0
    @test result.coreness[2:6] == zeros(5)
    @test result.quality == 1.0

    n = 5
    complete = ones(Float64, n, n) - Matrix{Float64}(I, n, n)
    result = random_walker_profiling(complete)
    @test result.coreness ≈ collect(0:(n - 1)) ./ (n - 1)
    @test result.quality ≈ 0.0 atol=eps(Float64)

    # The paper-exact mode samples both the initial minimum-strength node and
    # final ties. Explicit RNGs make the stochastic profile reproducible.
    paper_1 = random_walker_profiling(
        complete; tie_break=:paper, rng=MersenneTwister(17),
    )
    paper_2 = random_walker_profiling(
        complete; tie_break=:paper, rng=MersenneTwister(17),
    )
    paper_3 = random_walker_profiling(
        complete; tie_break=:paper, rng=MersenneTwister(18),
    )
    @test paper_1.coreness == paper_2.coreness
    @test paper_1.coreness != paper_3.coreness
    @test sort(paper_1.coreness) ≈ collect(0:(n - 1)) ./ (n - 1)
    @test paper_1.quality ≈ result.quality atol=eps(Float64)

    ensemble_serial = rossa_profile_ensemble(
        complete; n_runs=20, rng=MersenneTwister(29), threaded=false,
    )
    ensemble_threaded = rossa_profile_ensemble(
        complete; n_runs=20, rng=MersenneTwister(29), threaded=true,
    )
    @test ensemble_serial == ensemble_threaded
    @test ensemble_serial isa RossaEnsembleResult
    @test size(ensemble_serial.replicate_coreness) == (n, 20)
    @test ensemble_serial.centralizations ≈ zeros(20) atol=eps(Float64)
    @test ensemble_serial.unique_profiles > 1
    @test all(ensemble_serial.lower_coreness .<= ensemble_serial.mean_coreness .<=
              ensemble_serial.upper_coreness)
    compact_ensemble = rossa_profile_ensemble(
        complete;
        n_runs=3,
        rng=MersenneTwister(2),
        store_replicates=false,
    )
    @test compact_ensemble.replicate_coreness === nothing
    @test length(compact_ensemble.centralizations) == 3

    star_ensemble = rossa_profile_ensemble(
        star; n_runs=20, rng=MersenneTwister(31), threaded=true,
    )
    @test star_ensemble.mean_coreness == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test star_ensemble.centralizations == ones(20)
    @test star_ensemble.rank_stability == 1.0

    # A cycle supplies a symmetric tied fixture whose paper-mode node order is
    # unstable. Its summaries remain seed-reproducible and bounded.
    cycle = adjacency_to_matrix([(node, mod1(node + 1, 6)) for node in 1:6], 6)
    cycle_first = rossa_profile_ensemble(
        cycle; n_runs=30, rng=MersenneTwister(32), threaded=false,
    )
    cycle_second = rossa_profile_ensemble(
        cycle; n_runs=30, rng=MersenneTwister(32), threaded=true,
    )
    @test cycle_first == cycle_second
    @test cycle_first.unique_profiles > 1
    @test 0.0 <= cycle_first.rank_stability <= 1.0

    # A connected noisy planted graph checks that interval and centralization
    # summaries remain finite outside analytic symmetric cases.
    noisy = adjacency_to_matrix(
        vcat(
            [(1, node) for node in 2:10],
            [(2, node) for node in 3:7],
            [(6, 7), (8, 9), (9, 10)],
        ),
        10,
    )
    noisy_ensemble = rossa_profile_ensemble(
        noisy; n_runs=24, rng=MersenneTwister(33), store_replicates=false,
    )
    @test all(isfinite, noisy_ensemble.mean_coreness)
    @test all(isfinite, noisy_ensemble.centralizations)
    @test all(noisy_ensemble.lower_coreness .<= noisy_ensemble.upper_coreness)
    @test all(0.0 .<= noisy_ensemble.lower_coreness .<= 1.0)
    @test all(0.0 .<= noisy_ensemble.upper_coreness .<= 1.0)
    @test_throws ArgumentError rossa_profile_ensemble(complete; n_runs=0)
    @test_throws ArgumentError rossa_profile_ensemble(complete; interval=(0.9, 0.1))
    @test_throws ArgumentError rossa_profile_ensemble(
        complete; n_runs=10, memory_limit_bytes=8,
    )
    @test_throws ArgumentError rossa_profile_ensemble(
        complete; memory_limit_bytes=-1,
    )

    # Compatibility simulation keywords no longer affect this deterministic method.
    @test random_walker_profiling(star; n_walks=1, walk_length=1).coreness ==
          random_walker_profiling(star; n_walks=10_000, walk_length=100).coreness

    disconnected = adjacency_to_matrix([(1, 2), (3, 4)], 4)
    @test_throws ArgumentError random_walker_profiling(disconnected)

    asymmetric = copy(star)
    asymmetric[1, 2] = 0.0
    @test_throws ArgumentError random_walker_profiling(asymmetric)

    # Weighted, directed, and strongly connected. A direct stationary solve
    # gives pi = [39, 40, 40, 21] / 140. Equation (5), with deterministic
    # tie-breaking, inserts nodes [4, 3, 1, 2].
    directed = [
        0.0 2.0 0.0 1.0
        0.0 0.0 3.0 0.0
        4.0 0.0 0.0 1.0
        1.0 2.0 0.0 0.0
    ]
    transition = directed ./ sum(directed, dims=2)
    stationary_system = Matrix{Float64}(I, 4, 4) - transpose(transition)
    stationary_system[end, :] .= 1.0
    stationary = stationary_system \ [0.0, 0.0, 0.0, 1.0]
    @test stationary ≈ [39.0, 40.0, 40.0, 21.0] ./ 140
    @test transpose(transition) * stationary ≈ stationary

    lazy_stationary = rossa_stationary_distribution(directed; method=:lazy)
    linear_stationary = rossa_stationary_distribution(directed; method=:linear)
    automatic_stationary = rossa_stationary_distribution(directed; method=:auto)
    sparse_stationary = rossa_stationary_distribution(sparse(directed); method=:auto)
    @test lazy_stationary isa StationaryDistributionResult
    @test lazy_stationary.distribution ≈ stationary atol=2e-12
    @test linear_stationary.distribution ≈ stationary atol=2e-12
    @test automatic_stationary.method == :linear
    @test sparse_stationary.method == :lazy
    @test lazy_stationary.residual <= 1e-12
    @test linear_stationary.residual <= 1e-12
    @test lazy_stationary.iterations > 0

    directed_result = random_walker_profiling(directed)
    @test directed_result.coreness ≈ [3 / 5, 1, 8 / 61, 0] atol=2e-12
    @test directed_result.quality ≈ 82 / 305 atol=2e-12
    @test random_walker_profiling(sparse(directed)).coreness ≈
          directed_result.coreness
    @test random_walker_profiling(Float32.(directed)).coreness ≈
          directed_result.coreness
    @test random_walker_profiling(
        directed; stationary_method=:linear,
    ).coreness ≈ directed_result.coreness atol=2e-12
    @test random_walker_profiling(
        directed; stationary_distribution=stationary,
    ).coreness ≈ directed_result.coreness atol=2e-12

    directed_cycle = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
    @test random_walker_profiling(directed_cycle).coreness ≈ [0.0, 0.5, 1.0]

    not_strongly_connected = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    @test_throws ArgumentError random_walker_profiling(not_strongly_connected)
    @test_throws ArgumentError random_walker_profiling(star; tie_break=:unknown)
    @test_throws ArgumentError random_walker_profiling(star; stationary_tol=-1)
    @test_throws ArgumentError random_walker_profiling(star; stationary_max_iter=0)
    @test_throws ErrorException random_walker_profiling(
        directed; stationary_tol=0, stationary_max_iter=1,
    )
    @test_throws ArgumentError random_walker_profiling(
        directed; stationary_method=:invalid,
    )
    @test_throws DimensionMismatch random_walker_profiling(
        directed; stationary_distribution=[0.5, 0.5],
    )
    @test_throws ArgumentError random_walker_profiling(
        directed; stationary_distribution=fill(0.25, 4),
    )
    @test_throws ArgumentError random_walker_profiling(
        star; stationary_distribution=fill(1 / 6, 6),
    )

    # Two dense pairs linked by weak directed bridges form a slowly mixing but
    # strongly connected chain. The direct solver remains accurate when a tiny
    # lazy-iteration budget is deliberately insufficient.
    epsilon = 1e-6
    slow = [
        0.0 1.0 0.0 0.0
        1.0 0.0 epsilon 0.0
        0.0 0.0 0.0 2.0
        epsilon 0.0 2.0 0.0
    ]
    @test_throws ErrorException rossa_stationary_distribution(
        slow; method=:lazy, max_iter=2,
    )
    slow_linear = rossa_stationary_distribution(slow; method=:linear)
    @test slow_linear.residual <= 1e-12

    signed = copy(star)
    signed[1, 2] = signed[2, 1] = -1.0
    @test_throws ArgumentError random_walker_profiling(signed)
end
