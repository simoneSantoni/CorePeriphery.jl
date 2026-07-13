using CorePeriphery
using Random
using SparseArrays
using Test

function planted_pair_graph()
    A = zeros(8, 8)
    for offset in (0, 4)
        core = offset + 1
        for node in (offset + 2):(offset + 4)
            A[core, node] = A[node, core] = 1.0
        end
        A[offset + 1, offset + 2] = A[offset + 2, offset + 1] = 2.0
    end
    return A
end

@testset "Cached joint objective and move deltas" begin
    A = planted_pair_graph()
    graph = CorePeriphery._joint_build_graph(A)
    labels = [1, 1, 1, 2, 2, 2, 3, 3]
    core = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    aggregates = CorePeriphery._joint_pair_aggregates(
        labels, core, graph.degrees, 3,
    )
    normalization = graph.total_weight
    density = graph.total_weight / (size(A, 1) * (size(A, 1) - 1) / 2)

    for null_model in (:er, :configuration)
        context = if null_model === :er
            CorePeriphery._JointERContext(normalization, density)
        else
            CorePeriphery._JointConfigurationContext(
                normalization, 1.0 / (2.0 * graph.total_weight), graph.degrees,
            )
        end
        residual(i, j) = begin
            expected = null_model === :er ? density :
                graph.degrees[i] * graph.degrees[j] / (2.0 * graph.total_weight)
            (A[min(i, j), max(i, j)] - expected) / normalization
        end

        cached_quality = CorePeriphery._joint_objective(
            graph, labels, core, context, aggregates,
        )
        brute_quality = CorePeriphery._joint_objective(labels, core, residual)
        @test cached_quality ≈ brute_quality atol=1e-14

        for node in eachindex(labels)
            observed_pair = zeros(3)
            observed_core = zeros(3)
            for index in eachindex(graph.neighbors[node])
                neighbor = graph.neighbors[node][index]
                pair = labels[neighbor]
                weight = graph.weights[node][index]
                observed_pair[pair] += weight
                core[neighbor] == 1.0 && (observed_core[pair] += weight)
            end
            old_pair = labels[node]
            old_core = core[node]
            old_contribution = CorePeriphery._joint_pair_contribution(
                context, node, old_pair, old_core, old_pair, old_core,
                aggregates, observed_pair, observed_core,
            )
            for new_pair in 1:3, new_core in (0.0, 1.0)
                cached_delta = CorePeriphery._joint_pair_contribution(
                    context, node, new_pair, new_core, old_pair, old_core,
                    aggregates, observed_pair, observed_core,
                ) - old_contribution
                brute_delta = CorePeriphery._joint_move_delta(
                    node, new_pair, new_core, labels, core, residual,
                )
                @test cached_delta ≈ brute_delta atol=1e-14
            end
        end
    end

    updated = deepcopy(aggregates)
    updated_labels = copy(labels)
    updated_core = copy(core)
    node = 1
    CorePeriphery._joint_update_aggregates!(
        updated, node, labels[node], core[node], 2, 0.0, graph.degrees[node],
    )
    updated_labels[node] = 2
    updated_core[node] = 0.0
    rebuilt = CorePeriphery._joint_pair_aggregates(
        updated_labels, updated_core, graph.degrees, 3,
    )
    @test updated.counts == rebuilt.counts
    @test updated.core_counts == rebuilt.core_counts
    @test updated.degree_sums ≈ rebuilt.degree_sums
    @test updated.core_degree_sums ≈ rebuilt.core_degree_sums
    @test updated.degree_square_sums ≈ rebuilt.degree_square_sums
    @test updated.core_degree_square_sums ≈ rebuilt.core_degree_square_sums
end

@testset "Joint Kojaku-Masuda multi-pair engine" begin
    A = planted_pair_graph()
    result = CorePeriphery._multiple_cp_pairs_joint(
        A; max_pairs=2, min_pair_size=2, max_iter=100, n_runs=12,
        rng=MersenneTwister(41),
    )

    @test result.n_pairs == 2
    @test length(unique(result.pair_labels[1:4])) == 1
    @test length(unique(result.pair_labels[5:8])) == 1
    @test result.pair_labels[1] != result.pair_labels[5]
    @test all(diff(result.objective_history) .>= -1e-10)
    @test all(pair -> count(==(pair), result.pair_labels) >= 2,
              unique(result.pair_labels))
    total_weight = sum(A) / 2
    density = total_weight / (8 * 7 / 2)
    manual_quality = sum(
        (A[i, j] - density) / total_weight *
        max(result.coreness[i], result.coreness[j]) *
        (result.pair_labels[i] == result.pair_labels[j])
        for i in 1:7 for j in (i + 1):8
    )
    @test result.quality ≈ manual_quality

    capped = CorePeriphery._multiple_cp_pairs_joint(
        A; max_pairs=1, min_pair_size=2, n_runs=3,
        rng=MersenneTwister(2),
    )
    @test capped.n_pairs == 1

    repeated = CorePeriphery._multiple_cp_pairs_joint(
        A; max_pairs=2, min_pair_size=2, max_iter=100, n_runs=12,
        rng=MersenneTwister(41),
    )
    @test result.pair_labels == repeated.pair_labels
    @test result.coreness == repeated.coreness
    @test result.quality == repeated.quality

    configuration = CorePeriphery._multiple_cp_pairs_joint(
        A; null_model=:configuration, max_pairs=2, min_pair_size=2,
        n_runs=12, rng=MersenneTwister(41),
    )
    @test configuration.null_model == :configuration
    @test isfinite(configuration.quality)
    @test configuration.quality != result.quality

    sparse_result = CorePeriphery._multiple_cp_pairs_joint(
        sparse(A); max_pairs=2, min_pair_size=2, n_runs=12,
        rng=MersenneTwister(41),
    )
    @test sparse_result.pair_labels == result.pair_labels
    @test sparse_result.coreness == result.coreness
    @test sparse_result.quality ≈ result.quality

    selected = multiple_cp_pairs(
        A;
        max_pairs=2,
        min_pair_size=2,
        n_runs=20,
        pair_selection=:penalized,
        pair_penalty=0.5,
        rng=MersenneTwister(41),
    )
    @test selected.n_pairs == 2
    @test selected.pair_selection == :penalized
    @test selected.candidate_pair_counts == [1, 2]
    @test length(selected.candidate_qualities) == 2
    @test length(selected.selection_scores) == 2

    single_n = 20
    single_edges = [
        (i, j) for i in 1:(single_n - 1) for j in (i + 1):single_n if i <= 5
    ]
    single = adjacency_to_matrix(single_edges, single_n)
    single_selected = multiple_cp_pairs(
        single;
        max_pairs=4,
        min_pair_size=2,
        n_runs=40,
        pair_selection=:penalized,
        pair_penalty=1.0,
        rng=MersenneTwister(3),
    )
    @test single_selected.n_pairs == 1
    @test single_selected.candidate_pair_counts == collect(1:4)
    @test argmax(single_selected.selection_scores) == 1

    noisy_single = adjacency_to_matrix(
        vcat(single_edges, [(6, 7), (8, 9), (10, 11), (12, 13)]),
        single_n,
    )
    noisy_selected = multiple_cp_pairs(
        noisy_single;
        max_pairs=4,
        min_pair_size=2,
        n_runs=40,
        pair_selection=:penalized,
        pair_penalty=1.0,
        rng=MersenneTwister(4),
    )
    noisy_repeated = multiple_cp_pairs(
        noisy_single;
        max_pairs=4,
        min_pair_size=2,
        n_runs=40,
        pair_selection=:penalized,
        pair_penalty=1.0,
        rng=MersenneTwister(4),
    )
    @test noisy_selected.n_pairs == 1
    @test noisy_selected == noisy_repeated

    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(A; null_model=:invalid)
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(A; max_pairs=0)
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(A; min_pair_size=9)
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(
        A; pair_selection=:invalid,
    )
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(
        A; pair_penalty=-1,
    )
    @test_throws DimensionMismatch CorePeriphery._multiple_cp_pairs_joint(zeros(2, 3))

    asymmetric = copy(A)
    asymmetric[2, 1] = 0.0
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(asymmetric)
    self_loop = copy(A)
    self_loop[1, 1] = 1.0
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(self_loop)
    negative_upper = copy(A)
    negative_upper[1, 2] = -1.0
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(negative_upper)
    negative_lower = copy(A)
    negative_lower[2, 1] = -1.0
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(negative_lower)
    nonfinite_upper = copy(A)
    nonfinite_upper[1, 2] = Inf
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(nonfinite_upper)
    nonfinite_lower = copy(A)
    nonfinite_lower[2, 1] = NaN
    @test_throws ArgumentError CorePeriphery._multiple_cp_pairs_joint(nonfinite_lower)
end
