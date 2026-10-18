using CorePeriphery
using Logging
using Random
using SparseArrays
using Test

const _CPC = CorePeriphery

function two_cliques(k)
    A = zeros(2k, 2k)
    for offset in (0, k), i in 1:(k - 1), j in (i + 1):k
        A[offset + i, offset + j] = A[offset + j, offset + i] = 1.0
    end
    return A
end

# K4 on nodes 1:4, pendant 5 attached to node 1, isolate 6.
function clique_pendant_isolate()
    A = zeros(6, 6)
    for i in 1:3, j in (i + 1):4
        A[i, j] = A[j, i] = 1.0
    end
    A[1, 5] = A[5, 1] = 1.0
    return A
end

const owned_detectors = (
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    lip_discrete,
    label_switching_cp,
    (A; kwargs...) -> rombach_continuous(A; n_runs=2, rng=MersenneTwister(1), kwargs...),
    (A; kwargs...) -> rombach_core_score(
        A; alphas=0.5:0.5:1.0, betas=0.5:0.5:1.0, n_runs=1, rng=MersenneTwister(1), kwargs...),
    minres_symmetric,
    minres_svd,
    minres_svd_directed,
    (A; kwargs...) -> multiple_cp_pairs(A; n_runs=4, rng=MersenneTwister(1), kwargs...),
)

@testset "Component labels and isolates" begin
    A = clique_pendant_isolate()
    @test _CPC._component_labels(A) == [1, 1, 1, 1, 1, 2]
    @test _CPC._component_labels(two_cliques(4)) == [1, 1, 1, 1, 2, 2, 2, 2]
    @test _CPC._component_labels(sparse(two_cliques(3))) == [1, 1, 1, 2, 2, 2]
    @test _CPC._component_labels(zeros(3, 3)) == [1, 2, 3]
    directed = [0.0 1.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test _CPC._component_labels(directed) == [1, 1, 2]
    @test _CPC._isolates(A) == [6]
    @test _CPC._isolates(directed) == [3]
    @test _CPC._isolates(sparse(A)) == [6]
    @test isempty(_CPC._isolates(two_cliques(3)))
    coreness = [0.5, 0.2, 0.9]
    @test _CPC._enforce_isolates!(coreness, [1, 3], [2], [3]) == ([1], [2, 3])
    @test coreness == [0.5, 0.2, 0.0]
end

@testset "Two disjoint K4 warn, a connected graph does not" begin
    disconnected = two_cliques(4)
    connected = clique_pendant_isolate()[1:5, 1:5]
    for detector in owned_detectors
        @test_logs (:warn, r"connected components with at least two nodes") match_mode=:any detector(disconnected)
        @test_logs min_level=Logging.Warn detector(connected)
        # An isolate is a trivial component and does not trigger the warning.
        @test_logs min_level=Logging.Warn detector(clique_pendant_isolate())
    end
    @test (@test_logs (:warn, r"test_detector: the network has 2 connected components") _CPC._warn_disconnected(
        disconnected, :test_detector)) == 2
    @test _CPC._warn_disconnected(connected, :test_detector) == 1
    @test _CPC._warn_disconnected(zeros(4, 4), :test_detector) == 0

    # The warning carries `maxlog=1` keyed on the detector name, so a logger
    # that honours `maxlog` (the default console logger and `TestLogger`)
    # shows it at most once per detector for the lifetime of that logger,
    # however many disconnected inputs the detector sees; a second detector
    # has its own budget.
    logger = Test.TestLogger()
    Logging.with_logger(logger) do
        lip_discrete(disconnected)
        lip_discrete(two_cliques(3))
        lip_discrete(disconnected)
        borgatti_everett_discrete(disconnected)
        borgatti_everett_discrete(disconnected)
    end
    warnings = [record for record in logger.logs if record.level == Logging.Warn]
    @test length(warnings) == 2
    @test sort([record.id for record in warnings]) == [
        :disconnected_input_borgatti_everett_discrete, :disconnected_input_lip_discrete]
end

@testset "Isolates are never core" begin
    A = clique_pendant_isolate()
    for detector in owned_detectors
        result = detector(A)
        if result isa CPMultiResult
            @test result.pair_labels[6] == 0
            @test result.coreness[6] == 0.0
            @test all(>=(1), result.pair_labels[1:5])
            @test sort(unique(result.pair_labels[1:5])) == collect(1:result.n_pairs)
            @test result.pair_labels[1:4] == fill(result.pair_labels[1], 4)
            @test result.coreness[1] == 1.0
        else
            @test result.coreness[6] == 0.0
            @test 6 ∉ result.core_nodes
            @test 6 in result.periphery_nodes
            @test sort(vcat(result.core_nodes, result.periphery_nodes)) == 1:6
            @test length(result.coreness) == 6
            @test 1 in result.core_nodes
        end
        if result isa CPDirectedResult
            @test result.out_coreness[6] == 0.0
            @test result.in_coreness[6] == 0.0
        end
    end

    # Explicit size rules cannot pull an isolate into the core either.
    for detector in (borgatti_everett_continuous, minres_symmetric, minres_svd)
        sized = detector(A; partition=:core_size, core_size=6)
        @test 6 ∉ sized.core_nodes
        @test sized.core_nodes == collect(1:5)
    end
    @test 6 ∉ rombach_continuous(A; beta=0.0, n_runs=1).core_nodes

    # A network without edges has no core at all.
    empty = zeros(4, 4)
    for detector in (borgatti_everett_discrete, lip_discrete, borgatti_everett_continuous,
                     minres_symmetric, minres_svd,
                     (A; kwargs...) -> rombach_continuous(A; n_runs=1, kwargs...),
                     (A; kwargs...) -> rombach_core_score(
                         A; alphas=[0.5], betas=[0.5], n_runs=1, kwargs...))
        result = detector(empty)
        @test isempty(result.core_nodes)
        @test result.coreness == zeros(4)
    end
    all_isolates = multiple_cp_pairs(empty)
    @test all_isolates.pair_labels == zeros(Int, 4)
    @test all_isolates.n_pairs == 0
    @test all_isolates.coreness == zeros(4)

    # Detectors that drop isolates before fitting (Lip's degree prefix, Rombach,
    # multiple pairs) reproduce the fit of the non-isolated subnetwork. The
    # Borgatti-Everett correlations keep isolates as peripheral dyads of the
    # sociogram, so their fit may differ from the subnetwork's; the exhaustive
    # optimum nevertheless never puts an isolate in the core.
    core = A[1:5, 1:5]
    @test lip_discrete(A).coreness[1:5] == lip_discrete(core).coreness
    @test borgatti_everett_discrete(A; method=:exhaustive).coreness[6] == 0.0
    @test 6 ∉ borgatti_everett_discrete(A; method=:exhaustive, cp_block=:missing).core_nodes
    @test rombach_continuous(A; n_runs=2, rng=MersenneTwister(1)).coreness[1:5] ==
          rombach_continuous(core; n_runs=2, rng=MersenneTwister(1)).coreness
    multi = multiple_cp_pairs(A; n_runs=4, rng=MersenneTwister(1))
    multi_core = multiple_cp_pairs(core; n_runs=4, rng=MersenneTwister(1))
    @test multi.pair_labels[1:5] == multi_core.pair_labels
    @test multi.coreness[1:5] == multi_core.coreness
    @test multi.quality == multi_core.quality
    @test multi.candidate_qualities == multi_core.candidate_qualities
    @test multiple_cp_pairs(sparse(A); n_runs=4, rng=MersenneTwister(1)) == multi
end

@testset "Null samples skip the disconnected-input check" begin
    # A 40-node ring is connected, but ER null samples at its density (~0.05)
    # are almost always disconnected. They must not warn.
    n = 40
    ring = adjacency_to_matrix([(i, mod1(i + 1, n)) for i in 1:n], n)
    for storage in (identity, sparse)
        @test_logs min_level = Logging.Warn cp_significance(
            storage(ring), lip_discrete; null_model=:er, n_samples=30,
            rng=MersenneTwister(4), threaded=false)
        @test_logs min_level = Logging.Warn cp_significance(
            storage(ring), lip_discrete; null_model=:er, n_samples=30,
            rng=MersenneTwister(4), threaded=true)
    end

    # The observed network is still checked.
    @test_logs (:warn, r"connected components") match_mode = :any cp_significance(
        two_cliques(4), lip_discrete; null_model=:er, n_samples=10,
        rng=MersenneTwister(5), threaded=false)

    # Suppression is task-scoped: it is off outside the evaluation.
    @test !_CPC._suppressing_component_checks()
    @test _CPC._with_component_checks_suppressed(_CPC._suppressing_component_checks)
    @test !_CPC._suppressing_component_checks()

    # Statistics are unchanged by the suppression (same estimator, same seeds).
    serial = cp_significance(ring, lip_discrete; null_model=:er, n_samples=30,
        rng=MersenneTwister(4), threaded=false)
    @test isequal(serial, cp_significance(sparse(ring), lip_discrete; null_model=:er,
        n_samples=30, rng=MersenneTwister(4), threaded=false))
end
