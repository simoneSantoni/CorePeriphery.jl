using CorePeriphery
using Logging
using Random
using SparseArrays
using Test

@testset "Result types: equality, hashing, and display" begin
    A = adjacency_to_matrix(
        [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 6), (4, 7)],
        7,
    )
    directed = [
        0.0 2.0 0.0 1.0
        0.0 0.0 3.0 0.0
        4.0 0.0 0.0 1.0
        1.0 2.0 0.0 0.0
    ]

    results = Any[
        borgatti_everett_discrete(A),
        minres_svd_directed(directed),
        multiple_cp_pairs(A; n_runs=2, rng=MersenneTwister(1)),
        # Configuration-null samples of this small graph can be disconnected,
        # which triggers the detector's one-time warning; not under test here.
        Logging.with_logger(Logging.NullLogger()) do
            cp_significance(A, lip_discrete; n_samples=4, rng=MersenneTwister(2), n_swaps=10)
        end,
        rossa_stationary_distribution(directed),
        rossa_profile_ensemble(A; n_runs=3, rng=MersenneTwister(3)),
    ]
    repeats = Any[
        borgatti_everett_discrete(A),
        minres_svd_directed(directed),
        multiple_cp_pairs(A; n_runs=2, rng=MersenneTwister(1)),
        # Configuration-null samples of this small graph can be disconnected,
        # which triggers the detector's one-time warning; not under test here.
        Logging.with_logger(Logging.NullLogger()) do
            cp_significance(A, lip_discrete; n_samples=4, rng=MersenneTwister(2), n_swaps=10)
        end,
        rossa_stationary_distribution(directed),
        rossa_profile_ensemble(A; n_runs=3, rng=MersenneTwister(3)),
    ]
    expected_types = (
        CPResult, CPDirectedResult, CPMultiResult, CPSignificanceResult,
        StationaryDistributionResult, RossaEnsembleResult,
    )

    for (result, repeat, T) in zip(results, repeats, expected_types)
        @test result isa T
        @test result == repeat
        @test isequal(result, repeat)
        @test hash(result) == hash(repeat)
        @test length(Set([result, repeat])) == 1
        @test length(unique([result, repeat])) == 1
        @test haskey(Dict(result => 1), repeat)

        # A compact multi-line summary that never dumps whole vectors.
        text = sprint(show, MIME("text/plain"), result)
        @test occursin(string(nameof(T)), text)
        @test count(==('\n'), text) <= 3
        @test length(text) < 300
        @test !occursin("[", text)
    end

    # Different contents are neither equal nor (with overwhelming probability)
    # hash-equal.
    first = results[1]
    other = CPResult(first.coreness, first.core_nodes, first.periphery_nodes,
        first.quality + 1e-9, first.algorithm)
    @test other != first
    @test hash(other) != hash(first)
    renamed = CPResult(first.coreness, first.core_nodes, first.periphery_nodes,
        first.quality, "Other")
    @test renamed != first
    @test hash(renamed) != hash(first)

    # Equality follows field semantics: NaN quality is not == to itself but is
    # isequal, so hashing stays consistent.
    two_node = random_walker_profiling([0.0 1.0; 1.0 0.0])
    @test isnan(two_node.quality)
    @test !(two_node == two_node)
    @test isequal(two_node, two_node)
    @test hash(two_node) == hash(random_walker_profiling([0.0 1.0; 1.0 0.0]))

    # Null diagnostics are a concrete struct with the documented field names.
    diagnostics = results[4].null_diagnostics
    @test diagnostics isa CorePeriphery.NullDiagnostics
    @test isconcretetype(fieldtype(CPSignificanceResult, :null_diagnostics))
    @test propertynames(diagnostics) == (
        :requested_swaps, :accepted_swaps, :swap_attempts, :max_swap_attempts,
        :complete, :preserved,
    )
    @test diagnostics.complete
    @test diagnostics.preserved == :degree_sequence
    @test diagnostics == repeats[4].null_diagnostics
    @test hash(diagnostics) == hash(repeats[4].null_diagnostics)
    @test occursin("NullDiagnostics", sprint(show, MIME("text/plain"), diagnostics))

    # Compatibility constructor and manual construction still work.
    manual = CPResult([0.8, 0.9, 0.2, 0.1], [1, 2], [3, 4], 0.75, "Test")
    @test manual == CPResult([0.8, 0.9, 0.2, 0.1], [1, 2], [3, 4], 0.75, "Test")
    @test occursin("core = 2", sprint(show, MIME("text/plain"), manual))
    multi = CPMultiResult([1, 1, 2, 2], [1.0, 0.0, 1.0, 0.0], 2, 0.5, "Multi")
    @test multi.pair_selection == :objective
    @test occursin("pairs = 2", sprint(show, MIME("text/plain"), multi))

    # ConvergenceError is exported and prints its message.
    @test ConvergenceError <: Exception
    @test sprint(showerror, ConvergenceError("solver stalled")) ==
          "ConvergenceError: solver stalled"
end
