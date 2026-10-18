using CorePeriphery
using Random
using SparseArrays
using Test

@testset "Core/periphery partition rule" begin
    @testset "_partition helper" begin
        c = [0.9, 0.5, 0.5, 0.1, 0.5]
        core, periphery = CorePeriphery._partition(c)
        @test core == [1]
        @test periphery == [2, 3, 4, 5]
        @test CorePeriphery._partition(c; rule=:none) == (Int[], collect(1:5))
        @test CorePeriphery._partition(c; rule=:threshold, threshold=0.5) ==
              ([1], [2, 3, 4, 5])
        @test CorePeriphery._partition(c; rule=:threshold, threshold=0.4) ==
              ([1, 2, 3, 5], [4])
        @test CorePeriphery._partition(c; rule=:threshold, threshold=1) ==
              (Int[], collect(1:5))
        # Top-k with ties broken by lower index.
        @test CorePeriphery._partition(c; rule=:core_size, core_size=2) ==
              ([1, 2], [3, 4, 5])
        @test CorePeriphery._partition(c; rule=:core_size, core_size=3) ==
              ([1, 2, 3], [4, 5])
        @test CorePeriphery._partition(c; rule=:core_size, core_size=0) ==
              (Int[], collect(1:5))
        @test CorePeriphery._partition(c; rule=:core_size, core_size=5) ==
              (collect(1:5), Int[])
        @test CorePeriphery._partition([1, 0, 1]; rule=:core_size, core_size=1) ==
              ([1], [2, 3])
        @test CorePeriphery._partition([0.2, 0.2, 0.2]) == (Int[], [1, 2, 3])
        @test CorePeriphery._partition([0.3, 0.7]) == ([2], [1])
        @test CorePeriphery._partition(Float64[]) == (Int[], Int[])

        @test_throws ArgumentError CorePeriphery._partition(c; rule=:unknown)
        @test_throws ArgumentError CorePeriphery._partition(c; rule=:threshold)
        @test_throws ArgumentError CorePeriphery._partition(c; rule=:core_size)
        @test_throws ArgumentError CorePeriphery._partition(c; rule=:core_size, core_size=6)
        @test_throws ArgumentError CorePeriphery._partition(c; rule=:core_size, core_size=-1)
        @test_throws ArgumentError CorePeriphery._partition(c; threshold=0.5)
        @test_throws ArgumentError CorePeriphery._partition(c; core_size=2)
        @test_throws ArgumentError CorePeriphery._partition(c; rule=:threshold, threshold=NaN)
    end

    @testset "Star: median ties fall into the periphery" begin
        star = adjacency_to_matrix([(1, i) for i in 2:6], 6)
        for detector in (
            borgatti_everett_continuous,
            random_walker_profiling,
            minres_symmetric,
            minres_svd,
        )
            result = detector(star)
            @test result.core_nodes == [1]
            @test result.periphery_nodes == collect(2:6)
        end
        directed_star = minres_svd_directed(star)
        @test directed_star.core_nodes == [1]
        @test directed_star.periphery_nodes == collect(2:6)

        # Rombach's coreness is a permutation of n distinct template values, so
        # its median split is not a tie-aware hub detector; the default `:beta`
        # rule keeps the N - ⌊βN⌋ core ranks of the template, and the explicit
        # size rule recovers the hub as well.
        rombach = rombach_continuous(
            star; n_runs=2, rng=MersenneTwister(1), partition=:core_size, core_size=1)
        @test rombach.core_nodes == [1]
        @test rombach.periphery_nodes == collect(2:6)
        by_beta = rombach_continuous(star; beta=0.84, n_runs=2, rng=MersenneTwister(1))
        @test by_beta.core_nodes == [1]
        @test by_beta.periphery_nodes == collect(2:6)
    end

    @testset "partition keywords on every continuous detector" begin
        A = adjacency_to_matrix(
            [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 6), (4, 7)],
            7,
        )
        detectors = (
            borgatti_everett_continuous,
            random_walker_profiling,
            minres_symmetric,
            minres_svd,
            minres_svd_directed,
            # Rombach defaults to `partition=:beta`; the common rules below are
            # exercised through an explicit `:median` default.
            (A; partition=:median, kwargs...) -> rombach_continuous(
                A; n_runs=2, rng=MersenneTwister(1), partition=partition, kwargs...),
        )
        for detector in detectors
            default = detector(A)
            @test detector(A; partition=:median) == default
            @test length(default.core_nodes) + length(default.periphery_nodes) == 7
            @test sort(vcat(default.core_nodes, default.periphery_nodes)) == 1:7

            none = detector(A; partition=:none)
            @test isempty(none.core_nodes)
            @test none.periphery_nodes == collect(1:7)
            @test none.coreness == default.coreness
            @test none.quality == default.quality

            sized = detector(A; partition=:core_size, core_size=2)
            @test length(sized.core_nodes) == 2
            @test length(sized.periphery_nodes) == 5
            @test all(sized.coreness[sized.core_nodes] .>=
                      maximum(sized.coreness[sized.periphery_nodes]))

            cut = detector(A; partition=:threshold, threshold=0.5)
            @test cut.core_nodes == findall(>(0.5), cut.coreness)

            @test_throws ArgumentError detector(A; partition=:unknown)
            @test_throws ArgumentError detector(A; partition=:threshold)
            @test_throws ArgumentError detector(A; core_size=2)
            @test_throws ArgumentError detector(A; partition=:core_size, core_size=8)
        end

        # Sparse input and the ensemble forward the keywords unchanged.
        @test random_walker_profiling(sparse(A); partition=:core_size, core_size=3) ==
              random_walker_profiling(A; partition=:core_size, core_size=3)
        ensemble = rossa_profile_ensemble(
            A; n_runs=3, rng=MersenneTwister(2), partition=:none)
        @test ensemble isa RossaEnsembleResult
    end
end
