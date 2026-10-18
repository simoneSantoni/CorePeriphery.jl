using CorePeriphery
using LinearAlgebra
using Logging
using Random
using SparseArrays
using Statistics
using Test

# Independent oracle: the Pearson correlation between adjacency values and
# the ideal pattern, assembled as explicit dyad vectors. `:ones` is Borgatti &
# Everett (2000) eq. (2), max(c_i, c_j); `:missing` is their eq. (4), in which
# core-periphery dyads are missing data and only the core-core (pattern 1)
# and periphery-periphery (pattern 0) dyads are correlated.
function be_oracle_quality(A, c; directed=false, cp_block=:ones)
    n = length(c)
    a = Float64[]
    d = Float64[]
    for i in 1:n, j in 1:n
        i == j && continue
        (directed || i < j) || continue
        if cp_block === :missing
            c[i] == c[j] || continue
            push!(d, c[i])
        else
            push!(d, max(c[i], c[j]))
        end
        push!(a, A[i, j])
    end
    (length(a) >= 2 && var(a) > 0 && var(d) > 0) || return 0.0
    return cor(a, d)
end

function random_symmetric(rng, n, p)
    A = zeros(n, n)
    for i in 1:(n - 1), j in (i + 1):n
        rand(rng) < p && (A[i, j] = A[j, i] = rand(rng, (1.0, 2.0, 3.0)))
    end
    return A
end

function random_directed(rng, n, p)
    A = zeros(n, n)
    for i in 1:n, j in 1:n
        i != j && rand(rng) < p && (A[i, j] = rand(rng, (1.0, 2.0)))
    end
    return A
end

function binary_vectors(n, min_core, max_core)
    return (
        [Float64((mask >> (i - 1)) & 1) for i in 1:n]
        for mask in 0:(2^n - 1)
        if min_core <= count_ones(mask) <= max_core
    )
end

function graph_from_mask(n, mask)
    A = zeros(n, n)
    bit = 0
    for i in 1:(n - 1), j in (i + 1):n
        if (mask >> bit) & 1 == 1
            A[i, j] = A[j, i] = 1.0
        end
        bit += 1
    end
    return A
end

@testset "cp_block=:missing pattern and quality" begin
    Δ = ideal_cp_matrix([1, 1, 0, 0]; discrete=true, cp_block=:missing)
    @test Δ isa Matrix{Union{Missing,Float64}}
    @test isequal(Δ[1, 2], 1.0) && isequal(Δ[2, 1], 1.0)
    @test isequal(Δ[3, 4], 0.0)
    @test ismissing(Δ[1, 3]) && ismissing(Δ[4, 2])
    @test all(i -> isequal(Δ[i, i], 0.0), 1:4)
    @test ideal_cp_matrix([1, 1, 0, 0]; discrete=true) ==
          ideal_cp_matrix([1, 1, 0, 0]; discrete=true, cp_block=:ones)
    @test_throws ArgumentError ideal_cp_matrix([1, 0]; discrete=true, cp_block=:other)
    @test_throws ArgumentError ideal_cp_matrix([1, 0]; cp_block=:missing)
    @test_throws ArgumentError ideal_cp_matrix([1, 0.5]; discrete=true, cp_block=:missing)

    rng = MersenneTwister(11)
    A = random_symmetric(rng, 9, 0.5)
    D = random_directed(rng, 9, 0.4)
    for _ in 1:20
        c = Float64.(rand(rng, Bool, 9))
        for cp_block in (:ones, :missing)
            @test core_quality(A, c; discrete=true, cp_block=cp_block) ≈
                  be_oracle_quality(A, c; cp_block=cp_block) atol=1e-12
            @test core_quality(sparse(A), c; discrete=true, cp_block=cp_block) ≈
                  be_oracle_quality(A, c; cp_block=cp_block) atol=1e-12
            @test core_quality(D, c; discrete=true, directed=true, cp_block=cp_block) ≈
                  be_oracle_quality(D, c; directed=true, cp_block=cp_block) atol=1e-12
            @test core_quality(sparse(D), c; discrete=true, directed=true, cp_block=cp_block) ≈
                  be_oracle_quality(D, c; directed=true, cp_block=cp_block) atol=1e-12
        end
    end
    @test_throws ArgumentError core_quality(A, [1.0, 0.5, 0.0, 0, 0, 0, 0, 0, 0];
                                            discrete=true, cp_block=:missing)
    @test core_quality(A, ones(9); discrete=true, cp_block=:missing) == 0.0
end

@testset "Directed input for the Borgatti-Everett family" begin
    rng = MersenneTwister(5)
    D = random_directed(rng, 8, 0.45)
    A = random_symmetric(rng, 8, 0.5)
    @test !issymmetric(D)

    @test_throws ArgumentError borgatti_everett_discrete(D)
    @test_throws ArgumentError borgatti_everett_continuous(D)
    @test_throws ArgumentError label_switching_cp(D)

    for cp_block in (:ones, :missing)
        discrete = borgatti_everett_discrete(D; directed=true, cp_block=cp_block)
        @test all(x -> x == 0 || x == 1, discrete.coreness)
        @test discrete.quality ≈
              be_oracle_quality(D, discrete.coreness; directed=true, cp_block=cp_block)
        @test borgatti_everett_discrete(sparse(D); directed=true, cp_block=cp_block) == discrete
        # The ordered-dyad objective is transpose-invariant, and so is the fit.
        @test borgatti_everett_discrete(
            Matrix(transpose(D)); directed=true, cp_block=cp_block) == discrete

        switching = label_switching_cp(
            D; directed=true, cp_block=cp_block, rng=MersenneTwister(1), n_runs=3)
        @test switching.quality ≈
              be_oracle_quality(D, switching.coreness; directed=true, cp_block=cp_block)
        @test label_switching_cp(
            sparse(D); directed=true, cp_block=cp_block, rng=MersenneTwister(1), n_runs=3) ==
              switching

        # Symmetric input takes the same path with either flag.
        @test borgatti_everett_discrete(A; directed=true, cp_block=cp_block) ==
              borgatti_everett_discrete(A; cp_block=cp_block)
        @test label_switching_cp(A; directed=true, cp_block=cp_block, rng=MersenneTwister(2)) ==
              label_switching_cp(A; cp_block=cp_block, rng=MersenneTwister(2))
    end

    continuous = borgatti_everett_continuous(D; directed=true)
    @test continuous.quality ≈ core_quality(D, continuous.coreness; directed=true)
    # Continuous oracle: Pearson correlation of A[i, j] with c[i] c[j] over
    # all ordered off-diagonal dyads.
    ordered = [(i, j) for i in 1:8 for j in 1:8 if i != j]
    @test continuous.quality ≈ cor(
        [D[i, j] for (i, j) in ordered],
        [continuous.coreness[i] * continuous.coreness[j] for (i, j) in ordered])
    strengths = vec(sum(D, dims=2)) .+ vec(sum(D, dims=1))
    initial = (strengths .- minimum(strengths)) ./ (maximum(strengths) - minimum(strengths))
    @test continuous.quality >= core_quality(D, initial; directed=true) - 1e-12
    @test borgatti_everett_continuous(sparse(D); directed=true) == continuous
    @test borgatti_everett_continuous(Matrix(transpose(D)); directed=true) == continuous
    @test borgatti_everett_continuous(A; directed=true) == borgatti_everett_continuous(A)
end

@testset "n_runs and rng on borgatti_everett_discrete" begin
    rng = MersenneTwister(21)
    A = random_symmetric(rng, 12, 0.35)
    single = borgatti_everett_discrete(A)
    @test borgatti_everett_discrete(A; n_runs=1, rng=MersenneTwister(1)) == single
    @test borgatti_everett_discrete(A; n_runs=1, rng=MersenneTwister(2)) == single
    multi = borgatti_everett_discrete(A; n_runs=8, rng=MersenneTwister(3))
    @test multi == borgatti_everett_discrete(A; n_runs=8, rng=MersenneTwister(3))
    @test multi.quality >= single.quality
    @test multi.quality ≈ core_quality(A, multi.coreness; discrete=true)
    @test_throws ArgumentError borgatti_everett_discrete(A; n_runs=0)
    @test_throws ArgumentError borgatti_everett_discrete(A; cp_block=:other)
    @test_throws ArgumentError label_switching_cp(A; cp_block=:other)
end

@testset "Brute-force optimum on every graph with n ≤ 6" begin
    # Every labelled graph on 2..6 nodes (34,859 graphs) under both ideal
    # patterns. `method=:exhaustive` must attain the brute-force optimum on
    # all of them; the greedy flipping heuristics return locally optimal
    # partitions whose quality is reported exactly, and they reach the global
    # optimum on the large majority of graphs (the miss counts are asserted
    # below so a regression of the optimiser is caught).
    greedy_misses = 0
    switching_misses = 0
    cases = 0
    # Many of the enumerated graphs are disconnected; the one-time
    # disconnected-input warning is legitimate there and not under test, so
    # the enumeration runs under a null logger.
    Logging.with_logger(Logging.NullLogger()) do
    for n in 2:6
        n_edges = n * (n - 1) ÷ 2
        for mask in 0:(2^n_edges - 1)
            A = graph_from_mask(n, mask)
            isolates = findall(iszero, vec(sum(A, dims=2)))
            for cp_block in (:ones, :missing)
                cases += 1
                # Discrete Borgatti-Everett: core sizes 1:n-1. The optimum over
                # partitions that keep isolates peripheral coincides with the
                # unrestricted optimum, so the isolate guarantee costs nothing
                # (an edgeless graph has no admissible partition and scores 0).
                qualities = [be_oracle_quality(A, c; cp_block=cp_block)
                             for c in binary_vectors(n, 1, n - 1)]
                restricted = [be_oracle_quality(A, c; cp_block=cp_block)
                              for c in binary_vectors(n, 1, n - 1)
                              if all(iszero, c[isolates])]
                optimum = isempty(restricted) ? 0.0 : maximum(restricted)
                @test optimum ≈ maximum(qualities) atol=1e-12

                exact = borgatti_everett_discrete(A; cp_block=cp_block, method=:exhaustive)
                @test exact.quality ≈ optimum atol=1e-9
                @test exact.quality ≈
                      be_oracle_quality(A, exact.coreness; cp_block=cp_block) atol=1e-12
                @test all(iszero, exact.coreness[isolates])

                greedy = borgatti_everett_discrete(
                    A; cp_block=cp_block, n_runs=8, rng=MersenneTwister(mask))
                @test greedy.quality ≈
                      be_oracle_quality(A, greedy.coreness; cp_block=cp_block) atol=1e-12
                @test greedy.quality <= optimum + 1e-9
                isapprox(greedy.quality, optimum; atol=1e-9) || (greedy_misses += 1)

                # Label switching: core sizes 1:n-2 (needs n ≥ 3).
                n >= 3 || continue
                feasible = [be_oracle_quality(A, c; cp_block=cp_block)
                            for c in binary_vectors(n, 1, n - 2)
                            if all(iszero, c[isolates])]
                feasible_optimum = isempty(feasible) ? 0.0 : maximum(feasible)
                switching = label_switching_cp(
                    A; cp_block=cp_block, n_runs=8, rng=MersenneTwister(mask))
                @test switching.quality ≈
                      be_oracle_quality(A, switching.coreness; cp_block=cp_block) atol=1e-12
                @test switching.quality <= feasible_optimum + 1e-9
                isapprox(switching.quality, feasible_optimum; atol=1e-9) ||
                    (switching_misses += 1)
            end
        end
    end
    end
    @test cases == 2 * (2 + 8 + 64 + 1024 + 32768)
    # Measured on this fixed seed set: about 1.5% of cases for either greedy
    # optimiser with 8 starts (single-start greedy misses about 10%).
    @test greedy_misses <= 0.02 * cases
    @test switching_misses <= 0.02 * cases
    @test greedy_misses > 0  # the exhaustive method exists for a reason
end

@testset "Exhaustive method on larger inputs" begin
    rng = MersenneTwister(31)
    A = random_symmetric(rng, 11, 0.3)
    D = random_directed(rng, 9, 0.35)
    for cp_block in (:ones, :missing)
        exact = borgatti_everett_discrete(A; cp_block=cp_block, method=:exhaustive)
        @test exact.quality ≈ maximum(
            be_oracle_quality(A, c; cp_block=cp_block) for c in binary_vectors(11, 1, 10))
        @test exact.quality >= borgatti_everett_discrete(
            A; cp_block=cp_block, n_runs=8, rng=MersenneTwister(1)).quality - 1e-12
        @test borgatti_everett_discrete(sparse(A); cp_block=cp_block, method=:exhaustive) == exact
        exact_directed = borgatti_everett_discrete(
            D; cp_block=cp_block, directed=true, method=:exhaustive)
        @test exact_directed.quality ≈ maximum(
            be_oracle_quality(D, c; directed=true, cp_block=cp_block)
            for c in binary_vectors(9, 1, 8))
    end
    @test_throws ArgumentError borgatti_everett_discrete(A; method=:other)
    @test_throws ArgumentError borgatti_everett_discrete(zeros(25, 25); method=:exhaustive)
end
