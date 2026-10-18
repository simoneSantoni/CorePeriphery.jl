using CorePeriphery
using Random
using SparseArrays
using Test

const CP = CorePeriphery

# Exact bimodular surprise of de Jeude, Caldarelli and Squartini (2019),
# equation (4), by rational arithmetic on the full double sum:
# S = Σ_{i ≥ a} Σ_{j ≥ b} C(Vc, i) C(Vcp, j) C(Vpp, L - i - j) / C(V, L).
function exact_log_surprise(n::Int, L::Int, n_core::Int, a::Int, b::Int)
    V = n * (n - 1) ÷ 2
    (L == 0 || L == V) && return 0.0
    n_periphery = n - n_core
    Vc = n_core * (n_core - 1) ÷ 2
    Vcp = n_core * n_periphery
    Vpp = n_periphery * (n_periphery - 1) ÷ 2
    numerator = BigInt(0)
    for i in a:min(Vc, L), j in b:min(Vcp, L - i)
        L - i - j <= Vpp || continue
        numerator += binomial(BigInt(Vc), i) * binomial(BigInt(Vcp), j) *
                     binomial(BigInt(Vpp), L - i - j)
    end
    S = numerator // binomial(BigInt(V), L)
    return Float64(max(0, -log(BigFloat(S))))
end

# The 0.2.0 kernel: direct double loop over the tail accumulated with logaddexp.
function reference_log_surprise(n::Int, L::Int, n_core::Int, a::Int, b::Int)
    V = n * (n - 1) ÷ 2
    logfactorial = zeros(Float64, V + 1)
    for value in 2:V
        logfactorial[value + 1] = logfactorial[value] + log(value)
    end
    logchoose(total, selected) = (selected < 0 || selected > total) ? -Inf :
        logfactorial[total + 1] - logfactorial[selected + 1] -
        logfactorial[total - selected + 1]
    logaddexp(x, y) = x == -Inf ? y : y == -Inf ? x :
        max(x, y) + log1p(exp(min(x, y) - max(x, y)))
    (L == 0 || L == V) && return 0.0
    n_periphery = n - n_core
    Vc = n_core * (n_core - 1) ÷ 2
    Vcp = n_core * n_periphery
    Vpp = n_periphery * (n_periphery - 1) ÷ 2
    log_tail = -Inf
    for i in a:min(Vc, L)
        remaining = L - i
        for j in max(b, remaining - Vpp):min(Vcp, remaining)
            log_tail = logaddexp(log_tail,
                logchoose(Vc, i) + logchoose(Vcp, j) + logchoose(Vpp, remaining - j))
        end
    end
    return max(0.0, logchoose(V, L) - log_tail)
end

function block_counts(A, core::BitVector)
    n = size(A, 1)
    cc = cp = 0
    for j in 2:n, i in 1:(j - 1)
        A[i, j] > 0 || continue
        if core[i] && core[j]
            cc += 1
        elseif core[i] || core[j]
            cp += 1
        end
    end
    return cc, cp
end

function random_binary_graph(rng, n, p)
    A = zeros(Float64, n, n)
    for j in 2:n, i in 1:(j - 1)
        rand(rng) < p && (A[i, j] = A[j, i] = 1.0)
    end
    return A
end

@testset "Surprise kernel: exhaustive small-graph oracle" begin
    # Every partition of every graph in a seeded family with n ≤ 7, compared
    # with the exact rational evaluation of equation (4).
    rng = MersenneTwister(2019)
    checked = 0
    for n in 3:7, trial in 1:6
        A = random_binary_graph(rng, n, rand(rng, (0.2, 0.5, 0.8)))
        L = count(!iszero, A) ÷ 2
        ctx = CP._SurpriseContext(n, L)
        for mask in 0:(2^n - 1)
            core = BitVector(((mask >> (node - 1)) & 1) == 1 for node in 1:n)
            cc, cp = block_counts(A, core)
            value = CP._log_surprise(ctx, count(core), cc, cp)
            exact = exact_log_surprise(n, L, count(core), cc, cp)
            @test isapprox(value, exact; atol=1e-11, rtol=1e-12)
            checked += 1
        end
    end
    @test checked == 6 * sum(2^n for n in 3:7)
    # Star (hub plus five leaves) with the hub as the only core node.
    star = adjacency_to_matrix([(1, node) for node in 2:6], 6)
    ctx = CP._SurpriseContext(6, 5)
    @test CP._log_surprise(ctx, 1, 0, 5) ≈ exact_log_surprise(6, 5, 1, 0, 5) atol = 1e-12
    @test CP._log_surprise(ctx, 1, 0, 5) ≈ log(binomial(BigInt(15), 5)) atol = 1e-12
    # Empty and complete graphs carry no surprise for any partition.
    @test CP._log_surprise(CP._SurpriseContext(5, 0), 2, 0, 0) == 0.0
    @test CP._log_surprise(CP._SurpriseContext(5, 10), 2, 1, 6) == 0.0
end

@testset "Surprise kernel: agreement with the direct double sum" begin
    rng = MersenneTwister(4)
    worst = 0.0
    for trial in 1:200
        n = rand(rng, 5:45)
        V = n * (n - 1) ÷ 2
        L = rand(rng, 1:(V - 1))
        n_core = rand(rng, 0:n)
        n_periphery = n - n_core
        Vc = n_core * (n_core - 1) ÷ 2
        Vcp = n_core * n_periphery
        Vpp = n_periphery * (n_periphery - 1) ÷ 2
        a = rand(rng, max(0, L - Vcp - Vpp):min(Vc, L))
        b = rand(rng, max(0, L - a - Vpp):min(Vcp, L - a))
        value = CP._log_surprise(CP._SurpriseContext(n, L), n_core, a, b)
        reference = reference_log_surprise(n, L, n_core, a, b)
        @test isapprox(value, reference; rtol=1e-10, atol=1e-10)
        worst = max(worst, abs(value - reference))
    end
    @test worst < 1e-10

    # Statistics of partitions visited on actual graphs, dense and sparse.
    for trial in 1:5
        n = 30 + 5 * trial
        A = random_binary_graph(rng, n, 0.15)
        L = count(!iszero, A) ÷ 2
        ctx = CP._SurpriseContext(n, L)
        for _ in 1:20
            core = BitVector(rand(rng, Bool, n))
            cc, cp = block_counts(A, core)
            @test isapprox(CP._log_surprise(ctx, count(core), cc, cp),
                reference_log_surprise(n, L, count(core), cc, cp);
                rtol=1e-10, atol=1e-10)
        end
        dense = surprise_cp(A)
        sparse_result = surprise_cp(sparse(A))
        @test dense.coreness == sparse_result.coreness
        @test dense.quality == sparse_result.quality
        cc, cp = block_counts(A, BitVector(dense.coreness .== 1.0))
        @test dense.quality ≈ reference_log_surprise(n, L, length(dense.core_nodes), cc, cp) rtol = 1e-10
    end
end

@testset "Surprise kernel: log-factorial table and Stirling tail" begin
    lf = CP._LogFactorial(50)
    @test lf(0) == 0.0
    @test lf(1) == 0.0
    @test lf(20) ≈ log(factorial(big(20))) rtol = 1e-14
    @test CP._logchoose(lf, 10, 3) ≈ log(120) rtol = 1e-13
    @test CP._logchoose(lf, 10, 11) == -Inf
    @test CP._logchoose(lf, 10, -1) == -Inf
    # Values above the table use the Stirling series, continuous with the table.
    beyond = CP._LOGFACTORIAL_TABLE_CAP + 3
    @test lf(beyond) == CP._stirling_logfactorial(beyond)
    @test CP._stirling_logfactorial(beyond) ≈ Float64(log(factorial(big(beyond)))) rtol = 1e-14
    @test length(CP._LogFactorial(10 * CP._LOGFACTORIAL_TABLE_CAP).table) ==
          CP._LOGFACTORIAL_TABLE_CAP + 1
    # A network whose dyad count exceeds the table cap runs on the series and
    # agrees between storage formats.
    rng = MersenneTwister(11)
    n = 1600
    sources = Int[]
    targets = Int[]
    for node in 1:(n - 1)
        push!(sources, node)
        push!(targets, node + 1)
    end
    for _ in 1:(2 * n)
        i, j = rand(rng, 1:n, 2)
        i == j && continue
        push!(sources, min(i, j))
        push!(targets, max(i, j))
    end
    large = sparse([sources; targets], [targets; sources], 1.0, n, n)
    large.nzval .= 1.0
    @test n * (n - 1) ÷ 2 > CP._LOGFACTORIAL_TABLE_CAP
    result = surprise_cp(large; max_iter=3)
    @test isfinite(result.quality)
    @test result.quality >= 0
    @test surprise_cp(Matrix(large); max_iter=3).coreness == result.coreness
end
