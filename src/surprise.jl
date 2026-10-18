function _binary_block_counts(A::AbstractMatrix{<:Real}, c::AbstractVector{<:Real})
    n = length(c)
    core_core = 0
    core_periphery = 0
    periphery_periphery = 0
    @inbounds for j in 2:n, i in 1:(j - 1)
        A[i, j] > 0 || continue
        if !iszero(c[i]) && !iszero(c[j])
            core_core += 1
        elseif iszero(c[i]) && iszero(c[j])
            periphery_periphery += 1
        else
            core_periphery += 1
        end
    end
    return core_core, core_periphery, periphery_periphery
end

function _binary_block_counts(A::SparseMatrixCSC{<:Real}, c::AbstractVector{<:Real})
    core_core = 0
    core_periphery = 0
    periphery_periphery = 0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        i < j && values[index] > 0 || continue
        if !iszero(c[i]) && !iszero(c[j])
            core_core += 1
        elseif iszero(c[i]) && iszero(c[j])
            periphery_periphery += 1
        else
            core_periphery += 1
        end
    end
    return core_core, core_periphery, periphery_periphery
end

function _binary_neighbor_counts(A::AbstractMatrix{<:Real}, node::Int,
    c::AbstractVector{<:Real})
    core = 0
    periphery = 0
    @inbounds for neighbor in axes(A, 1)
        neighbor == node && continue
        A[neighbor, node] > 0 || continue
        if iszero(c[neighbor])
            periphery += 1
        else
            core += 1
        end
    end
    return core, periphery
end

function _binary_neighbor_counts(A::SparseMatrixCSC{<:Real}, node::Int,
    c::AbstractVector{<:Real})
    core = 0
    periphery = 0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, node)
        neighbor = rows[index]
        (neighbor != node && values[index] > 0) || continue
        if iszero(c[neighbor])
            periphery += 1
        else
            core += 1
        end
    end
    return core, periphery
end

# ---------------------------------------------------------------------------
# Log-factorials and the bivariate hypergeometric tail of de Jeude et al. (2019)
# ---------------------------------------------------------------------------

# Exact cumulative log-factorial table up to `_LOGFACTORIAL_TABLE_CAP` and the
# Stirling series beyond it. The table is sized by the number of dyads, so it
# stays below 8 MB (n ≲ 1450); larger networks switch to the series, whose
# truncation error (< 1 / (1680 N^7)) is far below double precision there.
const _LOGFACTORIAL_TABLE_CAP = 1 << 20

struct _LogFactorial
    table::Vector{Float64}
end

function _LogFactorial(max_value::Int)
    size = clamp(max_value, 1, _LOGFACTORIAL_TABLE_CAP)
    table = zeros(Float64, size + 1)
    @inbounds for value in 2:size
        table[value + 1] = table[value] + log(value)
    end
    return _LogFactorial(table)
end

@inline function (lf::_LogFactorial)(value::Int)
    value < length(lf.table) && return @inbounds lf.table[value + 1]
    return _stirling_logfactorial(value)
end

function _stirling_logfactorial(value::Int)
    x = Float64(value)
    inv_x = 1 / x
    inv_x2 = inv_x * inv_x
    return (x + 0.5) * log(x) - x + 0.9189385332046728 +
           inv_x * (1 / 12 - inv_x2 * (1 / 360 - inv_x2 * (1 / 1260 - inv_x2 / 1680)))
end

@inline function _logchoose(lf::_LogFactorial, total::Int, selected::Int)
    (selected < 0 || selected > total) && return -Inf
    return lf(total) - lf(selected) - lf(total - selected)
end

# Relative truncation tolerance of the tail sums. Terms of a log-concave
# sequence decay at least geometrically beyond its mode, so once the geometric
# bound on the remaining terms drops below this fraction of the accumulated
# sum the omitted mass is below double-precision resolution of the result.
const _SURPRISE_TOL = 1e-18

# log Σ_{j = lo}^{hi} C(K, j) C(M, d - j): a contiguous slice of the
# hypergeometric sum, evaluated from its mode outwards with O(1) ratio updates
# per term and geometric truncation. Returns -Inf for an empty slice.
function _log_hypergeometric_slice(lf::_LogFactorial, K::Int, M::Int, d::Int,
    lo::Int, hi::Int)
    lo = max(lo, d - M, 0)
    hi = min(hi, K, d)
    lo > hi && return -Inf
    mode = clamp(((d + 1) * (K + 1)) ÷ (K + M + 2), lo, hi)
    log_mode = _logchoose(lf, K, mode) + _logchoose(lf, M, d - mode)
    acc = 1.0
    term = 1.0
    j = mode
    @inbounds while j < hi
        ratio = ((K - j) / (j + 1)) * ((d - j) / (M - d + j + 1))
        term *= ratio
        j += 1
        acc += term
        ratio < 1 && term * ratio / (1 - ratio) <= _SURPRISE_TOL * acc && break
    end
    term = 1.0
    j = mode
    @inbounds while j > lo
        ratio = (j / (K - j + 1)) * ((M - d + j) / (d - j + 1))
        term *= ratio
        j -= 1
        acc += term
        ratio < 1 && term * ratio / (1 - ratio) <= _SURPRISE_TOL * acc && break
    end
    return log_mode + log(acc)
end

"""
    _log_bimodular_tail(lf, Vc, Vcp, Vpp, L, a, b)

Logarithm of the numerator of de Jeude et al. (2019) equation (4): the number
of ways to place `L` links on `Vc + Vcp + Vpp` dyads with at least `a` on the
`Vc` core-core dyads and at least `b` on the `Vcp` core-periphery dyads,

    Σ_{i ≥ a} Σ_{j ≥ b} C(Vc, i) C(Vcp, j) C(Vpp, L - i - j).

The inner sum over `j` is a hypergeometric slice (see
[`_log_hypergeometric_slice`](@ref)); the outer sequence
`g(i) = C(Vc, i) · slice(L - i)` is log-concave in `i` (a binomial times a
convolution of log-concave sequences), so its mode is found by bisection and
the sum is accumulated outwards with geometric truncation. Each query costs
`O((log L + w_i) · w_j)` with `w` the widths of the hypergeometric peaks,
instead of the `O(L²)` double loop of the direct sum, and no `O(n²)`
per-core-size tables are built.
"""
function _log_bimodular_tail(lf::_LogFactorial, Vc::Int, Vcp::Int, Vpp::Int,
    L::Int, a::Int, b::Int)
    g(i) = _logchoose(lf, Vc, i) + _log_hypergeometric_slice(lf, Vcp, Vpp, L - i, b, Vcp)
    lo = max(a, L - Vcp - Vpp, 0)
    hi = min(Vc, L - b)
    lo > hi && return -Inf

    left, right = lo, hi
    while left < right
        mid = (left + right) >> 1
        if g(mid) < g(mid + 1)
            left = mid + 1
        else
            right = mid
        end
    end
    mode = left
    log_mode = g(mode)
    isfinite(log_mode) || return -Inf

    acc = 1.0
    previous = log_mode
    i = mode
    while i < hi
        i += 1
        current = g(i)
        ratio = exp(current - previous)
        term = exp(current - log_mode)
        acc += term
        previous = current
        ratio < 1 && term * ratio / (1 - ratio) <= _SURPRISE_TOL * acc && break
    end
    previous = log_mode
    i = mode
    while i > lo
        i -= 1
        current = g(i)
        ratio = exp(current - previous)
        term = exp(current - log_mode)
        acc += term
        previous = current
        ratio < 1 && term * ratio / (1 - ratio) <= _SURPRISE_TOL * acc && break
    end
    return log_mode + log(acc)
end

# Network-level constants of the surprise statistic.
struct _SurpriseContext
    lf::_LogFactorial
    n::Int
    total_dyads::Int
    total_edges::Int
    log_denominator::Float64
end

function _SurpriseContext(n::Int, total_edges::Int)
    total_dyads = n * (n - 1) ÷ 2
    lf = _LogFactorial(total_dyads)
    return _SurpriseContext(lf, n, total_dyads, total_edges,
        _logchoose(lf, total_dyads, total_edges))
end

# Negative log of the bimodular surprise S (equation 4) for a partition with
# `n_core` core nodes, `observed_cc` core-core links and `observed_cp`
# core-periphery links. Zero when the graph is empty or complete.
function _log_surprise(ctx::_SurpriseContext, n_core::Int, observed_cc::Int,
    observed_cp::Int)
    (ctx.total_edges == 0 || ctx.total_edges == ctx.total_dyads) && return 0.0
    n_periphery = ctx.n - n_core
    capacity_cc = n_core * (n_core - 1) ÷ 2
    capacity_cp = n_core * n_periphery
    capacity_pp = n_periphery * (n_periphery - 1) ÷ 2
    log_tail = _log_bimodular_tail(ctx.lf, capacity_cc, capacity_cp, capacity_pp,
        ctx.total_edges, observed_cc, observed_cp)
    # An empty tail cannot occur for feasible block counts; treat it as no surprise.
    isfinite(log_tail) || return 0.0
    return max(0.0, ctx.log_denominator - log_tail)
end

"""
    surprise_cp(A; max_iter=100)

Greedily maximize the exact bivariate-hypergeometric joint-tail surprise of
de Jeude, Caldarelli and Squartini (2019, equation 4) for a simple, undirected
binary graph. The returned quality is the negative log tail probability.

The tail is evaluated from the peak of the hypergeometric terms outwards and
truncated once the geometrically bounded remainder is below double-precision
resolution, so every candidate flip costs roughly the width of the
hypergeometric peaks rather than the full `O(m²)` double sum. The optimizer is
single-start greedy node flipping from a mean-degree split.
"""
function surprise_cp(A::AbstractMatrix{<:Real}; max_iter::Int=100)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true, binary=true)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n == 1 && return CPResult([1.0], [1], Int[], 0.0, "Surprise CP")

    # Initialize based on degree
    degrees = vec(sum(A, dims=2))
    threshold = mean(degrees)
    c = Float64.(degrees .>= threshold)

    total_edges = Int(round(sum(degrees) / 2))
    ctx = _SurpriseContext(n, total_edges)
    log_surprise(n_core, observed_cc, observed_cp) =
        _log_surprise(ctx, n_core, observed_cc, observed_cp)

    n_core = count(x -> !iszero(x), c)
    edges_cc, edges_cp, edges_pp = _binary_block_counts(A, c)
    best_surprise = log_surprise(n_core, edges_cc, edges_cp)

    # Greedy optimization
    for _ in 1:max_iter
        improved = false

        for i in 1:n
            links_core, links_periphery = _binary_neighbor_counts(A, i, c)
            moving_to_periphery = !iszero(c[i])
            if moving_to_periphery
                candidate_n_core = n_core - 1
                candidate_cc = edges_cc - links_core
                candidate_cp = edges_cp + links_core - links_periphery
                candidate_pp = edges_pp + links_periphery
            else
                candidate_n_core = n_core + 1
                candidate_cc = edges_cc + links_core
                candidate_cp = edges_cp + links_periphery - links_core
                candidate_pp = edges_pp - links_periphery
            end
            s_new = log_surprise(candidate_n_core, candidate_cc, candidate_cp)
            if s_new > best_surprise
                best_surprise = s_new
                c[i] = 1.0 - c[i]
                n_core = candidate_n_core
                edges_cc = candidate_cc
                edges_cp = candidate_cp
                edges_pp = candidate_pp
                improved = true
            end
        end

        if !improved
            break
        end
    end

    best_c = c
    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, best_surprise, "Surprise CP")
end
