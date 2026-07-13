"""
    CorePeriphery

A Julia module implementing various core-periphery detection algorithms for network analysis.

Core-periphery structure is a mesoscale pattern in networks where nodes are divided into
a densely interconnected core and a sparsely connected periphery.

# Algorithms Implemented
- Borgatti-Everett continuous model (2000)
- Borgatti-Everett discrete model (2000)
- Lip's fast discrete algorithm (2011)
- Rombach's generalized model (2017)
- Spectral method (Cucuringu et al., 2016)
- Random walker profiling (Rossa et al., 2013)
- MINRES/SVD for asymmetric networks (Boyd et al., 2010)
- Multiple core-periphery pairs (Kojaku & Masuda, 2017)
- Surprise-based detection (Jeude et al., 2019)
- Label-switching algorithm (Yanchenko & Sengupta, 2025)

# References
- Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
- Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
- Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
- Cucuringu, M., et al. (2016). Detection of core-periphery structure using spectral methods.
- Rossa, F.D., et al. (2013). Profiling core-periphery network structure by random walkers.
- Boyd, J.P., et al. (2010). Computing core/periphery structures and permutation tests.
- Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks.
- Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise.
- Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection.
"""
module CorePeriphery

using LinearAlgebra
using Statistics
using Random
using SparseArrays

export
    # Data structures
    CPResult,
    CPMultiResult,
    CPDirectedResult,
    CPSignificanceResult,
    StationaryDistributionResult,
    RossaEnsembleResult,
    # Algorithms
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    lip_discrete,
    rombach_continuous,
    spectral_method,
    random_walker_profiling,
    rossa_stationary_distribution,
    rossa_profile_ensemble,
    minres_svd,
    minres_symmetric,
    minres_svd_directed,
    multiple_cp_pairs,
    multiple_cp_pairs_config,
    surprise_cp,
    label_switching_cp,
    cp_significance,
    # Utilities
    coreness_scores,
    core_quality,
    ideal_cp_matrix,
    adjacency_to_matrix

"""
    CPResult

Result structure for core-periphery detection algorithms.

# Fields
- `coreness::Vector{Float64}`: Coreness score for each node (higher = more core-like)
- `core_nodes::Vector{Int}`: Indices of nodes classified as core
- `periphery_nodes::Vector{Int}`: Indices of nodes classified as periphery
- `quality::Float64`: Quality score of the partition (algorithm-specific)
- `algorithm::String`: Name of the algorithm used
"""
struct CPResult
    coreness::Vector{Float64}
    core_nodes::Vector{Int}
    periphery_nodes::Vector{Int}
    quality::Float64
    algorithm::String
end

# Common adjacency-matrix contract. Self-loops are ignored by every objective
# and are therefore rejected rather than silently discarded.
function _validate_adjacency(A::AbstractMatrix{<:Real};
    symmetric::Bool=false,
    binary::Bool=false)
    n, m = size(A)
    n == m || throw(DimensionMismatch("adjacency matrix must be square, got $(n)×$(m)"))
    n > 0 || throw(ArgumentError("adjacency matrix must contain at least one node"))

    for j in axes(A, 2), i in axes(A, 1)
        a = A[i, j]
        isfinite(a) || throw(ArgumentError("adjacency matrix entries must be finite"))
        a >= zero(a) || throw(ArgumentError("adjacency matrix entries must be nonnegative"))
        if binary && !(iszero(a) || a == one(a))
            throw(ArgumentError("this algorithm requires a binary adjacency matrix"))
        end
    end
    all(iszero, diag(A)) || throw(ArgumentError("self-loops are not supported"))
    if symmetric && !issymmetric(A)
        throw(ArgumentError("this algorithm requires an undirected (symmetric) adjacency matrix"))
    end
    return n
end

function _validate_adjacency(A::SparseMatrixCSC{<:Real};
    symmetric::Bool=false,
    binary::Bool=false)
    n, m = size(A)
    n == m || throw(DimensionMismatch("adjacency matrix must be square, got $(n)×$(m)"))
    n > 0 || throw(ArgumentError("adjacency matrix must contain at least one node"))

    rows = rowvals(A)
    values = nonzeros(A)
    for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        a = values[index]
        isfinite(a) || throw(ArgumentError("adjacency matrix entries must be finite"))
        a >= zero(a) || throw(ArgumentError("adjacency matrix entries must be nonnegative"))
        i == j && !iszero(a) && throw(ArgumentError("self-loops are not supported"))
        if binary && !(iszero(a) || a == one(a))
            throw(ArgumentError("this algorithm requires a binary adjacency matrix"))
        end
    end
    if symmetric && !issymmetric(A)
        throw(ArgumentError("this algorithm requires an undirected (symmetric) adjacency matrix"))
    end
    return n
end

function _normalize_scores(x::AbstractVector{<:Real})
    isempty(x) && return Float64[]
    lo, hi = extrema(x)
    hi > lo || return fill(0.5, length(x))
    return (Float64.(x) .- lo) ./ (hi - lo)
end

"""
    adjacency_to_matrix(edges, n)

Convert edge list to adjacency matrix.

Edges are undirected, endpoints are one-based, and self-loops, negative weights,
and non-finite weights are rejected. If an undirected edge occurs more than once,
the last supplied weight wins.

# Arguments
- `edges`: Vector of tuples (i, j) or (i, j, weight)
- `n`: Number of nodes

# Returns
- Symmetric adjacency matrix
"""
function adjacency_to_matrix(edges::AbstractVector{<:Tuple{<:Integer,<:Integer}}, n::Integer)
    n >= 1 || throw(ArgumentError("n must be positive"))
    A = zeros(Float64, n, n)
    for (i, j) in edges
        1 <= i <= n && 1 <= j <= n || throw(BoundsError(A, (i, j)))
        i == j && throw(ArgumentError("self-loops are not supported"))
        w = 1.0
        A[i, j] = w
        A[j, i] = w
    end
    return A
end


function adjacency_to_matrix(
    edges::AbstractVector{<:Tuple{<:Integer,<:Integer,<:Real}}, n::Integer)
    n >= 1 || throw(ArgumentError("n must be positive"))
    A = zeros(Float64, n, n)
    for (i, j, weight) in edges
        1 <= i <= n && 1 <= j <= n || throw(BoundsError(A, (i, j)))
        i == j && throw(ArgumentError("self-loops are not supported"))
        isfinite(weight) && weight >= 0 ||
            throw(ArgumentError("edge weights must be finite and nonnegative"))
        w = Float64(weight)
        A[i, j] = w
        A[j, i] = w
    end
    return A
end

# Compatibility fallback for deliberately heterogeneous edge containers. The
# concrete two- and three-tuple methods above are the allocation-free fast path.
function adjacency_to_matrix(edges::AbstractVector, n::Integer)
    n >= 1 || throw(ArgumentError("n must be positive"))
    A = zeros(Float64, n, n)
    for edge in edges
        edge isa Tuple || throw(ArgumentError("each edge must be a tuple"))
        length(edge) in (2, 3) ||
            throw(ArgumentError("each edge must be (i, j) or (i, j, weight)"))
        i, j = edge[1], edge[2]
        i isa Integer && j isa Integer ||
            throw(ArgumentError("edge endpoints must be integers"))
        1 <= i <= n && 1 <= j <= n || throw(BoundsError(A, (i, j)))
        i == j && throw(ArgumentError("self-loops are not supported"))
        weight = length(edge) == 2 ? 1.0 : edge[3]
        weight isa Real && isfinite(weight) && weight >= 0 ||
            throw(ArgumentError("edge weights must be finite and nonnegative"))
        A[i, j] = A[j, i] = Float64(weight)
    end
    return A
end

"""
    ideal_cp_matrix(c)

Generate the ideal core-periphery pattern matrix for a given coreness vector.

For continuous model: Δ[i,j] = c[i] * c[j]
For discrete model: Δ[i,j] = max(c[i], c[j])

# Arguments
- `c::Vector{Float64}`: Coreness vector (values in [0,1])
- `discrete::Bool`: If true, use discrete ideal pattern

# Returns
- Ideal core-periphery matrix
"""
function ideal_cp_matrix(c::Vector{Float64}; discrete::Bool=false)
    n = length(c)
    Δ = zeros(Float64, n, n)

    if discrete
        for i in 1:n
            for j in 1:n
                if i != j
                    Δ[i, j] = max(c[i], c[j])
                end
            end
        end
    else
        for i in 1:n
            for j in 1:n
                if i != j
                    Δ[i, j] = c[i] * c[j]
                end
            end
        end
    end

    return Δ
end

#=
    core_quality(A, c; discrete=false, directed=false)

Compute quality (correlation) between adjacency matrix and ideal core-periphery pattern.

For undirected data, correlation is accumulated over the strict upper triangle.
With `directed=true`, all ordered off-diagonal dyads are used. No ideal matrix
or dyad vectors are materialized. If there are fewer than two dyads, or if
either dyad vector has zero variance, quality is defined as `0.0`.

# Arguments
- `A`: Adjacency matrix
- `c`: Coreness vector
- `discrete`: Use discrete ideal pattern
- `directed`: Use every ordered off-diagonal dyad

# Returns
- Pearson correlation coefficient
=#
function _adjacency_pattern_product(A::AbstractMatrix{<:Real},
    c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    n = length(c)
    sum_a = 0.0
    sum_a2 = 0.0
    sum_ad = 0.0
    if directed
        @inbounds for j in 1:n, i in 1:n
            i == j && continue
            a = Float64(A[i, j])
            d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                           Float64(c[i]) * Float64(c[j])
            sum_a += a
            sum_a2 += a * a
            sum_ad += a * d
        end
    else
        @inbounds for j in 2:n, i in 1:(j - 1)
            a = Float64(A[i, j])
            d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                           Float64(c[i]) * Float64(c[j])
            sum_a += a
            sum_a2 += a * a
            sum_ad += a * d
        end
    end
    return sum_a, sum_a2, sum_ad
end

function _adjacency_pattern_product(A::SparseMatrixCSC{<:Real},
    c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    sum_a = 0.0
    sum_a2 = 0.0
    sum_ad = 0.0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        (directed ? i != j : i < j) || continue
        a = Float64(values[index])
        d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                       Float64(c[i]) * Float64(c[j])
        sum_a += a
        sum_a2 += a * a
        sum_ad += a * d
    end
    return sum_a, sum_a2, sum_ad
end

function _pattern_moments(c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    n = length(c)
    count = directed ? n * (n - 1) : n * (n - 1) ÷ 2
    if discrete
        if all(value -> value == 0 || value == 1, c)
            n_periphery = Base.count(iszero, c)
            sum_d = Float64(n * (n - 1) ÷ 2 -
                            n_periphery * (n_periphery - 1) ÷ 2)
            directed && (sum_d *= 2)
            return count, sum_d, sum_d
        end
        ordered = sort!(Float64.(c))
        sum_d = 0.0
        sum_d2 = 0.0
        @inbounds for k in eachindex(ordered)
            multiplicity = k - 1
            value = ordered[k]
            sum_d += multiplicity * value
            sum_d2 += multiplicity * value * value
        end
        directed && (sum_d *= 2; sum_d2 *= 2)
        return count, sum_d, sum_d2
    end

    sum_c = 0.0
    sum_c2 = 0.0
    sum_c4 = 0.0
    @inbounds for value in c
        x = Float64(value)
        x2 = x * x
        sum_c += x
        sum_c2 += x2
        sum_c4 += x2 * x2
    end
    sum_d = (sum_c * sum_c - sum_c2) / 2
    sum_d2 = (sum_c2 * sum_c2 - sum_c4) / 2
    directed && (sum_d *= 2; sum_d2 *= 2)
    return count, sum_d, sum_d2
end

@inline function _pearson_from_sums(count::Int, sum_a::Float64, sum_a2::Float64,
    sum_d::Float64, sum_d2::Float64, sum_ad::Float64)
    count >= 2 || return 0.0
    inv_count = inv(Float64(count))
    variance_a = sum_a2 - sum_a * sum_a * inv_count
    variance_d = sum_d2 - sum_d * sum_d * inv_count
    (variance_a > 0.0 && variance_d > 0.0) || return 0.0
    covariance = sum_ad - sum_a * sum_d * inv_count
    return covariance / sqrt(variance_a * variance_d)
end

"""
    core_quality(A, c; discrete=false, directed=false)

Compute the Pearson correlation between adjacency values and the ideal
core-periphery pattern without materializing dyad vectors. Undirected inputs use
the strict upper triangle; `directed=true` uses all ordered off-diagonal dyads.
Returns `0.0` when either pattern has zero variance.
"""
function core_quality(A::AbstractMatrix{<:Real}, c::AbstractVector{<:Real};
    discrete::Bool=false,
    directed::Bool=false)
    n = length(c)
    size(A) == (n, n) ||
        throw(DimensionMismatch("A must be a square matrix with size ($(n), $(n))"))
    count, sum_d, sum_d2 = _pattern_moments(c, discrete, directed)
    sum_a, sum_a2, sum_ad = _adjacency_pattern_product(A, c, discrete, directed)
    return _pearson_from_sums(count, sum_a, sum_a2, sum_d, sum_d2, sum_ad)
end

@inline function _binary_quality(n::Int, n_periphery::Int, periphery_weight::Float64,
    sum_a::Float64, sum_a2::Float64)
    count = n * (n - 1) ÷ 2
    sum_d = Float64(count - n_periphery * (n_periphery - 1) ÷ 2)
    return _pearson_from_sums(
        count, sum_a, sum_a2, sum_d, sum_d, sum_a - periphery_weight)
end

function _links_to_periphery(A::AbstractMatrix{<:Real}, node::Int,
    c::AbstractVector{<:Real})
    total = 0.0
    @inbounds for neighbor in axes(A, 1)
        neighbor == node && continue
        iszero(c[neighbor]) && (total += Float64(A[neighbor, node]))
    end
    return total
end

function _links_to_periphery(A::SparseMatrixCSC{<:Real}, node::Int,
    c::AbstractVector{<:Real})
    total = 0.0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, node)
        neighbor = rows[index]
        neighbor == node && continue
        iszero(c[neighbor]) && (total += Float64(values[index]))
    end
    return total
end

function _axpy_column!(destination::AbstractVector{Float64},
    A::AbstractMatrix{<:Real}, column::Int, scale::Float64)
    @inbounds for row in axes(A, 1)
        destination[row] += scale * Float64(A[row, column])
    end
    return destination
end

function _axpy_column!(destination::AbstractVector{Float64},
    A::SparseMatrixCSC{<:Real}, column::Int, scale::Float64)
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, column)
        destination[rows[index]] += scale * Float64(values[index])
    end
    return destination
end

"""
    borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, step=0.05, init=nothing)

Borgatti-Everett continuous core-periphery model.

Finds coreness vector c that maximizes correlation with ideal pattern Δ[i,j] = c[i]*c[j].

# Arguments
- `A`: Adjacency matrix (n x n)
- `max_iter`: Maximum iterations
- `tol`: Convergence tolerance
- `step`: Coordinate-search resolution in (0, 1]
- `init`: Initial coreness vector (optional)

# Returns
- CPResult with continuous coreness scores

# Reference
Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
"""
function borgatti_everett_continuous(A::AbstractMatrix{<:Real};
    max_iter::Int=1000,
    tol::Float64=1e-6,
    step::Float64=0.05,
    init::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = _validate_adjacency(A; symmetric=true)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    tol >= 0 || throw(ArgumentError("tol must be nonnegative"))
    0 < step <= 1 || throw(ArgumentError("step must be in (0, 1]"))

    if init === nothing
        c = _normalize_scores(vec(sum(A, dims=2)))
    else
        length(init) == n || throw(DimensionMismatch("init must have length $n"))
        all(x -> isfinite(x) && 0 <= x <= 1, init) ||
            throw(ArgumentError("init values must be finite and in [0, 1]"))
        c = Float64.(init)
    end

    # Maintain the sufficient statistics of the Pearson objective. Candidate
    # coordinates are O(1); only accepted moves update A*c (O(n) dense or
    # O(degree) sparse).
    count = n * (n - 1) ÷ 2
    sum_a, sum_a2, sum_ad = _adjacency_pattern_product(A, c, false, false)
    sum_c = sum(c)
    sum_c2 = sum(abs2, c)
    sum_c4 = sum(x -> x^4, c)
    sum_d = (sum_c * sum_c - sum_c2) / 2
    sum_d2 = (sum_c2 * sum_c2 - sum_c4) / 2
    quality = _pearson_from_sums(count, sum_a, sum_a2, sum_d, sum_d2, sum_ad)
    neighbor_scores = Vector{Float64}(A * c)
    grid = collect(0.0:step:1.0)
    last(grid) < 1.0 && push!(grid, 1.0)
    for _ in 1:max_iter
        improved = false
        for i in 1:n
            old_ci = c[i]
            best_ci = old_ci
            best_quality = quality
            for candidate in grid
                candidate == old_ci && continue
                difference = candidate - old_ci
                candidate_sum_c = sum_c + difference
                candidate_sum_c2 = sum_c2 - old_ci^2 + candidate^2
                candidate_sum_c4 = sum_c4 - old_ci^4 + candidate^4
                candidate_sum_d =
                    (candidate_sum_c^2 - candidate_sum_c2) / 2
                candidate_sum_d2 =
                    (candidate_sum_c2^2 - candidate_sum_c4) / 2
                candidate_sum_ad = sum_ad + difference * neighbor_scores[i]
                candidate_quality = _pearson_from_sums(count, sum_a, sum_a2,
                    candidate_sum_d, candidate_sum_d2, candidate_sum_ad)
                if candidate_quality > best_quality + tol
                    best_quality = candidate_quality
                    best_ci = candidate
                end
            end
            if best_ci != old_ci
                difference = best_ci - old_ci
                sum_c += difference
                sum_c2 += best_ci^2 - old_ci^2
                sum_c4 += best_ci^4 - old_ci^4
                sum_ad += difference * neighbor_scores[i]
                _axpy_column!(neighbor_scores, A, i, difference)
                c[i] = best_ci
                improved = true
            end
            quality = best_quality
        end
        improved || break
    end

    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)
    quality = core_quality(A, c)

    return CPResult(c, core_nodes, periphery_nodes, quality, "Borgatti-Everett Continuous")
end

"""
    borgatti_everett_discrete(A; max_iter=1000, init=nothing)

Borgatti-Everett discrete core-periphery model.

Finds binary partition maximizing correlation with ideal discrete pattern.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum iterations for optimization
- `init`: Initial binary assignment (optional)

# Returns
- CPResult with binary coreness (0 or 1)

# Reference
Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
"""
function borgatti_everett_discrete(A::AbstractMatrix{<:Real};
    max_iter::Int=1000,
    init::Union{Nothing,AbstractVector{<:Real}}=nothing)
    n = _validate_adjacency(A; symmetric=true)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n == 1 && return CPResult([1.0], [1], Int[], 0.0, "Borgatti-Everett Discrete")

    if init === nothing
        degrees = vec(sum(A, dims=2))
        order = sortperm(degrees; rev=true, alg=Base.Sort.MergeSort)
        c = zeros(Float64, n)
        c[order[1:clamp(div(n, 2), 1, n - 1)]] .= 1.0
    else
        length(init) == n || throw(DimensionMismatch("init must have length $n"))
        all(x -> x == 0 || x == 1, init) ||
            throw(ArgumentError("init must be binary"))
        c = Float64.(init)
    end

    sum_a, sum_a2, sum_ad = _adjacency_pattern_product(A, c, true, false)
    n_periphery = count(iszero, c)
    periphery_weight = sum_a - sum_ad
    best_quality = _binary_quality(n, n_periphery, periphery_weight, sum_a, sum_a2)

    # Greedy optimization: try swapping each node
    improved = true
    iter = 0

    while improved && iter < max_iter
        improved = false
        iter += 1

        for i in 1:n
            links = _links_to_periphery(A, i, c)
            moving_to_periphery = !iszero(c[i])
            candidate_n_periphery = n_periphery + (moving_to_periphery ? 1 : -1)
            1 <= candidate_n_periphery < n || continue
            candidate_periphery_weight = periphery_weight +
                (moving_to_periphery ? links : -links)
            quality_new = _binary_quality(n, candidate_n_periphery,
                candidate_periphery_weight, sum_a, sum_a2)
            if quality_new > best_quality
                best_quality = quality_new
                c[i] = 1.0 - c[i]
                n_periphery = candidate_n_periphery
                periphery_weight = candidate_periphery_weight
                improved = true
            end
        end
    end

    best_c = c
    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)
    best_quality = core_quality(A, best_c; discrete=true)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Borgatti-Everett Discrete")
end

"""
    lip_discrete(A; max_iter=1000)

Lip's fast algorithm for discrete core-periphery bipartitioning.

For a binary, undirected network, this implements Lip's exact algorithm. Nodes are
ordered by decreasing degree, and the loss (missing core-core edges plus observed
periphery-periphery edges) is evaluated for every proper, nonempty degree-prefix
core. The smallest core is returned when several sizes attain the same minimum.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Deprecated compatibility keyword; ignored because the algorithm is exact

# Returns
- CPResult with binary partition

# Reference
Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
"""
function lip_discrete(A::AbstractMatrix{<:Real}; max_iter::Int=1000)
    n = _validate_adjacency(A; symmetric=true, binary=true)
    n >= 2 || throw(ArgumentError("Lip's bipartition requires at least two nodes"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))

    degrees = vec(sum(A, dims=2))
    sorted_idx = sortperm(degrees; rev=true, alg=Base.Sort.MergeSort)

    # Lip (2011), equations (6)--(7).
    z = sum(degrees) / 2
    best_z = Inf
    best_k = 1
    for k in 1:(n - 1)
        z += k - 1 - degrees[sorted_idx[k]]
        if z < best_z
            best_z = z
            best_k = k
        end
    end

    best_c = zeros(Float64, n)
    best_c[sorted_idx[1:best_k]] .= 1.0

    # Compute final quality as correlation
    quality = core_quality(A, best_c; discrete=true)

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, quality, "Lip Discrete")
end

"""
    rombach_continuous(A; alpha=0.5, beta=0.5, max_iter=1000, tol=1e-6,
                       n_runs=10, rng=Random.default_rng())

Rombach's generalized continuous core-periphery model.

Constructs the fixed Rombach coreness template controlled by `alpha` and
`beta`, then maximizes `sum(A[i,j] * c[i] * c[j] for i < j)` over permutations
of that template using greedy pair swaps and random restarts. `quality` is this
raw Rombach objective, rather than a Pearson correlation.

# Arguments
- `A`: Adjacency matrix
- `alpha`: Controls sharpness of core-periphery boundary
- `beta`: Controls size of core
- `max_iter`: Maximum iterations per run
- `tol`: Convergence tolerance
- `n_runs`: Number of random restarts
- `rng`: Random-number generator used for restart permutations

# Returns
- CPResult with continuous coreness scores

# Reference
Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
"""
function rombach_continuous(A::AbstractMatrix{<:Real};
    alpha::Float64=0.5,
    beta::Float64=0.5,
    max_iter::Int=1000,
    tol::Float64=1e-6,
    n_runs::Int=10,
    rng::AbstractRNG=Random.default_rng())
    n = _validate_adjacency(A; symmetric=true)
    0.0 <= alpha <= 1.0 || throw(ArgumentError("alpha must be in [0, 1]"))
    0.0 < beta < 1.0 || throw(ArgumentError("beta must be in (0, 1)"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    tol >= 0.0 || throw(ArgumentError("tol must be nonnegative"))

    # Fixed parameterized core-vector family, evaluated at ranks i/n.
    # Optimization only permutes these values among vertices.
    template = Vector{Float64}(undef, n)
    for i in 1:n
        x = i / n
        template[i] = if x <= beta
            x * (1.0 - alpha) / (2.0 * beta)
        else
            (x - beta) * (1.0 - alpha) / (2.0 * (1.0 - beta)) +
                (1.0 + alpha) / 2.0
        end
    end

    best_quality = -Inf
    best_c = copy(template)

    for _ in 1:n_runs
        c = template[randperm(rng, n)]

        # h[i] is the weighted coreness of all neighbors of i. The strict
        # upper triangle defines both h and the objective.
        h = zeros(Float64, n)
        quality = 0.0
        for i in 1:(n - 1)
            for j in (i + 1):n
                w = Float64(A[i, j])
                quality += w * c[i] * c[j]
                h[i] += w * c[j]
                h[j] += w * c[i]
            end
        end

        for _ in 1:max_iter
            swap_i = 0
            swap_j = 0
            best_delta = tol

            # The i-j term is unchanged by a swap. Cached neighbor scores make
            # every candidate swap delta constant-time.
            for i in 1:(n - 1)
                a = c[i]
                for j in (i + 1):n
                    b = c[j]
                    a == b && continue
                    w = Float64(A[i, j])
                    delta = (b - a) * ((h[i] - w * b) - (h[j] - w * a))
                    if delta > best_delta
                        best_delta = delta
                        swap_i = i
                        swap_j = j
                    end
                end
            end

            swap_i == 0 && break
            i, j = swap_i, swap_j
            a, b = c[i], c[j]
            difference = b - a
            w_ij = Float64(A[i, j])

            # An accepted swap updates all cached scores in linear time.
            for k in 1:n
                if k != i && k != j
                    w_ki = Float64(A[min(k, i), max(k, i)])
                    w_kj = Float64(A[min(k, j), max(k, j)])
                    h[k] += difference * (w_ki - w_kj)
                end
            end
            h[i] -= difference * w_ij
            h[j] += difference * w_ij
            c[i], c[j] = b, a
            quality += best_delta
        end

        if quality > best_quality
            best_quality = quality
            best_c = copy(c)
        end
    end

    threshold = median(best_c)
    core_nodes = findall(best_c .>= threshold)
    periphery_nodes = findall(best_c .< threshold)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Rombach Continuous")
end

#=
    random_walker_profiling(A; tie_break=:deterministic,
                            rng=Random.default_rng(),
                            stationary_method=:lazy,
                            stationary_distribution=nothing,
                            stationary_tol=1e-12,
                            stationary_max_iter=100_000,
                            n_walks=1000, walk_length=10)

Core-periphery profiling by random-walker persistence probabilities.

This implementation constructs nested peripheral sets, adding at each step the node
that minimizes the set's persistence probability. The persistence value at which a
node is inserted is its coreness. With `tie_break=:deterministic`, ties are resolved by
lower strength and then lower node index. Set `tie_break=:paper` to reproduce the
paper's random selection among the minimum-strength nodes tied for minimum
persistence. The returned `quality` is the paper's core-periphery centralization
(equation 8; undefined, and therefore `NaN`, for two-node graphs).

Inputs may be undirected connected networks or directed strongly connected networks,
with finite, nonnegative weights and no self-loops. Directed profiles use the
stationary distribution of the row-normalized random walk. A lazy power iteration
preserves that distribution while also converging for periodic strongly connected
networks.

# Arguments
- `A`: Adjacency matrix
- `tie_break`: `:deterministic` (default) or the stochastic, paper-exact `:paper`
- `rng`: Random-number generator used when `tie_break=:paper`
- `stationary_method`: Directed solver, `:lazy` (default), `:linear`, or `:auto`
- `stationary_distribution`: Optional validated stationary probability vector
- `stationary_tol`: L1 convergence tolerance for the directed stationary distribution
- `stationary_max_iter`: Maximum stationary iterations for directed inputs
- `n_walks`: Deprecated compatibility keyword; ignored
- `walk_length`: Deprecated compatibility keyword; ignored

# Returns
- CPResult with persistence-profile coreness

# Reference
Rossa, F.D., et al. (2013). Profiling core-periphery network structure by random walkers.
=#
function _positive_reachable(A::AbstractMatrix{<:Real}, start::Int)
    n = size(A, 1)
    seen = falses(n)
    seen[start] = true
    stack = [start]
    while !isempty(stack)
        node = pop!(stack)
        @inbounds for neighbor in axes(A, 1)
            if !seen[neighbor] && A[neighbor, node] > 0
                seen[neighbor] = true
                push!(stack, neighbor)
            end
        end
    end
    return seen
end

function _positive_reachable(A::SparseMatrixCSC{<:Real}, start::Int)
    n = size(A, 1)
    seen = falses(n)
    seen[start] = true
    stack = [start]
    rows = rowvals(A)
    values = nonzeros(A)
    while !isempty(stack)
        node = pop!(stack)
        @inbounds for index in nzrange(A, node)
            neighbor = rows[index]
            if !seen[neighbor] && values[index] > 0
                seen[neighbor] = true
                push!(stack, neighbor)
            end
        end
    end
    return seen
end

function _minimum_strength_node(strengths::AbstractVector{<:Real},
    tie_break::Symbol, rng::AbstractRNG)
    if tie_break === :deterministic
        return argmin(strengths)
    end

    minimum_strength = minimum(strengths)
    selected = 0
    tied = 0
    @inbounds for node in eachindex(strengths)
        if strengths[node] == minimum_strength
            tied += 1
            rand(rng, 1:tied) == 1 && (selected = node)
        end
    end
    return selected
end

function _transition_transpose_mul!(destination::Vector{Float64},
    A::AbstractMatrix{<:Real}, scaled_probability::Vector{Float64})
    @inbounds for column in axes(A, 2)
        total = 0.0
        for row in axes(A, 1)
            total += Float64(A[row, column]) * scaled_probability[row]
        end
        destination[column] = total
    end
    return destination
end

function _transition_transpose_mul!(destination::Vector{Float64},
    A::SparseMatrixCSC{<:Real}, scaled_probability::Vector{Float64})
    fill!(destination, 0.0)
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for column in axes(A, 2), index in nzrange(A, column)
        destination[column] += Float64(values[index]) * scaled_probability[rows[index]]
    end
    return destination
end

"""
    StationaryDistributionResult

Stationary random-walk distribution and numerical convergence diagnostics used by
directed Della Rossa profiling.
"""
struct StationaryDistributionResult
    distribution::Vector{Float64}
    residual::Float64
    iterations::Int
    converged::Bool
    method::Symbol
end

function _stationary_residual!(walked::Vector{Float64},
    scaled_probability::Vector{Float64}, A::AbstractMatrix{<:Real},
    probability::Vector{Float64}, out_strengths::Vector{Float64})
    @inbounds for node in eachindex(probability)
        scaled_probability[node] = probability[node] / out_strengths[node]
    end
    _transition_transpose_mul!(walked, A, scaled_probability)
    return sum(abs(walked[node] - probability[node]) for node in eachindex(probability))
end

function _stationary_lazy(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, tolerance::Float64, max_iter::Int)
    n = size(A, 1)
    probability = fill(inv(Float64(n)), n)
    candidate = similar(probability)
    walked = similar(probability)
    scaled_probability = similar(probability)

    for iteration in 1:max_iter
        @inbounds for node in 1:n
            scaled_probability[node] = probability[node] / out_strengths[node]
        end
        _transition_transpose_mul!(walked, A, scaled_probability)

        total = 0.0
        @inbounds for node in 1:n
            candidate[node] = 0.5 * (probability[node] + walked[node])
            total += candidate[node]
        end
        inverse_total = inv(total)
        distance = 0.0
        @inbounds for node in 1:n
            candidate[node] *= inverse_total
            distance += abs(candidate[node] - probability[node])
        end
        if distance <= 0.5 * tolerance
            residual = _stationary_residual!(
                walked, scaled_probability, A, candidate, out_strengths,
            )
            if residual <= tolerance
                return StationaryDistributionResult(
                    copy(candidate), residual, iteration, true, :lazy,
                )
            end
        end
        probability, candidate = candidate, probability
    end

    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    throw(ErrorException(
        "stationary method :lazy did not converge in $max_iter iterations " *
        "at tolerance $tolerance (final L1 residual $residual); increase " *
        "stationary_max_iter, relax stationary_tol, or use stationary_method=:linear",
    ))
end

function _stationary_linear(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, tolerance::Float64)
    n = size(A, 1)
    transition = Matrix{Float64}(A)
    @inbounds for column in 1:n, row in 1:n
        transition[row, column] /= out_strengths[row]
    end
    system = Matrix{Float64}(I, n, n) - transpose(transition)
    system[end, :] .= 1.0
    right_hand_side = zeros(Float64, n)
    right_hand_side[end] = 1.0
    probability = system \ right_hand_side
    minimum(probability) >= -sqrt(eps(Float64)) || throw(ErrorException(
        "stationary method :linear produced a negative probability"))
    probability .= max.(probability, 0.0)
    probability ./= sum(probability)
    walked = similar(probability)
    scaled_probability = similar(probability)
    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    residual <= tolerance || throw(ErrorException(
        "stationary method :linear failed tolerance $tolerance " *
        "(final L1 residual $residual)"))
    return StationaryDistributionResult(probability, residual, 1, true, :linear)
end

function _validated_stationary(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, probability_input::AbstractVector{<:Real},
    tolerance::Float64)
    n = size(A, 1)
    length(probability_input) == n || throw(DimensionMismatch(
        "stationary_distribution must have length $n"))
    probability = Float64.(probability_input)
    all(isfinite, probability) || throw(ArgumentError(
        "stationary_distribution entries must be finite"))
    all(>=(0.0), probability) || throw(ArgumentError(
        "stationary_distribution entries must be nonnegative"))
    total = sum(probability)
    total > 0.0 || throw(ArgumentError(
        "stationary_distribution must have positive mass"))
    probability ./= total
    walked = similar(probability)
    scaled_probability = similar(probability)
    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    residual <= tolerance || throw(ArgumentError(
        "stationary_distribution has L1 residual $residual, exceeding $tolerance"))
    return StationaryDistributionResult(probability, residual, 0, true, :supplied)
end

function _compute_stationary(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}; method::Symbol, tolerance::Float64,
    max_iter::Int, supplied::Union{Nothing,AbstractVector})
    supplied === nothing || return _validated_stationary(
        A, out_strengths, supplied, tolerance)
    selected = method === :auto ?
        ((A isa SparseMatrixCSC || size(A, 1) > 2_000) ? :lazy : :linear) : method
    selected === :lazy && return _stationary_lazy(
        A, out_strengths, tolerance, max_iter)
    selected === :linear && return _stationary_linear(A, out_strengths, tolerance)
    throw(ArgumentError("stationary_method must be :auto, :lazy, or :linear"))
end

"""
    rossa_stationary_distribution(A; method=:auto, tol=1e-12,
                                  max_iter=100_000)

Compute the stationary distribution used by directed Della Rossa profiling. The
returned `StationaryDistributionResult` includes the L1 balance residual, iteration
count, and selected method. `:linear` is a robust dense constrained solve; `:lazy`
uses stored-edge matrix-vector products and is suitable for large sparse graphs;
`:auto` selects between them. Input must be strongly connected, directed or
undirected, finite, nonnegative, and loop-free.
"""
function rossa_stationary_distribution(A::AbstractMatrix{<:Real};
    method::Symbol=:auto, tol::Real=1e-12, max_iter::Int=100_000)
    n = _validate_adjacency(A)
    n >= 2 || throw(ArgumentError("stationary distribution requires at least two nodes"))
    method in (:auto, :lazy, :linear) ||
        throw(ArgumentError("method must be :auto, :lazy, or :linear"))
    isfinite(tol) && tol >= 0 ||
        throw(ArgumentError("tol must be finite and nonnegative"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be positive"))
    reverse_A = A isa SparseMatrixCSC ? copy(transpose(A)) : transpose(A)
    all(_positive_reachable(A, 1)) && all(_positive_reachable(reverse_A, 1)) ||
        throw(ArgumentError("stationary distribution requires a strongly connected graph"))
    out_strengths = Float64.(vec(sum(A, dims=2)))
    return _compute_stationary(
        A, out_strengths;
        method=method,
        tolerance=Float64(tol),
        max_iter=max_iter,
        supplied=nothing,
    )
end


function _axpy_scaled_column!(destination::Vector{Float64},
    A::AbstractMatrix{<:Real}, column::Int, row_scale::Vector{Float64})
    @inbounds for row in axes(A, 1)
        destination[row] += Float64(A[row, column]) * row_scale[row]
    end
    return destination
end

function _axpy_scaled_column!(destination::Vector{Float64},
    A::SparseMatrixCSC{<:Real}, column::Int, row_scale::Vector{Float64})
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, column)
        row = rows[index]
        destination[row] += Float64(values[index]) * row_scale[row]
    end
    return destination
end

function _select_profile_candidate(in_periphery::BitVector,
    boundary::Vector{Float64}, mass::Float64, internal::Float64,
    node_mass::Vector{Float64}, strengths::Vector{Float64},
    tie_break::Symbol, rng::AbstractRNG)
    best_node = 0
    best_alpha = Inf
    best_strength = Inf
    tied = 0

    @inbounds for node in eachindex(in_periphery)
        in_periphery[node] && continue
        candidate_alpha = (internal + boundary[node]) / (mass + node_mass[node])
        if candidate_alpha < best_alpha ||
           (candidate_alpha == best_alpha && strengths[node] < best_strength)
            best_node = node
            best_alpha = candidate_alpha
            best_strength = strengths[node]
            tied = 1
        elseif candidate_alpha == best_alpha && strengths[node] == best_strength
            if tie_break === :paper
                tied += 1
                rand(rng, 1:tied) == 1 && (best_node = node)
            elseif node < best_node
                best_node = node
            end
        end
    end
    return best_node, best_alpha, boundary[best_node]
end

function _profile_result(c::Vector{Float64}, profile::Vector{Float64})
    n = length(c)
    quality = n == 2 ? NaN : 1 - (2 / (n - 2)) * sum(@view profile[1:(n - 1)])
    quality = isfinite(quality) ? clamp(quality, 0.0, 1.0) : quality
    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)
    return CPResult(c, core_nodes, periphery_nodes, quality, "Random Walker Profiling")
end

"""
    random_walker_profiling(A; tie_break=:deterministic,
                            rng=Random.default_rng(),
                            stationary_method=:lazy,
                            stationary_distribution=nothing,
                            stationary_tol=1e-12,
                            stationary_max_iter=100_000,
                            n_walks=1000, walk_length=10)

Compute the Della Rossa persistence profile of an undirected connected or directed
strongly connected network. The default deterministic mode replaces the paper's random
choices with lower node index; `tie_break=:paper` restores the published stochastic
selection rule using `rng`. The legacy simulation keywords are retained for
compatibility and ignored. The result quality is the paper's cp-centralization.
Directed inputs accept `stationary_method=:lazy`, `:linear`, or `:auto`, or a
validated `stationary_distribution` supplied by the caller.
"""
function random_walker_profiling(A::AbstractMatrix{<:Real};
    tie_break::Symbol=:deterministic,
    rng::AbstractRNG=Random.default_rng(),
    stationary_method::Symbol=:lazy,
    stationary_distribution::Union{Nothing,AbstractVector}=nothing,
    stationary_tol::Real=1e-12,
    stationary_max_iter::Int=100_000,
    n_walks::Int=1000,
    walk_length::Int=10)
    n = _validate_adjacency(A)
    n >= 2 || throw(ArgumentError("random-walker profiling requires at least two nodes"))
    tie_break in (:deterministic, :paper) ||
        throw(ArgumentError("tie_break must be :deterministic or :paper"))
    stationary_method in (:auto, :lazy, :linear) ||
        throw(ArgumentError("stationary_method must be :auto, :lazy, or :linear"))
    isfinite(stationary_tol) && stationary_tol >= 0 ||
        throw(ArgumentError("stationary_tol must be finite and nonnegative"))
    stationary_max_iter >= 1 ||
        throw(ArgumentError("stationary_max_iter must be positive"))
    n_walks >= 0 || throw(ArgumentError("n_walks must be nonnegative"))
    walk_length >= 0 || throw(ArgumentError("walk_length must be nonnegative"))

    symmetric = issymmetric(A)
    symmetric && stationary_distribution !== nothing && throw(ArgumentError(
        "stationary_distribution is only used for directed inputs"))
    reverse_A = symmetric ? A :
        (A isa SparseMatrixCSC ? copy(transpose(A)) : transpose(A))

    # `_positive_reachable` follows positive entries down a column. Checking A and
    # its transpose therefore checks both directed reachability orientations.
    all(_positive_reachable(A, 1)) &&
        (symmetric || all(_positive_reachable(reverse_A, 1))) ||
        throw(ArgumentError(symmetric ?
            "random-walker profiling requires a connected graph" :
            "directed random-walker profiling requires a strongly connected graph"))

    in_periphery = falses(n)
    profile = zeros(Float64, n)
    c = zeros(Float64, n)
    boundary = zeros(Float64, n)

    if symmetric
        strengths = Float64.(vec(sum(A, dims=2)))
        first_node = _minimum_strength_node(strengths, tie_break, rng)
        in_periphery[first_node] = true
        mass = strengths[first_node]
        internal = 0.0
        _axpy_column!(boundary, A, first_node, 2.0)

        for k in 2:n
            best_node, best_alpha, best_boundary = _select_profile_candidate(
                in_periphery, boundary, mass, internal, strengths, strengths,
                tie_break, rng,
            )
            in_periphery[best_node] = true
            internal += best_boundary
            mass += strengths[best_node]
            profile[k] = clamp(best_alpha, 0.0, 1.0)
            c[best_node] = profile[k]
            _axpy_column!(boundary, A, best_node, 2.0)
        end
    else
        out_strengths = Float64.(vec(sum(A, dims=2)))
        in_strengths = Float64.(vec(sum(A, dims=1)))
        strengths = out_strengths .+ in_strengths
        stationary_result = _compute_stationary(
            A, out_strengths;
            method=stationary_method,
            tolerance=Float64(stationary_tol),
            max_iter=stationary_max_iter,
            supplied=stationary_distribution,
        )
        probability = stationary_result.distribution
        flow_scale = probability ./ out_strengths

        first_node = _minimum_strength_node(strengths, tie_break, rng)
        in_periphery[first_node] = true
        mass = probability[first_node]
        internal = 0.0
        _axpy_scaled_column!(boundary, A, first_node, flow_scale)
        _axpy_column!(boundary, reverse_A, first_node, flow_scale[first_node])

        for k in 2:n
            best_node, best_alpha, best_boundary = _select_profile_candidate(
                in_periphery, boundary, mass, internal, probability, strengths,
                tie_break, rng,
            )
            in_periphery[best_node] = true
            internal += best_boundary
            mass += probability[best_node]
            profile[k] = clamp(best_alpha, 0.0, 1.0)
            c[best_node] = profile[k]
            _axpy_scaled_column!(boundary, A, best_node, flow_scale)
            _axpy_column!(boundary, reverse_A, best_node, flow_scale[best_node])
        end
    end

    return _profile_result(c, profile)
end

"""
    RossaEnsembleResult

Uncertainty summary from repeated paper-mode Della Rossa profiles. Optional raw
replicate coreness occupies `O(n * n_runs)` memory; set `store_replicates=false`
to retain only summaries after fitting.
"""
struct RossaEnsembleResult
    mean_coreness::Vector{Float64}
    std_coreness::Vector{Float64}
    lower_coreness::Vector{Float64}
    upper_coreness::Vector{Float64}
    centralizations::Vector{Float64}
    mean_centralization::Float64
    std_centralization::Float64
    rank_stability::Float64
    unique_profiles::Int
    replicate_coreness::Union{Nothing,Matrix{Float64}}
    seeds::Vector{UInt}
end

Base.:(==)(left::RossaEnsembleResult, right::RossaEnsembleResult) =
    left.mean_coreness == right.mean_coreness &&
    left.std_coreness == right.std_coreness &&
    left.lower_coreness == right.lower_coreness &&
    left.upper_coreness == right.upper_coreness &&
    left.centralizations == right.centralizations &&
    left.mean_centralization == right.mean_centralization &&
    left.std_centralization == right.std_centralization &&
    left.rank_stability == right.rank_stability &&
    left.unique_profiles == right.unique_profiles &&
    left.replicate_coreness == right.replicate_coreness &&
    left.seeds == right.seeds

function _average_ranks(values::AbstractVector{<:Real})
    n = length(values)
    order = sortperm(values; alg=Base.Sort.MergeSort)
    ranks = zeros(Float64, n)
    first = 1
    while first <= n
        last = first
        while last < n && values[order[last + 1]] == values[order[first]]
            last += 1
        end
        rank = (first + last) / 2
        @inbounds for position in first:last
            ranks[order[position]] = rank
        end
        first = last + 1
    end
    return ranks
end

function _rank_correlation(left::AbstractVector{<:Real}, right::AbstractVector{<:Real})
    left_ranks = _average_ranks(left)
    right_ranks = _average_ranks(right)
    left_centered = left_ranks .- mean(left_ranks)
    right_centered = right_ranks .- mean(right_ranks)
    denominator = norm(left_centered) * norm(right_centered)
    return iszero(denominator) ? 1.0 : dot(left_centered, right_centered) / denominator
end

"""
    rossa_profile_ensemble(A; n_runs=100, rng=Random.default_rng(),
                           interval=(0.025, 0.975), threaded=false,
                           store_replicates=true,
                           memory_limit_bytes=512 * 1024^2, kwargs...)

Fit repeated `tie_break=:paper` Rossa profiles and summarize node coreness,
cp-centralization, rank stability relative to the ensemble mean, and the number of
unique profiles. Child seeds are drawn before fitting, so fixed seeds produce the
same result in serial and threaded execution. Remaining keywords are forwarded to
`random_walker_profiling`, except `tie_break` and `rng`.

The working profile matrix requires `8 * n * n_runs` bytes even when raw replicates
are omitted from the returned compact result. `memory_limit_bytes` checks that
allocation before fitting so large ensembles fail predictably instead of exhausting
process memory.
"""
function rossa_profile_ensemble(A::AbstractMatrix{<:Real};
    n_runs::Int=100,
    rng::AbstractRNG=Random.default_rng(),
    interval::Tuple{<:Real,<:Real}=(0.025, 0.975),
    threaded::Bool=false,
    store_replicates::Bool=true,
    memory_limit_bytes::Integer=512 * 1024^2,
    kwargs...)
    n = size(A, 1)
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    memory_limit_bytes >= 0 ||
        throw(ArgumentError("memory_limit_bytes must be nonnegative"))
    required_profile_bytes = Base.checked_mul(Base.checked_mul(8, n), n_runs)
    required_profile_bytes <= memory_limit_bytes || throw(ArgumentError(
        "Rossa ensemble profile workspace requires $required_profile_bytes bytes, " *
        "exceeding memory_limit_bytes=$memory_limit_bytes; reduce n_runs or " *
        "raise the explicit limit",
    ))
    lower_probability, upper_probability = Float64.(interval)
    0.0 <= lower_probability <= upper_probability <= 1.0 ||
        throw(ArgumentError("interval probabilities must satisfy 0 <= lower <= upper <= 1"))
    seeds = rand(rng, UInt, n_runs)
    coreness = Matrix{Float64}(undef, n, n_runs)
    centralizations = Vector{Float64}(undef, n_runs)

    function fit_run!(run)
        result = random_walker_profiling(
            A;
            tie_break=:paper,
            rng=MersenneTwister(seeds[run]),
            kwargs...,
        )
        coreness[:, run] = result.coreness
        centralizations[run] = result.quality
        return nothing
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads :static for run in 1:n_runs
            fit_run!(run)
        end
    else
        for run in 1:n_runs
            fit_run!(run)
        end
    end

    mean_coreness = vec(mean(coreness, dims=2))
    std_coreness = vec(std(coreness, dims=2; corrected=false))
    lower_coreness = [quantile(@view(coreness[node, :]), lower_probability) for node in 1:n]
    upper_coreness = [quantile(@view(coreness[node, :]), upper_probability) for node in 1:n]
    rank_stability = mean(
        _rank_correlation(@view(coreness[:, run]), mean_coreness) for run in 1:n_runs)
    unique_profiles = length(Set(Tuple(@view(coreness[:, run])) for run in 1:n_runs))
    raw = store_replicates ? coreness : nothing
    return RossaEnsembleResult(
        mean_coreness,
        std_coreness,
        lower_coreness,
        upper_coreness,
        centralizations,
        mean(centralizations),
        std(centralizations; corrected=false),
        rank_stability,
        unique_profiles,
        raw,
        seeds,
    )
end

"""
    coreness_scores(result::CPResult)

Get coreness scores from a CPResult.

# Returns
- Vector of coreness scores (higher = more core-like)
"""
function coreness_scores(result::CPResult)
    return result.coreness
end

"""
    CPMultiResult

Result structure for multiple core-periphery pairs detection.

# Fields
- `pair_labels::Vector{Int}`: Pair assignment for each node (1, 2, ..., K)
- `coreness::Vector{Float64}`: Binary coreness within each pair
- `n_pairs::Int`: Number of detected core-periphery pairs
- `quality::Float64`: Quality score (Q^cp)
- `algorithm::String`: Name of the algorithm used
- `pair_selection::Symbol`: `:objective` or `:penalized` pair-count rule
- `candidate_pair_counts::Vector{Int}`: Evaluated model sizes
- `candidate_qualities::Vector{Float64}`: Raw best objective at each model size
- `selection_scores::Vector{Float64}`: Scores used to select the returned model
"""
struct CPMultiResult
    pair_labels::Vector{Int}
    coreness::Vector{Float64}
    n_pairs::Int
    quality::Float64
    algorithm::String
    pair_selection::Symbol
    candidate_pair_counts::Vector{Int}
    candidate_qualities::Vector{Float64}
    selection_scores::Vector{Float64}
end

Base.:(==)(left::CPMultiResult, right::CPMultiResult) =
    left.pair_labels == right.pair_labels &&
    left.coreness == right.coreness &&
    left.n_pairs == right.n_pairs &&
    left.quality == right.quality &&
    left.algorithm == right.algorithm &&
    left.pair_selection == right.pair_selection &&
    left.candidate_pair_counts == right.candidate_pair_counts &&
    left.candidate_qualities == right.candidate_qualities &&
    left.selection_scores == right.selection_scores

CPMultiResult(pair_labels::Vector{Int}, coreness::Vector{Float64}, n_pairs::Int,
    quality::Float64, algorithm::String) = CPMultiResult(
    pair_labels, coreness, n_pairs, quality, algorithm,
    :objective, [n_pairs], [quality], [quality],
)

#=
    surprise_cp(A; max_iter=100)

Surprise-based core-periphery detection.

Maximizes the negative log probability of observing at least as many core-core
and core-periphery edges under the multivariate hypergeometric null model.
This method is defined for simple, undirected, binary graphs.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum optimization iterations

# Returns
- CPResult with binary coreness

# Reference
Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise.
=#
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

"""
    surprise_cp(A; max_iter=100)

Greedily maximize the exact multivariate-hypergeometric joint-tail surprise for a
simple, undirected binary graph. The returned quality is the negative log tail
probability.
"""
function surprise_cp(A::AbstractMatrix{<:Real}; max_iter::Int=100)
    n = _validate_adjacency(A; symmetric=true, binary=true)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n == 1 && return CPResult([1.0], [1], Int[], 0.0, "Surprise CP")

    # Initialize based on degree
    degrees = vec(sum(A, dims=2))
    threshold = mean(degrees)
    c = Float64.(degrees .>= threshold)

    total_dyads = n * (n - 1) ÷ 2
    logfactorial = zeros(Float64, total_dyads + 1)
    @inbounds for value in 2:total_dyads
        logfactorial[value + 1] = logfactorial[value] + log(value)
    end
    logchoose(total, selected) = if selected < 0 || selected > total
        -Inf
    else
        logfactorial[total + 1] - logfactorial[selected + 1] -
            logfactorial[total - selected + 1]
    end

    function logaddexp(x, y)
        x == -Inf && return y
        y == -Inf && return x
        hi, lo = max(x, y), min(x, y)
        return hi + log1p(exp(lo - hi))
    end

    log_tables = Vector{Union{Nothing,NTuple{3,Vector{Float64}}}}(undef, n + 1)
    fill!(log_tables, nothing)
    function tables_for(n_core)
        cached = log_tables[n_core + 1]
        cached === nothing || return cached
        n_periphery = n - n_core
        capacities = (n_core * (n_core - 1) ÷ 2,
            n_core * n_periphery,
            n_periphery * (n_periphery - 1) ÷ 2)
        tables = ntuple(block ->
            [logchoose(capacities[block], selected)
             for selected in 0:capacities[block]], 3)
        log_tables[n_core + 1] = tables
        return tables
    end

    total_edges = Int(round(sum(degrees) / 2))
    log_denominator = logchoose(total_dyads, total_edges)

    # Exact joint upper tail for enriched core-core and core-periphery blocks.
    function log_surprise(n_core, observed_cc, observed_cp)
        (total_edges == 0 || total_edges == total_dyads) && return 0.0
        n_periphery = n - n_core
        capacity_cc = n_core * (n_core - 1) ÷ 2
        capacity_cp = n_core * n_periphery
        capacity_pp = n_periphery * (n_periphery - 1) ÷ 2
        log_cc, log_cp, log_pp = tables_for(n_core)

        log_tail = -Inf
        for l_cc in observed_cc:min(capacity_cc, total_edges)
            remaining = total_edges - l_cc
            lower_cp = max(observed_cp, remaining - capacity_pp)
            upper_cp = min(capacity_cp, remaining)
            for l_cp in lower_cp:upper_cp
                l_pp = remaining - l_cp
                term = log_cc[l_cc + 1] + log_cp[l_cp + 1] + log_pp[l_pp + 1]
                log_tail = logaddexp(log_tail, term)
            end
        end
        return max(0.0, log_denominator - log_tail)
    end

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

"""
    label_switching_cp(A; max_iter=100, rng=Random.default_rng(), n_runs=1)

Label-switching algorithm for core-periphery detection.

Uses greedy optimization of the discrete Pearson core quality with random
processing orders and optional random restarts.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum iterations
- `rng`: Random number generator used for node ordering and random restarts
- `n_runs`: Number of starts (the first is degree-based)

# Returns
- CPResult with binary coreness

# Reference
Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection.
"""
function label_switching_cp(A::AbstractMatrix{<:Real};
    max_iter::Int=100,
    rng::AbstractRNG=Random.default_rng(),
    n_runs::Int=1)
    n = _validate_adjacency(A; symmetric=true)
    n >= 3 || throw(ArgumentError("label_switching_cp requires at least three nodes"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))

    # A nonconstant discrete ideal pattern requires at least one core node and
    # two peripheral nodes. Start with the highest-degree nodes in the core.
    degrees = vec(sum(A, dims=2))
    degree_order = sortperm(degrees; rev=true)
    initial_core_size = clamp(div(n, 2), 1, n - 2)
    moment_c = zeros(Float64, n)
    sum_a, sum_a2, _ = _adjacency_pattern_product(A, moment_c, true, false)

    best_quality = -Inf
    best_c = zeros(Float64, n)

    for run in 1:n_runs
        c = zeros(Float64, n)
        if run == 1
            c[degree_order[1:initial_core_size]] .= 1.0
        else
            random_order = randperm(rng, n)
            random_core_size = rand(rng, 1:(n - 2))
            c[random_order[1:random_core_size]] .= 1.0
        end

        _, _, sum_ad = _adjacency_pattern_product(A, c, true, false)
        n_periphery = count(iszero, c)
        periphery_weight = sum_a - sum_ad
        current_quality =
            _binary_quality(n, n_periphery, periphery_weight, sum_a, sum_a2)
        isfinite(current_quality) || continue

        for _ in 1:max_iter
            improved = false

            for i in randperm(rng, n)
                links = _links_to_periphery(A, i, c)
                moving_to_periphery = !iszero(c[i])
                candidate_n_periphery =
                    n_periphery + (moving_to_periphery ? 1 : -1)
                2 <= candidate_n_periphery <= n - 1 || continue
                candidate_periphery_weight = periphery_weight +
                    (moving_to_periphery ? links : -links)
                candidate_quality = _binary_quality(n, candidate_n_periphery,
                    candidate_periphery_weight, sum_a, sum_a2)
                if isfinite(candidate_quality) && candidate_quality > current_quality
                    c[i] = 1.0 - c[i]
                    n_periphery = candidate_n_periphery
                    periphery_weight = candidate_periphery_weight
                    current_quality = candidate_quality
                    improved = true
                end
            end

            improved || break
        end

        if current_quality > best_quality
            best_quality = current_quality
            best_c = copy(c)
        end
    end

    isfinite(best_quality) || throw(ArgumentError(
        "label_switching_cp requires an adjacency matrix with nonzero variance"))

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)
    best_quality = core_quality(A, best_c; discrete=true)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Label Switching CP")
end

include("spectral.jl")
include("multipair.jl")
include("directed.jl")
include("significance.jl")

"""
    spectral_method(A; beta=0.1)

Apply Cucuringu et al.'s LowRank-Core rank-two reconstruction followed by
Find-Cut. The published method requires a simple undirected binary graph.
"""
function spectral_method(A::AbstractMatrix{<:Real}; beta::Real=0.1)
    result = _lowrank_core(A; beta=beta)
    coreness = _normalize_scores(result.coreness)
    return CPResult(
        coreness,
        result.core_nodes,
        result.periphery_nodes,
        result.quality,
        "Spectral Method",
    )
end

"""
    minres_svd(A; max_iter=1000, tol=1e-6)

Compatibility projection of `minres_svd_directed` onto `CPResult`. Use
`minres_svd_directed` to retain distinct out- and in-coreness vectors and the
least-squares residual.
"""
function minres_svd(A::AbstractMatrix{<:Real}; max_iter::Int=1000, tol::Real=1e-6)
    directed = minres_svd_directed(A; max_iter=max_iter, tol=tol)
    return CPResult(
        directed.coreness,
        directed.core_nodes,
        directed.periphery_nodes,
        directed.quality,
        "MINRES/SVD",
    )
end

"""
    multiple_cp_pairs(A; null_model=:er, max_pairs=10, min_pair_size=2,
                      max_iter=100, n_runs=10, pair_selection=:objective,
                      pair_penalty=0.5, rng=Random.default_rng())

Detect multiple core-periphery pairs by jointly switching each node's pair and
binary core/periphery role. `null_model` may be `:er` or `:configuration`.
The compatibility default selects the largest raw objective. Set
`pair_selection=:penalized` to subtract
`pair_penalty * (k - 1) * log(n) / n` from the best `k`-pair objective.
"""
function multiple_cp_pairs(
    A::AbstractMatrix{<:Real};
    null_model::Symbol=:er,
    max_pairs::Int=10,
    min_pair_size::Int=2,
    max_iter::Int=100,
    n_runs::Int=10,
    pair_selection::Symbol=:objective,
    pair_penalty::Real=0.5,
    rng::AbstractRNG=Random.default_rng(),
)
    _validate_adjacency(A; symmetric=true)
    result = _multiple_cp_pairs_joint(
        A;
        null_model=null_model,
        max_pairs=max_pairs,
        min_pair_size=min_pair_size,
        max_iter=max_iter,
        n_runs=n_runs,
        pair_selection=pair_selection,
        pair_penalty=pair_penalty,
        rng=rng,
    )
    return CPMultiResult(
        result.pair_labels,
        result.coreness,
        result.n_pairs,
        result.quality,
        "Multiple CP Pairs",
        result.pair_selection,
        result.candidate_pair_counts,
        result.candidate_qualities,
        result.selection_scores,
    )
end

"""
    multiple_cp_pairs_config(A; kwargs...)

Convenience wrapper for `multiple_cp_pairs(A; null_model=:configuration)`.
"""
multiple_cp_pairs_config(A::AbstractMatrix{<:Real}; kwargs...) =
    multiple_cp_pairs(A; null_model=:configuration, kwargs...)

end # module
