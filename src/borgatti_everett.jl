"""
    borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, step=0.05, init=nothing,
                                directed=false, partition=:median,
                                threshold=nothing, core_size=nothing)

Borgatti-Everett continuous core-periphery model.

Finds coreness vector c that maximizes correlation with ideal pattern Δ[i,j] = c[i]*c[j].

# Arguments
- `A`: Adjacency matrix (n x n)
- `max_iter`: Maximum iterations
- `tol`: Convergence tolerance
- `step`: Coordinate-search resolution in (0, 1]
- `init`: Initial coreness vector (optional)
- `directed`: With `false` (default) `A` must be symmetric and the correlation
  runs over the strict upper triangle. With `true`, asymmetric input is accepted
  and the correlation runs over all off-diagonal ordered dyads, as UCINET does
  for directed sociometric data; symmetric input gives the same fit either way.
- `partition`: Core/periphery split rule applied to the coreness vector:
  `:median` (default; nodes strictly above the median are core), `:none`
  (all nodes periphery), `:threshold` (nodes strictly above `threshold`), or
  `:core_size` (the `core_size` highest-coreness nodes, ties broken by lower index)
- `threshold`: Cutoff used by `partition=:threshold`
- `core_size`: Number of core nodes used by `partition=:core_size`

Isolates (nodes with zero strength) always receive coreness `0` and are never in
`core_nodes`. A warning is emitted (once per session) when `A` has more than one
connected component with at least two nodes.

# Returns
- CPResult with continuous coreness scores

# Reference
Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
"""
function borgatti_everett_continuous(A::AbstractMatrix{<:Real};
    max_iter::Int=1000,
    tol::Real=1e-6,
    step::Real=0.05,
    init::Union{Nothing,AbstractVector{<:Real}}=nothing,
    directed::Bool=false,
    partition::Symbol=:median,
    threshold::Union{Nothing,Real}=nothing,
    core_size::Union{Nothing,Integer}=nothing)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=!directed)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    0 < step <= 1 || throw(ArgumentError("step must be in (0, 1]"))
    tol = Float64(tol)
    step = Float64(step)
    # Symmetric input takes the (identical) undirected path even when
    # `directed=true` is requested.
    directed = directed && !_issymmetric(A)
    _warn_disconnected(A, :borgatti_everett_continuous)
    isolates = _isolates(A)

    if init === nothing
        strengths = Float64.(vec(sum(A, dims=2)))
        directed && (strengths .+= vec(sum(A, dims=1)))
        c = _normalize_scores(strengths)
    else
        length(init) == n || throw(DimensionMismatch("init must have length $n"))
        all(x -> isfinite(x) && 0 <= x <= 1, init) ||
            throw(ArgumentError("init values must be finite and in [0, 1]"))
        c = Float64.(init)
    end

    # Maintain the sufficient statistics of the Pearson objective. Candidate
    # coordinates are O(1); only accepted moves update A*c (O(n) dense or
    # O(degree) sparse). For directed input the pattern c[i]c[j] is symmetric,
    # so the pattern moments double and the neighbor score of node i is
    # (A*c)[i] + (A'*c)[i].
    count = directed ? n * (n - 1) : n * (n - 1) ÷ 2
    scale = directed ? 2.0 : 1.0
    sum_a, sum_a2, sum_ad = _adjacency_pattern_product(A, c, false, directed)
    sum_c = sum(c)
    sum_c2 = sum(abs2, c)
    sum_c4 = sum(x -> x^4, c)
    sum_d = scale * (sum_c * sum_c - sum_c2) / 2
    sum_d2 = scale * (sum_c2 * sum_c2 - sum_c4) / 2
    quality = _pearson_from_sums(count, sum_a, sum_a2, sum_d, sum_d2, sum_ad)
    At = directed ? _transposed_adjacency(A) : A
    neighbor_scores = Vector{Float64}(A * c)
    directed && (neighbor_scores .+= Vector{Float64}(At * c))
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
                    scale * (candidate_sum_c^2 - candidate_sum_c2) / 2
                candidate_sum_d2 =
                    scale * (candidate_sum_c2^2 - candidate_sum_c4) / 2
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
                directed && _axpy_column!(neighbor_scores, At, i, difference)
                c[i] = best_ci
                improved = true
            end
            quality = best_quality
        end
        improved || break
    end

    c[isolates] .= 0.0
    core_nodes, periphery_nodes = _partition(
        c; rule=partition, threshold=threshold, core_size=core_size)
    core_nodes, periphery_nodes = _enforce_isolates!(c, core_nodes, periphery_nodes, isolates)
    quality = core_quality(A, c; directed=directed)

    return CPResult(c, core_nodes, periphery_nodes, quality, "Borgatti-Everett Continuous")
end

"""
    borgatti_everett_discrete(A; max_iter=1000, init=nothing, directed=false,
                              cp_block=:ones, n_runs=1, rng=Random.default_rng(),
                              method=:greedy)

Borgatti-Everett discrete core-periphery model.

Finds a binary partition maximizing the Pearson correlation between `A` and the
ideal discrete pattern by greedy single-node flipping.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum sweeps per start
- `init`: Initial binary assignment (optional; used by the first start)
- `directed`: With `false` (default) `A` must be symmetric and the correlation
  runs over the strict upper triangle. With `true`, asymmetric input is accepted
  and the correlation runs over all off-diagonal ordered dyads (UCINET's
  convention for directed data); symmetric input gives the same fit either way.
- `cp_block`: Ideal pattern for the core-periphery blocks. `:ones` (default) is
  Borgatti-Everett eq. (2), `Δ[i,j] = max(c[i], c[j])`, in which core-periphery
  ties are expected. `:missing` is their eq. (4), the model they recommend and
  UCINET's categorical default: core-periphery dyads are treated as missing
  data, so only the core-core and periphery-periphery blocks are correlated.
  See [`ideal_cp_matrix`](@ref).
- `n_runs`: Number of greedy starts. The first start is degree-based (or
  `init`); the others draw a random core of random size from `rng`. The best
  final correlation is returned. With `n_runs=1` the fit is deterministic.
- `rng`: Random-number generator for the restarts
- `method`: `:greedy` (default) is the single-node flipping heuristic above.
  `:exhaustive` enumerates every core set of size `1:n-1` (Gray-code order,
  `O(2^n · n)` time) and returns a global optimum of the correlation, the
  first one met in enumeration order among ties; it is available for
  `n ≤ 24` and ignores `n_runs`, `rng`, `init`, and `max_iter`. Greedy
  flipping, even with restarts, misses the optimum on a few percent of small
  graphs because the correlation surface has many local optima, so prefer
  `:exhaustive` for small sociograms.

Isolates (nodes with zero strength) always receive coreness `0` and are never in
`core_nodes`. A warning is emitted (once per session) when `A` has more than one
connected component with at least two nodes.

# Returns
- CPResult with binary coreness (0 or 1); `quality` is the correlation under the
  requested `cp_block` and dyad set

# Reference
Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
"""
function borgatti_everett_discrete(A::AbstractMatrix{<:Real};
    max_iter::Int=1000,
    init::Union{Nothing,AbstractVector{<:Real}}=nothing,
    directed::Bool=false,
    cp_block::Symbol=:ones,
    n_runs::Int=1,
    rng::AbstractRNG=Random.default_rng(),
    method::Symbol=:greedy)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=!directed)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    cp_block in (:ones, :missing) ||
        throw(ArgumentError("cp_block must be :ones or :missing"))
    method in (:greedy, :exhaustive) ||
        throw(ArgumentError("method must be :greedy or :exhaustive"))
    method === :exhaustive && n > 24 && throw(ArgumentError(
        "method=:exhaustive enumerates 2^n core sets and supports at most 24 nodes"))
    directed = directed && !_issymmetric(A)
    n == 1 && return CPResult([1.0], [1], Int[], 0.0, "Borgatti-Everett Discrete")
    _warn_disconnected(A, :borgatti_everett_discrete)
    isolates = _isolates(A)

    if method === :exhaustive
        At = directed ? _transposed_adjacency(A) : A
        best_c = _exhaustive_discrete(A, At, directed, cp_block)
        core_nodes = findall(best_c .== 1.0)
        periphery_nodes = findall(best_c .== 0.0)
        core_nodes, periphery_nodes = _enforce_isolates!(
            best_c, core_nodes, periphery_nodes, isolates)
        best_quality = core_quality(
            A, best_c; discrete=true, directed=directed, cp_block=cp_block)
        return CPResult(best_c, core_nodes, periphery_nodes, best_quality,
            "Borgatti-Everett Discrete")
    end

    if init === nothing
        strengths = Float64.(vec(sum(A, dims=2)))
        directed && (strengths .+= vec(sum(A, dims=1)))
        order = sortperm(strengths; rev=true, alg=Base.Sort.MergeSort)
        initial = zeros(Float64, n)
        initial[order[1:clamp(div(n, 2), 1, n - 1)]] .= 1.0
    else
        length(init) == n || throw(DimensionMismatch("init must have length $n"))
        all(x -> x == 0 || x == 1, init) ||
            throw(ArgumentError("init must be binary"))
        initial = Float64.(init)
    end

    At = directed ? _transposed_adjacency(A) : A
    best_quality = -Inf
    best_c = copy(initial)
    for run in 1:n_runs
        c = if run == 1
            copy(initial)
        else
            start = zeros(Float64, n)
            start[randperm(rng, n)[1:rand(rng, 1:(n - 1))]] .= 1.0
            start
        end
        quality = _greedy_discrete_flips!(c, A, At, directed, cp_block, max_iter, 1, n - 1)
        if quality > best_quality
            best_quality = quality
            best_c = c
        end
    end

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)
    core_nodes, periphery_nodes = _enforce_isolates!(
        best_c, core_nodes, periphery_nodes, isolates)
    best_quality = core_quality(A, best_c; discrete=true, directed=directed, cp_block=cp_block)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Borgatti-Everett Discrete")
end

"""
    _greedy_discrete_flips!(c, A, At, directed, cp_block, max_iter,
                            min_core, max_core; order=1:n, rng=nothing)

Greedy first-improvement single-node flipping of the discrete Pearson
objective starting from the binary vector `c` (modified in place). A flip is
accepted when it strictly increases the objective and keeps the core size in
`min_core:max_core`. Nodes are visited in `1:n` order, or in a fresh random
order each sweep when `rng` is given. Returns the final objective value.
"""
function _greedy_discrete_flips!(c::Vector{Float64}, A::AbstractMatrix{<:Real},
    At::AbstractMatrix{<:Real}, directed::Bool, cp_block::Symbol, max_iter::Int,
    min_core::Int, max_core::Int; rng::Union{Nothing,AbstractRNG}=nothing)
    n = length(c)
    state = _discrete_state(A, c, directed, cp_block)
    quality = _discrete_quality(state)
    for _ in 1:max_iter
        improved = false
        order = rng === nothing ? (1:n) : randperm(rng, n)
        for i in order
            to_core_role = iszero(c[i])
            candidate_core = state.n_core + (to_core_role ? 1 : -1)
            min_core <= candidate_core <= max_core || continue
            links = _node_block_links(A, At, i, c)
            candidate_quality = _flip_quality(state, to_core_role, links...)
            if candidate_quality > quality
                _apply_flip!(state, to_core_role, links...)
                c[i] = to_core_role ? 1.0 : 0.0
                quality = candidate_quality
                improved = true
            end
        end
        improved || break
    end
    return quality
end

"""
    _exhaustive_discrete(A, At, directed, cp_block)

Global maximiser of the discrete Pearson objective over every core set of size
`1:n-1`, visited in binary-reflected Gray-code order so that consecutive states
differ by one flip and the sufficient statistics update in `O(degree)`. The
first optimum met in enumeration order is returned among ties.
"""
function _exhaustive_discrete(A::AbstractMatrix{<:Real}, At::AbstractMatrix{<:Real},
    directed::Bool, cp_block::Symbol)
    n = size(A, 1)
    c = zeros(Float64, n)
    state = _discrete_state(A, c, directed, cp_block)
    best_quality = -Inf
    best_c = copy(c)
    for k in 1:(2^n - 1)
        i = trailing_zeros(k) + 1
        to_core_role = iszero(c[i])
        links = _node_block_links(A, At, i, c)
        _apply_flip!(state, to_core_role, links...)
        c[i] = to_core_role ? 1.0 : 0.0
        1 <= state.n_core <= n - 1 || continue
        quality = _discrete_quality(state)
        if quality > best_quality + 1e-12
            best_quality = quality
            copyto!(best_c, c)
        end
    end
    return best_c
end
