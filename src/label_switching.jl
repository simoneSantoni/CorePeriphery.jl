"""
    label_switching_cp(A; max_iter=100, rng=Random.default_rng(), n_runs=1,
                       directed=false, cp_block=:ones)

Label-switching algorithm for core-periphery detection.

Uses greedy optimization of the discrete Pearson core quality with random
processing orders and optional random restarts. It optimises the same objective
as [`borgatti_everett_discrete`](@ref) and is the stronger optimiser of the
two: the random sweep order and restarts escape more local optima.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum iterations
- `rng`: Random number generator used for node ordering and random restarts
- `n_runs`: Number of starts (the first is degree-based)
- `directed`: With `false` (default) `A` must be symmetric and the correlation
  runs over the strict upper triangle. With `true`, asymmetric input is accepted
  and the correlation runs over all off-diagonal ordered dyads; symmetric input
  gives the same fit either way.
- `cp_block`: `:ones` (default, Borgatti-Everett eq. 2) or `:missing`
  (Borgatti-Everett eq. 4, core-periphery dyads excluded from the correlation);
  see [`borgatti_everett_discrete`](@ref) and [`ideal_cp_matrix`](@ref).

Isolates (nodes with zero strength) always receive coreness `0` and are never in
`core_nodes`. A warning is emitted (once per session) when `A` has more than one
connected component with at least two nodes.

# Returns
- CPResult with binary coreness

# Reference
Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection.
"""
function label_switching_cp(A::AbstractMatrix{<:Real};
    max_iter::Int=100,
    rng::AbstractRNG=Random.default_rng(),
    n_runs::Int=1,
    directed::Bool=false,
    cp_block::Symbol=:ones)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=!directed)
    n >= 3 || throw(ArgumentError("label_switching_cp requires at least three nodes"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    cp_block in (:ones, :missing) ||
        throw(ArgumentError("cp_block must be :ones or :missing"))
    directed = directed && !_issymmetric(A)
    _warn_disconnected(A, :label_switching_cp)
    isolates = _isolates(A)

    # A nonconstant discrete ideal pattern requires at least one core node and
    # two peripheral nodes. Start with the highest-degree nodes in the core.
    strengths = Float64.(vec(sum(A, dims=2)))
    directed && (strengths .+= vec(sum(A, dims=1)))
    degree_order = sortperm(strengths; rev=true)
    initial_core_size = clamp(div(n, 2), 1, n - 2)
    At = directed ? _transposed_adjacency(A) : A

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

        current_quality = _greedy_discrete_flips!(
            c, A, At, directed, cp_block, max_iter, 1, n - 2; rng=rng)
        isfinite(current_quality) || continue

        if current_quality > best_quality
            best_quality = current_quality
            best_c = c
        end
    end

    isfinite(best_quality) || throw(ArgumentError(
        "label_switching_cp requires an adjacency matrix with nonzero variance"))

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)
    core_nodes, periphery_nodes = _enforce_isolates!(
        best_c, core_nodes, periphery_nodes, isolates)
    best_quality = core_quality(A, best_c; discrete=true, directed=directed, cp_block=cp_block)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Label Switching CP")
end
