"""
    _rombach_boundary(n, beta)

The index `⌊βN⌋` separating the periphery (ranks `1:⌊βN⌋`) from the core
(ranks `⌊βN⌋+1:N`) in Rombach et al. (2017) eq. (4.9). A tolerance of `1e-9`
guards the floor against floating-point error (`0.29 * 100 = 28.999…`).
"""
_rombach_boundary(n::Int, beta::Float64) = clamp(floor(Int, beta * n + 1e-9), 0, n)

"""
    _rombach_template(n, alpha, beta)

Transition-function values `C*_i(α, β) = g_{α,β}(i)` of Rombach et al. (2017)
eq. (4.9), for ranks `i = 1:n`, before normalisation:

    g(i) = i (1 - α) / (2 ⌊βN⌋)                              for i ≤ ⌊βN⌋
    g(i) = (i - ⌊βN⌋)(1 - α) / (2 (N - ⌊βN⌋)) + (1 + α) / 2  for i > ⌊βN⌋

When `β = 0` only the second case applies and when `β = 1` only the first.
The values are nondecreasing in `i` and the last one equals `1`.
"""
function _rombach_template(n::Int, alpha::Float64, beta::Float64)
    boundary = _rombach_boundary(n, beta)
    template = Vector{Float64}(undef, n)
    for i in 1:n
        template[i] = if i <= boundary
            i * (1.0 - alpha) / (2.0 * boundary)
        else
            (i - boundary) * (1.0 - alpha) / (2.0 * (n - boundary)) +
                (1.0 + alpha) / 2.0
        end
    end
    return template
end

"""
    _rombach_normalized_template(n, alpha, beta)

The eq. (4.9) template scaled to sum to `1`, the normalisation the paper
imposes on the core vector. At `alpha = beta = 1` the template is identically
zero (every node is a peripheral node with zero score) and is returned as is.
"""
function _rombach_normalized_template(n::Int, alpha::Float64, beta::Float64)
    template = _rombach_template(n, alpha, beta)
    total = sum(template)
    total > 0.0 && (template ./= total)
    return template
end

"""
    rombach_continuous(A; alpha=0.5, beta=0.5, max_iter=1000, tol=1e-6,
                       n_runs=10, rng=Random.default_rng(),
                       partition=:beta, threshold=nothing, core_size=nothing)

Rombach's generalized continuous core-periphery model at a single point
`(alpha, beta)` of the parameter plane (Rombach et al. 2017, §4.2).

The core vector `C` is a shuffle of the transition-function template
`C*(α, β)` of eq. (4.9), normalised so that its entries sum to `1` as in the
paper, chosen to maximise the core quality of eq. (4.10),
`R_{α,β} = Σ_{i,j} A[i,j] C[i] C[j]` (the full double sum, so each undirected
edge counts twice). The search uses greedy pair swaps and random restarts.
`quality` is `R_{α,β}` of the returned vector, not a Pearson correlation;
`coreness` is `C` (its entries sum to `1`; multiply by `sum(C*)` to recover the
raw template values, which lie in `[0, 1]`). For the parameter-averaged
aggregate core score of eq. (4.11) see [`rombach_core_score`](@ref).

# Arguments
- `A`: Adjacency matrix (symmetric)
- `alpha ∈ [0, 1]`: Sharpness of the core-periphery boundary, i.e. the size of
  the score jump between the highest-scoring peripheral node and the
  lowest-scoring core node (`alpha = 1` is a discrete classification)
- `beta ∈ [0, 1]`: Fraction of nodes in the *periphery*: the first `⌊βN⌋`
  template ranks are peripheral and the remaining `N - ⌊βN⌋` are core, so the
  core shrinks from `N` nodes at `beta = 0` to none at `beta = 1`
- `max_iter`: Maximum swap iterations per run
- `tol`: Minimum objective gain for a swap to be accepted
- `n_runs`: Number of random restarts
- `rng`: Random-number generator used for restart permutations
- `partition`: Core/periphery split rule applied to the coreness vector:
  `:beta` (default; the `N - ⌊βN⌋` nodes holding the core ranks of the
  template, i.e. the core implied by `beta`), `:median` (nodes strictly above
  the median are core), `:none` (all nodes periphery), `:threshold` (nodes
  strictly above `threshold`), or `:core_size` (the `core_size`
  highest-coreness nodes, ties broken by lower index)
- `threshold`: Cutoff used by `partition=:threshold`
- `core_size`: Number of core nodes used by `partition=:core_size`

Isolates (nodes with zero strength) are removed before fitting: the template is
built for the `N'` non-isolated nodes, isolates receive coreness `0`, and they
are never in `core_nodes`. A warning is emitted (once per session) when `A`
has more than one connected component with at least two nodes.

# Returns
- CPResult with continuous coreness scores

# Reference
Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
"""
function rombach_continuous(A::AbstractMatrix{<:Real};
    alpha::Real=0.5,
    beta::Real=0.5,
    max_iter::Int=1000,
    tol::Real=1e-6,
    n_runs::Int=10,
    rng::AbstractRNG=Random.default_rng(),
    partition::Symbol=:beta,
    threshold::Union{Nothing,Real}=nothing,
    core_size::Union{Nothing,Integer}=nothing)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true)
    0 <= alpha <= 1 || throw(ArgumentError("alpha must be in [0, 1]"))
    0 <= beta <= 1 || throw(ArgumentError("beta must be in [0, 1]"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    partition in (:beta, :median, :none, :threshold, :core_size) ||
        throw(ArgumentError("partition must be :beta, :median, :none, :threshold, or :core_size"))
    alpha = Float64(alpha)
    beta = Float64(beta)
    tol = Float64(tol)
    _warn_disconnected(A, :rombach_continuous)

    isolates = _isolates(A)
    keep = setdiff(1:n, isolates)
    n_kept = length(keep)
    coreness = zeros(Float64, n)
    quality = 0.0
    if n_kept > 0
        sub = n_kept == n ? A : A[keep, keep]
        template = _rombach_normalized_template(n_kept, alpha, beta)
        c, quality = _rombach_optimize(sub, template, max_iter, tol, n_runs, rng)
        coreness[keep] = c
    end

    if partition === :beta
        threshold === nothing && core_size === nothing || throw(ArgumentError(
            "threshold and core_size are not used with partition=:beta"))
        core_nodes, periphery_nodes = _partition(
            coreness; rule=:core_size, core_size=n_kept - _rombach_boundary(n_kept, beta))
    else
        core_nodes, periphery_nodes = _partition(
            coreness; rule=partition, threshold=threshold, core_size=core_size)
    end
    core_nodes, periphery_nodes = _enforce_isolates!(
        coreness, core_nodes, periphery_nodes, isolates)

    return CPResult(coreness, core_nodes, periphery_nodes, quality, "Rombach Continuous")
end

"""
    _rombach_optimize(A, template, max_iter, tol, n_runs, rng)

Maximise `Σ_{i,j} A[i,j] c[i] c[j]` over shuffles `c` of `template` by greedy
best-improvement pair swaps from `n_runs` random permutations. Returns the
best `(c, objective)`; the objective is the full double sum of eq. (4.10).
"""
function _rombach_optimize(A::AbstractMatrix{<:Real}, template::Vector{Float64},
    max_iter::Int, tol::Float64, n_runs::Int, rng::AbstractRNG)
    n = length(template)
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

    # The swap loop accumulates the strict upper triangle; eq. (4.10) sums
    # over all ordered pairs, which for symmetric A counts each edge twice.
    return best_c, 2.0 * best_quality
end

"""
    rombach_core_score(A; alphas=0.1:0.1:1.0, betas=0.1:0.1:1.0, max_iter=1000,
                       tol=1e-6, n_runs=3, rng=Random.default_rng(),
                       partition=:median, threshold=nothing, core_size=nothing)

Aggregate core score of Rombach et al. (2017) eq. (4.11),

    CS(i) = Z Σ_γ C_i(γ) R_γ,    γ = (α, β) ∈ alphas × betas,

where `C(γ)` is the optimal normalised shuffle of the eq. (4.9) template found
by [`rombach_continuous`](@ref) at `γ`, `R_γ` is its core quality (eq. 4.10),
and `Z` is chosen so that `max_k CS(k) = 1`. This is the paper's recommended
deliverable: a continuous centrality-like coreness that averages over the
parameter plane instead of fixing one `(α, β)`. The paper averages over the
grid `α, β ∈ 0.01:0.01:1` (10,000 points); the default `10 × 10` grid keeps
the cost at 100 optimisations and can be refined through `alphas`/`betas`.

# Arguments
- `A`: Adjacency matrix (symmetric)
- `alphas`, `betas`: Grid values in `[0, 1]` (see [`rombach_continuous`](@ref)
  for their meaning; `beta` is the periphery fraction)
- `max_iter`, `tol`, `n_runs`, `rng`: Optimiser settings passed to each
  `(α, β)` fit, in grid order (`alphas` outer, `betas` inner) so results are
  reproducible for a seeded `rng`
- `partition`, `threshold`, `core_size`: Split rule applied to the aggregate
  score, as in [`borgatti_everett_continuous`](@ref); default `:median`

Isolates receive score `0` and are never in `core_nodes`. The result's
`quality` is the mean core quality `R_γ` over the grid; `algorithm` is
`"Rombach aggregate core score"`.

# Reference
Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
"""
function rombach_core_score(A::AbstractMatrix{<:Real};
    alphas::AbstractVector{<:Real}=0.1:0.1:1.0,
    betas::AbstractVector{<:Real}=0.1:0.1:1.0,
    max_iter::Int=1000,
    tol::Real=1e-6,
    n_runs::Int=3,
    rng::AbstractRNG=Random.default_rng(),
    partition::Symbol=:median,
    threshold::Union{Nothing,Real}=nothing,
    core_size::Union{Nothing,Integer}=nothing)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true)
    isempty(alphas) && throw(ArgumentError("alphas must not be empty"))
    isempty(betas) && throw(ArgumentError("betas must not be empty"))
    all(a -> 0 <= a <= 1, alphas) || throw(ArgumentError("alphas must lie in [0, 1]"))
    all(b -> 0 <= b <= 1, betas) || throw(ArgumentError("betas must lie in [0, 1]"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    tol = Float64(tol)
    _warn_disconnected(A, :rombach_core_score)

    isolates = _isolates(A)
    keep = setdiff(1:n, isolates)
    n_kept = length(keep)
    score = zeros(Float64, n)
    total_quality = 0.0
    if n_kept > 0
        sub = n_kept == n ? A : A[keep, keep]
        aggregate = zeros(Float64, n_kept)
        for alpha in alphas, beta in betas
            template = _rombach_normalized_template(n_kept, Float64(alpha), Float64(beta))
            c, quality = _rombach_optimize(sub, template, max_iter, tol, n_runs, rng)
            aggregate .+= c .* quality
            total_quality += quality
        end
        largest = maximum(aggregate)
        largest > 0.0 && (aggregate ./= largest)
        score[keep] = aggregate
    end
    mean_quality = total_quality / (length(alphas) * length(betas))

    core_nodes, periphery_nodes = _partition(
        score; rule=partition, threshold=threshold, core_size=core_size)
    core_nodes, periphery_nodes = _enforce_isolates!(
        score, core_nodes, periphery_nodes, isolates)

    return CPResult(score, core_nodes, periphery_nodes, mean_quality,
        "Rombach aggregate core score")
end
