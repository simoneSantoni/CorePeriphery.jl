struct _JointGraph
    neighbors::Vector{Vector{Int}}
    weights::Vector{Vector{Float64}}
    degrees::Vector{Float64}
    total_weight::Float64
end

abstract type _JointNullContext end

struct _JointERContext <: _JointNullContext
    normalization::Float64
    density::Float64
end

struct _JointConfigurationContext <: _JointNullContext
    normalization::Float64
    inv_double_weight::Float64
    degrees::Vector{Float64}
end

mutable struct _JointPairAggregates
    counts::Vector{Int}
    core_counts::Vector{Int}
    degree_sums::Vector{Float64}
    core_degree_sums::Vector{Float64}
    degree_square_sums::Vector{Float64}
    core_degree_square_sums::Vector{Float64}
end

"""
    _multiple_cp_pairs_joint(A; null_model=:er, max_pairs=10,
        min_pair_size=2, max_iter=100, n_runs=10,
        rng=Random.default_rng())

Internal joint label-switching engine for the Kojaku--Masuda multi-pair
objective. A state assigns every node both a pair label and a binary
core/periphery label. For each node, a sweep jointly considers core and
periphery membership in its current pair and all pairs represented among its
neighbors. Only strictly improving moves are accepted.

The objective is

`sum((A[i,j] - P[i,j]) * max(x[i], x[j]) * (label[i] == label[j]) for i < j) / m`,

where `P` is either an Erdős--Rényi density null (`:er`, Kojaku & Masuda 2017
eq. 3 divided by `m`) or the configuration null `degree[i] * degree[j] / (2m)`
(`:configuration`). The configuration objective additionally carries the
`i = j` null term of Kojaku & Masuda (2018) eq. (13),
`- sum(degree[i]^2 for i in core) / (4 m^2)`, so that it equals eq. (13)
exactly (and reduces to modularity when every node is core). Input must be a finite,
nonnegative, symmetric adjacency matrix with a zero diagonal; the objective is
then evaluated once per dyad over its upper triangle. Pair-size feasibility is
preserved throughout each run, and final labels are consecutive. The returned
named tuple is compatible with the fields of `CPMultiResult` and additionally
contains `objective_history`.
"""
function _multiple_cp_pairs_joint(A::AbstractMatrix{<:Real};
    null_model::Symbol=:er,
    max_pairs::Int=10,
    min_pair_size::Int=2,
    max_iter::Int=100,
    n_runs::Int=10,
    pair_selection::Symbol=:objective,
    pair_penalty::Real=0.5,
    rng::AbstractRNG=Random.default_rng())
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true)
    null_model in (:er, :configuration) ||
        throw(ArgumentError("null_model must be :er or :configuration"))
    max_pairs >= 1 || throw(ArgumentError("max_pairs must be positive"))
    min_pair_size >= 1 || throw(ArgumentError("min_pair_size must be positive"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    pair_selection in (:objective, :penalized) ||
        throw(ArgumentError("pair_selection must be :objective or :penalized"))
    isfinite(pair_penalty) && pair_penalty >= 0 ||
        throw(ArgumentError("pair_penalty must be finite and nonnegative"))

    min_pair_size <= n ||
        throw(ArgumentError("min_pair_size cannot exceed the number of nodes"))

    graph = _joint_build_graph(A)
    degrees = graph.degrees
    total_weight = graph.total_weight
    max_feasible_pairs = min(max_pairs, fld(n, min_pair_size))
    max_feasible_pairs = max(max_feasible_pairs, 1)
    normalization = total_weight > 0.0 ? total_weight : 1.0
    density = n > 1 ? total_weight / (n * (n - 1) / 2) : 0.0

    context = if null_model === :er
        _JointERContext(normalization, density)
    else
        inv_double_weight = total_weight > 0.0 ? 1.0 / (2.0 * total_weight) : 0.0
        _JointConfigurationContext(normalization, inv_double_weight, degrees)
    end
    # Strict comparison avoids the all-core trap when a large peripheral block
    # shares the median degree. An all-periphery regular graph receives one
    # deterministic seed core so role moves remain feasible.
    initial_core = Float64.(degrees .> _joint_median(degrees))
    all(iszero, initial_core) && (initial_core[argmax(degrees)] = 1.0)

    candidate_qualities = fill(-Inf, max_feasible_pairs)
    candidate_labels = [Int[] for _ in 1:max_feasible_pairs]
    candidate_core = [Float64[] for _ in 1:max_feasible_pairs]
    candidate_histories = [Float64[] for _ in 1:max_feasible_pairs]

    for run in 1:n_runs
        # Cycle through feasible pair counts so enough restarts compare simpler
        # and richer models. Randomness controls the balanced assignment.
        n_initial_pairs = 1 + mod(run - 1, max_feasible_pairs)
        order = randperm(rng, n)
        labels = Vector{Int}(undef, n)
        for (position, node) in enumerate(order)
            labels[node] = 1 + mod(position - 1, n_initial_pairs)
        end
        core = copy(initial_core)
        aggregates = _joint_pair_aggregates(labels, core, degrees, n_initial_pairs)
        quality = _joint_objective(graph, labels, core, context, aggregates)
        history = Float64[quality]
        observed_pair = zeros(Float64, n_initial_pairs)
        observed_core = zeros(Float64, n_initial_pairs)
        candidate_mask = falses(n_initial_pairs)
        candidate_pairs = Vector{Int}(undef, n_initial_pairs)

        for _ in 1:max_iter
            improved = false
            for node in randperm(rng, n)
                old_pair = labels[node]
                old_core = core[node]
                n_candidate_pairs = 1
                candidate_pairs[1] = old_pair
                candidate_mask[old_pair] = true
                for index in eachindex(graph.neighbors[node])
                    neighbor = graph.neighbors[node][index]
                    weight = graph.weights[node][index]
                    pair = labels[neighbor]
                    observed_pair[pair] += weight
                    if core[neighbor] == 1.0
                        observed_core[pair] += weight
                    end
                    if !candidate_mask[pair]
                        n_candidate_pairs += 1
                        candidate_pairs[n_candidate_pairs] = pair
                        candidate_mask[pair] = true
                    end
                end
                sort!(@view candidate_pairs[1:n_candidate_pairs])
                old_contribution = _joint_pair_contribution(
                    context, node, old_pair, old_core, old_pair, old_core,
                    aggregates, observed_pair, observed_core,
                )

                best_delta = 0.0
                best_pair = old_pair
                best_core_value = old_core
                for candidate_index in 1:n_candidate_pairs
                    candidate_pair = candidate_pairs[candidate_index]
                    if candidate_pair != old_pair &&
                            aggregates.counts[old_pair] <= min_pair_size
                        continue
                    end
                    for candidate_core in (0.0, 1.0)
                        candidate_pair == old_pair && candidate_core == old_core && continue
                        new_contribution = _joint_pair_contribution(
                            context, node, candidate_pair, candidate_core,
                            old_pair, old_core, aggregates, observed_pair, observed_core,
                        )
                        delta = new_contribution - old_contribution
                        if delta > best_delta + 1e-12
                            best_delta = delta
                            best_pair = candidate_pair
                            best_core_value = candidate_core
                        end
                    end
                end
                for candidate_index in 1:n_candidate_pairs
                    pair = candidate_pairs[candidate_index]
                    candidate_mask[pair] = false
                    observed_pair[pair] = 0.0
                    observed_core[pair] = 0.0
                end

                if best_delta > 0.0
                    _joint_update_aggregates!(
                        aggregates, node, old_pair, old_core,
                        best_pair, best_core_value, degrees[node],
                    )
                    labels[node] = best_pair
                    core[node] = best_core_value
                    quality += best_delta
                    push!(history, quality)
                    improved = true
                end
            end
            improved || break
        end

        labels = _joint_relabel(labels)
        aggregates = _joint_pair_aggregates(labels, core, degrees, n_initial_pairs)
        quality = _joint_objective(graph, labels, core, context, aggregates)
        history[end] = quality
        pair_count = length(unique(labels))
        if quality > candidate_qualities[pair_count]
            candidate_qualities[pair_count] = quality
            candidate_labels[pair_count] = copy(labels)
            candidate_core[pair_count] = copy(core)
            candidate_histories[pair_count] = copy(history)
        end
    end

    evaluated_pairs = findall(isfinite, candidate_qualities)
    selection_scores = copy(candidate_qualities)
    if pair_selection === :penalized
        penalty_unit = Float64(pair_penalty) * log(max(n, 2)) / n
        for pair_count in evaluated_pairs
            selection_scores[pair_count] -= penalty_unit * (pair_count - 1)
        end
    end
    best_pair_count = evaluated_pairs[1]
    for pair_count in @view evaluated_pairs[2:end]
        if selection_scores[pair_count] > selection_scores[best_pair_count]
            best_pair_count = pair_count
        end
    end
    best_labels = candidate_labels[best_pair_count]
    best_core = candidate_core[best_pair_count]
    best_quality = candidate_qualities[best_pair_count]
    best_history = candidate_histories[best_pair_count]

    return (
        pair_labels=best_labels,
        coreness=best_core,
        n_pairs=length(unique(best_labels)),
        quality=best_quality,
        algorithm="Multiple CP Pairs (Joint Label Switching)",
        objective_history=best_history,
        null_model=null_model,
        pair_selection=pair_selection,
        candidate_pair_counts=evaluated_pairs,
        candidate_qualities=candidate_qualities[evaluated_pairs],
        selection_scores=selection_scores[evaluated_pairs],
    )
end

function _joint_build_graph(A::AbstractMatrix{<:Real})
    n = size(A, 1)
    neighbors = [Int[] for _ in 1:n]
    weights = [Float64[] for _ in 1:n]
    degrees = zeros(Float64, n)
    total_weight = 0.0
    for i in 1:(n - 1), j in (i + 1):n
        weight = Float64(A[i, j])
        isfinite(weight) || throw(ArgumentError("A must contain only finite weights"))
        weight >= 0.0 || throw(ArgumentError("A must contain nonnegative weights"))
        if !iszero(weight)
            push!(neighbors[i], j)
            push!(neighbors[j], i)
            push!(weights[i], weight)
            push!(weights[j], weight)
            degrees[i] += weight
            degrees[j] += weight
            total_weight += weight
        end
    end
    return _JointGraph(neighbors, weights, degrees, total_weight)
end

function _joint_build_graph(A::SparseMatrixCSC{<:Real,<:Integer})
    n = size(A, 1)
    neighbors = [Int[] for _ in 1:n]
    weights = [Float64[] for _ in 1:n]
    degrees = zeros(Float64, n)
    total_weight = 0.0
    rows, columns, values = findnz(A)
    for index in eachindex(values)
        i, j = rows[index], columns[index]
        i < j || continue
        weight = Float64(values[index])
        isfinite(weight) || throw(ArgumentError("A must contain only finite weights"))
        weight >= 0.0 || throw(ArgumentError("A must contain nonnegative weights"))
        iszero(weight) && continue
        push!(neighbors[i], j)
        push!(neighbors[j], i)
        push!(weights[i], weight)
        push!(weights[j], weight)
        degrees[i] += weight
        degrees[j] += weight
        total_weight += weight
    end
    return _JointGraph(neighbors, weights, degrees, total_weight)
end

function _joint_pair_aggregates(labels, core, degrees, n_pairs)
    aggregates = _JointPairAggregates(
        zeros(Int, n_pairs), zeros(Int, n_pairs),
        zeros(Float64, n_pairs), zeros(Float64, n_pairs),
        zeros(Float64, n_pairs), zeros(Float64, n_pairs),
    )
    for node in eachindex(labels)
        pair = labels[node]
        degree = degrees[node]
        degree_square = degree * degree
        aggregates.counts[pair] += 1
        aggregates.degree_sums[pair] += degree
        aggregates.degree_square_sums[pair] += degree_square
        if core[node] == 1.0
            aggregates.core_counts[pair] += 1
            aggregates.core_degree_sums[pair] += degree
            aggregates.core_degree_square_sums[pair] += degree_square
        end
    end
    return aggregates
end

function _joint_pair_contribution(
    context::_JointNullContext, node, pair, role, old_pair, old_role,
    aggregates, observed_pair, observed_core,
)
    same_pair = pair == old_pair
    if role == 1.0
        observed = observed_pair[pair]
        count = aggregates.counts[pair] - same_pair
        degree_sum = aggregates.degree_sums[pair]
        if same_pair
            degree_sum -= _joint_node_degree(context, node)
        end
        return _joint_residual_sum(context, node, observed, count, degree_sum) +
            _joint_diagonal_term(context, node)
    else
        observed = observed_core[pair]
        remove_old_core = same_pair && old_role == 1.0
        count = aggregates.core_counts[pair] - remove_old_core
        degree_sum = aggregates.core_degree_sums[pair]
        if remove_old_core
            degree_sum -= _joint_node_degree(context, node)
        end
    end
    return _joint_residual_sum(context, node, observed, count, degree_sum)
end

_joint_node_degree(context::_JointERContext, node) = 0.0
_joint_node_degree(context::_JointConfigurationContext, node) = context.degrees[node]

function _joint_residual_sum(context::_JointERContext, node, observed, count, degree_sum)
    return (observed - context.density * count) / context.normalization
end

function _joint_residual_sum(
    context::_JointConfigurationContext, node, observed, count, degree_sum,
)
    expected = context.degrees[node] * context.inv_double_weight * degree_sum
    return (observed - expected) / context.normalization
end

# Kojaku & Masuda (2018) eq. (13) sums over all (i, j) including i = j. With a
# zero diagonal the i = j null term -d_i^2 / (2M) survives for every core node
# (x_i + x_i - x_i^2 = x_i), and in the per-dyad normalisation used here
# (1 / M per unordered dyad, i.e. 1 / (2M) per ordered pair) it contributes
# -d_i^2 / (4M^2) to the objective. It depends only on the node's own role, so
# it is part of the node's contribution in every move delta. The ER objective
# (Kojaku & Masuda 2017 eq. 3) is strictly lower-triangular and has no such term.
_joint_diagonal_term(::_JointERContext, node) = 0.0
function _joint_diagonal_term(context::_JointConfigurationContext, node)
    degree = context.degrees[node]
    return -degree * degree * context.inv_double_weight / 2.0 / context.normalization
end

function _joint_update_aggregates!(
    aggregates, node, old_pair, old_core, new_pair, new_core, degree,
)
    degree_square = degree * degree
    if old_pair != new_pair
        aggregates.counts[old_pair] -= 1
        aggregates.degree_sums[old_pair] -= degree
        aggregates.degree_square_sums[old_pair] -= degree_square
        aggregates.counts[new_pair] += 1
        aggregates.degree_sums[new_pair] += degree
        aggregates.degree_square_sums[new_pair] += degree_square
        if old_core == 1.0
            aggregates.core_counts[old_pair] -= 1
            aggregates.core_degree_sums[old_pair] -= degree
            aggregates.core_degree_square_sums[old_pair] -= degree_square
        end
        if new_core == 1.0
            aggregates.core_counts[new_pair] += 1
            aggregates.core_degree_sums[new_pair] += degree
            aggregates.core_degree_square_sums[new_pair] += degree_square
        end
    elseif old_core != new_core
        direction = new_core == 1.0 ? 1 : -1
        aggregates.core_counts[old_pair] += direction
        aggregates.core_degree_sums[old_pair] += direction * degree
        aggregates.core_degree_square_sums[old_pair] += direction * degree_square
    end
    return aggregates
end

function _joint_objective(graph::_JointGraph, labels, core, context, aggregates)
    observed = 0.0
    for node in eachindex(labels)
        for index in eachindex(graph.neighbors[node])
            neighbor = graph.neighbors[node][index]
            neighbor > node || continue
            if labels[node] == labels[neighbor]
                observed += graph.weights[node][index] * max(core[node], core[neighbor])
            end
        end
    end
    return (observed - _joint_expected(context, aggregates)) / context.normalization
end

function _joint_expected(context::_JointERContext, aggregates)
    eligible_dyads = 0
    for pair in eachindex(aggregates.counts)
        count = aggregates.counts[pair]
        core_count = aggregates.core_counts[pair]
        eligible_dyads += core_count * (core_count - 1) ÷ 2
        eligible_dyads += core_count * (count - core_count)
    end
    return context.density * eligible_dyads
end

function _joint_expected(context::_JointConfigurationContext, aggregates)
    expected = 0.0
    for pair in eachindex(aggregates.counts)
        degree_sum = aggregates.degree_sums[pair]
        degree_square_sum = aggregates.degree_square_sums[pair]
        periphery_degree_sum = degree_sum - aggregates.core_degree_sums[pair]
        periphery_degree_square_sum =
            degree_square_sum - aggregates.core_degree_square_sums[pair]
        eligible_degree_products = (
            degree_sum * degree_sum - degree_square_sum -
            periphery_degree_sum * periphery_degree_sum + periphery_degree_square_sum
        ) / 2.0
        # Off-diagonal eligible dyads plus the i = j null term of every core
        # node (eq. 13 of Kojaku & Masuda 2018; see `_joint_diagonal_term`).
        expected += context.inv_double_weight *
            (eligible_degree_products + aggregates.core_degree_square_sums[pair] / 2.0)
    end
    return expected
end

function _joint_relabel(labels)
    mapping = Dict{Int,Int}()
    next_label = 0
    return [get!(mapping, label) do
        next_label += 1
        next_label
    end for label in labels]
end

function _joint_median(values)
    ordered = sort(values)
    n = length(ordered)
    isodd(n) ? ordered[(n + 1) ÷ 2] : (ordered[n ÷ 2] + ordered[n ÷ 2 + 1]) / 2
end

"""
    multiple_cp_pairs(A; null_model=:er, max_pairs=10, min_pair_size=2,
                      max_iter=100, n_runs=10, pair_selection=:penalized,
                      pair_penalty=0.5, rng=Random.default_rng())

Detect multiple core-periphery pairs (Kojaku & Masuda 2017, 2018) by jointly
switching each node's pair and binary core/periphery role.

# Keyword arguments
- `null_model`: `:er` (default) uses the Erdős–Rényi null of Kojaku & Masuda
  (2017) eq. (3), `p = m / binom(n, 2)`; `:configuration` uses the
  degree-preserving null of Kojaku & Masuda (2018) eq. (13),
  `P[i,j] = d[i] d[j] / (2m)`, implemented exactly (including the paper's
  `i = j` null term). **Kojaku & Masuda (2018) recommend the configuration
  null**: under the ER null any network with a heterogeneous degree
  distribution exhibits core-periphery pairs, because high-degree nodes are
  "core" relative to a uniform density even when the structure is fully
  explained by the degrees; the configuration null only rewards
  core-periphery structure beyond what the degree sequence implies. `:er`
  remains the default for continuity with 0.1 and with the 2017 paper, so pass
  `null_model=:configuration` (or use [`multiple_cp_pairs_config`](@ref)) for
  substantive analyses.
- `max_pairs`: Largest number of pairs considered; any value from `1` to `n`
  is accepted (the effective cap is `n ÷ min_pair_size`).
- `min_pair_size`: Smallest pair allowed.
- `max_iter`, `n_runs`: Sweeps per run and number of random restarts; restarts
  cycle through `1:max_pairs` initial pair counts.
- `pair_selection`: Rule choosing the number of pairs among the best `k`-pair
  solutions found. `:penalized` (default) maximises the objective minus
  `pair_penalty * (k - 1) * log(n) / n`; `:objective` maximises the raw
  objective, which is nondecreasing in `k` in practice and therefore tends to
  over-partition (and to vary with the seed) on a single planted pair.
- `pair_penalty`: Penalty strength for `:penalized`.
- `rng`: Random-number generator.

Isolates (nodes with zero strength) are removed before fitting and reported
with pair label `0` and coreness `0`; `n_pairs` counts the pairs of the
non-isolated nodes, whose labels are `1:n_pairs`. The null densities are those
of the non-isolated subnetwork. A warning is emitted (once per session) when
`A` has more than one connected component with at least two nodes.
"""
function multiple_cp_pairs(
    A::AbstractMatrix{<:Real};
    null_model::Symbol=:er,
    max_pairs::Int=10,
    min_pair_size::Int=2,
    max_iter::Int=100,
    n_runs::Int=10,
    pair_selection::Symbol=:penalized,
    pair_penalty::Real=0.5,
    rng::AbstractRNG=Random.default_rng(),
)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true)
    _warn_disconnected(A, :multiple_cp_pairs)
    isolates = _isolates(A)
    keep = setdiff(1:n, isolates)
    if isempty(keep)
        # Nothing to fit: every node is an isolate.
        return CPMultiResult(
            zeros(Int, n), zeros(Float64, n), 0, 0.0, "Multiple CP Pairs",
            pair_selection, Int[], Float64[], Float64[],
        )
    end
    sub = length(keep) == n ? A : A[keep, keep]
    result = _multiple_cp_pairs_joint(
        sub;
        null_model=null_model,
        max_pairs=max_pairs,
        min_pair_size=min_pair_size,
        max_iter=max_iter,
        n_runs=n_runs,
        pair_selection=pair_selection,
        pair_penalty=pair_penalty,
        rng=rng,
    )
    pair_labels = zeros(Int, n)
    pair_labels[keep] = result.pair_labels
    coreness = zeros(Float64, n)
    coreness[keep] = result.coreness
    return CPMultiResult(
        pair_labels,
        coreness,
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

Convenience wrapper for `multiple_cp_pairs(A; null_model=:configuration)`, the
degree-preserving null that Kojaku & Masuda (2018) recommend. All other
keywords are forwarded to [`multiple_cp_pairs`](@ref); passing `null_model`
raises an `ArgumentError` because the wrapper fixes it.
"""
function multiple_cp_pairs_config(A::AbstractMatrix{<:Real}; kwargs...)
    haskey(kwargs, :null_model) && throw(ArgumentError(
        "multiple_cp_pairs_config fixes null_model=:configuration; " *
        "call multiple_cp_pairs to choose another null"))
    return multiple_cp_pairs(A; kwargs..., null_model=:configuration)
end
