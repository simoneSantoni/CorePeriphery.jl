"""
    CorePeriphery

Implementation of core-periphery detection algorithms for networks.

This module implements the KM-config algorithm from:
Kojaku, S., & Masuda, N. (2018). Core-periphery structure requires something else
in the network. New Journal of Physics, 20(4), 043012.

The algorithm detects multiple core-periphery pairs using the configuration model
as the null model, controlling for node degree effects.
"""
module CorePeriphery

using Graphs
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using StatsBase
using Distributions
using KernelDensity

export KMConfig, detect_core_periphery, quality, pair_quality
export CorePeripheryResult, is_significant

# =============================================================================
# Data Structures
# =============================================================================

"""
    CorePeripheryPair

Represents a single core-periphery pair detected in the network.

# Fields
- `core_nodes::Vector{Int}`: Indices of core nodes
- `periphery_nodes::Vector{Int}`: Indices of peripheral nodes
- `quality::Float64`: Quality score q of this pair
- `is_significant::Bool`: Whether this pair is statistically significant
- `p_value::Float64`: P-value from statistical test
"""
struct CorePeripheryPair
    core_nodes::Vector{Int}
    periphery_nodes::Vector{Int}
    quality::Float64
    is_significant::Bool
    p_value::Float64
end

"""
    CorePeripheryResult

Complete result from core-periphery detection.

# Fields
- `pairs::Vector{CorePeripheryPair}`: Detected core-periphery pairs
- `node_assignments::Vector{Int}`: Core-periphery pair index for each node (0 = residual)
- `node_roles::Vector{Int}`: Role of each node (1 = core, 0 = periphery)
- `total_quality::Float64`: Total quality Q_config^cp
- `residual_nodes::Vector{Int}`: Nodes not in any significant pair
"""
struct CorePeripheryResult
    pairs::Vector{CorePeripheryPair}
    node_assignments::Vector{Int}
    node_roles::Vector{Int}
    total_quality::Float64
    residual_nodes::Vector{Int}
end

"""
    KMConfig

Parameters for the KM-config algorithm.

# Fields
- `n_runs::Int`: Number of optimization runs (default: 10)
- `n_random_networks::Int`: Number of randomized networks for significance testing (default: 500)
- `significance_level::Float64`: Target significance level α' (default: 0.05)
- `seed::Union{Int,Nothing}`: Random seed for reproducibility
"""
struct KMConfig
    n_runs::Int
    n_random_networks::Int
    significance_level::Float64
    seed::Union{Int,Nothing}

    function KMConfig(;
        n_runs::Int = 10,
        n_random_networks::Int = 500,
        significance_level::Float64 = 0.05,
        seed::Union{Int,Nothing} = nothing
    )
        new(n_runs, n_random_networks, significance_level, seed)
    end
end

# =============================================================================
# Core Algorithm
# =============================================================================

"""
    detect_core_periphery(g::AbstractGraph; config::KMConfig = KMConfig())

Detect core-periphery pairs in an undirected network using the KM-config algorithm.

# Arguments
- `g::AbstractGraph`: Input graph (must be undirected, no self-loops)
- `config::KMConfig`: Algorithm parameters

# Returns
- `CorePeripheryResult`: Complete detection results including pairs and statistics

# Algorithm
Uses a label switching heuristic to maximize Q_config^cp, which measures the
similarity between the network and an idealized core-periphery structure,
relative to the configuration model.

# Reference
Kojaku & Masuda (2018), Eq. 13
"""
function detect_core_periphery(g::AbstractGraph; config::KMConfig = KMConfig())
    if config.seed !== nothing
        Random.seed!(config.seed)
    end

    N = nv(g)
    M = ne(g)

    if M == 0
        return CorePeripheryResult(
            CorePeripheryPair[],
            zeros(Int, N),
            zeros(Int, N),
            0.0,
            collect(1:N)
        )
    end

    # Get adjacency and degrees
    A = adjacency_matrix(g)
    degrees = degree(g)

    # Run optimization multiple times
    best_c = zeros(Int, N)
    best_x = zeros(Int, N)
    best_Q = -Inf

    for run in 1:config.n_runs
        c, x, Q = _optimize_labels(A, degrees, N, M)
        if Q > best_Q
            best_Q = Q
            best_c = copy(c)
            best_x = copy(x)
        end
    end

    # Relabel pairs consecutively
    best_c, n_pairs = _relabel_consecutive(best_c)

    # Calculate quality for each pair
    pair_qualities = zeros(n_pairs)
    pair_sizes = zeros(Int, n_pairs)
    for i in 1:N
        if best_c[i] > 0
            pair_sizes[best_c[i]] += 1
        end
    end

    for pair_idx in 1:n_pairs
        pair_qualities[pair_idx] = _calculate_pair_quality(
            A, degrees, best_c, best_x, pair_idx, N, M
        )
    end

    # Statistical significance testing
    is_significant_vec, p_values = _statistical_test(
        g, A, degrees, best_c, best_x, pair_qualities, pair_sizes, n_pairs, config
    )

    # Build result
    pairs = CorePeripheryPair[]
    for pair_idx in 1:n_pairs
        core_nodes = Int[]
        periphery_nodes = Int[]
        for i in 1:N
            if best_c[i] == pair_idx
                if best_x[i] == 1
                    push!(core_nodes, i)
                else
                    push!(periphery_nodes, i)
                end
            end
        end

        # Check if bipartite-like (fewer intra-core edges than expected)
        is_bipartite_like = _check_bipartite_like(A, degrees, core_nodes, M)

        push!(pairs, CorePeripheryPair(
            core_nodes,
            periphery_nodes,
            pair_qualities[pair_idx],
            is_significant_vec[pair_idx] && !is_bipartite_like,
            p_values[pair_idx]
        ))
    end

    # Identify residual nodes
    residual_nodes = Int[]
    final_assignments = copy(best_c)
    for i in 1:N
        if best_c[i] == 0 || (best_c[i] > 0 && !pairs[best_c[i]].is_significant)
            push!(residual_nodes, i)
            if best_c[i] > 0 && !pairs[best_c[i]].is_significant
                final_assignments[i] = 0
            end
        end
    end

    return CorePeripheryResult(
        pairs,
        final_assignments,
        best_x,
        best_Q,
        residual_nodes
    )
end

"""
    _optimize_labels(A, degrees, N, M)

Run label switching optimization to find core-periphery structure.

# Reference
Kojaku & Masuda (2018), Section 3.1.3
"""
function _optimize_labels(A::AbstractMatrix, degrees::Vector{Int}, N::Int, M::Int)
    # Initialize: each node in its own pair as core
    c = collect(1:N)
    x = ones(Int, N)

    # Precompute degree sums for each block
    # Block (c, x) has nodes with pair index c and role x
    D_cx = Dict{Tuple{Int,Int}, Float64}()
    for i in 1:N
        key = (c[i], x[i])
        D_cx[key] = get(D_cx, key, 0.0) + degrees[i]
    end

    # Precompute edge counts between node and blocks
    # d_tilde[i][(c,x)] = number of edges from node i to block (c,x)
    d_tilde = [Dict{Tuple{Int,Int}, Int}() for _ in 1:N]
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in 1:N
        for k in nzrange(A, i)
            j = rows[k]
            if i != j  # Skip self-loops
                key = (c[j], x[j])
                d_tilde[i][key] = get(d_tilde[i], key, 0) + 1
            end
        end
    end

    # Iterate until convergence
    changed = true
    max_iterations = 100 * N
    iteration = 0

    while changed && iteration < max_iterations
        changed = false
        iteration += 1

        # Random order of nodes
        node_order = randperm(N)

        for i in node_order
            old_c = c[i]
            old_x = x[i]
            old_key = (old_c, old_x)

            # Find neighboring pairs
            neighbor_pairs = Set{Int}()
            for k in nzrange(A, i)
                j = rows[k]
                if i != j
                    push!(neighbor_pairs, c[j])
                end
            end

            if isempty(neighbor_pairs)
                continue
            end

            best_delta = 0.0
            best_new_c = old_c
            best_new_x = old_x

            # Try each neighboring pair
            for new_c in neighbor_pairs
                for new_x in [0, 1]
                    if new_c == old_c && new_x == old_x
                        continue
                    end

                    new_key = (new_c, new_x)

                    # Calculate delta Q using Equation 15
                    delta = _calculate_delta_Q(
                        i, old_c, old_x, new_c, new_x,
                        degrees, D_cx, d_tilde, M
                    )

                    if delta > best_delta
                        best_delta = delta
                        best_new_c = new_c
                        best_new_x = new_x
                    end
                end
            end

            # Update if improvement found
            if best_delta > 1e-10
                changed = true

                # Update D_cx
                D_cx[old_key] = get(D_cx, old_key, 0.0) - degrees[i]
                if D_cx[old_key] ≈ 0.0
                    delete!(D_cx, old_key)
                end
                new_key = (best_new_c, best_new_x)
                D_cx[new_key] = get(D_cx, new_key, 0.0) + degrees[i]

                # Update d_tilde for neighbors
                for k in nzrange(A, i)
                    j = rows[k]
                    if i != j
                        # Remove from old block count
                        d_tilde[j][old_key] = get(d_tilde[j], old_key, 0) - 1
                        if d_tilde[j][old_key] == 0
                            delete!(d_tilde[j], old_key)
                        end
                        # Add to new block count
                        d_tilde[j][new_key] = get(d_tilde[j], new_key, 0) + 1
                    end
                end

                c[i] = best_new_c
                x[i] = best_new_x
            end
        end
    end

    # Calculate final Q
    Q = _calculate_Q(A, degrees, c, x, N, M)

    return c, x, Q
end

"""
    _calculate_delta_Q(i, old_c, old_x, new_c, new_x, degrees, D_cx, d_tilde, M)

Calculate the change in Q when moving node i from (old_c, old_x) to (new_c, new_x).

# Reference
Kojaku & Masuda (2018), Equation 15
"""
function _calculate_delta_Q(
    i::Int, old_c::Int, old_x::Int, new_c::Int, new_x::Int,
    degrees::Vector{Int}, D_cx::Dict{Tuple{Int,Int}, Float64},
    d_tilde::Vector{Dict{Tuple{Int,Int}, Int}}, M::Int
)
    d_i = degrees[i]

    # Get block degree sums
    D_new_c1 = get(D_cx, (new_c, 1), 0.0)
    D_new_c0 = get(D_cx, (new_c, 0), 0.0)
    D_old_c1 = get(D_cx, (old_c, 1), 0.0)
    D_old_c0 = get(D_cx, (old_c, 0), 0.0)

    # Get edge counts from i to blocks
    d_i_new_c1 = get(d_tilde[i], (new_c, 1), 0)
    d_i_new_c0 = get(d_tilde[i], (new_c, 0), 0)
    d_i_old_c1 = get(d_tilde[i], (old_c, 1), 0)
    d_i_old_c0 = get(d_tilde[i], (old_c, 0), 0)

    # First term: contribution from new assignment
    term1 = (1.0 / M) * (
        d_i_new_c1 + new_x * d_i_new_c0 -
        d_i * (D_new_c1 + new_x * D_new_c0) / (2 * M) -
        (d_i^2 / (4 * M)) * (new_x^2 - 2 * (new_x + new_x - new_x * new_x))
    )

    # Simplified using paper's formula (Eq. 15)
    delta = (1.0 / M) * (
        d_i_new_c1 + new_x * d_i_new_c0 -
        d_i * (D_new_c1 + new_x * D_new_c0) / (2 * M) +
        (d_i^2 / (4 * M)) * new_x
    ) - (1.0 / M) * (
        d_i_old_c1 + old_x * d_i_old_c0 -
        d_i * (D_old_c1 + old_x * D_old_c0) / (2 * M) +
        (d_i^2 / (4 * M)) * old_x
    )

    return delta
end

"""
    _calculate_Q(A, degrees, c, x, N, M)

Calculate the total quality Q_config^cp.

# Reference
Kojaku & Masuda (2018), Equation 13
"""
function _calculate_Q(A::AbstractMatrix, degrees::Vector{Int},
                      c::Vector{Int}, x::Vector{Int}, N::Int, M::Int)
    Q = 0.0
    two_M = 2 * M

    rows = rowvals(A)
    vals = nonzeros(A)

    for i in 1:N
        for k in nzrange(A, i)
            j = rows[k]
            if j >= i  # Count each edge once
                if c[i] == c[j]  # Same pair
                    A_ij = vals[k]
                    null_term = (degrees[i] * degrees[j]) / two_M
                    cp_term = x[i] + x[j] - x[i] * x[j]
                    Q += (A_ij - null_term) * cp_term
                end
            end
        end
    end

    return Q / M
end

"""
    _calculate_pair_quality(A, degrees, c, x, pair_idx, N, M)

Calculate the quality q of a specific core-periphery pair.

# Reference
Kojaku & Masuda (2018), Equation 16
"""
function _calculate_pair_quality(A::AbstractMatrix, degrees::Vector{Int},
                                  c::Vector{Int}, x::Vector{Int},
                                  pair_idx::Int, N::Int, M::Int)
    q = 0.0
    two_M = 2 * M

    rows = rowvals(A)
    vals = nonzeros(A)

    for i in 1:N
        if c[i] != pair_idx
            continue
        end
        for k in nzrange(A, i)
            j = rows[k]
            if j >= i && c[j] == pair_idx
                A_ij = vals[k]
                null_term = (degrees[i] * degrees[j]) / two_M
                cp_term = x[i] + x[j] - x[i] * x[j]
                q += (A_ij - null_term) * cp_term
            end
        end
    end

    return q / M
end

"""
    _relabel_consecutive(c)

Relabel pair indices to be consecutive starting from 1.
"""
function _relabel_consecutive(c::Vector{Int})
    unique_pairs = sort(unique(c))
    mapping = Dict{Int, Int}()
    new_idx = 0
    for old_idx in unique_pairs
        if old_idx > 0
            new_idx += 1
            mapping[old_idx] = new_idx
        end
    end

    new_c = [c[i] > 0 ? mapping[c[i]] : 0 for i in eachindex(c)]
    return new_c, new_idx
end

"""
    _check_bipartite_like(A, degrees, core_nodes, M)

Check if a pair is bipartite-like (fewer intra-core edges than expected).
"""
function _check_bipartite_like(A::AbstractMatrix, degrees::Vector{Int},
                                core_nodes::Vector{Int}, M::Int)
    if isempty(core_nodes)
        return false
    end

    # Count actual intra-core edges
    actual_edges = 0
    for i in core_nodes
        for j in core_nodes
            if j > i
                actual_edges += A[i, j]
            end
        end
    end

    # Expected under configuration model
    D_core = sum(degrees[i] for i in core_nodes)
    expected = D_core^2 / (4 * M)

    return actual_edges < expected
end

# =============================================================================
# Statistical Testing
# =============================================================================

"""
    _statistical_test(g, A, degrees, c, x, pair_qualities, pair_sizes, n_pairs, config)

Perform statistical significance testing using configuration model randomization.

# Reference
Kojaku & Masuda (2018), Section 3.1.4
"""
function _statistical_test(g::AbstractGraph, A::AbstractMatrix, degrees::Vector{Int},
                           c::Vector{Int}, x::Vector{Int},
                           pair_qualities::Vector{Float64}, pair_sizes::Vector{Int},
                           n_pairs::Int, config::KMConfig)

    N = nv(g)
    M = ne(g)

    if n_pairs == 0
        return Bool[], Float64[]
    end

    # Generate samples from randomized networks
    random_qualities = Float64[]
    random_sizes = Int[]

    for _ in 1:config.n_random_networks
        # Generate configuration model random graph
        g_rand = _configuration_model_sample(degrees)

        if ne(g_rand) > 0
            A_rand = adjacency_matrix(g_rand)
            degrees_rand = degree(g_rand)

            # Detect core-periphery in randomized network
            c_rand, x_rand, _ = _optimize_labels(A_rand, degrees_rand, N, ne(g_rand))
            c_rand, n_pairs_rand = _relabel_consecutive(c_rand)

            # Record qualities and sizes
            for pair_idx in 1:n_pairs_rand
                q = _calculate_pair_quality(A_rand, degrees_rand, c_rand, x_rand,
                                           pair_idx, N, ne(g_rand))
                n = count(i -> c_rand[i] == pair_idx, 1:N)
                push!(random_qualities, q)
                push!(random_sizes, n)
            end
        end
    end

    # Šidák correction for multiple comparisons
    alpha_corrected = 1 - (1 - config.significance_level)^(1/n_pairs)

    # Test each pair
    is_significant = zeros(Bool, n_pairs)
    p_values = ones(n_pairs)

    if isempty(random_qualities)
        return is_significant, p_values
    end

    # Use kernel density estimation
    for pair_idx in 1:n_pairs
        q = pair_qualities[pair_idx]
        n = pair_sizes[pair_idx]

        p_value = _estimate_p_value(q, n, random_qualities, random_sizes)
        p_values[pair_idx] = p_value
        is_significant[pair_idx] = p_value < alpha_corrected
    end

    return is_significant, p_values
end

"""
    _estimate_p_value(q, n, random_qualities, random_sizes)

Estimate p-value using Gaussian kernel density estimation.

# Reference
Kojaku & Masuda (2018), Appendix C, Equation C8
"""
function _estimate_p_value(q::Float64, n::Int,
                           random_qualities::Vector{Float64},
                           random_sizes::Vector{Int})
    S = length(random_qualities)
    if S == 0
        return 1.0
    end

    # Calculate statistics
    σ_q = std(random_qualities)
    σ_n = std(Float64.(random_sizes))

    if σ_q ≈ 0 || σ_n ≈ 0
        # Fallback: simple proportion
        count_larger = count(i -> random_qualities[i] >= q, 1:S)
        return count_larger / S
    end

    # Correlation coefficient
    mean_q = mean(random_qualities)
    mean_n = mean(Float64.(random_sizes))
    γ = cov(random_qualities, Float64.(random_sizes)) / (σ_q * σ_n)
    γ = clamp(γ, -0.99, 0.99)  # Avoid numerical issues

    # Bandwidth (Scott's rule)
    h = S^(-1/6)

    # Calculate p-value using Eq. C8
    numerator = 0.0
    denominator = 0.0

    for s in 1:S
        q_s = random_qualities[s]
        n_s = random_sizes[s]

        exp_term = exp(-(n - n_s)^2 / (2 * σ_n^2 * h^2))

        # Φ argument
        phi_arg = (σ_n * (q - q_s) - γ * σ_q * (n - n_s)) /
                  (σ_n * σ_q * h * sqrt(1 - γ^2))

        # Normal CDF
        phi_val = cdf(Normal(), phi_arg)

        numerator += exp_term * phi_val
        denominator += exp_term
    end

    if denominator ≈ 0
        return 1.0
    end

    # P(q̂ >= q | n) = 1 - P(q̂ < q | n)
    p_value = 1.0 - numerator / denominator

    return clamp(p_value, 0.0, 1.0)
end

"""
    _configuration_model_sample(degrees)

Generate a random graph from the configuration model.
"""
function _configuration_model_sample(degrees::Vector{Int})
    N = length(degrees)

    # Create stubs
    stubs = Int[]
    for i in 1:N
        for _ in 1:degrees[i]
            push!(stubs, i)
        end
    end

    # Shuffle and pair
    shuffle!(stubs)

    # Create graph
    g = SimpleGraph(N)

    for k in 1:2:length(stubs)-1
        i = stubs[k]
        j = stubs[k+1]
        if i != j  # No self-loops
            add_edge!(g, i, j)  # Multi-edges become single edges
        end
    end

    return g
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    quality(result::CorePeripheryResult)

Get the total quality Q_config^cp of the detection result.
"""
quality(result::CorePeripheryResult) = result.total_quality

"""
    pair_quality(pair::CorePeripheryPair)

Get the quality q of a specific core-periphery pair.
"""
pair_quality(pair::CorePeripheryPair) = pair.quality

"""
    is_significant(pair::CorePeripheryPair)

Check if a core-periphery pair is statistically significant.
"""
is_significant(pair::CorePeripheryPair) = pair.is_significant

"""
    n_significant_pairs(result::CorePeripheryResult)

Count the number of significant core-periphery pairs.
"""
n_significant_pairs(result::CorePeripheryResult) =
    count(p -> p.is_significant, result.pairs)

"""
    get_core_nodes(result::CorePeripheryResult, pair_idx::Int)

Get core nodes of a specific pair (1-indexed).
"""
function get_core_nodes(result::CorePeripheryResult, pair_idx::Int)
    if pair_idx < 1 || pair_idx > length(result.pairs)
        return Int[]
    end
    return result.pairs[pair_idx].core_nodes
end

"""
    get_periphery_nodes(result::CorePeripheryResult, pair_idx::Int)

Get periphery nodes of a specific pair (1-indexed).
"""
function get_periphery_nodes(result::CorePeripheryResult, pair_idx::Int)
    if pair_idx < 1 || pair_idx > length(result.pairs)
        return Int[]
    end
    return result.pairs[pair_idx].periphery_nodes
end

"""
    summarize(result::CorePeripheryResult)

Print a summary of the core-periphery detection results.
"""
function summarize(result::CorePeripheryResult)
    println("Core-Periphery Detection Results (KM-config)")
    println("=" ^ 50)
    println("Total quality Q: $(round(result.total_quality, digits=4))")
    println("Number of pairs detected: $(length(result.pairs))")
    println("Number of significant pairs: $(n_significant_pairs(result))")
    println("Number of residual nodes: $(length(result.residual_nodes))")
    println()

    for (i, pair) in enumerate(result.pairs)
        status = pair.is_significant ? "SIGNIFICANT" : "not significant"
        println("Pair $i ($status):")
        println("  - Core nodes: $(length(pair.core_nodes))")
        println("  - Periphery nodes: $(length(pair.periphery_nodes))")
        println("  - Quality q: $(round(pair.quality, digits=4))")
        println("  - p-value: $(round(pair.p_value, digits=4))")
        println()
    end
end

end # module
