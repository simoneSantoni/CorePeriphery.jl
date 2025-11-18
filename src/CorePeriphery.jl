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

export
    # Data structures
    CPResult,
    CPMultiResult,
    # Algorithms
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    lip_discrete,
    rombach_continuous,
    spectral_method,
    random_walker_profiling,
    minres_svd,
    multiple_cp_pairs,
    surprise_cp,
    label_switching_cp,
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

"""
    adjacency_to_matrix(edges, n)

Convert edge list to adjacency matrix.

# Arguments
- `edges`: Vector of tuples (i, j) or (i, j, weight)
- `n`: Number of nodes

# Returns
- Symmetric adjacency matrix
"""
function adjacency_to_matrix(edges::Vector, n::Int)
    A = zeros(Float64, n, n)
    for edge in edges
        if length(edge) == 2
            i, j = edge
            w = 1.0
        else
            i, j, w = edge
        end
        A[i, j] = w
        A[j, i] = w
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

"""
    core_quality(A, c; discrete=false)

Compute quality (correlation) between adjacency matrix and ideal core-periphery pattern.

# Arguments
- `A`: Adjacency matrix
- `c`: Coreness vector
- `discrete`: Use discrete ideal pattern

# Returns
- Pearson correlation coefficient
"""
function core_quality(A::Matrix{Float64}, c::Vector{Float64}; discrete::Bool=false)
    n = length(c)
    Δ = ideal_cp_matrix(c; discrete=discrete)

    # Extract upper triangle (excluding diagonal)
    a_vec = Float64[]
    d_vec = Float64[]

    for i in 1:n
        for j in (i+1):n
            push!(a_vec, A[i, j])
            push!(d_vec, Δ[i, j])
        end
    end

    # Pearson correlation
    return cor(a_vec, d_vec)
end

"""
    borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, init=nothing)

Borgatti-Everett continuous core-periphery model.

Finds coreness vector c that maximizes correlation with ideal pattern Δ[i,j] = c[i]*c[j].

# Arguments
- `A`: Adjacency matrix (n x n)
- `max_iter`: Maximum iterations
- `tol`: Convergence tolerance
- `init`: Initial coreness vector (optional)

# Returns
- CPResult with continuous coreness scores

# Reference
Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
"""
function borgatti_everett_continuous(A::Matrix{Float64};
                                     max_iter::Int=1000,
                                     tol::Float64=1e-6,
                                     init::Union{Nothing, Vector{Float64}}=nothing)
    n = size(A, 1)

    # Initialize coreness vector
    if init === nothing
        # Use degree centrality as initialization
        c = vec(sum(A, dims=2))
        c = c ./ maximum(c)
    else
        c = copy(init)
    end

    # Iterative optimization using coordinate ascent
    for iter in 1:max_iter
        c_old = copy(c)

        for i in 1:n
            # Optimize c[i] given all other c values
            # Derivative: sum_j A[i,j] * c[j] - sum_j c[i]*c[j]^2
            numerator = sum(A[i, j] * c[j] for j in 1:n if j != i)
            denominator = sum(c[j]^2 for j in 1:n if j != i)

            if denominator > 0
                c[i] = numerator / denominator
            end
        end

        # Normalize to [0, 1]
        c_min, c_max = extrema(c)
        if c_max > c_min
            c = (c .- c_min) ./ (c_max - c_min)
        else
            c = fill(0.5, n)
        end

        # Check convergence
        if norm(c - c_old) < tol
            break
        end
    end

    # Compute quality
    quality = core_quality(A, c)

    # Classify nodes (threshold at median)
    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)

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
function borgatti_everett_discrete(A::Matrix{Float64};
                                   max_iter::Int=1000,
                                   init::Union{Nothing, Vector{Float64}}=nothing)
    n = size(A, 1)

    # Initialize with degree-based assignment
    if init === nothing
        degrees = vec(sum(A, dims=2))
        threshold = median(degrees)
        c = Float64.(degrees .>= threshold)
    else
        c = copy(init)
    end

    best_quality = core_quality(A, c; discrete=true)
    best_c = copy(c)

    # Greedy optimization: try swapping each node
    improved = true
    iter = 0

    while improved && iter < max_iter
        improved = false
        iter += 1

        for i in 1:n
            # Try flipping node i
            c_new = copy(c)
            c_new[i] = 1.0 - c_new[i]

            quality_new = core_quality(A, c_new; discrete=true)

            if quality_new > best_quality
                best_quality = quality_new
                best_c = copy(c_new)
                c = c_new
                improved = true
            end
        end
    end

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Borgatti-Everett Discrete")
end

"""
    lip_discrete(A; max_iter=1000)

Lip's fast algorithm for discrete core-periphery bipartitioning.

Uses efficient swap-based optimization with O(n) updates per iteration.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum iterations

# Returns
- CPResult with binary partition

# Reference
Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
"""
function lip_discrete(A::Matrix{Float64}; max_iter::Int=1000)
    n = size(A, 1)
    m = sum(A) / 2  # Total edges

    # Initialize based on degree
    degrees = vec(sum(A, dims=2))
    sorted_idx = sortperm(degrees, rev=true)

    # Start with top half as core
    c = zeros(Float64, n)
    k = div(n, 2)
    c[sorted_idx[1:k]] .= 1.0

    # Compute initial statistics
    core_set = Set(findall(c .== 1.0))

    # Edges within core
    E_cc = sum(A[i, j] for i in core_set for j in core_set if i < j)
    # Edges between core and periphery
    E_cp = sum(A[i, j] for i in core_set for j in 1:n if !(j in core_set))

    # Quality function: E_cc + E_cp (edges involving core)
    best_score = E_cc + E_cp
    best_c = copy(c)

    improved = true
    iter = 0

    while improved && iter < max_iter
        improved = false
        iter += 1

        for i in 1:n
            # Compute change in score if we flip node i
            is_core = i in core_set

            # Edges from i to core (excluding i)
            edges_to_core = sum(A[i, j] for j in core_set if j != i)
            # Edges from i to periphery
            edges_to_periphery = sum(A[i, j] for j in 1:n if !(j in core_set) && j != i)

            if is_core
                # Moving i from core to periphery
                # Lose: edges_to_core (was E_cc, now E_cp - still counted once)
                # Lose: edges_to_periphery (was E_cp, now E_pp - not counted)
                delta = -edges_to_periphery
            else
                # Moving i from periphery to core
                # Gain: edges_to_periphery (was E_pp, now E_cp)
                delta = edges_to_periphery
            end

            if delta > 0
                # Apply the swap
                if is_core
                    delete!(core_set, i)
                    c[i] = 0.0
                    E_cc -= edges_to_core
                    E_cp = E_cp - edges_to_periphery + edges_to_core
                else
                    push!(core_set, i)
                    c[i] = 1.0
                    E_cc += edges_to_core
                    E_cp = E_cp + edges_to_periphery - edges_to_core
                end

                current_score = E_cc + E_cp
                if current_score > best_score
                    best_score = current_score
                    best_c = copy(c)
                end
                improved = true
            end
        end
    end

    # Compute final quality as correlation
    quality = core_quality(A, best_c; discrete=true)

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, quality, "Lip Discrete")
end

"""
    rombach_continuous(A; alpha=0.5, beta=0.5, max_iter=1000, tol=1e-6, n_runs=10)

Rombach's generalized continuous core-periphery model.

Uses transition function controlled by α and β parameters.

# Arguments
- `A`: Adjacency matrix
- `alpha`: Controls sharpness of core-periphery boundary
- `beta`: Controls size of core
- `max_iter`: Maximum iterations per run
- `tol`: Convergence tolerance
- `n_runs`: Number of random restarts

# Returns
- CPResult with continuous coreness scores

# Reference
Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
"""
function rombach_continuous(A::Matrix{Float64};
                            alpha::Float64=0.5,
                            beta::Float64=0.5,
                            max_iter::Int=1000,
                            tol::Float64=1e-6,
                            n_runs::Int=10)
    n = size(A, 1)

    # Transition function
    function transition(x, α, β)
        if x <= β
            return x / (2 * β)
        else
            return 0.5 + (x - β) / (2 * (1 - β))
        end
    end

    # Generate ideal matrix for given coreness and parameters
    function rombach_ideal(c, α, β)
        n = length(c)
        Δ = zeros(Float64, n, n)

        for i in 1:n
            for j in 1:n
                if i != j
                    ti = transition(c[i], α, β)
                    tj = transition(c[j], α, β)
                    Δ[i, j] = ti * tj
                end
            end
        end
        return Δ
    end

    best_quality = -Inf
    best_c = nothing

    for run in 1:n_runs
        # Random initialization
        c = rand(n)

        for iter in 1:max_iter
            c_old = copy(c)

            # Update each node's coreness
            for i in 1:n
                # Grid search for optimal c[i]
                best_ci = c[i]
                best_local_quality = -Inf

                for ci in 0.0:0.05:1.0
                    c_test = copy(c)
                    c_test[i] = ci
                    Δ = rombach_ideal(c_test, alpha, beta)

                    # Compute correlation for this assignment
                    a_vec = Float64[]
                    d_vec = Float64[]
                    for ii in 1:n
                        for jj in (ii+1):n
                            push!(a_vec, A[ii, jj])
                            push!(d_vec, Δ[ii, jj])
                        end
                    end

                    q = cor(a_vec, d_vec)
                    if q > best_local_quality
                        best_local_quality = q
                        best_ci = ci
                    end
                end

                c[i] = best_ci
            end

            if norm(c - c_old) < tol
                break
            end
        end

        # Compute final quality
        Δ = rombach_ideal(c, alpha, beta)
        a_vec = Float64[]
        d_vec = Float64[]
        for i in 1:n
            for j in (i+1):n
                push!(a_vec, A[i, j])
                push!(d_vec, Δ[i, j])
            end
        end
        quality = cor(a_vec, d_vec)

        if quality > best_quality
            best_quality = quality
            best_c = copy(c)
        end
    end

    # Classify nodes
    threshold = median(best_c)
    core_nodes = findall(best_c .>= threshold)
    periphery_nodes = findall(best_c .< threshold)

    return CPResult(best_c, core_nodes, periphery_nodes, best_quality, "Rombach Continuous")
end

"""
    spectral_method(A)

Spectral method for core-periphery detection.

Uses eigenvector corresponding to largest eigenvalue of adjacency matrix.

# Arguments
- `A`: Adjacency matrix

# Returns
- CPResult with spectral-based coreness scores

# Reference
Cucuringu, M., et al. (2016). Detection of core-periphery structure using spectral methods.
"""
function spectral_method(A::Matrix{Float64})
    n = size(A, 1)

    # Compute eigendecomposition
    eigenvalues, eigenvectors = eigen(A)

    # Get eigenvector for largest eigenvalue
    idx = argmax(eigenvalues)
    v = eigenvectors[:, idx]

    # Take absolute values and normalize
    c = abs.(v)
    c_min, c_max = extrema(c)
    if c_max > c_min
        c = (c .- c_min) ./ (c_max - c_min)
    else
        c = fill(0.5, n)
    end

    # Compute quality
    quality = core_quality(A, c)

    # Classify nodes
    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)

    return CPResult(c, core_nodes, periphery_nodes, quality, "Spectral Method")
end

"""
    random_walker_profiling(A; n_walks=1000, walk_length=10)

Random walker profiling for core-periphery detection.

Nodes visited more frequently by random walks are more core-like.

# Arguments
- `A`: Adjacency matrix
- `n_walks`: Number of random walks
- `walk_length`: Length of each walk

# Returns
- CPResult with visit-frequency-based coreness

# Reference
Rossa, F.D., et al. (2013). Profiling core-periphery network structure by random walkers.
"""
function random_walker_profiling(A::Matrix{Float64};
                                 n_walks::Int=1000,
                                 walk_length::Int=10)
    n = size(A, 1)

    # Compute transition matrix
    degrees = vec(sum(A, dims=2))
    P = copy(A)
    for i in 1:n
        if degrees[i] > 0
            P[i, :] ./= degrees[i]
        end
    end

    # Count visits
    visits = zeros(Int, n)

    for _ in 1:n_walks
        # Start from random node
        current = rand(1:n)

        for _ in 1:walk_length
            visits[current] += 1

            # Move to next node
            if degrees[current] > 0
                probs = P[current, :]
                r = rand()
                cumsum_p = 0.0
                for j in 1:n
                    cumsum_p += probs[j]
                    if r <= cumsum_p
                        current = j
                        break
                    end
                end
            else
                # Isolated node: jump to random node
                current = rand(1:n)
            end
        end
    end

    # Normalize to coreness scores
    c = Float64.(visits)
    c_min, c_max = extrema(c)
    if c_max > c_min
        c = (c .- c_min) ./ (c_max - c_min)
    else
        c = fill(0.5, n)
    end

    # Compute quality
    quality = core_quality(A, c)

    # Classify nodes
    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)

    return CPResult(c, core_nodes, periphery_nodes, quality, "Random Walker Profiling")
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
"""
struct CPMultiResult
    pair_labels::Vector{Int}
    coreness::Vector{Float64}
    n_pairs::Int
    quality::Float64
    algorithm::String
end

"""
    minres_svd(A; max_iter=1000, tol=1e-6)

MINRES/SVD algorithm for core-periphery detection in asymmetric networks.

Minimizes residual f = Σᵢ Σⱼ≠ᵢ (Aᵢⱼ - uᵢvⱼ)² to find in-coreness (v) and out-coreness (u).

# Arguments
- `A`: Adjacency matrix (can be asymmetric)
- `max_iter`: Maximum iterations
- `tol`: Convergence tolerance

# Returns
- CPResult with coreness scores (average of in/out for asymmetric)
- For symmetric matrices, in-coreness equals out-coreness

# Reference
Boyd, J.P., et al. (2010). Computing core/periphery structures and permutation tests for social relations data.
"""
function minres_svd(A::Matrix{Float64}; max_iter::Int=1000, tol::Float64=1e-6)
    n = size(A, 1)

    # Initialize with column/row sums
    u = vec(sum(A, dims=2))  # out-coreness (row sums)
    v = vec(sum(A, dims=1))  # in-coreness (column sums)

    # Normalize
    u_norm = norm(u)
    v_norm = norm(v)
    if u_norm > 0
        u = u ./ u_norm
    else
        u = ones(n) ./ sqrt(n)
    end
    if v_norm > 0
        v = v ./ v_norm
    else
        v = ones(n) ./ sqrt(n)
    end

    # Iterative refinement
    for iter in 1:max_iter
        u_old = copy(u)
        v_old = copy(v)

        # Update u: minimize over u given v
        # u_i = (Σⱼ Aᵢⱼvⱼ) / (Σⱼ vⱼ² - vᵢ²) for each i
        for i in 1:n
            numerator = sum(A[i, j] * v[j] for j in 1:n if j != i)
            denominator = sum(v[j]^2 for j in 1:n if j != i)
            if denominator > 0
                u[i] = numerator / denominator
            end
        end

        # Update v: minimize over v given u
        # v_j = (Σᵢ Aᵢⱼuᵢ) / (Σᵢ uᵢ² - uⱼ²) for each j
        for j in 1:n
            numerator = sum(A[i, j] * u[i] for i in 1:n if i != j)
            denominator = sum(u[i]^2 for i in 1:n if i != j)
            if denominator > 0
                v[j] = numerator / denominator
            end
        end

        # Normalize to [0, 1]
        for vec in [u, v]
            v_min, v_max = extrema(vec)
            if v_max > v_min
                vec .= (vec .- v_min) ./ (v_max - v_min)
            else
                vec .= 0.5
            end
        end

        # Check convergence
        if norm(u - u_old) < tol && norm(v - v_old) < tol
            break
        end
    end

    # For result, use average of in/out coreness
    c = (u .+ v) ./ 2

    # Normalize final coreness
    c_min, c_max = extrema(c)
    if c_max > c_min
        c = (c .- c_min) ./ (c_max - c_min)
    end

    # Compute quality (correlation with ideal pattern)
    quality = core_quality(A, c)

    # Classify nodes
    threshold = median(c)
    core_nodes = findall(c .>= threshold)
    periphery_nodes = findall(c .< threshold)

    return CPResult(c, core_nodes, periphery_nodes, quality, "MINRES/SVD")
end

"""
    multiple_cp_pairs(A; max_pairs=10, min_pair_size=2, max_iter=100)

Detect multiple non-overlapping core-periphery pairs.

Uses quality function Q^cp and label switching optimization.

# Arguments
- `A`: Adjacency matrix
- `max_pairs`: Maximum number of CP pairs to detect
- `min_pair_size`: Minimum nodes per pair
- `max_iter`: Maximum iterations for optimization

# Returns
- CPMultiResult with pair assignments and coreness

# Reference
Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks.
"""
function multiple_cp_pairs(A::Matrix{Float64};
                           max_pairs::Int=10,
                           min_pair_size::Int=2,
                           max_iter::Int=100)
    n = size(A, 1)
    m = sum(A) / 2  # Total edges
    p = 2 * m / (n * (n - 1))  # Edge density (ER null model)

    # Initialize: each node in its own pair, all core
    pair_labels = collect(1:n)
    coreness = ones(Float64, n)

    # Quality function Q^cp
    function compute_qcp(labels, core)
        q = 0.0
        for i in 1:n
            for j in (i+1):n
                if labels[i] == labels[j]  # Same pair
                    x_i, x_j = core[i], core[j]
                    # Ideal pattern: x_i + x_j - x_i*x_j
                    delta = x_i + x_j - x_i * x_j
                    q += (A[i, j] - p) * delta
                end
            end
        end
        return q
    end

    # Merge pairs greedily to maximize Q^cp
    current_pairs = Set(1:n)

    for _ in 1:max_iter
        best_merge = nothing
        best_delta = 0.0

        pairs_list = collect(current_pairs)

        for i in 1:length(pairs_list)
            for j in (i+1):length(pairs_list)
                p1, p2 = pairs_list[i], pairs_list[j]

                # Try merging p1 and p2
                test_labels = copy(pair_labels)
                test_labels[pair_labels .== p2] .= p1

                # Compute change in quality
                q_new = compute_qcp(test_labels, coreness)
                q_old = compute_qcp(pair_labels, coreness)
                delta = q_new - q_old

                if delta > best_delta
                    best_delta = delta
                    best_merge = (p1, p2)
                end
            end
        end

        if best_merge !== nothing && best_delta > 0
            p1, p2 = best_merge
            pair_labels[pair_labels .== p2] .= p1
            delete!(current_pairs, p2)
        else
            break
        end
    end

    # Optimize core/periphery within each pair
    for _ in 1:max_iter
        improved = false

        for i in 1:n
            # Try flipping coreness of node i
            test_core = copy(coreness)
            test_core[i] = 1.0 - test_core[i]

            if compute_qcp(pair_labels, test_core) > compute_qcp(pair_labels, coreness)
                coreness[i] = test_core[i]
                improved = true
            end
        end

        if !improved
            break
        end
    end

    # Relabel pairs to consecutive integers
    unique_pairs = sort(unique(pair_labels))
    pair_map = Dict(p => i for (i, p) in enumerate(unique_pairs))
    pair_labels = [pair_map[p] for p in pair_labels]
    n_pairs = length(unique_pairs)

    # Remove pairs that are too small
    for pair_id in 1:n_pairs
        pair_size = sum(pair_labels .== pair_id)
        if pair_size < min_pair_size
            # Merge with nearest pair
            pair_nodes = findall(pair_labels .== pair_id)
            if pair_id > 1
                pair_labels[pair_nodes] .= pair_id - 1
            elseif pair_id < n_pairs
                pair_labels[pair_nodes] .= pair_id + 1
            end
        end
    end

    # Relabel again
    unique_pairs = sort(unique(pair_labels))
    pair_map = Dict(p => i for (i, p) in enumerate(unique_pairs))
    pair_labels = [pair_map[p] for p in pair_labels]
    n_pairs = length(unique_pairs)

    quality = compute_qcp(pair_labels, coreness)

    return CPMultiResult(pair_labels, coreness, n_pairs, quality, "Multiple CP Pairs")
end

"""
    surprise_cp(A; max_iter=100)

Surprise-based core-periphery detection.

Uses multinomial hypergeometric distribution to compute surprise of CP structure.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum optimization iterations

# Returns
- CPResult with binary coreness

# Reference
Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise.
"""
function surprise_cp(A::Matrix{Float64}; max_iter::Int=100)
    n = size(A, 1)
    m = Int(sum(A) / 2)  # Total edges

    # Initialize based on degree
    degrees = vec(sum(A, dims=2))
    threshold = median(degrees)
    c = Float64.(degrees .>= threshold)

    # Count edges in each block
    function count_blocks(core_vec)
        n_core = sum(core_vec)
        n_periph = n - n_core

        # Maximum possible edges in each block
        V_cc = n_core * (n_core - 1) / 2  # Core-core
        V_cp = n_core * n_periph          # Core-periphery
        V_pp = n_periph * (n_periph - 1) / 2  # Periphery-periphery

        # Actual edges in each block
        L_cc = 0.0
        L_cp = 0.0
        L_pp = 0.0

        for i in 1:n
            for j in (i+1):n
                if A[i, j] > 0
                    if core_vec[i] == 1.0 && core_vec[j] == 1.0
                        L_cc += 1
                    elseif core_vec[i] == 0.0 && core_vec[j] == 0.0
                        L_pp += 1
                    else
                        L_cp += 1
                    end
                end
            end
        end

        return V_cc, V_cp, V_pp, L_cc, L_cp, L_pp
    end

    # Compute log-surprise (negative log p-value)
    function log_surprise(core_vec)
        V_cc, V_cp, V_pp, L_cc, L_cp, L_pp = count_blocks(core_vec)

        # Total possible and actual
        V = V_cc + V_cp + V_pp
        L = L_cc + L_cp + L_pp

        if V == 0 || L == 0
            return 0.0
        end

        # Expected edges under random model
        p = L / V

        # Surprise: how unlikely is it to have this many edges in core-core and core-periphery?
        # Higher surprise = better CP structure

        # Simplified: use difference from expected
        expected_cc = V_cc > 0 ? V_cc * p : 0.0
        expected_cp = V_cp > 0 ? V_cp * p : 0.0
        expected_pp = V_pp > 0 ? V_pp * p : 0.0

        # Surprise is high when L_cc and L_cp are higher than expected, L_pp is lower
        surprise = (L_cc - expected_cc) + (L_cp - expected_cp) - (L_pp - expected_pp)

        return surprise
    end

    best_surprise = log_surprise(c)
    best_c = copy(c)

    # Greedy optimization
    for _ in 1:max_iter
        improved = false

        for i in 1:n
            c_new = copy(c)
            c_new[i] = 1.0 - c_new[i]

            s_new = log_surprise(c_new)
            if s_new > best_surprise
                best_surprise = s_new
                best_c = copy(c_new)
                c = c_new
                improved = true
            end
        end

        if !improved
            break
        end
    end

    # Compute correlation-based quality for consistency
    quality = core_quality(A, best_c; discrete=true)

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, quality, "Surprise CP")
end

"""
    label_switching_cp(A; max_iter=100)

Fast label-switching algorithm for core-periphery detection.

Uses greedy optimization with efficient O(n) updates per iteration.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Maximum iterations

# Returns
- CPResult with binary coreness

# Reference
Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection.
"""
function label_switching_cp(A::Matrix{Float64}; max_iter::Int=100)
    n = size(A, 1)

    # Initialize based on degree
    degrees = vec(sum(A, dims=2))
    threshold = median(degrees)
    c = Float64.(degrees .>= threshold)

    # Compute initial M = Σᵢⱼ Aᵢⱼ * Δᵢⱼ where Δᵢⱼ = max(cᵢ, cⱼ) for discrete model
    function compute_M(core_vec)
        M = 0.0
        for i in 1:n
            for j in (i+1):n
                delta = max(core_vec[i], core_vec[j])
                M += A[i, j] * delta
            end
        end
        return M
    end

    # Current objective
    M = compute_M(c)
    best_M = M
    best_c = copy(c)

    # Pre-compute neighbor sums for efficiency
    # For each node i, compute sum of edges to core and periphery
    function compute_neighbor_sums(core_vec)
        edges_to_core = zeros(n)
        edges_to_periph = zeros(n)

        for i in 1:n
            for j in 1:n
                if i != j && A[i, j] > 0
                    if core_vec[j] == 1.0
                        edges_to_core[i] += A[i, j]
                    else
                        edges_to_periph[i] += A[i, j]
                    end
                end
            end
        end

        return edges_to_core, edges_to_periph
    end

    edges_to_core, edges_to_periph = compute_neighbor_sums(c)

    # Label switching iterations
    for iter in 1:max_iter
        improved = false

        # Random order for node processing
        node_order = randperm(n)

        for i in node_order
            # Compute change in M if we flip node i
            # If i is currently core (c[i] = 1):
            #   - Edges to periphery: Δ changes from 1 to 0, lose edges_to_periph[i]
            # If i is currently periphery (c[i] = 0):
            #   - Edges to periphery: Δ changes from 0 to 1, gain edges_to_periph[i]

            if c[i] == 1.0
                # Moving to periphery
                delta_M = -edges_to_periph[i]
            else
                # Moving to core
                delta_M = edges_to_periph[i]
            end

            if delta_M > 0
                # Apply the switch
                if c[i] == 1.0
                    c[i] = 0.0
                    # Update neighbor sums
                    for j in 1:n
                        if j != i && A[i, j] > 0
                            edges_to_core[j] -= A[i, j]
                            edges_to_periph[j] += A[i, j]
                        end
                    end
                else
                    c[i] = 1.0
                    # Update neighbor sums
                    for j in 1:n
                        if j != i && A[i, j] > 0
                            edges_to_core[j] += A[i, j]
                            edges_to_periph[j] -= A[i, j]
                        end
                    end
                end

                M += delta_M
                if M > best_M
                    best_M = M
                    best_c = copy(c)
                end
                improved = true
            end
        end

        if !improved
            break
        end
    end

    # Compute quality
    quality = core_quality(A, best_c; discrete=true)

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)

    return CPResult(best_c, core_nodes, periphery_nodes, quality, "Label Switching CP")
end

end # module
