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

# References
- Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
- Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
- Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
- Cucuringu, M., et al. (2016). Detection of core-periphery structure using spectral methods.
- Rossa, F.D., et al. (2013). Profiling core-periphery network structure by random walkers.
"""
module CorePeriphery

using LinearAlgebra
using Statistics
using Random

export
    # Data structures
    CPResult,
    # Algorithms
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    lip_discrete,
    rombach_continuous,
    spectral_method,
    random_walker_profiling,
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

end # module
