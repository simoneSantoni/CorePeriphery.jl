using Arpack

"""
    _lowrank_core(A; beta=0.1)

Apply Cucuringu et al.'s LowRank-Core algorithm to a simple undirected graph.

The two eigenpairs of `A` with largest eigenvalue magnitude define a rank-two
reconstruction. Entries strictly larger than `0.5` are set to one and all other
entries to zero; row sums of that thresholded matrix are the LowRank-Core scores.
Find-Cut orders vertices by decreasing score and maximizes the CP-density objective
over cuts for which both parts contain at least `ceil(beta * n)` vertices.

`beta` is the Find-Cut boundary fraction, not a known core fraction. The return value
is a named tuple with raw LowRank-Core scores in `coreness`, node-index vectors
`core_nodes` and `periphery_nodes`, and the unscaled CP-density in `quality`.

Sparse matrices use ARPACK to compute only the two eigenpairs required by the
method. Rank-two threshold scores and Find-Cut statistics are streamed without
materializing a dense reconstruction or a permuted adjacency matrix.
"""
function _lowrank_core(A::AbstractMatrix{<:Real}; beta::Real=0.1)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("adjacency matrix must be square"))
    n >= 2 || throw(ArgumentError("LowRank-Core requires at least two vertices"))
    0 <= beta <= 0.5 || throw(ArgumentError("beta must be in [0, 0.5]"))
    all(isfinite, A) || throw(ArgumentError("adjacency matrix entries must be finite"))
    issymmetric(A) || throw(ArgumentError("LowRank-Core requires a symmetric matrix"))
    all(iszero, diag(A)) || throw(ArgumentError("LowRank-Core does not allow self-loops"))
    all(x -> x == 0 || x == 1, A) ||
        throw(ArgumentError("published LowRank-Core requires a binary adjacency matrix"))

    boundary = max(1, ceil(Int, beta * n))
    boundary <= n - boundary ||
        throw(ArgumentError("beta leaves no feasible nonempty core/periphery cut"))

    values, vectors = _lowrank_eigenpairs(A)
    scores = _thresholded_rank_two_scores(values, vectors)
    order = sortperm(scores; rev=true, alg=Base.Sort.MergeSort)

    best_k, quality = _find_cut(A, order, boundary)
    core_nodes = sort!(order[1:best_k])
    periphery_nodes = sort!(order[(best_k + 1):n])

    return (
        coreness=scores,
        core_nodes=core_nodes,
        periphery_nodes=periphery_nodes,
        quality=quality,
    )
end

function _lowrank_eigenpairs(A::AbstractMatrix{<:Real})
    decomposition = eigen(Symmetric(Matrix{Float64}(A)))
    selected = partialsortperm(abs.(decomposition.values), 1:2; rev=true)
    return decomposition.values[selected], decomposition.vectors[:, selected]
end

function _lowrank_eigenpairs(A::SparseMatrixCSC{<:Real,<:Integer})
    n = size(A, 1)
    n <= 3 && return _lowrank_eigenpairs(Matrix{Float64}(A))

    float_A = SparseMatrixCSC{Float64,Int}(A)
    values, vectors, converged, _, _, _ = eigs(
        Symmetric(float_A);
        nev=2,
        which=:LM,
        tol=0.0,
        maxiter=max(300, 20 * n),
        ritzvec=true,
    )
    converged == 2 || throw(ErrorException(
        "ARPACK converged $converged of the two required LowRank-Core eigenpairs"))
    selected = sortperm(abs.(values); rev=true, alg=Base.Sort.MergeSort)
    return Float64.(real(values[selected])), Float64.(real(vectors[:, selected]))
end

function _thresholded_rank_two_scores(values, vectors)
    n = size(vectors, 1)
    scores = zeros(Float64, n)
    for i in 1:n
        diagonal = values[1] * vectors[i, 1]^2 + values[2] * vectors[i, 2]^2
        diagonal > 0.5 && (scores[i] += 1.0)
        for j in (i + 1):n
            reconstructed = values[1] * vectors[i, 1] * vectors[j, 1] +
                            values[2] * vectors[i, 2] * vectors[j, 2]
            if reconstructed > 0.5
                scores[i] += 1.0
                scores[j] += 1.0
            end
        end
    end
    return scores
end

"""
Find the score-prefix cut maximizing equation (4.4) with gamma = 0.
"""
function _find_cut(A::AbstractMatrix{<:Real}, order::Vector{Int}, boundary::Int)
    n = length(order)
    strengths, edges_to_earlier, total_edges = _ordered_cut_statistics(A, order)

    core_edges = 0.0
    cross_edges = 0.0
    best_quality = -Inf
    best_k = boundary

    for k in 1:(n - 1)
        to_core = edges_to_earlier[k]
        to_periphery = strengths[k] - to_core
        core_edges += to_core
        cross_edges += to_periphery - to_core

        if boundary <= k <= n - boundary
            periphery_edges = total_edges - core_edges - cross_edges
            core_volume = k * (k - 1) / 2
            cross_volume = k * (n - k)
            periphery_volume = (n - k) * (n - k - 1) / 2

            core_density = iszero(core_volume) ? 0.0 : core_edges / core_volume
            cross_density = cross_edges / cross_volume
            periphery_density =
                iszero(periphery_volume) ? 0.0 : periphery_edges / periphery_volume
            quality = core_density + cross_density - periphery_density

            # Stable score sorting plus strict improvement makes the smallest
            # maximizing core deterministic.
            if quality > best_quality
                best_quality = quality
                best_k = k
            end
        end
    end

    return best_k, best_quality
end


function _ordered_cut_statistics(A::AbstractMatrix{<:Real}, order::Vector{Int})
    n = length(order)
    strengths = Float64.(vec(sum(A, dims=2)))[order]
    edges_to_earlier = zeros(Float64, n)
    total_edges = 0.0
    for earlier in 1:(n - 1), later in (earlier + 1):n
        weight = Float64(A[order[earlier], order[later]])
        edges_to_earlier[later] += weight
        total_edges += weight
    end
    return strengths, edges_to_earlier, total_edges
end

function _ordered_cut_statistics(
    A::SparseMatrixCSC{<:Real,<:Integer}, order::Vector{Int})
    n = length(order)
    position = Vector{Int}(undef, n)
    for (rank, node) in enumerate(order)
        position[node] = rank
    end

    strengths = Float64.(vec(sum(A, dims=2)))[order]
    edges_to_earlier = zeros(Float64, n)
    total_edges = 0.0
    rows, columns, values = findnz(A)
    for index in eachindex(values)
        source = rows[index]
        target = columns[index]
        source < target || continue
        weight = Float64(values[index])
        later = max(position[source], position[target])
        edges_to_earlier[later] += weight
        total_edges += weight
    end
    return strengths, edges_to_earlier, total_edges
end
