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
    A = _unwrap_adjacency(A)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("adjacency matrix must be square"))
    n >= 2 || throw(ArgumentError("LowRank-Core requires at least two vertices"))
    0 <= beta <= 0.5 || throw(ArgumentError("beta must be in [0, 0.5]"))
    _allfinite(A) || throw(ArgumentError("adjacency matrix entries must be finite"))
    _issymmetric(A) || throw(ArgumentError("LowRank-Core requires a symmetric matrix"))
    all(iszero, diag(A)) || throw(ArgumentError("LowRank-Core does not allow self-loops"))
    _isbinary(A) ||
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

_allfinite(A::AbstractMatrix{<:Real}) = all(isfinite, A)
_allfinite(A::SparseMatrixCSC{<:Real}) = all(isfinite, nonzeros(A))

function _lowrank_eigenpairs(A::AbstractMatrix{<:Real})
    decomposition = eigen(Symmetric(Matrix{Float64}(A)))
    selected = partialsortperm(abs.(decomposition.values), 1:2; rev=true)
    return decomposition.values[selected], decomposition.vectors[:, selected]
end

# ARPACK settings for the two largest-magnitude eigenpairs. A finite tolerance
# stops the implicitly restarted Lanczos iteration once the Ritz residuals are
# below `tol · |λ|`; `tol = 0` asks for machine precision and can stall on
# clustered spectra (a 2000-node ring needed ~14,000 iterations). A Krylov
# basis of `_LOWRANK_NCV` vectors resolves clustered extremal eigenvalues in far
# fewer restarts than the minimal `2·nev + 1`.
const _LOWRANK_TOL = 1e-10
const _LOWRANK_NCV = 48

# Arpack.jl is not reentrant (it keeps workspace in global Fortran state), so
# concurrent sparse LowRank-Core calls, for example null samples evaluated by
# `cp_significance(threaded=true)`, are serialized through this lock.
const _ARPACK_LOCK = ReentrantLock()

function _lowrank_eigenpairs(A::SparseMatrixCSC{<:Real,<:Integer};
    tol::Float64=_LOWRANK_TOL,
    ncv::Int=_LOWRANK_NCV,
    maxiter::Union{Nothing,Int}=nothing)
    n = size(A, 1)
    n <= 3 && return _lowrank_eigenpairs(Matrix{Float64}(A))

    float_A = SparseMatrixCSC{Float64,Int}(A)
    ncv = clamp(ncv, 3, n)
    maxiter = maxiter === nothing ? max(300, 20 * n) : maxiter
    values, vectors, converged, iterations = try
        lock(_ARPACK_LOCK) do
            eigs(
                Symmetric(float_A);
                nev=2,
                which=:LM,
                tol=tol,
                ncv=ncv,
                maxiter=maxiter,
                ritzvec=true,
            )
        end
    catch error
        error isa InterruptException && rethrow()
        throw(ConvergenceError(
            "ARPACK failed while computing the two largest-magnitude eigenpairs " *
            "for LowRank-Core (n=$n, tol=$tol, ncv=$ncv, maxiter=$maxiter): " *
            sprint(showerror, error)))
    end
    converged == 2 || throw(ConvergenceError(
        "ARPACK converged $converged of the two required LowRank-Core eigenpairs " *
        "after $iterations restarts (n=$n, tol=$tol, ncv=$ncv, " *
        "maxiter=$maxiter); the extremal eigenvalues may be clustered, so " *
        "consider a dense input (LAPACK) for this network"))
    selected = sortperm(abs.(values); rev=true, alg=Base.Sort.MergeSort)
    return Float64.(real(values[selected])), Float64.(real(vectors[:, selected]))
end

# Row sums of the thresholded rank-two reconstruction
# Â_ij = λ₁ u_i u_j + λ₂ v_i v_j > 1/2 (the diagonal included, as in the oracle
# convention of the test suite). Counting the pairs above the threshold is a
# half-plane range-counting problem in the plane of (u_j, v_j), for which no
# exact sub-quadratic method is practical without geometric range-search
# structures, so the pair loop stays O(|active|²). It is pruned exactly: row i
# can only exceed the threshold for some j if its best case
# max_j (λ₁ u_i u_j + λ₂ v_i v_j), bounded by the extreme coordinates, exceeds
# 1/2, and a pair needs both rows active. On sparse core-periphery networks
# the Perron vector concentrates on the core, so |active| ≪ n.
function _thresholded_rank_two_scores(values, vectors)
    n = size(vectors, 1)
    λ1 = Float64(values[1])
    λ2 = Float64(values[2])
    u = vectors[:, 1]
    v = vectors[:, 2]
    scores = zeros(Float64, n)
    u_min, u_max = extrema(u)
    v_min, v_max = extrema(v)
    active = Int[]
    for i in 1:n
        a = λ1 * u[i]
        b = λ2 * v[i]
        bound = max(a * u_min, a * u_max) + max(b * v_min, b * v_max)
        bound > 0.5 && push!(active, i)
        diagonal = a * u[i] + b * v[i]
        diagonal > 0.5 && (scores[i] += 1.0)
    end
    @inbounds for (position, i) in enumerate(active)
        a = λ1 * u[i]
        b = λ2 * v[i]
        for later in (position + 1):length(active)
            j = active[later]
            if a * u[j] + b * v[j] > 0.5
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
