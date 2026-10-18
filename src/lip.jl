"""
    lip_discrete(A; max_iter=1000)

Lip's fast algorithm for discrete core-periphery bipartitioning.

For a binary, undirected network, this implements Lip's exact algorithm. Nodes are
ordered by decreasing degree, and the loss (missing core-core edges plus observed
periphery-periphery edges) is evaluated for every proper, nonempty degree-prefix
core. The smallest core is returned when several sizes attain the same minimum.

# Arguments
- `A`: Adjacency matrix
- `max_iter`: Deprecated compatibility keyword; ignored because the algorithm is exact.
  Passing it emits a deprecation warning (shown with `--depwarn=yes`).

# Returns
- CPResult with binary partition

# Reference
Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
"""
function lip_discrete(A::AbstractMatrix{<:Real}; max_iter::Union{Nothing,Int}=nothing)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=true, binary=true)
    n >= 2 || throw(ArgumentError("Lip's bipartition requires at least two nodes"))
    if max_iter !== nothing
        max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
        Base.depwarn(
            "the `max_iter` keyword of `lip_discrete` is deprecated and ignored: " *
            "the algorithm is exact and needs no iteration budget",
            :lip_discrete,
        )
    end

    _warn_disconnected(A, :lip_discrete)

    degrees = vec(sum(A, dims=2))
    # Input is validated symmetric and nonnegative, so isolates are exactly the
    # zero-degree nodes; reusing the degree pass avoids a second O(n^2) scan.
    isolates = findall(iszero, degrees)
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

    core_nodes = findall(best_c .== 1.0)
    periphery_nodes = findall(best_c .== 0.0)
    # Isolates sort last, so they enter the degree-prefix core only when the
    # network has no edges at all; they are peripheral by definition.
    core_nodes, periphery_nodes = _enforce_isolates!(
        best_c, core_nodes, periphery_nodes, isolates)

    # Compute final quality as correlation
    quality = core_quality(A, best_c; discrete=true)

    return CPResult(best_c, core_nodes, periphery_nodes, quality, "Lip Discrete")
end
