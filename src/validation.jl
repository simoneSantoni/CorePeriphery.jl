# Common adjacency-matrix contract. Self-loops are ignored by every objective
# and are therefore rejected rather than silently discarded.

"""
    _unwrap_adjacency(A)

Return a plain `SparseMatrixCSC` when `A` is a `Symmetric`, `Hermitian`,
`Transpose`, or `Adjoint` wrapper around one, so that the sparse kernels are
dispatched instead of the generic dense `AbstractMatrix` fallbacks. Every other
matrix is returned unchanged.
"""
_unwrap_adjacency(A::AbstractMatrix{<:Real}) = A
_unwrap_adjacency(A::Symmetric{<:Real,<:SparseMatrixCSC}) = sparse(A)
_unwrap_adjacency(A::Hermitian{<:Real,<:SparseMatrixCSC}) = sparse(A)
_unwrap_adjacency(A::Transpose{<:Real,<:SparseMatrixCSC}) = copy(A)
_unwrap_adjacency(A::Adjoint{<:Real,<:SparseMatrixCSC}) = copy(A)

# Symmetry is tested with an absolute tolerance scaled by the largest entry so
# that floating-point noise (for example from `(A + A') / 2`) does not make an
# undirected network look directed, while genuinely directed input still fails.
_symmetry_atol(largest::Real) = 1e-12 * max(1.0, Float64(largest))

function _issymmetric(A::AbstractMatrix{<:Real}, atol::Float64)
    n, m = size(A)
    n == m || return false
    @inbounds for j in 2:n, i in 1:(j - 1)
        abs(Float64(A[i, j]) - Float64(A[j, i])) <= atol || return false
    end
    return true
end

function _issymmetric(A::SparseMatrixCSC{<:Real}, atol::Float64)
    n, m = size(A)
    n == m || return false
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        i == j && continue
        abs(Float64(values[index]) - Float64(A[j, i])) <= atol || return false
    end
    return true
end

"""
    _issymmetric(A)

Symmetry test used by every undirected-only detector and by the detectors that
choose between undirected and directed treatment. Entries are compared with the
absolute tolerance `1e-12 * max(1, maximum(abs, A))`.
"""
function _issymmetric(A::AbstractMatrix{<:Real})
    largest = isempty(A) ? 0.0 : Float64(maximum(abs, A))
    return _issymmetric(A, _symmetry_atol(largest))
end

function _issymmetric(A::SparseMatrixCSC{<:Real})
    values = nonzeros(A)
    largest = isempty(values) ? 0.0 : Float64(maximum(abs, values))
    return _issymmetric(A, _symmetry_atol(largest))
end

function _validate_adjacency(A::AbstractMatrix{<:Real};
    symmetric::Bool=false,
    binary::Bool=false)
    n, m = size(A)
    n == m || throw(DimensionMismatch("adjacency matrix must be square, got $(n)×$(m)"))
    n > 0 || throw(ArgumentError("adjacency matrix must contain at least one node"))

    largest = 0.0
    for j in axes(A, 2), i in axes(A, 1)
        a = A[i, j]
        isfinite(a) || throw(ArgumentError("adjacency matrix entries must be finite"))
        a >= zero(a) || throw(ArgumentError("adjacency matrix entries must be nonnegative"))
        if binary && !(iszero(a) || a == one(a))
            throw(ArgumentError("this algorithm requires a binary adjacency matrix"))
        end
        largest = max(largest, Float64(a))
    end
    all(iszero, diag(A)) || throw(ArgumentError("self-loops are not supported"))
    if symmetric && !_issymmetric(A, _symmetry_atol(largest))
        throw(ArgumentError("this algorithm requires an undirected (symmetric) adjacency matrix"))
    end
    return n
end

function _validate_adjacency(A::SparseMatrixCSC{<:Real};
    symmetric::Bool=false,
    binary::Bool=false)
    n, m = size(A)
    n == m || throw(DimensionMismatch("adjacency matrix must be square, got $(n)×$(m)"))
    n > 0 || throw(ArgumentError("adjacency matrix must contain at least one node"))

    rows = rowvals(A)
    values = nonzeros(A)
    largest = 0.0
    for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        a = values[index]
        isfinite(a) || throw(ArgumentError("adjacency matrix entries must be finite"))
        a >= zero(a) || throw(ArgumentError("adjacency matrix entries must be nonnegative"))
        i == j && !iszero(a) && throw(ArgumentError("self-loops are not supported"))
        if binary && !(iszero(a) || a == one(a))
            throw(ArgumentError("this algorithm requires a binary adjacency matrix"))
        end
        largest = max(largest, Float64(a))
    end
    if symmetric && !_issymmetric(A, _symmetry_atol(largest))
        throw(ArgumentError("this algorithm requires an undirected (symmetric) adjacency matrix"))
    end
    return n
end

function _normalize_scores(x::AbstractVector{<:Real})
    isempty(x) && return Float64[]
    lo, hi = extrema(x)
    hi > lo || return fill(0.5, length(x))
    return (Float64.(x) .- lo) ./ (hi - lo)
end

# Connectivity and isolate helpers shared by the detectors. Every detector
# treats the network as undirected for these purposes (weak connectivity for
# directed input): a node is an isolate when it has neither in- nor out-weight.

function _neighbor_lists(A::AbstractMatrix{<:Real})
    n = size(A, 1)
    neighbors = [Int[] for _ in 1:n]
    @inbounds for j in 1:n, i in 1:n
        i == j && continue
        if !iszero(A[i, j])
            push!(neighbors[i], j)
            push!(neighbors[j], i)
        end
    end
    return neighbors
end

function _neighbor_lists(A::SparseMatrixCSC{<:Real})
    n = size(A, 1)
    neighbors = [Int[] for _ in 1:n]
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in 1:n, index in nzrange(A, j)
        i = rows[index]
        i == j && continue
        if !iszero(values[index])
            push!(neighbors[i], j)
            push!(neighbors[j], i)
        end
    end
    return neighbors
end

"""
    _component_labels(A)

Label the weakly connected components of the network with adjacency matrix `A`
(an edge `i -- j` exists when `A[i, j]` or `A[j, i]` is nonzero). Returns a
`Vector{Int}` of consecutive component labels in order of first appearance; an
isolate forms its own component.
"""
function _component_labels(A::AbstractMatrix{<:Real})
    n = size(A, 1)
    neighbors = _neighbor_lists(A)
    labels = zeros(Int, n)
    next_label = 0
    stack = Int[]
    for start in 1:n
        labels[start] == 0 || continue
        next_label += 1
        labels[start] = next_label
        push!(stack, start)
        while !isempty(stack)
            node = pop!(stack)
            for neighbor in neighbors[node]
                if labels[neighbor] == 0
                    labels[neighbor] = next_label
                    push!(stack, neighbor)
                end
            end
        end
    end
    return labels
end

"""
    _isolates(A)

Sorted indices of the nodes with zero in- and out-strength.
"""
function _isolates(A::AbstractMatrix{<:Real})
    # One column-major pass: a node is isolated iff its row and column are zero.
    n = size(A, 1)
    touched = falses(n)
    @inbounds for j in 1:n, i in 1:n
        if !iszero(A[i, j])
            touched[i] = true
            touched[j] = true
        end
    end
    return findall(!, touched)
end

function _isolates(A::SparseMatrixCSC{<:Real})
    n = size(A, 1)
    touched = falses(n)
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in 1:n, k in nzrange(A, j)
        if !iszero(values[k])
            touched[rows[k]] = true
            touched[j] = true
        end
    end
    return findall(!, touched)
end

const _NULL_SAMPLE_KEY = :CorePeriphery_null_sample_evaluation

"""
    _suppressing_component_checks()

True while the current task is evaluating a `cp_significance` null sample (see
[`_with_component_checks_suppressed`](@ref)). Task-local, so concurrent user
calls on other tasks are unaffected.
"""
_suppressing_component_checks() =
    get(task_local_storage(), _NULL_SAMPLE_KEY, false)::Bool

"""
    _with_component_checks_suppressed(f)

Run `f()` with the disconnected-input check of `_warn_disconnected` disabled for
the current task only. Isolate handling (`_enforce_isolates!`) is unaffected, so
null and observed statistics are computed by the same estimator.
"""
_with_component_checks_suppressed(f) = task_local_storage(f, _NULL_SAMPLE_KEY, true)

"""
    _warn_disconnected(A, detector::Symbol)

Emit a warning (once per `detector` per session) when `A` contains more than
one non-trivial connected component, i.e. more than one component with at
least two nodes. Core-periphery objectives are defined for a single connected
network; on a union of components the optimiser may assemble a core from
several components or label a whole component as periphery.

The check is skipped (returning `0`) while `_suppressing_component_checks()` is
true, i.e. inside `cp_significance` null-sample evaluations: random null graphs
are routinely disconnected, the warning would describe the null rather than the
user's network, and the component search would add `O(n^2)` work per sample on
dense input.
"""
function _warn_disconnected(A::AbstractMatrix{<:Real}, detector::Symbol)
    _suppressing_component_checks() && return 0
    labels = _component_labels(A)
    sizes = zeros(Int, maximum(labels; init=0))
    for label in labels
        sizes[label] += 1
    end
    nontrivial = count(>=(2), sizes)
    if nontrivial > 1
        @warn string(
            detector, ": the network has ", nontrivial,
            " connected components with at least two nodes; core-periphery ",
            "objectives assume a single connected network, so consider analysing ",
            "each component separately",
        ) maxlog = 1 _id = Symbol(:disconnected_input_, detector)
    end
    return nontrivial
end

"""
    _enforce_isolates!(coreness, core_nodes, periphery_nodes, isolates)

Give every isolate coreness `0` and move it to the periphery. Returns the
corrected `(core_nodes, periphery_nodes)` pair (sorted).
"""
function _enforce_isolates!(coreness::AbstractVector{Float64},
    core_nodes::Vector{Int}, periphery_nodes::Vector{Int},
    isolates::AbstractVector{<:Integer})
    isempty(isolates) && return core_nodes, periphery_nodes
    coreness[isolates] .= 0.0
    in_core = falses(length(coreness))
    in_core[core_nodes] .= true
    in_core[isolates] .= false
    return findall(in_core), findall(!, in_core)
end
