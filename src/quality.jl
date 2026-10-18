"""
    adjacency_to_matrix(edges, n; directed=false, combine=nothing)
    adjacency_to_matrix(edges; directed=false, combine=nothing, sort_labels=false)

Convert an edge list to a dense `Matrix{Float64}` adjacency matrix.

Each edge is a tuple `(i, j)` (weight `1`) or `(i, j, weight)`. Self-loops,
negative weights, and non-finite weights are rejected.

- With `n` given, endpoints must be one-based integers in `1:n` and the
  matrix is returned.
- Without `n`, endpoints may be any node identifiers (strings, symbols,
  integers, ...). The method returns `(A, labels)`, where `labels[k]` is the
  identifier of node `k`. Labels are listed in order of first appearance in
  `edges` (first endpoint before second), or sorted with `sort_labels=true`.

`directed=false` (default) stores every edge in both `A[i, j]` and `A[j, i]`;
`(i, j)` and `(j, i)` then denote the same dyad. `directed=true` stores only
`A[i, j]`, so `(i, j)` and `(j, i)` are distinct dyads.

A dyad that occurs more than once with the same weight is stored once. A dyad
that occurs with *conflicting* weights raises an `ArgumentError` unless
`combine` names how to resolve it: `:sum`, `:max`, `:min`, or `:last` (the
last supplied weight wins). With `combine=:sum`, repeated identical weights
are summed as well.
"""
function adjacency_to_matrix(edges::AbstractVector{<:Tuple{<:Integer,<:Integer}}, n::Integer;
    directed::Bool=false, combine::Union{Nothing,Symbol}=nothing)
    return _edges_to_matrix(edges, n, directed, _check_combine(combine))
end

function adjacency_to_matrix(
    edges::AbstractVector{<:Tuple{<:Integer,<:Integer,<:Real}}, n::Integer;
    directed::Bool=false, combine::Union{Nothing,Symbol}=nothing)
    return _edges_to_matrix(edges, n, directed, _check_combine(combine))
end

# Compatibility fallback for deliberately heterogeneous edge containers. The
# concrete two- and three-tuple methods above are the allocation-free fast path.
function adjacency_to_matrix(edges::AbstractVector, n::Integer;
    directed::Bool=false, combine::Union{Nothing,Symbol}=nothing)
    return _edges_to_matrix(edges, n, directed, _check_combine(combine))
end

function adjacency_to_matrix(edges::AbstractVector;
    directed::Bool=false, combine::Union{Nothing,Symbol}=nothing,
    sort_labels::Bool=false)
    combine = _check_combine(combine)
    labels = Any[]
    index = Dict{Any,Int}()
    for edge in edges
        edge isa Tuple || throw(ArgumentError("each edge must be a tuple"))
        length(edge) in (2, 3) ||
            throw(ArgumentError("each edge must be (i, j) or (i, j, weight)"))
        for endpoint in (edge[1], edge[2])
            if !haskey(index, endpoint)
                push!(labels, endpoint)
                index[endpoint] = length(labels)
            end
        end
    end
    isempty(labels) && throw(ArgumentError("the edge list contains no nodes"))
    labels = [label for label in labels]  # narrow the element type
    if sort_labels
        labels = sort(labels)
        index = Dict{Any,Int}(label => k for (k, label) in enumerate(labels))
    end
    integer_edges = map(edges) do edge
        weight = length(edge) == 2 ? 1.0 : edge[3]
        (index[edge[1]], index[edge[2]], weight)
    end
    A = _edges_to_matrix(integer_edges, length(labels), directed, combine)
    return A, labels
end

function _check_combine(combine::Union{Nothing,Symbol})
    combine === nothing || combine in (:sum, :max, :min, :last) ||
        throw(ArgumentError("combine must be nothing, :sum, :max, :min, or :last"))
    return combine
end

function _edges_to_matrix(edges, n::Integer, directed::Bool,
    combine::Union{Nothing,Symbol})
    n >= 1 || throw(ArgumentError("n must be positive"))
    A = zeros(Float64, n, n)
    filled = falses(n, n)
    for edge in edges
        edge isa Tuple || throw(ArgumentError("each edge must be a tuple"))
        length(edge) in (2, 3) ||
            throw(ArgumentError("each edge must be (i, j) or (i, j, weight)"))
        i, j = edge[1], edge[2]
        i isa Integer && j isa Integer ||
            throw(ArgumentError("edge endpoints must be integers"))
        1 <= i <= n && 1 <= j <= n || throw(BoundsError(A, (i, j)))
        i == j && throw(ArgumentError("self-loops are not supported"))
        weight = length(edge) == 2 ? 1.0 : edge[3]
        weight isa Real && isfinite(weight) && weight >= 0 ||
            throw(ArgumentError("edge weights must be finite and nonnegative"))
        w = Float64(weight)
        if !directed && i > j
            i, j = j, i
        end
        if filled[i, j]
            previous = A[i, j]
            if combine === :sum
                w = previous + w
            elseif combine === :max
                w = max(previous, w)
            elseif combine === :min
                w = min(previous, w)
            elseif combine === nothing && previous != w
                throw(ArgumentError(string(
                    "dyad (", edge[1], ", ", edge[2], ") occurs more than once with ",
                    "conflicting weights ", previous, " and ", w,
                    "; pass combine=:sum, :max, :min, or :last, or directed=true ",
                    "if (i, j) and (j, i) are distinct dyads")))
            end
        end
        filled[i, j] = true
        A[i, j] = w
        directed || (A[j, i] = w)
    end
    return A
end

"""
    block_densities(A, result)
    block_densities(A, core_nodes)

Observed densities of the four blocks of the adjacency matrix induced by a
core/periphery partition: the standard reporting companion of a Borgatti-Everett
fit (UCINET's density table). Returns a `NamedTuple`
`(core_core, core_periphery, periphery_core, periphery_periphery)` whose entries
are the mean off-diagonal value of each block (the density for binary input,
the mean weight for weighted input). A block with no dyads has density `NaN`.

`core_periphery` averages `A[i, j]` over core `i` and periphery `j`, and
`periphery_core` over periphery `i` and core `j`; the two coincide for
symmetric `A`. `result` may be any result with `core_nodes` and
`periphery_nodes` fields (`CPResult`, `CPDirectedResult`).
"""
function block_densities(A::AbstractMatrix{<:Real},
    result::Union{CPResult,CPDirectedResult})
    return block_densities(A, result.core_nodes)
end

function block_densities(A::AbstractMatrix{<:Real}, core_nodes::AbstractVector{<:Integer})
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A)
    all(node -> 1 <= node <= n, core_nodes) ||
        throw(ArgumentError("core_nodes must be node indices in 1:$n"))
    in_core = falses(n)
    in_core[core_nodes] .= true
    core = findall(in_core)
    periphery = findall(!, in_core)
    n_core = length(core)
    n_periphery = length(periphery)
    core_core = _block_sum(A, core, core)
    core_periphery = _block_sum(A, core, periphery)
    periphery_core = _block_sum(A, periphery, core)
    periphery_periphery = _block_sum(A, periphery, periphery)
    density(total, dyads) = dyads > 0 ? total / dyads : NaN
    return (
        core_core=density(core_core, n_core * (n_core - 1)),
        core_periphery=density(core_periphery, n_core * n_periphery),
        periphery_core=density(periphery_core, n_periphery * n_core),
        periphery_periphery=density(periphery_periphery, n_periphery * (n_periphery - 1)),
    )
end

function _block_sum(A::AbstractMatrix{<:Real}, rows::Vector{Int}, columns::Vector{Int})
    total = 0.0
    @inbounds for j in columns, i in rows
        i == j && continue
        total += Float64(A[i, j])
    end
    return total
end

function _block_sum(A::SparseMatrixCSC{<:Real}, rows::Vector{Int}, columns::Vector{Int})
    in_rows = falses(size(A, 1))
    in_rows[rows] .= true
    total = 0.0
    row_indices = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in columns, index in nzrange(A, j)
        i = row_indices[index]
        (i != j && in_rows[i]) && (total += Float64(values[index]))
    end
    return total
end

"""
    ideal_cp_matrix(c; discrete=false, cp_block=:ones)

Generate the ideal core-periphery pattern matrix for a given coreness vector.

- Continuous model (`discrete=false`): `Δ[i,j] = c[i] * c[j]`.
- Discrete model, `cp_block=:ones`: `Δ[i,j] = max(c[i], c[j])`, Borgatti and
  Everett (2000) eq. (2): core-periphery dyads are expected to be present.
- Discrete model, `cp_block=:missing`: Borgatti and Everett (2000) eq. (4),
  the model they recommend. `Δ[i,j] = 1` for core-core dyads, `0` for
  periphery-periphery dyads, and `missing` for core-periphery dyads, which are
  excluded from the correlation. Requires a binary `c`; the result is then a
  `Matrix{Union{Missing,Float64}}`.

The diagonal is always `0`.

# Arguments
- `c`: Coreness vector (any real element type, values in [0,1])
- `discrete::Bool`: If true, use the discrete ideal pattern
- `cp_block::Symbol`: `:ones` or `:missing` (discrete model only)

# Returns
- Ideal core-periphery matrix (`Matrix{Float64}`, or
  `Matrix{Union{Missing,Float64}}` for `cp_block=:missing`)
"""
function ideal_cp_matrix(c::AbstractVector{<:Real}; discrete::Bool=false,
    cp_block::Symbol=:ones)
    _check_cp_block(cp_block, discrete, c)
    n = length(c)

    if discrete && cp_block === :missing
        Δ = Matrix{Union{Missing,Float64}}(undef, n, n)
        fill!(Δ, 0.0)
        for i in 1:n, j in 1:n
            i == j && continue
            Δ[i, j] = c[i] == c[j] ? Float64(c[i]) : missing
        end
        return Δ
    end

    Δ = zeros(Float64, n, n)
    if discrete
        for i in 1:n
            for j in 1:n
                if i != j
                    Δ[i, j] = max(Float64(c[i]), Float64(c[j]))
                end
            end
        end
    else
        for i in 1:n
            for j in 1:n
                if i != j
                    Δ[i, j] = Float64(c[i]) * Float64(c[j])
                end
            end
        end
    end

    return Δ
end

function _adjacency_pattern_product(A::AbstractMatrix{<:Real},
    c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    n = length(c)
    sum_a = 0.0
    sum_a2 = 0.0
    sum_ad = 0.0
    if directed
        @inbounds for j in 1:n, i in 1:n
            i == j && continue
            a = Float64(A[i, j])
            d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                           Float64(c[i]) * Float64(c[j])
            sum_a += a
            sum_a2 += a * a
            sum_ad += a * d
        end
    else
        @inbounds for j in 2:n, i in 1:(j - 1)
            a = Float64(A[i, j])
            d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                           Float64(c[i]) * Float64(c[j])
            sum_a += a
            sum_a2 += a * a
            sum_ad += a * d
        end
    end
    return sum_a, sum_a2, sum_ad
end

function _adjacency_pattern_product(A::SparseMatrixCSC{<:Real},
    c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    sum_a = 0.0
    sum_a2 = 0.0
    sum_ad = 0.0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        (directed ? i != j : i < j) || continue
        a = Float64(values[index])
        d = discrete ? max(Float64(c[i]), Float64(c[j])) :
                       Float64(c[i]) * Float64(c[j])
        sum_a += a
        sum_a2 += a * a
        sum_ad += a * d
    end
    return sum_a, sum_a2, sum_ad
end

function _pattern_moments(c::AbstractVector{<:Real}, discrete::Bool, directed::Bool)
    n = length(c)
    count = directed ? n * (n - 1) : n * (n - 1) ÷ 2
    if discrete
        if all(value -> value == 0 || value == 1, c)
            n_periphery = Base.count(iszero, c)
            sum_d = Float64(n * (n - 1) ÷ 2 -
                            n_periphery * (n_periphery - 1) ÷ 2)
            directed && (sum_d *= 2)
            return count, sum_d, sum_d
        end
        ordered = sort!(Float64.(c))
        sum_d = 0.0
        sum_d2 = 0.0
        @inbounds for k in eachindex(ordered)
            multiplicity = k - 1
            value = ordered[k]
            sum_d += multiplicity * value
            sum_d2 += multiplicity * value * value
        end
        directed && (sum_d *= 2; sum_d2 *= 2)
        return count, sum_d, sum_d2
    end

    sum_c = 0.0
    sum_c2 = 0.0
    sum_c4 = 0.0
    @inbounds for value in c
        x = Float64(value)
        x2 = x * x
        sum_c += x
        sum_c2 += x2
        sum_c4 += x2 * x2
    end
    sum_d = (sum_c * sum_c - sum_c2) / 2
    sum_d2 = (sum_c2 * sum_c2 - sum_c4) / 2
    directed && (sum_d *= 2; sum_d2 *= 2)
    return count, sum_d, sum_d2
end

@inline function _pearson_from_sums(count::Int, sum_a::Float64, sum_a2::Float64,
    sum_d::Float64, sum_d2::Float64, sum_ad::Float64)
    count >= 2 || return 0.0
    inv_count = inv(Float64(count))
    variance_a = sum_a2 - sum_a * sum_a * inv_count
    variance_d = sum_d2 - sum_d * sum_d * inv_count
    (variance_a > 0.0 && variance_d > 0.0) || return 0.0
    covariance = sum_ad - sum_a * sum_d * inv_count
    return covariance / sqrt(variance_a * variance_d)
end

"""
    core_quality(A, c; discrete=false, directed=false, cp_block=:ones)

Compute the Pearson correlation between adjacency values and the ideal
core-periphery pattern (see [`ideal_cp_matrix`](@ref)) without materializing
dyad vectors. Undirected inputs use the strict upper triangle; `directed=true`
uses all ordered off-diagonal dyads. With `discrete=true` and
`cp_block=:missing` (Borgatti-Everett eq. 4), only core-core and
periphery-periphery dyads enter the correlation and `c` must be binary.
Returns `0.0` when either pattern has zero variance.
"""
function core_quality(A::AbstractMatrix{<:Real}, c::AbstractVector{<:Real};
    discrete::Bool=false,
    directed::Bool=false,
    cp_block::Symbol=:ones)
    A = _unwrap_adjacency(A)
    n = length(c)
    size(A) == (n, n) ||
        throw(DimensionMismatch("A must be a square matrix with size ($(n), $(n))"))
    _check_cp_block(cp_block, discrete, c)
    if discrete && cp_block === :missing
        state = _discrete_state(A, c, directed, cp_block)
        return _discrete_quality(state)
    end
    count, sum_d, sum_d2 = _pattern_moments(c, discrete, directed)
    sum_a, sum_a2, sum_ad = _adjacency_pattern_product(A, c, discrete, directed)
    return _pearson_from_sums(count, sum_a, sum_a2, sum_d, sum_d2, sum_ad)
end

function _check_cp_block(cp_block::Symbol, discrete::Bool, c::AbstractVector{<:Real})
    cp_block in (:ones, :missing) ||
        throw(ArgumentError("cp_block must be :ones or :missing"))
    if cp_block === :missing
        discrete || throw(ArgumentError(
            "cp_block=:missing is defined for the discrete model only (discrete=true)"))
        all(value -> value == 0 || value == 1, c) ||
            throw(ArgumentError("cp_block=:missing requires a binary coreness vector"))
    end
    return cp_block
end

# Sufficient statistics of the discrete Pearson objective for a binary core
# vector. They cover both ideal patterns (`cp_block=:ones`, max(c_i, c_j);
# `cp_block=:missing`, core-periphery dyads excluded) over either the strict
# upper triangle (undirected) or all ordered off-diagonal dyads (directed).
# Flipping one node's role updates them in O(degree) through
# `_node_block_links`, which is what the greedy optimisers rely on.
mutable struct _DiscreteState
    directed::Bool
    cp_block::Symbol
    n::Int
    n_core::Int
    total_weight::Float64
    total_weight2::Float64
    core_weight::Float64
    core_weight2::Float64
    periphery_weight::Float64
    periphery_weight2::Float64
end

_dyads(k::Int, directed::Bool) = directed ? k * (k - 1) : k * (k - 1) ÷ 2

function _discrete_state(A::AbstractMatrix{<:Real}, c::AbstractVector{<:Real},
    directed::Bool, cp_block::Symbol)
    n = length(c)
    total, total2, cc, cc2, pp, pp2 = _block_weight_sums(A, c, directed)
    return _DiscreteState(directed, cp_block, n, count(!iszero, c),
        total, total2, cc, cc2, pp, pp2)
end

function _block_weight_sums(A::AbstractMatrix{<:Real}, c::AbstractVector{<:Real},
    directed::Bool)
    n = length(c)
    total = total2 = cc = cc2 = pp = pp2 = 0.0
    @inbounds for j in 1:n, i in 1:n
        i == j && continue
        (directed || i < j) || continue
        a = Float64(A[i, j])
        iszero(a) && continue
        a2 = a * a
        total += a
        total2 += a2
        if c[i] == c[j]
            if iszero(c[i])
                pp += a
                pp2 += a2
            else
                cc += a
                cc2 += a2
            end
        end
    end
    return total, total2, cc, cc2, pp, pp2
end

function _block_weight_sums(A::SparseMatrixCSC{<:Real}, c::AbstractVector{<:Real},
    directed::Bool)
    total = total2 = cc = cc2 = pp = pp2 = 0.0
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for j in axes(A, 2), index in nzrange(A, j)
        i = rows[index]
        (directed ? i != j : i < j) || continue
        a = Float64(values[index])
        iszero(a) && continue
        a2 = a * a
        total += a
        total2 += a2
        if c[i] == c[j]
            if iszero(c[i])
                pp += a
                pp2 += a2
            else
                cc += a
                cc2 += a2
            end
        end
    end
    return total, total2, cc, cc2, pp, pp2
end

@inline function _discrete_quality(directed::Bool, cp_block::Symbol, n::Int, n_core::Int,
    total::Float64, total2::Float64, cc::Float64, cc2::Float64,
    pp::Float64, pp2::Float64)
    n_periphery = n - n_core
    if cp_block === :ones
        count = _dyads(n, directed)
        sum_d = Float64(count - _dyads(n_periphery, directed))
        return _pearson_from_sums(count, total, total2, sum_d, sum_d, total - pp)
    end
    count = _dyads(n_core, directed) + _dyads(n_periphery, directed)
    sum_d = Float64(_dyads(n_core, directed))
    return _pearson_from_sums(count, cc + pp, cc2 + pp2, sum_d, sum_d, cc)
end

_discrete_quality(s::_DiscreteState) = _discrete_quality(
    s.directed, s.cp_block, s.n, s.n_core, s.total_weight, s.total_weight2,
    s.core_weight, s.core_weight2, s.periphery_weight, s.periphery_weight2)

# Weight (and squared weight) of the dyads joining `node` to the current core
# and to the current periphery. `At` is `A` for undirected input and the
# transpose of `A` for directed input, in which case both orientations count.
function _node_block_links(A::AbstractMatrix{<:Real}, At::AbstractMatrix{<:Real},
    node::Int, c::AbstractVector{<:Real})
    to_core = to_core2 = to_periphery = to_periphery2 = 0.0
    @inbounds for neighbor in axes(A, 1)
        neighbor == node && continue
        a = Float64(A[neighbor, node])
        if At !== A
            b = Float64(A[node, neighbor])
            if iszero(c[neighbor])
                to_periphery += a + b
                to_periphery2 += a * a + b * b
            else
                to_core += a + b
                to_core2 += a * a + b * b
            end
        else
            if iszero(c[neighbor])
                to_periphery += a
                to_periphery2 += a * a
            else
                to_core += a
                to_core2 += a * a
            end
        end
    end
    return to_core, to_core2, to_periphery, to_periphery2
end

function _node_block_links(A::SparseMatrixCSC{<:Real}, At::SparseMatrixCSC{<:Real},
    node::Int, c::AbstractVector{<:Real})
    to_core = to_core2 = to_periphery = to_periphery2 = 0.0
    for matrix in (At === A ? (A,) : (A, At))
        rows = rowvals(matrix)
        values = nonzeros(matrix)
        @inbounds for index in nzrange(matrix, node)
            neighbor = rows[index]
            neighbor == node && continue
            a = Float64(values[index])
            if iszero(c[neighbor])
                to_periphery += a
                to_periphery2 += a * a
            else
                to_core += a
                to_core2 += a * a
            end
        end
    end
    return to_core, to_core2, to_periphery, to_periphery2
end

# Objective value after flipping a node whose block links are given; `to_core`
# means the node moves from the periphery into the core.
@inline function _flip_quality(s::_DiscreteState, to_core_role::Bool,
    to_core::Float64, to_core2::Float64, to_periphery::Float64, to_periphery2::Float64)
    sign = to_core_role ? 1.0 : -1.0
    return _discrete_quality(s.directed, s.cp_block, s.n,
        s.n_core + (to_core_role ? 1 : -1),
        s.total_weight, s.total_weight2,
        s.core_weight + sign * to_core, s.core_weight2 + sign * to_core2,
        s.periphery_weight - sign * to_periphery,
        s.periphery_weight2 - sign * to_periphery2)
end

@inline function _apply_flip!(s::_DiscreteState, to_core_role::Bool,
    to_core::Float64, to_core2::Float64, to_periphery::Float64, to_periphery2::Float64)
    sign = to_core_role ? 1.0 : -1.0
    s.n_core += to_core_role ? 1 : -1
    s.core_weight += sign * to_core
    s.core_weight2 += sign * to_core2
    s.periphery_weight -= sign * to_periphery
    s.periphery_weight2 -= sign * to_periphery2
    return s
end

_transposed_adjacency(A::AbstractMatrix{<:Real}) = Matrix{Float64}(transpose(A))
_transposed_adjacency(A::SparseMatrixCSC{<:Real}) = sparse(transpose(A))

function _axpy_column!(destination::AbstractVector{Float64},
    A::AbstractMatrix{<:Real}, column::Int, scale::Float64)
    @inbounds for row in axes(A, 1)
        destination[row] += scale * Float64(A[row, column])
    end
    return destination
end

function _axpy_column!(destination::AbstractVector{Float64},
    A::SparseMatrixCSC{<:Real}, column::Int, scale::Float64)
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, column)
        destination[rows[index]] += scale * Float64(values[index])
    end
    return destination
end

"""
    _partition(coreness; rule=:median, threshold=nothing, core_size=nothing)

Split nodes into sorted `core_nodes` and `periphery_nodes` index vectors from a
coreness vector.

- `:median`: nodes with coreness strictly greater than the median are core, so
  nodes tied at the median are periphery (a star's leaves never join its hub).
- `:none`: empty core, every node in the periphery.
- `:threshold`: nodes with coreness strictly greater than `threshold` are core.
- `:core_size`: the `core_size` nodes with the highest coreness are core; ties
  are broken by lower node index.

`threshold` may only be given with `rule=:threshold`, and `core_size` only with
`rule=:core_size`.
"""
function _partition(coreness::AbstractVector{<:Real};
    rule::Symbol=:median,
    threshold::Union{Nothing,Real}=nothing,
    core_size::Union{Nothing,Integer}=nothing)
    n = length(coreness)
    rule in (:median, :none, :threshold, :core_size) ||
        throw(ArgumentError("partition must be :median, :none, :threshold, or :core_size"))
    threshold === nothing || rule === :threshold ||
        throw(ArgumentError("threshold is only used with partition=:threshold"))
    core_size === nothing || rule === :core_size ||
        throw(ArgumentError("core_size is only used with partition=:core_size"))

    core_nodes = if rule === :median
        n == 0 ? Int[] : findall(>(median(coreness)), coreness)
    elseif rule === :none
        Int[]
    elseif rule === :threshold
        threshold === nothing &&
            throw(ArgumentError("partition=:threshold requires the threshold keyword"))
        isfinite(threshold) || throw(ArgumentError("threshold must be finite"))
        findall(>(threshold), coreness)
    else
        core_size === nothing &&
            throw(ArgumentError("partition=:core_size requires the core_size keyword"))
        0 <= core_size <= n ||
            throw(ArgumentError("core_size must be between 0 and the number of nodes"))
        order = sortperm(coreness; rev=true, alg=Base.Sort.MergeSort)
        sort!(order[1:core_size])
    end
    in_core = falses(n)
    in_core[core_nodes] .= true
    periphery_nodes = findall(!, in_core)
    return core_nodes, periphery_nodes
end
