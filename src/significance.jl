# Null-model samplers and the Monte Carlo significance test.
#
# Every null model has an in-place sampler that draws into a `_NullWorkspace`.
# A workspace holds either a dense `Matrix{Float64}` or a
# `SparseMatrixCSC{Float64,Int}`, chosen from the storage of the observed
# network, so that detectors evaluate null samples with the same kernels
# (dense or sparse) they use for the observed network. Dense and sparse
# samplers consume the RNG stream identically and therefore produce the same
# graphs for the same seed.

"""
Diagnostics from one degree-preserving configuration-null sample.
"""
struct ConfigurationSwapStats
    requested::Int
    accepted::Int
    attempts::Int
    max_attempts::Int
end

_swap_complete(stats::ConfigurationSwapStats) = stats.accepted == stats.requested

function _configuration_attempt_budget(n_swaps::Int,
    max_attempts::Union{Nothing,Int})
    max_attempts === nothing && return max(1_000, 100 * n_swaps)
    max_attempts >= 0 || throw(ArgumentError("max_swap_attempts must be nonnegative"))
    return max_attempts
end

function _handle_swap_shortfall(stats::ConfigurationSwapStats, policy::Symbol)
    _swap_complete(stats) && return nothing
    message = "configuration null accepted $(stats.accepted) of " *
        "$(stats.requested) requested swaps in $(stats.attempts) attempts " *
        "(budget $(stats.max_attempts)). The degree sequence may admit no " *
        "valid double-edge swap at all (stars and other threshold graphs have a " *
        "unique simple realization), in which case a degree-preserving null " *
        "cannot be sampled: use null_model=:er, raise max_swap_attempts, or pass " *
        "swap_shortfall=:accept (or :warn) to keep incomplete samples"
    if policy === :error
        throw(ConvergenceError(message))
    elseif policy === :warn
        @warn message
    elseif policy !== :accept
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    end
    return nothing
end

_isbinary(A::AbstractMatrix{<:Real}) = all(a -> iszero(a) || a == one(a), A)
_isbinary(A::SparseMatrixCSC{<:Real}) = all(a -> iszero(a) || a == one(a), nonzeros(A))

function _sample_er!(sampled::Matrix{Float64}, density::Float64, rng::AbstractRNG)
    n = size(sampled, 1)
    fill!(sampled, 0.0)
    for j in 2:n, i in 1:(j - 1)
        if rand(rng) < density
            sampled[i, j] = sampled[j, i] = 1.0
        end
    end
    return sampled
end

function _sample_er(A::AbstractMatrix{<:Real}, rng::AbstractRNG)
    n = _validate_adjacency(A; symmetric=true, binary=true)
    possible = n * (n - 1) ÷ 2
    edge_count = length(_edge_list(A))
    density = possible == 0 ? 0.0 : edge_count / possible
    return _sample_er!(zeros(Float64, n, n), density, rng)
end

function _edge_list(A::AbstractMatrix{<:Real})
    n = size(A, 1)
    return [(i, j) for i in 1:(n - 1) for j in (i + 1):n if !iszero(A[i, j])]
end


function _edge_list(A::SparseMatrixCSC{<:Real,<:Integer})
    rows, columns, values = findnz(A)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, length(values) ÷ 2)
    for index in eachindex(values)
        i, j = rows[index], columns[index]
        i < j && !iszero(values[index]) && push!(edges, (i, j))
    end
    # Match the dense lexicographic edge order so an explicit RNG produces the
    # same switch chain independently of storage representation.
    return sort!(edges)
end

function _arc_list(A::AbstractMatrix{<:Real})
    n = size(A, 1)
    return [(i, j) for i in 1:n for j in 1:n if i != j && !iszero(A[i, j])]
end

function _arc_list(A::SparseMatrixCSC{<:Real,<:Integer})
    rows, columns, values = findnz(A)
    arcs = Tuple{Int,Int}[]
    sizehint!(arcs, length(values))
    for index in eachindex(values)
        i, j = rows[index], columns[index]
        i != j && !iszero(values[index]) && push!(arcs, (i, j))
    end
    return sort!(arcs)
end

# Edge-existence state consulted by the switching chains. The dense state is
# the sampled adjacency matrix itself; the sparse state is a hash set of
# endpoint pairs, so a configuration sample never materializes an n×n matrix.
struct _DenseEdgeState
    matrix::Matrix{Float64}
end

struct _SetEdgeState
    edges::Set{Tuple{Int,Int}}
end

@inline _has_edge(state::_DenseEdgeState, u::Int, v::Int) = !iszero(state.matrix[u, v])
@inline _has_edge(state::_SetEdgeState, u::Int, v::Int) = (u, v) in state.edges

@inline function _set_edge!(state::_DenseEdgeState, u::Int, v::Int, present::Bool,
    symmetric::Bool)
    value = present ? 1.0 : 0.0
    state.matrix[u, v] = value
    symmetric && (state.matrix[v, u] = value)
    return nothing
end

@inline function _set_edge!(state::_SetEdgeState, u::Int, v::Int, present::Bool,
    ::Bool)
    present ? push!(state.edges, (u, v)) : delete!(state.edges, (u, v))
    return nothing
end

# Simple-graph double-edge switch chain (undirected). `edges` holds the current
# edge list as (min, max) pairs and is updated in place. Rejections count as
# attempts, so the chain is a lazy Markov chain with the uniform stationary
# distribution over simple graphs with the given degree sequence. Returns the
# accepted and attempted counts. The RNG is consumed identically by the dense
# and sparse states.
function _undirected_switch_chain!(state, edges::Vector{Tuple{Int,Int}},
    rng::AbstractRNG, n_swaps::Int, attempt_budget::Int)
    accepted = 0
    attempts = 0
    while length(edges) >= 2 && accepted < n_swaps && attempts < attempt_budget
        attempts += 1
        first_index = rand(rng, eachindex(edges))
        second_index = rand(rng, eachindex(edges))
        first_index == second_index && continue
        a, b = edges[first_index]
        c, d = edges[second_index]
        (a == c || a == d || b == c || b == d) && continue
        rand(rng, Bool) && ((c, d) = (d, c))
        u, v = minmax(a, d)
        x, y = minmax(c, b)
        (u == v || x == y || (u == x && v == y)) && continue
        (_has_edge(state, u, v) || _has_edge(state, x, y)) && continue

        _set_edge!(state, minmax(a, b)..., false, true)
        _set_edge!(state, minmax(c, d)..., false, true)
        _set_edge!(state, u, v, true, true)
        _set_edge!(state, x, y, true, true)
        edges[first_index] = (u, v)
        edges[second_index] = (x, y)
        accepted += 1
    end
    return accepted, attempts
end

# Directed arc switch chain preserving in- and out-degrees of a simple digraph.
function _directed_switch_chain!(state, arcs::Vector{Tuple{Int,Int}},
    rng::AbstractRNG, n_swaps::Int, attempt_budget::Int)
    accepted = 0
    attempts = 0
    while length(arcs) >= 2 && accepted < n_swaps && attempts < attempt_budget
        attempts += 1
        first_index = rand(rng, eachindex(arcs))
        second_index = rand(rng, eachindex(arcs))
        first_index == second_index && continue
        a, b = arcs[first_index]
        c, d = arcs[second_index]
        (a == c || b == d || a == d || c == b) && continue
        (_has_edge(state, a, d) || _has_edge(state, c, b)) && continue

        _set_edge!(state, a, b, false, false)
        _set_edge!(state, c, d, false, false)
        _set_edge!(state, a, d, true, false)
        _set_edge!(state, c, b, true, false)
        arcs[first_index] = (a, d)
        arcs[second_index] = (c, b)
        accepted += 1
    end
    return accepted, attempts
end

function _sample_directed_configuration!(
    sampled::Matrix{Float64},
    A::AbstractMatrix{<:Real},
    arcs::Vector{Tuple{Int,Int}},
    base_arcs::Vector{Tuple{Int,Int}},
    rng::AbstractRNG,
    n_swaps::Int;
    max_attempts::Union{Nothing,Int}=nothing,
    on_shortfall::Symbol=:error,
    stats_ref=nothing,
)
    n_swaps >= 0 || throw(ArgumentError("n_swaps must be nonnegative"))
    on_shortfall in (:error, :warn, :accept) ||
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    attempt_budget = _configuration_attempt_budget(n_swaps, max_attempts)
    copyto!(sampled, A)
    resize!(arcs, length(base_arcs))
    copyto!(arcs, base_arcs)
    accepted, attempts = _directed_switch_chain!(
        _DenseEdgeState(sampled), arcs, rng, n_swaps, attempt_budget)
    stats = ConfigurationSwapStats(n_swaps, accepted, attempts, attempt_budget)
    stats_ref === nothing || (stats_ref[] = stats)
    _handle_swap_shortfall(stats, on_shortfall)
    return sampled
end

function _sample_directed_configuration(A::AbstractMatrix{<:Real}, rng::AbstractRNG;
    n_swaps::Union{Nothing,Int}=nothing,
    max_swap_attempts::Union{Nothing,Int}=nothing,
    swap_shortfall::Symbol=:error,
    return_stats::Bool=false)
    _validate_adjacency(A; binary=true)
    base_arcs = _arc_list(A)
    swaps = n_swaps === nothing ? max(10, 10 * length(base_arcs)) : n_swaps
    sampled = Matrix{Float64}(undef, size(A)...)
    stats = Ref{ConfigurationSwapStats}()
    _sample_directed_configuration!(
        sampled, A, copy(base_arcs), base_arcs, rng, swaps;
        max_attempts=max_swap_attempts,
        on_shortfall=swap_shortfall,
        stats_ref=stats,
    )
    return return_stats ? (network=sampled, stats=stats[]) : sampled
end

function _sample_weight_permutation!(sampled::Matrix{Float64},
    base_edges::Vector{Tuple{Int,Int}}, base_weights::Vector{Float64},
    weight_buffer::Vector{Float64}, rng::AbstractRNG; symmetric::Bool)
    fill!(sampled, 0.0)
    resize!(weight_buffer, length(base_weights))
    copyto!(weight_buffer, base_weights)
    shuffle!(rng, weight_buffer)
    @inbounds for index in eachindex(base_edges, weight_buffer)
        source, target = base_edges[index]
        weight = weight_buffer[index]
        sampled[source, target] = weight
        symmetric && (sampled[target, source] = weight)
    end
    return sampled
end

function _sample_weight_permutation(A::AbstractMatrix{<:Real}, rng::AbstractRNG)
    n = _validate_adjacency(A)
    symmetric = _issymmetric(A)
    edges = symmetric ? _edge_list(A) : _arc_list(A)
    weights = Float64[A[source, target] for (source, target) in edges]
    sampled = zeros(Float64, n, n)
    return _sample_weight_permutation!(
        sampled, edges, weights, similar(weights), rng; symmetric=symmetric)
end


function _sample_configuration!(
    sampled::Matrix{Float64},
    A::AbstractMatrix{<:Real},
    edges::Vector{Tuple{Int,Int}},
    base_edges::Vector{Tuple{Int,Int}},
    rng::AbstractRNG,
    n_swaps::Int,
    ;
    max_attempts::Union{Nothing,Int}=nothing,
    on_shortfall::Symbol=:error,
    stats_ref=nothing,
)
    n_swaps >= 0 || throw(ArgumentError("n_swaps must be nonnegative"))
    on_shortfall in (:error, :warn, :accept) ||
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    attempt_budget = _configuration_attempt_budget(n_swaps, max_attempts)
    copyto!(sampled, A)
    resize!(edges, length(base_edges))
    copyto!(edges, base_edges)
    accepted, attempts = _undirected_switch_chain!(
        _DenseEdgeState(sampled), edges, rng, n_swaps, attempt_budget)
    stats = ConfigurationSwapStats(n_swaps, accepted, attempts, attempt_budget)
    stats_ref === nothing || (stats_ref[] = stats)
    _handle_swap_shortfall(stats, on_shortfall)
    return sampled
end

function _sample_configuration(
    A::AbstractMatrix{<:Real},
    rng::AbstractRNG;
    n_swaps::Union{Nothing,Int}=nothing,
    max_swap_attempts::Union{Nothing,Int}=nothing,
    swap_shortfall::Symbol=:error,
    return_stats::Bool=false,
)
    _validate_adjacency(A; symmetric=true, binary=true)
    base_edges = _edge_list(A)
    swaps = n_swaps === nothing ? max(10, 10 * length(base_edges)) : n_swaps
    swaps >= 0 || throw(ArgumentError("n_swaps must be nonnegative"))
    sampled = Matrix{Float64}(undef, size(A)...)
    stats = Ref{ConfigurationSwapStats}()
    _sample_configuration!(
        sampled, A, copy(base_edges), base_edges, rng, swaps;
        max_attempts=max_swap_attempts,
        on_shortfall=swap_shortfall,
        stats_ref=stats,
    )
    return return_stats ? (network=sampled, stats=stats[]) : sampled
end

# ---------------------------------------------------------------------------
# Per-worker workspaces
# ---------------------------------------------------------------------------

# Everything one worker needs to draw null samples: the output network (dense
# or sparse), the switching-chain edge list, a weight buffer, the sparse
# edge-existence set, and counting-sort scratch for CSC assembly.
struct _NullWorkspace{M<:AbstractMatrix{Float64}}
    network::M
    edges::Vector{Tuple{Int,Int}}
    weights::Vector{Float64}
    edge_set::Set{Tuple{Int,Int}}
    row_count::Vector{Int}
    order::Vector{Int}
end

function _dense_workspace(n::Int, base_edges::Vector{Tuple{Int,Int}},
    base_weights::Vector{Float64})
    return _NullWorkspace(
        Matrix{Float64}(undef, n, n), copy(base_edges), copy(base_weights),
        Set{Tuple{Int,Int}}(), Int[], Int[])
end

function _sparse_workspace(n::Int, base_edges::Vector{Tuple{Int,Int}},
    base_weights::Vector{Float64})
    network = SparseMatrixCSC{Float64,Int}(n, n, ones(Int, n + 1), Int[], Float64[])
    edge_set = Set{Tuple{Int,Int}}()
    sizehint!(edge_set, max(1, length(base_edges)))
    return _NullWorkspace(
        network, copy(base_edges), copy(base_weights), edge_set,
        Vector{Int}(undef, n + 1), Int[])
end

# Fill the workspace's CSC matrix in place with the entries `(s, t, w)` of
# `edges`/`weights` (unit weights when `weights === nothing`), mirrored when
# `symmetric`. A two-pass counting sort (rows, then columns) leaves the row
# indices sorted within every column, so the result is a valid CSC matrix in
# O(entries + n) time without allocation beyond buffer growth.
function _assemble_csc!(ws::_NullWorkspace{SparseMatrixCSC{Float64,Int}},
    edges::Vector{Tuple{Int,Int}}, weights::Union{Nothing,Vector{Float64}},
    symmetric::Bool)
    network = ws.network
    n = size(network, 1)
    entries = symmetric ? 2 * length(edges) : length(edges)
    colptr = SparseArrays.getcolptr(network)
    rowval = rowvals(network)
    nzval = nonzeros(network)
    resize!(rowval, entries)
    resize!(nzval, entries)
    row_start = ws.row_count
    order = ws.order
    resize!(order, entries)

    # Entry k refers to edge e = (k + 1) ÷ 2 (symmetric) or e = k, oriented
    # (s, t) for odd k and (t, s) for even k when mirrored.
    @inline entry_row(k) = symmetric ?
        (isodd(k) ? edges[(k + 1) >> 1][1] : edges[k >> 1][2]) : edges[k][1]
    @inline entry_col(k) = symmetric ?
        (isodd(k) ? edges[(k + 1) >> 1][2] : edges[k >> 1][1]) : edges[k][2]
    @inline entry_weight(k) = weights === nothing ? 1.0 :
        (symmetric ? weights[(k + 1) >> 1] : weights[k])

    # Pass 1: order entries by row.
    fill!(row_start, 0)
    @inbounds for k in 1:entries
        row_start[entry_row(k) + 1] += 1
    end
    row_start[1] = 1
    @inbounds for r in 1:n
        row_start[r + 1] += row_start[r]
    end
    @inbounds for k in 1:entries
        r = entry_row(k)
        order[row_start[r]] = k
        row_start[r] += 1
    end

    # Pass 2: column pointers, then scatter in row order.
    fill!(colptr, 0)
    @inbounds for k in 1:entries
        colptr[entry_col(k) + 1] += 1
    end
    colptr[1] = 1
    @inbounds for c in 1:n
        colptr[c + 1] += colptr[c]
    end
    next = row_start
    copyto!(next, colptr)
    @inbounds for position in 1:entries
        k = order[position]
        c = entry_col(k)
        slot = next[c]
        rowval[slot] = entry_row(k)
        nzval[slot] = entry_weight(k)
        next[c] = slot + 1
    end
    return network
end

# Parameters of one null model, shared read-only by all workers.
struct _NullSpec{T<:AbstractMatrix{<:Real}}
    null_model::Symbol
    A::T
    base_edges::Vector{Tuple{Int,Int}}
    base_weights::Vector{Float64}
    density::Float64
    swaps::Int
    attempt_budget::Int
    swap_shortfall::Symbol
    symmetric::Bool
end

# Draw one null sample into the workspace. Returns the switching statistics for
# configuration nulls and `nothing` otherwise. Dense and sparse methods consume
# the RNG identically and produce the same graph.
function _draw_null!(ws::_NullWorkspace{Matrix{Float64}}, spec::_NullSpec,
    rng::AbstractRNG)
    if spec.null_model === :er
        _sample_er!(ws.network, spec.density, rng)
        return nothing
    elseif spec.null_model === :configuration
        stats = Ref{ConfigurationSwapStats}()
        _sample_configuration!(
            ws.network, spec.A, ws.edges, spec.base_edges, rng, spec.swaps;
            max_attempts=spec.attempt_budget,
            on_shortfall=spec.swap_shortfall,
            stats_ref=stats,
        )
        return stats[]
    elseif spec.null_model === :directed_configuration
        stats = Ref{ConfigurationSwapStats}()
        _sample_directed_configuration!(
            ws.network, spec.A, ws.edges, spec.base_edges, rng, spec.swaps;
            max_attempts=spec.attempt_budget,
            on_shortfall=spec.swap_shortfall,
            stats_ref=stats,
        )
        return stats[]
    else
        _sample_weight_permutation!(
            ws.network, spec.base_edges, spec.base_weights, ws.weights, rng;
            symmetric=spec.symmetric,
        )
        return nothing
    end
end

function _reset_chain!(ws::_NullWorkspace, base_edges::Vector{Tuple{Int,Int}})
    resize!(ws.edges, length(base_edges))
    copyto!(ws.edges, base_edges)
    empty!(ws.edge_set)
    for edge in base_edges
        push!(ws.edge_set, edge)
    end
    return ws
end

function _draw_null!(ws::_NullWorkspace{SparseMatrixCSC{Float64,Int}},
    spec::_NullSpec, rng::AbstractRNG)
    n = size(ws.network, 1)
    if spec.null_model === :er
        edges = ws.edges
        empty!(edges)
        for j in 2:n, i in 1:(j - 1)
            rand(rng) < spec.density && push!(edges, (i, j))
        end
        _assemble_csc!(ws, edges, nothing, true)
        return nothing
    elseif spec.null_model === :configuration
        _reset_chain!(ws, spec.base_edges)
        accepted, attempts = _undirected_switch_chain!(
            _SetEdgeState(ws.edge_set), ws.edges, rng, spec.swaps, spec.attempt_budget)
        stats = ConfigurationSwapStats(spec.swaps, accepted, attempts, spec.attempt_budget)
        _handle_swap_shortfall(stats, spec.swap_shortfall)
        _assemble_csc!(ws, ws.edges, nothing, true)
        return stats
    elseif spec.null_model === :directed_configuration
        _reset_chain!(ws, spec.base_edges)
        accepted, attempts = _directed_switch_chain!(
            _SetEdgeState(ws.edge_set), ws.edges, rng, spec.swaps, spec.attempt_budget)
        stats = ConfigurationSwapStats(spec.swaps, accepted, attempts, spec.attempt_budget)
        _handle_swap_shortfall(stats, spec.swap_shortfall)
        _assemble_csc!(ws, ws.edges, nothing, false)
        return stats
    else
        resize!(ws.weights, length(spec.base_weights))
        copyto!(ws.weights, spec.base_weights)
        shuffle!(rng, ws.weights)
        _assemble_csc!(ws, spec.base_edges, ws.weights, spec.symmetric)
        return nothing
    end
end

# ---------------------------------------------------------------------------
# Sample scheduling
# ---------------------------------------------------------------------------

# Run `evaluate!(sample, workspace)` for every sample index. Each worker owns
# the workspace it creates, every sample index is delivered exactly once, and
# `evaluate!` writes only to sample-indexed slots, so the result is independent
# of the schedule. All names assigned inside task bodies are declared `local`:
# a bare assignment inside a closure to a name that also exists in the
# enclosing function would make it a shared, boxed function variable (the
# 0.2.0 greedy race).
function _run_null_samples!(evaluate!, make_workspace, n_samples::Int,
    schedule::Symbol)
    if schedule === :serial
        local serial_workspace = make_workspace()
        for sample in 1:n_samples
            evaluate!(sample, serial_workspace)
        end
    elseif schedule === :greedy
        workers = min(Threads.nthreads(), n_samples)
        channel = Channel{Int}(workers; spawn=true) do ch
            for sample in 1:n_samples
                put!(ch, sample)
            end
        end
        tasks = map(1:workers) do _
            Threads.@spawn begin
                local task_workspace = make_workspace()
                for sample in channel
                    evaluate!(sample, task_workspace)
                end
            end
        end
        try
            foreach(wait, tasks)
        finally
            close(channel)
        end
    elseif schedule === :static
        workers = min(Threads.nthreads(), n_samples)
        workspaces = [make_workspace() for _ in 1:workers]
        Threads.@threads :static for worker in 1:workers
            for sample in worker:workers:n_samples
                evaluate!(sample, workspaces[worker])
            end
        end
    else
        throw(ArgumentError("unknown schedule $schedule"))
    end
    return nothing
end

"""
    cp_significance(A, detector; null_model=:configuration, n_samples=999,
                    alpha=0.05, rng=Random.default_rng(),
                    detector_kwargs=NamedTuple(), pass_rng=true, n_swaps=nothing,
                    max_swap_attempts=nothing, swap_shortfall=:error,
                    threaded=false, thread_schedule=:auto)

Test whether a detector's observed quality exceeds qualities obtained from
random graphs. `detector` must return a result with a finite `quality` field.
The `:configuration` null uses undirected degree-preserving double-edge swaps;
`:directed_configuration` preserves binary in- and out-degrees through directed
arc switches; `:weight_permutation` preserves topology and the edge-weight
multiset; and `:er` samples an undirected Erdős-Rényi graph with observed density.

The returned [`CPSignificanceResult`](@ref) carries the plus-one Monte Carlo
`pvalue` together with effect sizes: `null_mean`, `null_std` (sample standard
deviation, `NaN` for a single sample), and `zscore`, the observed quality minus
the null mean in null standard deviations (`NaN` when the null distribution is
degenerate).

Null samples are drawn in the storage format of `A`: a `SparseMatrixCSC`
input is sampled directly into compressed sparse column form, so detectors use
their sparse kernels and memory stays `O(workers · (n + m))` instead of
`O(workers · n²)`. Dense and sparse inputs of the same network give identical
results for the same `rng`.

Configuration samples must complete all requested swaps by default; a shortfall
raises [`ConvergenceError`](@ref). Use `max_swap_attempts` to set the proposal
budget and `swap_shortfall=:warn` or `:accept` to retain incomplete samples
explicitly. Per-sample accepted swaps and attempt counts are returned in
`null_diagnostics` (a [`NullDiagnostics`](@ref)).

With `pass_rng=true` (the default) a fresh child RNG, derived from `rng`, is
passed as the `rng` keyword to every detector call whose method accepts one,
so stochastic detectors produce reproducible null distributions. Detectors
without an `rng` keyword are called without it; for detectors defined outside
CorePeriphery a warning is emitted once, because their null distribution is
then not reproducible from `rng`. Pass `pass_rng=false` to never forward it.

With `threaded=true`, null samples are distributed over Julia threads. Sample
seeds are generated serially before execution and every worker owns its
sampling workspace, so serial and threaded runs are bit-for-bit identical for
deterministic detectors and for stochastic detectors receiving `rng`.
Detectors must be thread-safe. `thread_schedule=:auto` (and `:greedy`) uses a
work queue with one task per thread, which balances irregular detector
runtimes; `:static` assigns samples round-robin to one task per thread and
cannot be used from inside another threaded region.
"""
function cp_significance(
    A::AbstractMatrix{<:Real},
    detector::F;
    null_model::Symbol=:configuration,
    n_samples::Int=999,
    alpha::Real=0.05,
    rng::AbstractRNG=Random.default_rng(),
    detector_kwargs::NamedTuple=NamedTuple(),
    pass_rng::Bool=true,
    n_swaps::Union{Nothing,Int}=nothing,
    max_swap_attempts::Union{Nothing,Int}=nothing,
    swap_shortfall::Symbol=:error,
    threaded::Bool=false,
    thread_schedule::Symbol=:auto,
) where {F}
    null_model in (:configuration, :directed_configuration, :weight_permutation, :er) ||
        throw(ArgumentError(
            "null_model must be :configuration, :directed_configuration, " *
            ":weight_permutation, or :er"))
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A; symmetric=null_model in (:configuration, :er))
    if null_model !== :weight_permutation && !_isbinary(A)
        preserved = null_model === :er ? "the edge count" :
            null_model === :configuration ? "the degree sequence" :
            "the in- and out-degree sequences"
        throw(ArgumentError(
            "null_model=$(repr(null_model)) requires a binary adjacency matrix " *
            "because it preserves $preserved of a simple graph; for a weighted " *
            "network use null_model=:weight_permutation, which keeps the topology " *
            "and permutes the observed weights, or binarize the network first"))
    end
    n_samples >= 1 || throw(ArgumentError("n_samples must be positive"))
    0 < alpha < 1 || throw(ArgumentError("alpha must be in (0, 1)"))
    alpha_value = Float64(alpha)
    swap_shortfall in (:error, :warn, :accept) ||
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    thread_schedule in (:auto, :static, :greedy) ||
        throw(ArgumentError("thread_schedule must be :auto, :static, or :greedy"))

    forward_rng = pass_rng && hasmethod(detector, Tuple{typeof(A)}, (:rng,))
    if pass_rng && !forward_rng && parentmodule(typeof(detector)) !== @__MODULE__
        @warn string(
            "detector $(detector) does not accept an `rng` keyword, so ",
            "pass_rng=true cannot seed it; if the detector is stochastic its null ",
            "distribution is not reproducible from `rng`. Pass pass_rng=false to ",
            "silence this warning.") maxlog = 1 _id = Symbol(
            :cp_significance_pass_rng_, nameof(typeof(detector)))
    end

    function detect_quality(network, child_rng)
        local result = try
            forward_rng ?
                detector(network; detector_kwargs..., rng=child_rng) :
                detector(network; detector_kwargs...)
        catch error
            if error isa ArgumentError || error isa DimensionMismatch ||
                    error isa MethodError
                local message = sprint(showerror, error)
                throw(ArgumentError(
                    "detector is incompatible with null_model=$(repr(null_model)): $message",
                ))
            end
            rethrow()
        end
        local quality = Float64(result.quality)
        isfinite(quality) || throw(ArgumentError("detector returned non-finite quality"))
        return quality
    end

    # Draw every seed before any threaded work. Each sample owns one RNG stream,
    # making results independent of scheduling and thread count.
    seeds = rand(rng, UInt, n_samples + 1)
    observed_quality = detect_quality(A, MersenneTwister(seeds[1]))
    null_qualities = Vector{Float64}(undef, n_samples)

    symmetric_input = _issymmetric(A)
    observed_edges = symmetric_input ? _edge_list(A) : _arc_list(A)
    base_edges = null_model === :er ? Tuple{Int,Int}[] : observed_edges
    base_weights = null_model === :weight_permutation ?
        Float64[A[source, target] for (source, target) in base_edges] : Float64[]
    swaps = if null_model in (:configuration, :directed_configuration)
        n_swaps === nothing ? max(10, 10 * length(base_edges)) : n_swaps
    else
        0
    end
    swaps >= 0 || throw(ArgumentError("n_swaps must be nonnegative"))
    attempt_budget = _configuration_attempt_budget(swaps, max_swap_attempts)
    possible = n * (n - 1) ÷ 2
    density = possible == 0 ? 0.0 : length(observed_edges) / possible
    accepted_swaps = zeros(Int, n_samples)
    swap_attempts = zeros(Int, n_samples)
    spec = _NullSpec(null_model, A, base_edges, base_weights, density, swaps,
        attempt_budget, swap_shortfall, symmetric_input)
    make_workspace = A isa SparseMatrixCSC ?
        (() -> _sparse_workspace(n, base_edges, base_weights)) :
        (() -> _dense_workspace(n, base_edges, base_weights))

    function evaluate_sample!(sample, workspace)
        local child_rng = MersenneTwister(seeds[sample + 1])
        local stats = _draw_null!(workspace, spec, child_rng)
        if stats !== nothing
            accepted_swaps[sample] = stats.accepted
            swap_attempts[sample] = stats.attempts
        end
        null_qualities[sample] = _with_component_checks_suppressed() do
            detect_quality(workspace.network, child_rng)
        end
        return nothing
    end

    schedule = if threaded && Threads.nthreads() > 1
        thread_schedule === :static ? :static : :greedy
    else
        :serial
    end
    _run_null_samples!(evaluate_sample!, make_workspace, n_samples, schedule)

    pvalue = (1 + count(>=(observed_quality), null_qualities)) / (n_samples + 1)
    null_mean = mean(null_qualities)
    null_std = n_samples >= 2 ? std(null_qualities; mean=null_mean) : NaN
    zscore = isfinite(null_std) && null_std > 0 ?
        (observed_quality - null_mean) / null_std : NaN
    diagnostics = NullDiagnostics(
        swaps,
        accepted_swaps,
        swap_attempts,
        attempt_budget,
        null_model in (:configuration, :directed_configuration) ?
            all(==(swaps), accepted_swaps) : true,
        null_model === :configuration ? :degree_sequence :
            null_model === :directed_configuration ? :in_out_degree_sequences :
            null_model === :weight_permutation ? :topology_and_weight_multiset :
            :expected_density,
    )
    return CPSignificanceResult(
        observed_quality,
        null_qualities,
        pvalue,
        pvalue <= alpha_value,
        alpha_value,
        null_model,
        diagnostics,
        null_mean,
        null_std,
        zscore,
    )
end
