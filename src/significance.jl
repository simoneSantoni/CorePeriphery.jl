"""
    CPSignificanceResult

Monte Carlo significance result for a detected core-periphery structure.
`pvalue` uses the conservative plus-one estimator
`(1 + count(null_quality >= observed_quality)) / (n_samples + 1)`.
`null_diagnostics` records the preserved constraint and switching completion.
"""
struct CPSignificanceResult
    observed_quality::Float64
    null_qualities::Vector{Float64}
    pvalue::Float64
    significant::Bool
    alpha::Float64
    null_model::Symbol
    null_diagnostics::NamedTuple
end

Base.:(==)(left::CPSignificanceResult, right::CPSignificanceResult) =
    left.observed_quality == right.observed_quality &&
    left.null_qualities == right.null_qualities &&
    left.pvalue == right.pvalue &&
    left.significant == right.significant &&
    left.alpha == right.alpha &&
    left.null_model == right.null_model &&
    left.null_diagnostics == right.null_diagnostics

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
        "$(stats.requested) requested swaps in $(stats.attempts) attempts"
    if policy === :error
        throw(ErrorException(message))
    elseif policy === :warn
        @warn message
    elseif policy !== :accept
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    end
    return nothing
end

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
        (!iszero(sampled[a, d]) || !iszero(sampled[c, b])) && continue

        sampled[a, b] = 0.0
        sampled[c, d] = 0.0
        sampled[a, d] = 1.0
        sampled[c, b] = 1.0
        arcs[first_index] = (a, d)
        arcs[second_index] = (c, b)
        accepted += 1
    end
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
    symmetric = issymmetric(A)
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
        (!iszero(sampled[u, v]) || !iszero(sampled[x, y])) && continue

        sampled[a, b] = sampled[b, a] = 0.0
        sampled[c, d] = sampled[d, c] = 0.0
        sampled[u, v] = sampled[v, u] = 1.0
        sampled[x, y] = sampled[y, x] = 1.0
        edges[first_index] = (u, v)
        edges[second_index] = (x, y)
        accepted += 1
    end
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

"""
    cp_significance(A, detector; null_model=:configuration, n_samples=100,
                    alpha=0.05, rng=Random.default_rng(),
                    detector_kwargs=NamedTuple(), pass_rng=false, n_swaps=nothing,
                    max_swap_attempts=nothing, swap_shortfall=:error,
                    threaded=false, thread_schedule=:auto)

Test whether a detector's observed quality exceeds qualities obtained from
random graphs. `detector` must return a result with a finite `quality` field.
The `:configuration` null uses undirected degree-preserving double-edge swaps;
`:directed_configuration` preserves binary in- and out-degrees through directed
arc switches; `:weight_permutation` preserves topology and the edge-weight
multiset; and `:er` samples an undirected Erdős-Rényi graph with observed density.

Configuration samples must complete all requested swaps by default. Use
`max_swap_attempts` to set the proposal budget and `swap_shortfall=:warn` or
`:accept` to retain incomplete samples explicitly. Per-sample accepted swaps and
attempt counts are returned in `null_diagnostics`.

Set `pass_rng=true` for stochastic detectors accepting an `rng` keyword. A
fresh child RNG is then passed to every detection call.

With `threaded=true`, null samples are distributed using `Threads.@threads`.
Sample seeds are generated serially before execution, so serial and threaded
runs are bit-for-bit identical for deterministic detectors and stochastic
detectors using `pass_rng=true`. Detectors must otherwise be thread-safe.

`thread_schedule=:auto` uses a benchmark-informed Julia 1.12 policy: greedy
task-local scheduling for ER samples and for sufficiently large configuration
runs, otherwise static logical workers. Older supported Julia versions always
use the static path. Use `:static` or `:greedy` to override the policy.
"""
function cp_significance(
    A::AbstractMatrix{<:Real},
    detector::F;
    null_model::Symbol=:configuration,
    n_samples::Int=100,
    alpha::Float64=0.05,
    rng::AbstractRNG=Random.default_rng(),
    detector_kwargs::NamedTuple=NamedTuple(),
    pass_rng::Bool=false,
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
    n = if null_model in (:configuration, :er)
        _validate_adjacency(A; symmetric=true, binary=true)
    elseif null_model === :directed_configuration
        _validate_adjacency(A; binary=true)
    else
        _validate_adjacency(A)
    end
    n_samples >= 1 || throw(ArgumentError("n_samples must be positive"))
    0.0 < alpha < 1.0 || throw(ArgumentError("alpha must be in (0, 1)"))
    swap_shortfall in (:error, :warn, :accept) ||
        throw(ArgumentError("swap_shortfall must be :error, :warn, or :accept"))
    thread_schedule in (:auto, :static, :greedy) ||
        throw(ArgumentError("thread_schedule must be :auto, :static, or :greedy"))
    if thread_schedule === :greedy && VERSION < v"1.12"
        throw(ArgumentError("thread_schedule=:greedy requires Julia 1.12 or newer"))
    end

    function detect_quality(network, child_rng)
        result = try
            pass_rng ?
                detector(network; detector_kwargs..., rng=child_rng) :
                detector(network; detector_kwargs...)
        catch error
            if error isa ArgumentError || error isa DimensionMismatch ||
                    error isa MethodError
                message = sprint(showerror, error)
                throw(ArgumentError(
                    "detector is incompatible with null_model=$(repr(null_model)): $message",
                ))
            end
            rethrow()
        end
        quality = Float64(result.quality)
        isfinite(quality) || throw(ArgumentError("detector returned non-finite quality"))
        return quality
    end

    # Draw every seed before any threaded work. Each sample owns one RNG stream,
    # making results independent of scheduling and thread count.
    seeds = rand(rng, UInt, n_samples + 1)
    observed_rng = MersenneTwister(seeds[1])
    observed_quality = detect_quality(A, observed_rng)
    null_qualities = Vector{Float64}(undef, n_samples)

    symmetric_input = issymmetric(A)
    observed_edges = symmetric_input ? _edge_list(A) : _arc_list(A)
    base_edges = null_model in (:configuration, :directed_configuration,
                                :weight_permutation) ?
        observed_edges : Tuple{Int,Int}[]
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

    function evaluate_sample!(sample, matrix_buffer, edge_buffer, weight_buffer)
        child_rng = MersenneTwister(seeds[sample + 1])
        randomized = if null_model === :er
            _sample_er!(matrix_buffer, density, child_rng)
        elseif null_model === :configuration
            stats = Ref{ConfigurationSwapStats}()
            _sample_configuration!(
                matrix_buffer, A, edge_buffer, base_edges, child_rng, swaps;
                max_attempts=attempt_budget,
                on_shortfall=swap_shortfall,
                stats_ref=stats,
            )
            accepted_swaps[sample] = stats[].accepted
            swap_attempts[sample] = stats[].attempts
            matrix_buffer
        elseif null_model === :directed_configuration
            stats = Ref{ConfigurationSwapStats}()
            _sample_directed_configuration!(
                matrix_buffer, A, edge_buffer, base_edges, child_rng, swaps;
                max_attempts=attempt_budget,
                on_shortfall=swap_shortfall,
                stats_ref=stats,
            )
            accepted_swaps[sample] = stats[].accepted
            swap_attempts[sample] = stats[].attempts
            matrix_buffer
        else
            _sample_weight_permutation!(
                matrix_buffer, base_edges, base_weights, weight_buffer, child_rng;
                symmetric=symmetric_input,
            )
        end
        null_qualities[sample] = detect_quality(randomized, child_rng)
        return nothing
    end

    if threaded && Threads.nthreads() > 1
        auto_greedy = VERSION >= v"1.12" &&
            (null_model === :er || n_samples >= 4 * Threads.nthreads())
        schedule = thread_schedule === :auto ?
            (auto_greedy ? :greedy : :static) : thread_schedule
        if schedule === :greedy
            @static if VERSION >= v"1.12"
                # Greedy scheduling balances irregular detector/null-model
                # runtimes. Every scheduler task lazily owns one reusable buffer,
                # so task migration cannot create buffer races.
                task_buffers = Base.OncePerTask{
                    Tuple{Matrix{Float64},Vector{Tuple{Int,Int}},Vector{Float64}}
                }() do
                    (Matrix{Float64}(undef, n, n), copy(base_edges), copy(base_weights))
                end
                Threads.@threads :greedy for sample in 1:n_samples
                    matrix_buffer, edge_buffer, weight_buffer = task_buffers()
                    evaluate_sample!(sample, matrix_buffer, edge_buffer, weight_buffer)
                end
            else
                error("unreachable greedy scheduling branch")
            end
        else
            # Index buffers by logical worker, not `threadid()`: thread IDs need
            # not be contiguous, and static workers exclusively own a buffer.
            buffer_count = min(Threads.nthreads(), n_samples)
            matrix_buffers =
                [Matrix{Float64}(undef, n, n) for _ in 1:buffer_count]
            edge_buffers = [copy(base_edges) for _ in 1:buffer_count]
            weight_buffers = [copy(base_weights) for _ in 1:buffer_count]
            Threads.@threads :static for worker in 1:buffer_count
                for sample in worker:buffer_count:n_samples
                    evaluate_sample!(
                        sample, matrix_buffers[worker], edge_buffers[worker],
                        weight_buffers[worker])
                end
            end
        end
    else
        matrix_buffer = Matrix{Float64}(undef, n, n)
        edge_buffer = copy(base_edges)
        weight_buffer = copy(base_weights)
        for sample in 1:n_samples
            evaluate_sample!(sample, matrix_buffer, edge_buffer, weight_buffer)
        end
    end

    pvalue = (1 + count(>=(observed_quality), null_qualities)) / (n_samples + 1)
    diagnostics = (
        requested_swaps=swaps,
        accepted_swaps=accepted_swaps,
        swap_attempts=swap_attempts,
        max_swap_attempts=attempt_budget,
        complete=null_model in (:configuration, :directed_configuration) ?
            all(==(swaps), accepted_swaps) : true,
        preserved=null_model === :configuration ? :degree_sequence :
            null_model === :directed_configuration ? :in_out_degree_sequences :
            null_model === :weight_permutation ? :topology_and_weight_multiset :
            :expected_density,
    )
    return CPSignificanceResult(
        observed_quality,
        null_qualities,
        pvalue,
        pvalue <= alpha,
        alpha,
        null_model,
        diagnostics,
    )
end
