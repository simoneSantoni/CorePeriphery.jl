function _positive_reachable(A::AbstractMatrix{<:Real}, start::Int)
    n = size(A, 1)
    seen = falses(n)
    seen[start] = true
    stack = [start]
    while !isempty(stack)
        node = pop!(stack)
        @inbounds for neighbor in axes(A, 1)
            if !seen[neighbor] && A[neighbor, node] > 0
                seen[neighbor] = true
                push!(stack, neighbor)
            end
        end
    end
    return seen
end

function _positive_reachable(A::SparseMatrixCSC{<:Real}, start::Int)
    n = size(A, 1)
    seen = falses(n)
    seen[start] = true
    stack = [start]
    rows = rowvals(A)
    values = nonzeros(A)
    while !isempty(stack)
        node = pop!(stack)
        @inbounds for index in nzrange(A, node)
            neighbor = rows[index]
            if !seen[neighbor] && values[index] > 0
                seen[neighbor] = true
                push!(stack, neighbor)
            end
        end
    end
    return seen
end

function _minimum_strength_node(strengths::AbstractVector{<:Real},
    tie_break::Symbol, rng::AbstractRNG)
    if tie_break === :deterministic
        return argmin(strengths)
    end

    minimum_strength = minimum(strengths)
    selected = 0
    tied = 0
    @inbounds for node in eachindex(strengths)
        if strengths[node] == minimum_strength
            tied += 1
            rand(rng, 1:tied) == 1 && (selected = node)
        end
    end
    return selected
end

function _transition_transpose_mul!(destination::Vector{Float64},
    A::AbstractMatrix{<:Real}, scaled_probability::Vector{Float64})
    @inbounds for column in axes(A, 2)
        total = 0.0
        for row in axes(A, 1)
            total += Float64(A[row, column]) * scaled_probability[row]
        end
        destination[column] = total
    end
    return destination
end

function _transition_transpose_mul!(destination::Vector{Float64},
    A::SparseMatrixCSC{<:Real}, scaled_probability::Vector{Float64})
    fill!(destination, 0.0)
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for column in axes(A, 2), index in nzrange(A, column)
        destination[column] += Float64(values[index]) * scaled_probability[rows[index]]
    end
    return destination
end

function _stationary_residual!(walked::Vector{Float64},
    scaled_probability::Vector{Float64}, A::AbstractMatrix{<:Real},
    probability::Vector{Float64}, out_strengths::Vector{Float64})
    @inbounds for node in eachindex(probability)
        scaled_probability[node] = probability[node] / out_strengths[node]
    end
    _transition_transpose_mul!(walked, A, scaled_probability)
    return sum(abs(walked[node] - probability[node]) for node in eachindex(probability))
end

function _stationary_lazy(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, tolerance::Float64, max_iter::Int)
    n = size(A, 1)
    probability = fill(inv(Float64(n)), n)
    candidate = similar(probability)
    walked = similar(probability)
    scaled_probability = similar(probability)

    for iteration in 1:max_iter
        @inbounds for node in 1:n
            scaled_probability[node] = probability[node] / out_strengths[node]
        end
        _transition_transpose_mul!(walked, A, scaled_probability)

        total = 0.0
        @inbounds for node in 1:n
            candidate[node] = 0.5 * (probability[node] + walked[node])
            total += candidate[node]
        end
        inverse_total = inv(total)
        distance = 0.0
        @inbounds for node in 1:n
            candidate[node] *= inverse_total
            distance += abs(candidate[node] - probability[node])
        end
        if distance <= 0.5 * tolerance
            residual = _stationary_residual!(
                walked, scaled_probability, A, candidate, out_strengths,
            )
            if residual <= tolerance
                return StationaryDistributionResult(
                    copy(candidate), residual, iteration, true, :lazy,
                )
            end
        end
        probability, candidate = candidate, probability
    end

    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    throw(ConvergenceError(
        "stationary method :lazy did not converge in $max_iter iterations " *
        "at tolerance $tolerance (final L1 residual $residual); increase " *
        "stationary_max_iter, relax stationary_tol, or use stationary_method=:linear",
    ))
end

function _stationary_linear(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, tolerance::Float64)
    n = size(A, 1)
    transition = Matrix{Float64}(A)
    @inbounds for column in 1:n, row in 1:n
        transition[row, column] /= out_strengths[row]
    end
    system = Matrix{Float64}(I, n, n) - transpose(transition)
    system[end, :] .= 1.0
    right_hand_side = zeros(Float64, n)
    right_hand_side[end] = 1.0
    probability = system \ right_hand_side
    minimum(probability) >= -sqrt(eps(Float64)) || throw(ConvergenceError(
        "stationary method :linear produced a negative probability"))
    probability .= max.(probability, 0.0)
    probability ./= sum(probability)
    walked = similar(probability)
    scaled_probability = similar(probability)
    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    residual <= tolerance || throw(ConvergenceError(
        "stationary method :linear failed tolerance $tolerance " *
        "(final L1 residual $residual)"))
    return StationaryDistributionResult(probability, residual, 1, true, :linear)
end

function _validated_stationary(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}, probability_input::AbstractVector{<:Real},
    tolerance::Float64)
    n = size(A, 1)
    length(probability_input) == n || throw(DimensionMismatch(
        "stationary_distribution must have length $n"))
    probability = Float64.(probability_input)
    all(isfinite, probability) || throw(ArgumentError(
        "stationary_distribution entries must be finite"))
    all(>=(0.0), probability) || throw(ArgumentError(
        "stationary_distribution entries must be nonnegative"))
    total = sum(probability)
    total > 0.0 || throw(ArgumentError(
        "stationary_distribution must have positive mass"))
    probability ./= total
    walked = similar(probability)
    scaled_probability = similar(probability)
    residual = _stationary_residual!(
        walked, scaled_probability, A, probability, out_strengths,
    )
    residual <= tolerance || throw(ArgumentError(
        "stationary_distribution has L1 residual $residual, exceeding $tolerance"))
    return StationaryDistributionResult(probability, residual, 0, true, :supplied)
end

function _compute_stationary(A::AbstractMatrix{<:Real},
    out_strengths::Vector{Float64}; method::Symbol, tolerance::Float64,
    max_iter::Int, supplied::Union{Nothing,AbstractVector})
    supplied === nothing || return _validated_stationary(
        A, out_strengths, supplied, tolerance)
    selected = method === :auto ?
        ((A isa SparseMatrixCSC || size(A, 1) > 2_000) ? :lazy : :linear) : method
    selected === :lazy && return _stationary_lazy(
        A, out_strengths, tolerance, max_iter)
    selected === :linear && return _stationary_linear(A, out_strengths, tolerance)
    throw(ArgumentError("stationary_method must be :auto, :lazy, or :linear"))
end

"""
    rossa_stationary_distribution(A; method=:auto, tol=1e-12,
                                  max_iter=100_000)

Compute the stationary distribution used by directed Della Rossa profiling. The
returned `StationaryDistributionResult` includes the L1 balance residual, iteration
count, and selected method. `:linear` is a robust dense constrained solve; `:lazy`
uses stored-edge matrix-vector products and is suitable for large sparse graphs;
`:auto` selects between them. Input must be strongly connected, directed or
undirected, finite, nonnegative, and loop-free. A solver that does not reach
`tol` raises [`ConvergenceError`](@ref).
"""
function rossa_stationary_distribution(A::AbstractMatrix{<:Real};
    method::Symbol=:auto, tol::Real=1e-12, max_iter::Int=100_000)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A)
    n >= 2 || throw(ArgumentError("stationary distribution requires at least two nodes"))
    method in (:auto, :lazy, :linear) ||
        throw(ArgumentError("method must be :auto, :lazy, or :linear"))
    isfinite(tol) && tol >= 0 ||
        throw(ArgumentError("tol must be finite and nonnegative"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be positive"))
    reverse_A = A isa SparseMatrixCSC ? copy(transpose(A)) : transpose(A)
    all(_positive_reachable(A, 1)) && all(_positive_reachable(reverse_A, 1)) ||
        throw(ArgumentError("stationary distribution requires a strongly connected graph"))
    out_strengths = Float64.(vec(sum(A, dims=2)))
    return _compute_stationary(
        A, out_strengths;
        method=method,
        tolerance=Float64(tol),
        max_iter=max_iter,
        supplied=nothing,
    )
end


function _axpy_scaled_column!(destination::Vector{Float64},
    A::AbstractMatrix{<:Real}, column::Int, row_scale::Vector{Float64})
    @inbounds for row in axes(A, 1)
        destination[row] += Float64(A[row, column]) * row_scale[row]
    end
    return destination
end

function _axpy_scaled_column!(destination::Vector{Float64},
    A::SparseMatrixCSC{<:Real}, column::Int, row_scale::Vector{Float64})
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for index in nzrange(A, column)
        row = rows[index]
        destination[row] += Float64(values[index]) * row_scale[row]
    end
    return destination
end

function _select_profile_candidate(in_periphery::BitVector,
    boundary::Vector{Float64}, mass::Float64, internal::Float64,
    node_mass::Vector{Float64}, strengths::Vector{Float64},
    tie_break::Symbol, rng::AbstractRNG)
    best_node = 0
    best_alpha = Inf
    best_strength = Inf
    tied = 0

    @inbounds for node in eachindex(in_periphery)
        in_periphery[node] && continue
        candidate_alpha = (internal + boundary[node]) / (mass + node_mass[node])
        if candidate_alpha < best_alpha ||
           (candidate_alpha == best_alpha && strengths[node] < best_strength)
            best_node = node
            best_alpha = candidate_alpha
            best_strength = strengths[node]
            tied = 1
        elseif candidate_alpha == best_alpha && strengths[node] == best_strength
            if tie_break === :paper
                tied += 1
                rand(rng, 1:tied) == 1 && (best_node = node)
            elseif node < best_node
                best_node = node
            end
        end
    end
    return best_node, best_alpha, boundary[best_node]
end

function _profile_result(c::Vector{Float64}, profile::Vector{Float64};
    partition::Symbol, threshold::Union{Nothing,Real},
    core_size::Union{Nothing,Integer})
    n = length(c)
    quality = n == 2 ? NaN : 1 - (2 / (n - 2)) * sum(@view profile[1:(n - 1)])
    quality = isfinite(quality) ? clamp(quality, 0.0, 1.0) : quality
    core_nodes, periphery_nodes = _partition(
        c; rule=partition, threshold=threshold, core_size=core_size)
    return CPResult(c, core_nodes, periphery_nodes, quality, "Random Walker Profiling")
end

"""
    random_walker_profiling(A; tie_break=:deterministic,
                            rng=Random.default_rng(),
                            stationary_method=:lazy,
                            stationary_distribution=nothing,
                            stationary_tol=1e-12,
                            stationary_max_iter=100_000,
                            partition=:median, threshold=nothing,
                            core_size=nothing)

Compute the Della Rossa persistence profile of an undirected connected or directed
strongly connected network. The default deterministic mode replaces the paper's random
choices with lower node index; `tie_break=:paper` restores the published stochastic
selection rule using `rng`. The result quality is the paper's cp-centralization.
Directed inputs accept `stationary_method=:lazy`, `:linear`, or `:auto`, or a
validated `stationary_distribution` supplied by the caller; a stationary solver
that does not converge raises [`ConvergenceError`](@ref).

# Keyword arguments
- `tie_break`: `:deterministic` (default) or the stochastic, paper-exact `:paper`
- `rng`: Random-number generator used when `tie_break=:paper`
- `stationary_method`: Directed solver, `:lazy` (default), `:linear`, or `:auto`
- `stationary_distribution`: Optional validated stationary probability vector
- `stationary_tol`: L1 convergence tolerance for the directed stationary distribution
- `stationary_max_iter`: Maximum stationary iterations for directed inputs
- `partition`: Core/periphery split rule applied to the coreness vector:
  `:median` (default; nodes strictly above the median are core), `:none`
  (all nodes periphery), `:threshold` (nodes strictly above `threshold`), or
  `:core_size` (the `core_size` highest-coreness nodes, ties broken by lower index)
- `threshold`: Cutoff used by `partition=:threshold`
- `core_size`: Number of core nodes used by `partition=:core_size`
- `n_walks`, `walk_length`: Deprecated legacy simulation keywords. They are
  accepted and ignored (profiling is analytic); passing them emits a deprecation
  warning (shown with `--depwarn=yes`).
"""
function random_walker_profiling(A::AbstractMatrix{<:Real};
    tie_break::Symbol=:deterministic,
    rng::AbstractRNG=Random.default_rng(),
    stationary_method::Symbol=:lazy,
    stationary_distribution::Union{Nothing,AbstractVector}=nothing,
    stationary_tol::Real=1e-12,
    stationary_max_iter::Int=100_000,
    partition::Symbol=:median,
    threshold::Union{Nothing,Real}=nothing,
    core_size::Union{Nothing,Integer}=nothing,
    n_walks::Union{Nothing,Int}=nothing,
    walk_length::Union{Nothing,Int}=nothing)
    A = _unwrap_adjacency(A)
    n = _validate_adjacency(A)
    n >= 2 || throw(ArgumentError("random-walker profiling requires at least two nodes"))
    tie_break in (:deterministic, :paper) ||
        throw(ArgumentError("tie_break must be :deterministic or :paper"))
    stationary_method in (:auto, :lazy, :linear) ||
        throw(ArgumentError("stationary_method must be :auto, :lazy, or :linear"))
    isfinite(stationary_tol) && stationary_tol >= 0 ||
        throw(ArgumentError("stationary_tol must be finite and nonnegative"))
    stationary_max_iter >= 1 ||
        throw(ArgumentError("stationary_max_iter must be positive"))
    if n_walks !== nothing
        n_walks >= 0 || throw(ArgumentError("n_walks must be nonnegative"))
        Base.depwarn(
            "the `n_walks` keyword of `random_walker_profiling` is deprecated and " *
            "ignored: the persistence profile is computed analytically",
            :random_walker_profiling,
        )
    end
    if walk_length !== nothing
        walk_length >= 0 || throw(ArgumentError("walk_length must be nonnegative"))
        Base.depwarn(
            "the `walk_length` keyword of `random_walker_profiling` is deprecated and " *
            "ignored: the persistence profile is computed analytically",
            :random_walker_profiling,
        )
    end

    symmetric = _issymmetric(A)
    symmetric && stationary_distribution !== nothing && throw(ArgumentError(
        "stationary_distribution is only used for directed inputs"))
    reverse_A = symmetric ? A :
        (A isa SparseMatrixCSC ? copy(transpose(A)) : transpose(A))

    # `_positive_reachable` follows positive entries down a column. Checking A and
    # its transpose therefore checks both directed reachability orientations.
    all(_positive_reachable(A, 1)) &&
        (symmetric || all(_positive_reachable(reverse_A, 1))) ||
        throw(ArgumentError(symmetric ?
            "random-walker profiling requires a connected graph" :
            "directed random-walker profiling requires a strongly connected graph"))

    in_periphery = falses(n)
    profile = zeros(Float64, n)
    c = zeros(Float64, n)
    boundary = zeros(Float64, n)

    if symmetric
        strengths = Float64.(vec(sum(A, dims=2)))
        first_node = _minimum_strength_node(strengths, tie_break, rng)
        in_periphery[first_node] = true
        mass = strengths[first_node]
        internal = 0.0
        _axpy_column!(boundary, A, first_node, 2.0)

        for k in 2:n
            best_node, best_alpha, best_boundary = _select_profile_candidate(
                in_periphery, boundary, mass, internal, strengths, strengths,
                tie_break, rng,
            )
            in_periphery[best_node] = true
            internal += best_boundary
            mass += strengths[best_node]
            profile[k] = clamp(best_alpha, 0.0, 1.0)
            c[best_node] = profile[k]
            _axpy_column!(boundary, A, best_node, 2.0)
        end
    else
        out_strengths = Float64.(vec(sum(A, dims=2)))
        in_strengths = Float64.(vec(sum(A, dims=1)))
        strengths = out_strengths .+ in_strengths
        stationary_result = _compute_stationary(
            A, out_strengths;
            method=stationary_method,
            tolerance=Float64(stationary_tol),
            max_iter=stationary_max_iter,
            supplied=stationary_distribution,
        )
        probability = stationary_result.distribution
        flow_scale = probability ./ out_strengths

        first_node = _minimum_strength_node(strengths, tie_break, rng)
        in_periphery[first_node] = true
        mass = probability[first_node]
        internal = 0.0
        _axpy_scaled_column!(boundary, A, first_node, flow_scale)
        _axpy_column!(boundary, reverse_A, first_node, flow_scale[first_node])

        for k in 2:n
            best_node, best_alpha, best_boundary = _select_profile_candidate(
                in_periphery, boundary, mass, internal, probability, strengths,
                tie_break, rng,
            )
            in_periphery[best_node] = true
            internal += best_boundary
            mass += probability[best_node]
            profile[k] = clamp(best_alpha, 0.0, 1.0)
            c[best_node] = profile[k]
            _axpy_scaled_column!(boundary, A, best_node, flow_scale)
            _axpy_column!(boundary, reverse_A, best_node, flow_scale[best_node])
        end
    end

    return _profile_result(
        c, profile; partition=partition, threshold=threshold, core_size=core_size)
end

function _average_ranks(values::AbstractVector{<:Real})
    n = length(values)
    order = sortperm(values; alg=Base.Sort.MergeSort)
    ranks = zeros(Float64, n)
    first = 1
    while first <= n
        last = first
        while last < n && values[order[last + 1]] == values[order[first]]
            last += 1
        end
        rank = (first + last) / 2
        @inbounds for position in first:last
            ranks[order[position]] = rank
        end
        first = last + 1
    end
    return ranks
end

function _rank_correlation(left::AbstractVector{<:Real}, right::AbstractVector{<:Real})
    left_ranks = _average_ranks(left)
    right_ranks = _average_ranks(right)
    left_centered = left_ranks .- mean(left_ranks)
    right_centered = right_ranks .- mean(right_ranks)
    denominator = norm(left_centered) * norm(right_centered)
    return iszero(denominator) ? 1.0 : dot(left_centered, right_centered) / denominator
end

"""
    rossa_profile_ensemble(A; n_runs=100, rng=Random.default_rng(),
                           interval=(0.025, 0.975), threaded=false,
                           store_replicates=true,
                           memory_limit_bytes=512 * 1024^2, kwargs...)

Fit repeated `tie_break=:paper` Rossa profiles and summarize node coreness,
cp-centralization, rank stability relative to the ensemble mean, and the number of
unique profiles. Child seeds are drawn before fitting, so fixed seeds produce the
same result in serial and threaded execution. Remaining keywords are forwarded to
`random_walker_profiling`; passing `tie_break` (the ensemble always uses `:paper`)
raises `ArgumentError`, and `rng` seeds the ensemble rather than a single profile.

The working profile matrix requires `8 * n * n_runs` bytes even when raw replicates
are omitted from the returned compact result. `memory_limit_bytes` checks that
allocation before fitting so large ensembles fail predictably instead of exhausting
process memory.
"""
function rossa_profile_ensemble(A::AbstractMatrix{<:Real};
    n_runs::Int=100,
    rng::AbstractRNG=Random.default_rng(),
    interval::Tuple{<:Real,<:Real}=(0.025, 0.975),
    threaded::Bool=false,
    store_replicates::Bool=true,
    memory_limit_bytes::Integer=512 * 1024^2,
    kwargs...)
    haskey(kwargs, :tie_break) && throw(ArgumentError(
        "rossa_profile_ensemble always uses tie_break=:paper; " *
        "the tie_break keyword cannot be overridden"))
    haskey(kwargs, :rng) && throw(ArgumentError(
        "rng cannot be forwarded to individual profiles; " *
        "pass it as the ensemble's own rng keyword"))
    A = _unwrap_adjacency(A)
    n = size(A, 1)
    n_runs >= 1 || throw(ArgumentError("n_runs must be positive"))
    memory_limit_bytes >= 0 ||
        throw(ArgumentError("memory_limit_bytes must be nonnegative"))
    required_profile_bytes = Base.checked_mul(Base.checked_mul(8, n), n_runs)
    required_profile_bytes <= memory_limit_bytes || throw(ArgumentError(
        "Rossa ensemble profile workspace requires $required_profile_bytes bytes, " *
        "exceeding memory_limit_bytes=$memory_limit_bytes; reduce n_runs or " *
        "raise the explicit limit",
    ))
    lower_probability, upper_probability = Float64.(interval)
    0.0 <= lower_probability <= upper_probability <= 1.0 ||
        throw(ArgumentError("interval probabilities must satisfy 0 <= lower <= upper <= 1"))
    seeds = rand(rng, UInt, n_runs)
    coreness = Matrix{Float64}(undef, n, n_runs)
    centralizations = Vector{Float64}(undef, n_runs)

    function fit_run!(run)
        result = random_walker_profiling(
            A;
            tie_break=:paper,
            rng=MersenneTwister(seeds[run]),
            kwargs...,
        )
        coreness[:, run] = result.coreness
        centralizations[run] = result.quality
        return nothing
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads :static for run in 1:n_runs
            fit_run!(run)
        end
    else
        for run in 1:n_runs
            fit_run!(run)
        end
    end

    mean_coreness = vec(mean(coreness, dims=2))
    std_coreness = vec(std(coreness, dims=2; corrected=false))
    lower_coreness = [quantile(@view(coreness[node, :]), lower_probability) for node in 1:n]
    upper_coreness = [quantile(@view(coreness[node, :]), upper_probability) for node in 1:n]
    rank_stability = mean(
        _rank_correlation(@view(coreness[:, run]), mean_coreness) for run in 1:n_runs)
    unique_profiles = length(Set(Tuple(@view(coreness[:, run])) for run in 1:n_runs))
    raw = store_replicates ? coreness : nothing
    return RossaEnsembleResult(
        mean_coreness,
        std_coreness,
        lower_coreness,
        upper_coreness,
        centralizations,
        mean(centralizations),
        std(centralizations; corrected=false),
        rank_stability,
        unique_profiles,
        raw,
        seeds,
    )
end
