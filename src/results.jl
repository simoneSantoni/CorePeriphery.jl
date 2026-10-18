# Result containers shared by every detector, plus their equality, hashing,
# and compact display methods.

"""
    CPResult

Result structure for core-periphery detection algorithms.

# Fields
- `coreness::Vector{Float64}`: Coreness score for each node (higher = more
  core-like). Isolates (nodes with zero strength) always score `0`.
- `core_nodes::Vector{Int}`: Indices of nodes classified as core. For the
  continuous detectors this is the split produced by their `partition` keyword
  (a median rule by default), not an estimate of the model; isolates are never
  core.
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
    CPDirectedResult

Result of directed core-periphery detection. `out_coreness` describes how
strongly each node connects outward to core-like receivers, while
`in_coreness` describes how strongly it receives connections from core-like
senders. `coreness` is their arithmetic mean.
"""
struct CPDirectedResult
    out_coreness::Vector{Float64}
    in_coreness::Vector{Float64}
    coreness::Vector{Float64}
    core_nodes::Vector{Int}
    periphery_nodes::Vector{Int}
    residual::Float64
    quality::Float64
    algorithm::String
end

"""
    CPMultiResult

Result structure for multiple core-periphery pairs detection.

# Fields
- `pair_labels::Vector{Int}`: Pair assignment for each node (`1, 2, ..., K`).
  Isolates (nodes with zero strength) are not fitted and carry the label `0`;
  consumers that group nodes by pair should treat `0` as "no pair".
- `coreness::Vector{Float64}`: Binary coreness within each pair (`0` for isolates)
- `n_pairs::Int`: Number of detected core-periphery pairs `K` among the
  non-isolated nodes
- `quality::Float64`: Quality score (Q^cp)
- `algorithm::String`: Name of the algorithm used
- `pair_selection::Symbol`: `:objective` or `:penalized` pair-count rule
- `candidate_pair_counts::Vector{Int}`: Evaluated model sizes
- `candidate_qualities::Vector{Float64}`: Raw best objective at each model size
- `selection_scores::Vector{Float64}`: Scores used to select the returned model
"""
struct CPMultiResult
    pair_labels::Vector{Int}
    coreness::Vector{Float64}
    n_pairs::Int
    quality::Float64
    algorithm::String
    pair_selection::Symbol
    candidate_pair_counts::Vector{Int}
    candidate_qualities::Vector{Float64}
    selection_scores::Vector{Float64}
end

CPMultiResult(pair_labels::Vector{Int}, coreness::Vector{Float64}, n_pairs::Int,
    quality::Float64, algorithm::String) = CPMultiResult(
    pair_labels, coreness, n_pairs, quality, algorithm,
    :objective, [n_pairs], [quality], [quality],
)

"""
    NullDiagnostics

Per-run diagnostics of the null-model sampler used by [`cp_significance`](@ref).

# Fields
- `requested_swaps::Int`: Edge switches requested per configuration sample (0 for
  other null models)
- `accepted_swaps::Vector{Int}`: Accepted switches in each sample
- `swap_attempts::Vector{Int}`: Proposals made in each sample
- `max_swap_attempts::Int`: Proposal budget per sample
- `complete::Bool`: Whether every sample completed the requested switches
- `preserved::Symbol`: Constraint preserved by the null model
  (`:degree_sequence`, `:in_out_degree_sequences`,
  `:topology_and_weight_multiset`, or `:expected_density`)
"""
struct NullDiagnostics
    requested_swaps::Int
    accepted_swaps::Vector{Int}
    swap_attempts::Vector{Int}
    max_swap_attempts::Int
    complete::Bool
    preserved::Symbol
end

"""
    CPSignificanceResult

Monte Carlo significance result for a detected core-periphery structure.
`pvalue` uses the conservative plus-one estimator
`(1 + count(null_quality >= observed_quality)) / (n_samples + 1)`.
`null_diagnostics` is a [`NullDiagnostics`](@ref) recording the preserved
constraint and switching completion.

# Fields
- `observed_quality::Float64`: detector quality on the observed network
- `null_qualities::Vector{Float64}`: detector quality on every null sample
- `pvalue::Float64`: plus-one Monte Carlo p-value
- `significant::Bool`: `pvalue <= alpha`
- `alpha::Float64`: significance level used for `significant`
- `null_model::Symbol`: null model that generated the samples
- `null_diagnostics::NullDiagnostics`: sampler diagnostics
- `null_mean::Float64`: mean of `null_qualities`
- `null_std::Float64`: sample standard deviation of `null_qualities` (`NaN`
  with a single sample)
- `zscore::Float64`: `(observed_quality - null_mean) / null_std`, an effect
  size in null standard deviations; `NaN` when `null_std` is zero or undefined
"""
struct CPSignificanceResult
    observed_quality::Float64
    null_qualities::Vector{Float64}
    pvalue::Float64
    significant::Bool
    alpha::Float64
    null_model::Symbol
    null_diagnostics::NullDiagnostics
    null_mean::Float64
    null_std::Float64
    zscore::Float64
end

"""
    StationaryDistributionResult

Stationary random-walk distribution and numerical convergence diagnostics used by
directed Della Rossa profiling.
"""
struct StationaryDistributionResult
    distribution::Vector{Float64}
    residual::Float64
    iterations::Int
    converged::Bool
    method::Symbol
end

"""
    RossaEnsembleResult

Uncertainty summary from repeated paper-mode Della Rossa profiles. Optional raw
replicate coreness occupies `O(n * n_runs)` memory; set `store_replicates=false`
to retain only summaries after fitting.
"""
struct RossaEnsembleResult
    mean_coreness::Vector{Float64}
    std_coreness::Vector{Float64}
    lower_coreness::Vector{Float64}
    upper_coreness::Vector{Float64}
    centralizations::Vector{Float64}
    mean_centralization::Float64
    std_centralization::Float64
    rank_stability::Float64
    unique_profiles::Int
    replicate_coreness::Union{Nothing,Matrix{Float64}}
    seeds::Vector{UInt}
end

# Structural equality, `isequal`, and `hash` over every field, so results can be
# compared, deduplicated, and used as `Dict`/`Set` keys. `==` follows the field
# semantics (a `NaN` quality is not `==` to itself); `isequal` and `hash` agree.
for T in (
    :CPResult,
    :CPDirectedResult,
    :CPMultiResult,
    :NullDiagnostics,
    :CPSignificanceResult,
    :StationaryDistributionResult,
    :RossaEnsembleResult,
)
    @eval begin
        function Base.:(==)(left::$T, right::$T)
            return all(
                getfield(left, name) == getfield(right, name)
                for name in fieldnames($T)
            )
        end
        function Base.isequal(left::$T, right::$T)
            return all(
                isequal(getfield(left, name), getfield(right, name))
                for name in fieldnames($T)
            )
        end
        function Base.hash(result::$T, h::UInt)
            h = hash($(QuoteNode(T)), h)
            for name in fieldnames($T)
                h = hash(getfield(result, name), h)
            end
            return h
        end
    end
end

_show_number(x::Real) = isfinite(x) ? string(round(x; sigdigits=4)) : string(x)

function Base.show(io::IO, ::MIME"text/plain", result::CPResult)
    print(io, "CPResult: ", result.algorithm,
        "\n  nodes = ", length(result.coreness),
        ", core = ", length(result.core_nodes),
        ", periphery = ", length(result.periphery_nodes),
        "\n  quality = ", _show_number(result.quality))
end

function Base.show(io::IO, ::MIME"text/plain", result::CPDirectedResult)
    print(io, "CPDirectedResult: ", result.algorithm,
        "\n  nodes = ", length(result.coreness),
        ", core = ", length(result.core_nodes),
        ", periphery = ", length(result.periphery_nodes),
        "\n  quality = ", _show_number(result.quality),
        ", residual = ", _show_number(result.residual))
end

function Base.show(io::IO, ::MIME"text/plain", result::CPMultiResult)
    print(io, "CPMultiResult: ", result.algorithm,
        "\n  nodes = ", length(result.pair_labels),
        ", pairs = ", result.n_pairs,
        ", core = ", count(!iszero, result.coreness),
        "\n  quality = ", _show_number(result.quality),
        ", pair_selection = ", repr(result.pair_selection))
end

function Base.show(io::IO, ::MIME"text/plain", diagnostics::NullDiagnostics)
    print(io, "NullDiagnostics: preserved = ", repr(diagnostics.preserved),
        "\n  samples = ", length(diagnostics.accepted_swaps),
        ", requested_swaps = ", diagnostics.requested_swaps,
        ", max_swap_attempts = ", diagnostics.max_swap_attempts,
        "\n  complete = ", diagnostics.complete)
end

function Base.show(io::IO, ::MIME"text/plain", result::CPSignificanceResult)
    print(io, "CPSignificanceResult: null_model = ", repr(result.null_model),
        "\n  samples = ", length(result.null_qualities),
        ", observed quality = ", _show_number(result.observed_quality),
        "\n  null mean = ", _show_number(result.null_mean),
        " ± ", _show_number(result.null_std),
        ", z-score = ", _show_number(result.zscore),
        "\n  p-value = ", _show_number(result.pvalue),
        ", significant at alpha = ", _show_number(result.alpha), ": ",
        result.significant)
end

function Base.show(io::IO, ::MIME"text/plain", result::StationaryDistributionResult)
    print(io, "StationaryDistributionResult: method = ", repr(result.method),
        "\n  nodes = ", length(result.distribution),
        ", iterations = ", result.iterations,
        ", residual = ", _show_number(result.residual),
        ", converged = ", result.converged)
end

function Base.show(io::IO, ::MIME"text/plain", result::RossaEnsembleResult)
    print(io, "RossaEnsembleResult",
        "\n  nodes = ", length(result.mean_coreness),
        ", runs = ", length(result.centralizations),
        ", unique profiles = ", result.unique_profiles,
        "\n  centralization = ", _show_number(result.mean_centralization),
        " ± ", _show_number(result.std_centralization),
        ", rank stability = ", _show_number(result.rank_stability),
        ", replicates stored = ", result.replicate_coreness !== nothing)
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
