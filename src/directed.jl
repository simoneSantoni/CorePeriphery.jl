export CPDirectedResult, minres_svd_directed, minres_symmetric

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

function _minres_symmetric_factor(A::AbstractMatrix{<:Real}, n::Int,
    max_iter::Int, tolerance::Float64, weight::Float64)
    strengths = Float64.(vec(sum(A, dims=2)))
    largest_strength = maximum(strengths)
    factor = largest_strength > 0.0 ?
        sqrt.(strengths ./ largest_strength) : zeros(Float64, n)
    candidate = similar(factor)
    product = similar(factor)
    for _ in 1:max_iter
        mul!(product, A, factor)
        factor_squared = sum(abs2, factor)
        change_squared = 0.0
        @inbounds for node in 1:n
            denominator = factor_squared - abs2(factor[node])
            fixed_point = denominator > 0.0 ? product[node] / denominator : 0.0
            candidate[node] = (1.0 - weight) * factor[node] + weight * fixed_point
            change_squared += abs2(candidate[node] - factor[node])
        end
        factor, candidate = candidate, factor
        sqrt(change_squared) <= tolerance && break
    end
    return factor, product
end

"""
    minres_symmetric(A; max_iter=10_000, tol=1e-8, damping=0.5)

Fit the symmetric off-diagonal rank-one MINRES model
`A[i,j] ≈ w[i]w[j]`. Simultaneous damped fixed-point updates make the estimator
equivariant under node permutation and independent of row/column orientation.
The returned coreness is `w` normalized to `[0,1]`; quality is the fraction of
squared off-diagonal adjacency weight explained by the fitted model.

This is a distinct estimand from `minres_svd_directed`, which fits separate sender
and receiver factors.
"""
function minres_symmetric(A::AbstractMatrix{<:Real};
    max_iter::Int=10_000,
    tol::Real=1e-8,
    damping::Real=0.5)
    n = _validate_adjacency(A; symmetric=true)
    n >= 2 || throw(ArgumentError("symmetric MINRES requires at least two nodes"))
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    isfinite(tol) && tol >= 0 ||
        throw(ArgumentError("tol must be finite and nonnegative"))
    isfinite(damping) && 0 < damping <= 1 ||
        throw(ArgumentError("damping must be in (0, 1]"))

    factor, product = _minres_symmetric_factor(
        A, n, max_iter, Float64(tol), Float64(damping),
    )
    residual = _directed_residual!(product, A, factor, factor)
    baseline = Float64(sum(abs2, A))
    quality = iszero(baseline) ? 1.0 : clamp(1.0 - residual / baseline, 0.0, 1.0)
    coreness = _normalize_scores(factor)
    threshold = median(coreness)
    core_nodes = findall(coreness .>= threshold)
    periphery_nodes = findall(coreness .< threshold)
    return CPResult(
        coreness,
        core_nodes,
        periphery_nodes,
        quality,
        "MINRES Symmetric Rank One",
    )
end

function _balance_factors!(u::Vector{Float64}, v::Vector{Float64})
    u_norm = norm(u)
    v_norm = norm(v)
    if iszero(u_norm) || iszero(v_norm)
        fill!(u, 0.0)
        fill!(v, 0.0)
        return u, v
    end

    common_norm = sqrt(u_norm * v_norm)
    u .*= common_norm / u_norm
    v .*= common_norm / v_norm
    return u, v
end

function _directed_als(A::AbstractMatrix{<:Real}, u_initial::Vector{Float64},
    v_initial::Vector{Float64}, max_iter::Int, tol::Float64;
    out_first::Bool)
    n = length(u_initial)
    u = copy(u_initial)
    v = copy(v_initial)
    u_previous = similar(u)
    v_previous = similar(v)
    transposed_A = transpose(A)

    for _ in 1:max_iter
        copyto!(u_previous, u)
        copyto!(v_previous, v)
        if out_first
            v_squared = sum(abs2, v)
            mul!(u, A, v)
            for i in 1:n
                denominator = v_squared - abs2(v[i])
                u[i] = denominator > 0.0 ? u[i] / denominator : 0.0
            end

            u_squared = sum(abs2, u)
            mul!(v, transposed_A, u)
            for i in 1:n
                denominator = u_squared - abs2(u[i])
                v[i] = denominator > 0.0 ? v[i] / denominator : 0.0
            end
        else
            u_squared = sum(abs2, u)
            mul!(v, transposed_A, u)
            for i in 1:n
                denominator = u_squared - abs2(u[i])
                v[i] = denominator > 0.0 ? v[i] / denominator : 0.0
            end

            v_squared = sum(abs2, v)
            mul!(u, A, v)
            for i in 1:n
                denominator = v_squared - abs2(v[i])
                u[i] = denominator > 0.0 ? u[i] / denominator : 0.0
            end
        end
        _balance_factors!(u, v)

        u_change = 0.0
        v_change = 0.0
        for i in 1:n
            u_change += abs2(u[i] - u_previous[i])
            v_change += abs2(v[i] - v_previous[i])
        end
        sqrt(max(u_change, v_change)) <= tol && break
    end
    return u, v
end

function _directed_residual!(work::Vector{Float64}, A::AbstractMatrix{<:Real},
    u::Vector{Float64}, v::Vector{Float64})
    mul!(work, A, v)
    baseline = Float64(sum(abs2, A))
    cross_term = dot(u, work)
    model_squared = sum(abs2, u) * sum(abs2, v)
    for i in eachindex(u, v)
        model_squared -= abs2(u[i] * v[i])
    end
    return max(0.0, baseline - 2.0 * cross_term + model_squared)
end

"""
    minres_svd_directed(A; max_iter=1000, tol=1e-6)

Fit the directed rank-one model `A[i,j] ≈ u[i]v[j]` over off-diagonal dyads
by paired alternating least-squares searches. Matrix-vector products and work
vectors are reused across iterations, including for sparse matrices. The returned out- and in-coreness
scores are independently normalized to `[0, 1]`; `residual` is the unscaled
least-squares residual and `quality` is the fraction of squared adjacency
weight explained by the fit.

The two possible update orders are transpose-dual and their balanced factors are
averaged, making the fit equivariant under transposition (which swaps out- and
in-coreness) and under applying the same node permutation to rows and columns.
"""
function minres_svd_directed(A::AbstractMatrix{<:Real};
    max_iter::Int=1000,
    tol::Real=1e-6)
    n = _validate_adjacency(A)
    max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))

    # Row and column strengths are transpose-dual, permutation-equivariant
    # initial values. Coupled balancing changes neither their outer product nor
    # the symmetry between the two factors.
    u = Float64.(vec(sum(A, dims=2)))
    v = Float64.(vec(sum(A, dims=1)))
    _balance_factors!(u, v)

    tolerance = Float64(tol)
    u_out_first, v_out_first = _directed_als(
        A, u, v, max_iter, tolerance; out_first=true)
    u_in_first, v_in_first = _directed_als(
        A, u, v, max_iter, tolerance; out_first=false)

    u = (u_out_first .+ u_in_first) ./ 2
    v = (v_out_first .+ v_in_first) ./ 2
    _balance_factors!(u, v)

    # Evaluate ||A - uv'||² over off-diagonal dyads without materializing uv'.
    product = similar(u)
    residual = _directed_residual!(product, A, u, v)
    baseline = Float64(sum(abs2, A))
    quality = iszero(baseline) ? 1.0 : 1.0 - residual / baseline

    out_coreness = _normalize_scores(u)
    in_coreness = _normalize_scores(v)
    coreness = (out_coreness .+ in_coreness) ./ 2
    threshold = median(coreness)
    core_nodes = findall(coreness .>= threshold)
    periphery_nodes = findall(coreness .< threshold)

    return CPDirectedResult(
        out_coreness,
        in_coreness,
        coreness,
        core_nodes,
        periphery_nodes,
        residual,
        quality,
        "MINRES/SVD Directed",
    )
end
