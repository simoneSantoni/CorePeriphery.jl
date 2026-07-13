# Algorithm and Result Contracts

CorePeriphery.jl uses common result containers where useful, but the numerical meaning
of a score or quality value remains algorithm-specific. In particular, raw `quality`
values must not be used to rank different algorithms unless they optimize the same
objective on the same scale.

## Detector contracts

| Detector | Input contract | `coreness` meaning | Partition rule | `quality` meaning | Randomness | Main computational cost |
|---|---|---|---|---|---|---|
| `borgatti_everett_continuous` | Nonnegative, loop-free, undirected weighted | Optimized continuous product-pattern coordinate | Median score | Pearson correlation with `c*c'` over undirected dyads | None | Cached coordinate sweeps |
| `borgatti_everett_discrete` | Nonnegative, loop-free, undirected weighted | Binary discrete role | Returned binary state | Pearson correlation with `max(c_i,c_j)` | None | Greedy label sweeps |
| `lip_discrete` | Simple undirected binary | Binary degree-prefix role | Exact minimizing prefix; smallest tied prefix | Negative published Lip loss | None | Sort plus linear prefix scan, `O(m+n log n)` |
| `rombach_continuous` | Nonnegative, loop-free, undirected weighted | Fixed Rombach rank template assigned to nodes | Median template value | Weighted template-product objective | Explicit RNG across restarts | Pair-swap sweeps over template permutations |
| `spectral_method` | Simple undirected binary | Normalized thresholded rank-two reconstruction row score | Published Find-Cut optimum, not a score threshold | Core plus cross density minus periphery density | None | Two eigenpairs plus `O(n^2)` streamed reconstruction |
| `random_walker_profiling` | Connected undirected or strongly connected directed, nonnegative, loop-free | Persistence when the node enters the growing periphery | Median persistence for compatibility | Della Rossa equation-(8) cp-centralization | Node-index default or explicit paper RNG | `O(n^2)` greedy scan plus directed stationary solve |
| `minres_symmetric` | Nonnegative, loop-free, undirected weighted | Normalized symmetric rank-one factor | Median score | Fraction of squared off-diagonal weight explained by `w*w'` | None | Matrix-vector fixed-point iterations |
| `minres_svd_directed` | Nonnegative, loop-free, directed or undirected weighted | Separate normalized sender, receiver, and mean factors | Median combined score | Fraction of squared off-diagonal weight explained by `u*v'` | None | Paired alternating least-squares matrix-vector iterations |
| `minres_svd` | Same as directed MINRES | Compatibility projection of mean sender/receiver score | Directed MINRES partition | Directed MINRES explained squared weight | None | Same as `minres_svd_directed` |
| `multiple_cp_pairs` | Nonnegative, loop-free, undirected weighted | Binary role within an explicit pair | Joint pair and role state | Kojaku--Masuda residual objective under ER or configuration null | Explicit RNG across restarts | Local adjacency-list moves across pair counts |
| `multiple_cp_pairs_config` | Same as multiple pairs | Same | Same | Configuration-null specialization | Explicit RNG | Same |
| `surprise_cp` | Simple undirected binary | Binary role | Returned optimized state | Negative log multivariate-hypergeometric joint-tail probability | None | Cached block-count label sweeps |
| `label_switching_cp` | Nonnegative, loop-free, undirected weighted | Binary role | Returned optimized state | Discrete Borgatti--Everett Pearson correlation | Explicit RNG and optional restarts | Cached label-switching sweeps |

`rossa_profile_ensemble` is an uncertainty wrapper rather than another estimand. It
summarizes repeated paper-mode Rossa profiles, centralization, rank stability, and raw
replicates when requested. Its preflight `memory_limit_bytes` bounds the `O(n*n_runs)`
profile workspace (512 MiB by default); `store_replicates=false` omits that workspace
from the returned result after summaries are formed. `rossa_stationary_distribution` is a numerical helper that
returns the distribution, L1 balance residual, iteration count, and selected solver.

`cp_significance` is also a wrapper. It compares one detector quality with a null
distribution using a conservative plus-one Monte Carlo p-value. Its supported nulls
are undirected ER, undirected degree-preserving configuration, directed in/out-degree
configuration, and topology-preserving weight permutation.

The switching null follows the prescribed-degree switching family discussed by Milo
et al. (2003), including its important caveat that a finite switch chain has no general
mixing-time guarantee. The weight-permutation null is the fixed-topology weight
reshuffling design used to remove weight--topology association. It preserves the
global edge-weight multiset, **not** individual node strengths; use it only when that
conditional null is the scientific question.

## Result containers

### `CPResult`

- `coreness` is the method-specific node score described above.
- `core_nodes` and `periphery_nodes` are the method-specific partition. They are not
  necessarily obtained by thresholding `coreness`; LowRank-Core is the key example.
- `quality` is method-specific and is only comparable within a common objective.
- `algorithm` identifies the estimator.

### `CPDirectedResult`

`out_coreness` is the sender factor, `in_coreness` is the receiver factor, and
`coreness` is their arithmetic mean after separate normalization. `residual` is the
unscaled off-diagonal least-squares residual.

### `CPMultiResult`

In addition to pair labels and binary within-pair roles, this result records
`pair_selection`, every evaluated `candidate_pair_counts`, raw
`candidate_qualities`, and the `selection_scores` used to choose model complexity.

### `CPSignificanceResult`

`null_diagnostics` states the constraint preserved by the selected null model. For
switching nulls it also reports requested and accepted swaps, attempts, the attempt
budget, and whether every sample completed. Incomplete switching is an error by
default; `swap_shortfall=:warn` or `:accept` is an explicit opt-in.

## Directionality guide

Use directed Rossa when directionality means stationary flow and persistence inside a
node subset. The graph must be strongly connected. Use `minres_svd_directed` when
directionality means distinct sender and receiver roles in an asymmetric rank-one
model. Use `minres_symmetric` for the symmetric `w*w'` estimand; do not compare it
directly with the directed `u*v'` model as though they were the same optimizer.

## Reproducible examples

```jldoctest contracts
julia> using CorePeriphery, Random, SparseArrays

julia> A = adjacency_to_matrix([(1, 2), (1, 3), (1, 4), (2, 3)], 4);

julia> dense = lip_discrete(A); sparse_result = lip_discrete(sparse(A));

julia> @assert dense.coreness == sparse_result.coreness

julia> paper = random_walker_profiling(A; tie_break=:paper, rng=MersenneTwister(7));

julia> @assert isfinite(paper.quality)
```

```jldoctest directed-contracts
julia> using CorePeriphery

julia> D = [0.0 2 0 1; 0 0 3 0; 4 0 0 1; 1 2 0 0];

julia> stationary = rossa_stationary_distribution(D; method=:linear);

julia> @assert stationary.residual <= 1e-12

julia> directed = minres_svd_directed(D);

julia> @assert length(directed.out_coreness) == 4
```

## Edge cases and compatibility

- Rossa cp-centralization is undefined for a two-node graph and is returned as `NaN`.
- Rossa keywords `n_walks` and `walk_length` are deprecated compatibility arguments and
  do not alter the analytic persistence result.
- Deterministic Rossa final ties use node index and are therefore label-dependent. Use
  paper mode or the ensemble wrapper when tied order is scientifically meaningful.
- Sparse storage preserves estimator semantics but is not guaranteed to reduce latency
  on small graphs.
- Julia 1.10 is the supported compatibility floor; Julia 1.12 is the primary optimized
  and documentation runtime.

## Null-model references

- Milo, R., Kashtan, N., Itzkovitz, S., Newman, M. E. J., and Alon, U. (2003).
  *On the uniform generation of random graphs with prescribed degree sequences*.
  [arXiv:cond-mat/0312028](https://arxiv.org/abs/cond-mat/0312028).
- Opsahl, T., Colizza, V., Panzarasa, P., and Ramasco, J. J. (2008).
  *Prominence and control: the weighted rich-club effect*.
  [Physical Review Letters 101, 168702](https://doi.org/10.1103/PhysRevLett.101.168702).
