# Algorithm and Result Contracts

CorePeriphery.jl uses common result containers where useful, but the numerical meaning
of a score or quality value remains algorithm-specific. In particular, raw `quality`
values must not be used to rank different algorithms unless they optimize the same
objective on the same scale.

## Detector contracts

| Detector | Input contract | `coreness` meaning | Partition rule | `quality` meaning | Randomness | Main computational cost |
|---|---|---|---|---|---|---|
| `borgatti_everett_continuous` | Nonnegative, loop-free, undirected weighted; asymmetric input with `directed=true` | Optimized continuous product-pattern coordinate in `[0, 1]` | `partition` keyword, default `:median` (strictly above the median) | Pearson correlation with `c*c'` over undirected dyads (all ordered off-diagonal dyads when directed) | None | Cached coordinate sweeps |
| `borgatti_everett_discrete` | Nonnegative, loop-free, undirected weighted; asymmetric input with `directed=true` | Binary discrete role | Returned binary state | Pearson correlation with the `cp_block` pattern: `:ones` is `max(c_i,c_j)` (BE eq. 2), `:missing` correlates only core-core and periphery-periphery dyads (BE eq. 4) | None with `n_runs=1`; explicit RNG for restarts; `method=:exhaustive` (`n ≤ 24`) is a deterministic global optimum | Greedy label sweeps, or `O(2^n · n)` Gray-code enumeration |
| `lip_discrete` | Simple undirected binary | Binary degree-prefix role | Exact minimizing prefix; smallest tied prefix | Discrete Pearson correlation of the returned partition (the Lip loss is minimised but not returned) | None | Sort plus linear prefix scan, `O(m+n log n)` |
| `rombach_continuous` | Nonnegative, loop-free, undirected weighted | Shuffle of the Rombach eq. (4.9) template `C*(α, β)` with `⌊βN⌋` peripheral ranks, normalised to sum to `1` | `partition` keyword, default `:beta` (the `N − ⌊βN⌋` core ranks of the template) | Eq. (4.10) core quality `R = Σ_{i,j} A_ij C_i C_j` (full double sum) | Explicit RNG across restarts | Pair-swap sweeps over template permutations |
| `rombach_core_score` | Same as `rombach_continuous` | Eq. (4.11) aggregate core score over the `alphas × betas` grid, scaled so the maximum is `1` | `partition` keyword, default `:median` | Mean `R_γ` over the grid | Explicit RNG, consumed in grid order | One `rombach_continuous` optimisation per grid point |
| `spectral_method` | Simple undirected binary | Normalized thresholded rank-two reconstruction row score | Published Find-Cut optimum, not a score threshold | Core plus cross density minus periphery density | None | Two eigenpairs (LAPACK dense, ARPACK sparse) plus a thresholded reconstruction pruned to the rows that can exceed the threshold (quadratic in their number) |
| `random_walker_profiling` | Connected undirected or strongly connected directed, nonnegative, loop-free | Persistence when the node enters the growing periphery | `partition` keyword, default `:median` | Della Rossa equation-(8) cp-centralization | Node-index default or explicit paper RNG | `O(n^2)` greedy scan plus directed stationary solve |
| `minres_symmetric` | Nonnegative, loop-free, undirected weighted | Normalized symmetric rank-one factor | `partition` keyword, default `:median` | Fraction of squared off-diagonal weight explained by `w*w'` | None | Matrix-vector fixed-point iterations |
| `minres_svd_directed` | Nonnegative, loop-free, directed or undirected weighted | Separate normalized sender, receiver, and mean factors | `partition` keyword, default `:median` | Fraction of squared off-diagonal weight explained by `u*v'` | None | Paired alternating least-squares matrix-vector iterations |
| `minres_svd` | Same as directed MINRES | Compatibility projection of mean sender/receiver score | Directed MINRES partition | Directed MINRES explained squared weight | None | Same as `minres_svd_directed` |
| `multiple_cp_pairs` | Nonnegative, loop-free, undirected weighted | Binary role within an explicit pair (`0` for isolates, which get pair label `0`) | Joint pair and role state; pair count chosen by `pair_selection` (default `:penalized`) | Kojaku--Masuda (2017) eq. (3) / `m` under the ER null, or Kojaku--Masuda (2018) eq. (13) exactly under the configuration null | Explicit RNG across restarts | Local adjacency-list moves across pair counts |
| `multiple_cp_pairs_config` | Same as multiple pairs | Same | Same | Configuration-null specialization (eq. 13) | Explicit RNG | Same |
| `surprise_cp` | Simple undirected binary | Binary role | Returned optimized state | Negative log multivariate-hypergeometric joint-tail probability | None | Label sweeps with cached block counts; each tail evaluation costs roughly the product of the two hypergeometric peak widths |
| `label_switching_cp` | Nonnegative, loop-free, undirected weighted; asymmetric input with `directed=true` | Binary role | Returned optimized state | Discrete Borgatti--Everett Pearson correlation under the chosen `cp_block` | Explicit RNG and optional restarts | Cached label-switching sweeps |

Row notes:

- `rombach_continuous`: `beta` is the fraction of nodes in the **periphery** (paper
  semantics): the first `⌊βN⌋` template ranks are peripheral, so the core shrinks from
  `N` nodes at `beta = 0` to none at `beta = 1`. `alpha` sets the size of the score jump
  between the highest peripheral and lowest core rank. `coreness` sums to `1` (multiply
  by the template sum to recover raw `[0, 1]` template values). `quality` is the eq.
  (4.10) objective of one `(α, β)` point; it depends on the template as well as on the
  data, so do not select `(α, β)` by maximising it. The paper's deliverable is the
  grid-averaged `rombach_core_score`.
- `lip_discrete`: the optimiser minimises the published Lip loss, but the returned
  `quality` is the discrete Borgatti--Everett correlation of that partition, so it can
  be negative on inputs without core-periphery structure.
- `spectral_method`: `beta` is the paper's lower bound `b/n` on the size of either side
  of the Find-Cut, not a prescribed core fraction. Sparse input uses ARPACK
  (`tol = 1e-10`, `ncv = 48`, `maxiter = max(300, 20n)`) and raises
  `ConvergenceError` when the two eigenpairs do not converge; dense input uses LAPACK.
- `borgatti_everett_discrete`: `cp_block=:missing` is the model Borgatti and Everett
  recommend and UCINET's categorical default; `:ones` remains the package default for
  continuity. `method=:exhaustive` returns a global optimum (first met in enumeration
  order among ties) and ignores `n_runs`, `rng`, `init`, and `max_iter`.
- Continuous detectors: `core_nodes` is the split produced by the `partition` keyword,
  not an estimate of the model. `:median` (default) puts nodes **strictly above** the
  median coreness in the core, so nodes tied at the median (a star's leaves, for
  instance) stay peripheral; `:none` returns an empty core; `:threshold` uses
  `threshold` (strictly above); `:core_size` takes the `core_size` highest scores,
  ties broken by lower node index. `rombach_continuous` additionally offers `:beta`,
  its default. Borgatti and Everett (2000) themselves do not partition in the
  continuous model, so report the scores rather than the split when the split is not
  an explicit analytic choice.
- Isolates and components: every detector except `surprise_cp` and `spectral_method`
  gives isolates (zero in- and out-strength) coreness `0` and never places them in
  `core_nodes` (`multiple_cp_pairs` gives them pair label `0`); an edgeless network
  returns all-zero coreness and an empty core. The same detectors emit one warning per
  detector per session when the input has more than one connected component with at
  least two nodes; `random_walker_profiling` instead requires connected input.

`rossa_profile_ensemble` is an uncertainty wrapper rather than another estimand. It
summarizes repeated paper-mode Rossa profiles, centralization, rank stability, and raw
replicates when requested. Its preflight `memory_limit_bytes` bounds the `O(n*n_runs)`
profile workspace (512 MiB by default); `store_replicates=false` omits that workspace
from the returned result after summaries are formed. Passing `tie_break` or `rng`
through its forwarded keywords raises `ArgumentError` (the ensemble always uses
`:paper` and seeds the profiles from its own `rng`). `rossa_stationary_distribution`
is a numerical helper that returns the distribution, L1 balance residual, iteration
count, and selected solver.

`cp_significance` is also a wrapper. It compares one detector quality with a null
distribution using a conservative plus-one Monte Carlo p-value over `n_samples=999`
samples by default, and reports `null_mean`, `null_std`, and `zscore` as effect sizes.
Its supported nulls are undirected ER, undirected degree-preserving configuration,
directed in/out-degree configuration, and topology-preserving weight permutation. The
first three require a binary network; non-binary input raises an error naming
`:weight_permutation` as the weighted alternative. Null samples are drawn in the
storage format of the input (sparse input gives sparse samples) and dense and sparse
inputs produce identical results for the same `rng`. With `pass_rng=true` (default)
each detector call that accepts an `rng` keyword receives a child RNG derived from
`rng`. With `threaded=true` the samples are distributed over a work queue with one task
per thread (`thread_schedule=:auto` or `:greedy`) or round-robin (`:static`); every
schedule is bit-for-bit identical to the serial run.

The switching null follows the prescribed-degree switching family discussed by Milo
et al. (2003), including its important caveat that a finite switch chain has no general
mixing-time guarantee. The weight-permutation null is the fixed-topology weight
reshuffling design used to remove weight--topology association. It preserves the
global edge-weight multiset, **not** individual node strengths; use it only when that
conditional null is the scientific question.

## Result containers

All result types (and `CorePeriphery.NullDiagnostics`) define structural `==`,
`isequal`, and `hash` over every field, and a compact `show`. `==` follows the field
semantics, so a result whose `quality` or `zscore` is `NaN` is not `==` to itself;
compare such results with `isequal`.

### `CPResult`

- `coreness` is the method-specific node score described above. Isolates score `0`.
- `core_nodes` and `periphery_nodes` are the method-specific partition. They are not
  necessarily obtained by thresholding `coreness`; LowRank-Core is the key example,
  and for the continuous detectors they follow the `partition` keyword. Isolates are
  never in `core_nodes`.
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
Isolates carry pair label `0` and coreness `0`; `n_pairs` counts the pairs among the
non-isolated nodes, whose labels are `1:n_pairs`.

### `CPSignificanceResult`

`pvalue` is the plus-one Monte Carlo estimate `(1 + #{null ≥ observed}) / (n_samples + 1)`.
`null_mean`, `null_std` (sample standard deviation, `NaN` with a single sample), and
`zscore` (`(observed − null_mean) / null_std`, `NaN` when the null is degenerate) are
effect sizes to report alongside the p-value. `null_diagnostics` is a
`CorePeriphery.NullDiagnostics` stating the constraint preserved by the selected null
model. For switching nulls it also reports requested and accepted swaps, attempts, the
attempt budget, and whether every sample completed. Incomplete switching raises
`ConvergenceError` by default; `swap_shortfall=:warn` or `:accept` is an explicit
opt-in.

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

- Input validation failures raise `ArgumentError` or `DimensionMismatch`. Numerical
  procedures that do not reach the requested accuracy raise `ConvergenceError`: an
  ARPACK eigensolve in sparse `spectral_method`, a Rossa stationary solver that
  exhausts its iteration budget or tolerance, and a configuration-null switching chain
  that cannot complete the requested swaps under `swap_shortfall=:error`.
- Symmetry is tested with an absolute tolerance of `1e-12 · max(1, max|A|)`, so
  floating-point noise such as `(A + A') / 2` does not make an undirected network look
  directed. `Symmetric`, `Hermitian`, `Transpose`, and `Adjoint` wrappers around a
  sparse matrix are unwrapped to `SparseMatrixCSC` at entry so the sparse kernels are
  used.
- Rossa cp-centralization is undefined for a two-node graph and is returned as `NaN`.
- Rossa keywords `n_walks` and `walk_length` and Lip's `max_iter` are deprecated
  compatibility arguments: they do not alter the result and emit `Base.depwarn`
  (visible with `--depwarn=yes`).
- Deterministic Rossa final ties use node index and are therefore label-dependent. Use
  paper mode or the ensemble wrapper when tied order is scientifically meaningful.
- Numeric keywords accept any `Real` (`tol=1//1000`, `alpha=1f-2`) and
  `ideal_cp_matrix` accepts any real vector.
- `adjacency_to_matrix` rejects a dyad that occurs more than once with conflicting
  weights unless `combine=:sum|:max|:min|:last` is given; `directed=true` keeps `(i, j)`
  and `(j, i)` distinct. The method without `n` accepts arbitrary node identifiers and
  returns `(A, labels)`.
- Sparse storage preserves estimator semantics but is not guaranteed to reduce latency
  on small graphs.
- Julia 1.10 is the supported compatibility floor; Julia 1.13 is the primary
  development runtime and CI tests Julia 1.10, 1.12, and 1.13. The benchmark snapshots
  on this site were measured on Julia 1.10 and 1.12.

## Null-model references

- Milo, R., Kashtan, N., Itzkovitz, S., Newman, M. E. J., and Alon, U. (2003).
  *On the uniform generation of random graphs with prescribed degree sequences*.
  [arXiv:cond-mat/0312028](https://arxiv.org/abs/cond-mat/0312028).
- Opsahl, T., Colizza, V., Panzarasa, P., and Ramasco, J. J. (2008).
  *Prominence and control: the weighted rich-club effect*.
  [Physical Review Letters 101, 168702](https://doi.org/10.1103/PhysRevLett.101.168702).
