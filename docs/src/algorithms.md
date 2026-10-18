# Algorithms

CorePeriphery.jl implements ten published core-periphery detection methods (the MINRES method has distinct symmetric and directed estimators) plus Monte Carlo significance testing. This page provides detailed descriptions of each method, their theoretical foundations, and guidance on when to use them. Full citations are in the [References](@ref algorithm-references) section at the end of the page.

## Algorithm Categories

The algorithms can be grouped into three categories:

1. **Continuous methods**: Assign each node a coreness score in [0, 1]
2. **Discrete methods**: Classify nodes as either core (1) or periphery (0)
3. **Multi-pair methods**: Detect multiple non-overlapping core-periphery pairs

## Borgatti-Everett Continuous Model

```julia
result = borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, step=0.05, init=nothing,
                                     directed=false, partition=:median)
```

The classic continuous core-periphery model introduced by Borgatti and Everett (2000). It finds a coreness vector ``c`` that maximizes the correlation between the adjacency matrix ``A`` and an ideal core-periphery pattern:

```math
\Delta_{ij} = c_i \times c_j
```

**Algorithm**: Bounded grid coordinate ascent on the stated Pearson objective.
Only objective-improving updates are accepted.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)
- `tol`: Convergence tolerance (default: 1e-6)
- `step`: Coordinate-search resolution (default: 0.05)
- `init`: Initial coreness vector (default: degree centrality)
- `directed`: Accept asymmetric input and correlate over all ordered off-diagonal
  dyads, as UCINET does for directed sociometric data (default: `false`)
- `partition`, `threshold`, `core_size`: Split rule for `core_nodes`: `:median`
  (default, strictly above the median), `:none`, `:threshold`, or `:core_size`

Borgatti and Everett do not partition in the continuous model; `core_nodes` is a
reporting convenience. Report the scores, the fit, and the block densities
(`block_densities`) rather than the split alone.

**When to use**: When you need continuous coreness scores and want a well-established baseline method.

## Borgatti-Everett Discrete Model

```julia
result = borgatti_everett_discrete(A; max_iter=1000, init=nothing, directed=false,
                                   cp_block=:ones, n_runs=1, rng=Random.default_rng(),
                                   method=:greedy)
```

The discrete version of the Borgatti-Everett model. Finds a binary partition maximizing correlation with the ideal discrete pattern. With `cp_block=:ones` (default) the pattern is their equation (2),

```math
\Delta_{ij} = \max(c_i, c_j)
```

where ``c_i \in \{0, 1\}``, so core-periphery ties are expected to be present. With
`cp_block=:missing` the pattern is their equation (4), the model Borgatti and Everett
recommend and UCINET's categorical default: core-periphery dyads are treated as missing
data and only the core-core (pattern 1) and periphery-periphery (pattern 0) dyads enter
the correlation.

**Algorithm**: Greedy single-node flipping from a degree-based start (`method=:greedy`),
optionally with `n_runs` random restarts, or Gray-code enumeration of every core set
(`method=:exhaustive`, `n ≤ 24`), which returns a global optimum. Greedy flipping misses
the optimum on a few percent of small graphs, so prefer `:exhaustive` for small
sociograms.

**Parameters**:

- `max_iter`: Maximum sweeps per start (default: 1000)
- `init`: Initial binary assignment (default: degree-based)
- `directed`: Accept asymmetric input and correlate over all ordered off-diagonal dyads
- `cp_block`: `:ones` (equation 2) or `:missing` (equation 4)
- `n_runs`, `rng`: Random restarts and their generator (`n_runs=1` is deterministic)
- `method`: `:greedy` or `:exhaustive`

**When to use**: When you need a binary classification and want to maximize correlation with the ideal discrete pattern.

## Lip's Fast Discrete Algorithm

```julia
result = lip_discrete(A)
```

The exact binary algorithm introduced by Lip (2011). It minimizes missing
core-core edges plus observed periphery-periphery edges over degree-sorted
prefix partitions.

**Algorithm**: Stable degree sort followed by a linear scan of all nontrivial
prefix sizes. Input must be simple, undirected, and binary.

**Quality**: The returned `quality` is the discrete Borgatti-Everett Pearson
correlation of the returned partition. The Lip loss is minimised internally but is not
returned, so `quality` can be negative even though the loss is at its minimum (two
disjoint cliques are an example).

**When to use**: When you need an exact, fast discrete classification of a simple
binary graph. Lip minimises an error count rather than the correlation, so its
partition can differ from `borgatti_everett_discrete` and `label_switching_cp`, which
maximise the correlation directly.

## Rombach's Generalized Model

```julia
result = rombach_continuous(A; alpha=0.5, beta=0.5, n_runs=10, rng=Random.default_rng(),
                            partition=:beta)
score = rombach_core_score(A; alphas=0.1:0.1:1.0, betas=0.1:0.1:1.0, n_runs=3)
```

A generalized continuous model by Rombach et al. (2017) with parameters controlling the core-periphery structure:

- ``\alpha``: The size of the score jump between the highest-scoring peripheral node and
  the lowest-scoring core node (``\alpha = 1`` is a discrete classification)
- ``\beta``: The fraction of nodes in the **periphery** (paper semantics): the first
  ``\lfloor \beta N \rfloor`` template ranks are peripheral, so the core shrinks from
  ``N`` nodes at ``\beta = 0`` to none at ``\beta = 1``

The parameters define the transition-function template of the paper's equation (4.9),
implemented exactly with the ``\lfloor \beta N \rfloor`` boundary and normalised so the
core vector sums to one, as in the paper. The algorithm assigns a permutation of that
template to nodes.

**Algorithm**: Seedable random restarts with objective-improving pair swaps,
maximizing the equation-(4.10) core quality ``R_{\alpha,\beta} = \sum_{i,j} A_{ij} C_i C_j``
(the full double sum, so each undirected edge counts twice).

The returned `coreness` is the optimised normalised template for one
``(\alpha, \beta)`` pair, and the default `partition=:beta` reports the
``N - \lfloor \beta N \rfloor`` nodes holding the template's core ranks as `core_nodes`,
so `beta` sets the core size directly. The `quality` is ``R_{\alpha,\beta}`` of the
returned vector; it depends on the template as well as on the data, so it must not be
used to select ``(\alpha, \beta)``. The paper's recommended deliverable is the aggregate
core score of equation (4.11), ``CS(i) = Z \sum_\gamma C_i(\gamma) R_\gamma`` over a
grid of ``\gamma = (\alpha, \beta)`` values, which `rombach_core_score` provides (the
paper's ``100 \times 100`` grid is reduced to ``10 \times 10`` by default; refine with
`alphas`/`betas`).

**Parameters**:

- `alpha`: Boundary sharpness, 0 to 1 (default: 0.5)
- `beta`: Periphery fraction, 0 to 1 (default: 0.5); larger `beta`, smaller core
- `n_runs`: Number of random restarts (default: 10; 3 per grid point in `rombach_core_score`)
- `partition`, `threshold`, `core_size`: Split rule; `:beta` (default), `:median`,
  `:none`, `:threshold`, or `:core_size` (`rombach_core_score` defaults to `:median`)

**When to use**: When you want to tune the expected core-periphery structure, explore how results change with different boundary definitions, or report the parameter-averaged core score.

## Spectral Method

```julia
result = spectral_method(A)
```

Implements Cucuringu et al.'s LowRank-Core algorithm: reconstruct the adjacency
matrix from its two eigenpairs of largest magnitude, threshold the rank-two
matrix, rank nodes by its row sums, and apply Find-Cut.

**Algorithm**: Rank-two symmetric eigendecomposition plus the published
core/core-periphery/periphery density cut objective. The returned `CPResult`
contains the LowRank-Core row-score coreness values, while `core_nodes` and
`periphery_nodes` come from the best Find-Cut prefix rather than a threshold on
those scores. Input is binary, undirected, and loop-free.

**Parameters**: `beta` is the paper's lower bound ``b/n`` on the fraction of nodes allowed on either side of the Find-Cut (default 0.1), not the prescribed core fraction ``\beta`` of the paper's known-core-size variant, which is not implemented.

**Convergence**: Dense input uses a full LAPACK eigendecomposition. Sparse input
computes only the two largest-magnitude eigenpairs with ARPACK (`tol = 1e-10`,
`ncv = 48` Krylov vectors, `maxiter = max(300, 20n)` restarts; calls are serialised
through a lock because ARPACK is not reentrant). If fewer than two eigenpairs converge,
a `ConvergenceError` reports the settings; clustered extremal eigenvalues (a long ring,
for instance) are the usual cause, and a dense matrix is the fallback. The thresholded
reconstruction is pruned by a per-row bound, so only rows that can exceed the ``1/2``
threshold are paired; on sparse core-periphery networks that is a small fraction of
``n``, but the worst case stays quadratic in the number of active rows.

**When to use**: Fast method for continuous coreness scores. Works well when the core-periphery structure aligns with the network's spectral properties. No greedy optimization is needed; the sparse path is an iterative eigensolve.

## Random Walker Profiling

```julia
result = random_walker_profiling(A)
```

Implements Della Rossa et al.'s persistence profile. Starting from a
minimum-strength node, it grows a nested peripheral set by selecting the node
that minimizes the new set's persistence probability.

**Algorithm**: Greedy persistence profiling for undirected and directed networks.
The default `tie_break=:deterministic` replaces the paper's random choices with lower
node index while retaining its persistence and secondary strength criteria. Use
`tie_break=:paper` and an explicit `rng` to sample exactly according to the published
tie hierarchy. `quality` is the paper's equation (8) cp-centralization.

Input must be connected when undirected, strongly connected when directed, loop-free,
and nonnegative. Directed inputs use total in-plus-out strength for the tie hierarchy
and the stationary distribution of the row-normalized random walk for persistence.
`stationary_tol` and `stationary_max_iter` control its sparse-friendly lazy iteration.
Use `stationary_method=:linear` for a robust dense constrained solve or `:auto` to
select by storage and size. `rossa_stationary_distribution` exposes convergence
residual and iteration metadata, while `rossa_profile_ensemble` summarizes uncertainty
across paper-mode tie realizations.

**When to use**: When you want a dynamics-based definition of coreness. Useful for understanding network flow and accessibility.

## MINRES/SVD Method

```julia
result = minres_svd(A; max_iter=1000, tol=1e-6)
```

Minimizes the residual between the adjacency matrix and an outer product model (Boyd et al., 2010):

```math
\min_{u,v} \sum_{i \neq j} (A_{ij} - u_i v_j)^2
```

This yields separate in-coreness (``v``) and out-coreness (``u``) vectors.
`minres_svd_directed` retains both, their combined scores, the residual, and the
explained squared-weight quality. `minres_svd` returns a combined compatibility result.

**Algorithm**: Alternating least squares optimization.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)
- `tol`: Convergence tolerance (default: 1e-6)
- `partition`, `threshold`, `core_size`: Split rule for `core_nodes` (default `:median`)

**When to use**: For directed (asymmetric) outer-product structure. The current
`CPResult` reports the average of in- and out-coreness; directed quality uses
all ordered off-diagonal dyads.

For a symmetric rank-one model, use:

```julia
result = minres_symmetric(A; max_iter=10_000, tol=1e-8, damping=0.5, partition=:median)
```

This fits ``A_{ij}\approx w_iw_j`` and is a different estimand from the directed
``u_iv_j`` model.

## Multiple Core-Periphery Pairs

```julia
result = multiple_cp_pairs(A; max_pairs=10, min_pair_size=2, max_iter=100,
                           null_model=:er, pair_selection=:penalized)
```

Detects multiple non-overlapping core-periphery pairs using the ``Q^{cp}`` quality function (Kojaku & Masuda, 2017, 2018).

**Algorithm**:

1. Create balanced seeded initial pair assignments
2. Jointly propose each node as core or periphery in its own or a neighboring pair
3. Accept only objective-improving moves across multiple restarts

Returns a `CPMultiResult` with pair labels, within-pair coreness, the chosen
pair count, and the candidate qualities inspected during model selection. Isolates
are removed before fitting and reported with pair label `0`.

**Objectives**: Under `null_model=:er` the objective is Kojaku & Masuda (2017)
equation (3) divided by ``m``. Under `:configuration` it is Kojaku & Masuda (2018)
equation (13) exactly, including the paper's ``i = j`` null term
``-d_i^2 / (2m)`` for every core node (so the objective reduces to modularity when every
node is core, as the paper states); both are checked against literal transcriptions of
the equations in the test suite. Earlier development versions omitted the ``i = j``
term, which lowered the configuration objective by ``\sum_{i \in \text{core}} d_i^2 / (4m^2)``.

**Package variants**: The objectives and move gains follow the papers, but
the search and model-selection procedures are package variants rather than the
published heuristics. The optimiser seeds balanced initial assignments for each
candidate pair count from 1 to `max_pairs` instead of growing pairs from singleton
cores, and the pair count is chosen either by the raw objective (`:objective`) or by an
unpublished penalty `pair_penalty * (k - 1) * log(n) / n` (`:penalized`, the default,
because the raw objective is nondecreasing in the pair count in practice and
over-partitions a single planted pair). The per-pair statistical test of Kojaku &
Masuda (2018), which discards pairs that are not significant under the configuration
model, is not implemented. Kojaku & Masuda (2018) recommend the configuration null
because the ER null finds core-periphery pairs in any network with heterogeneous
degrees; `:er` remains the default for continuity with the 2017 paper, and
`multiple_cp_pairs_config` selects the configuration null directly.

**Parameters**:

- `max_pairs`: Maximum number of pairs (default: 10)
- `min_pair_size`: Minimum nodes per pair (default: 2)
- `max_iter`: Optimization iterations (default: 100)
- `n_runs`: Seeded restarts (default: 10)
- `null_model`: `:er` (default) or `:configuration` (the only nulls accepted here; `:directed_configuration` and `:weight_permutation` exist only in `cp_significance`)
- `rng`: Explicit random-number generator
- `pair_selection`: `:penalized` (default) for complexity control or `:objective` for the raw objective
- `pair_penalty`: Multiplier for ``(k-1)\log(n)/n`` (default: 0.5)

**When to use**: When you expect the network to contain multiple distinct core-periphery structures rather than a single global one.

Use `multiple_cp_pairs_config(A)` as a convenience wrapper for the degree-preserving
configuration null.

## Significance Testing

```julia
result = cp_significance(A, lip_discrete; null_model=:configuration, n_samples=999)
```

Compares one detector's observed `quality` against a Monte Carlo null
distribution of `n_samples` (default 999) random networks. The wrapper accepts any
detector that returns a finite `quality` field and supports four null models:

- `:er` for undirected Erdős-Rényi sampling
- `:configuration` for undirected degree-preserving edge switching
- `:directed_configuration` for directed in/out-degree-preserving arc switching
- `:weight_permutation` for topology-preserving weight shuffling

The first three preserve properties of a simple graph and therefore require binary
input; a weighted network raises an error that names `:weight_permutation` as the
alternative. The result reports the plus-one p-value and the effect sizes
`null_mean`, `null_std`, and `zscore`.

`pass_rng=true` (default) forwards a fresh child RNG, derived from `rng`, to every
detector call whose method accepts an `rng` keyword, so stochastic detectors give
reproducible null distributions; detectors defined outside the package that lack an
`rng` keyword trigger a one-time warning. Null samples are drawn in the storage format
of `A`, so a sparse input is sampled directly into sparse form and evaluated with the
sparse kernels; dense and sparse inputs give identical results. With `threaded=true`,
the null samples are distributed across threads from a precomputed seed stream:
`thread_schedule=:auto` (the default) and `:greedy` use a work queue with one task per
thread, `:static` assigns samples round-robin, and every schedule reproduces the
serial result bit for bit on every supported Julia version.

Configuration-style nulls require the requested swap budget to complete by
default; a shortfall raises `ConvergenceError` with an explanation (stars and other
threshold graphs admit no double-edge swap at all). Use `max_swap_attempts` to raise
the proposal budget, or set `swap_shortfall=:warn` or `:accept` if an incomplete chain
is scientifically acceptable and you want that fact recorded in `null_diagnostics`.

### Which null model answers which question

The null model decides what "significant" means, and the choice can flip the answer:

- `:er` preserves only the edge count. A small p-value says the network is more
  core-periphery-like than a random graph of the same density.
- `:configuration` preserves every node's degree. A small p-value says the
  core-periphery signal goes beyond degree heterogeneity. A few high-degree nodes
  attached to many low-degree nodes already produce a high core-periphery quality
  under `:er`; the configuration null removes that explanation (Kojaku & Masuda,
  2018).
- `:directed_configuration` preserves in- and out-degrees but not reciprocity.
- `:weight_permutation` preserves the topology and the multiset of weights, not node
  strengths; it asks whether the weights rather than the edges carry the structure.

A network with a planted core can be significant under `:er` and not under
`:configuration`. That is the expected outcome when the planted structure is mostly a
degree effect, not a contradiction. Report both p-values and name the null model
behind each one. The [Tutorial](@ref) contains a worked example. The analytic
inference framework of Yanchenko, Sengupta & Mukherjee (2025, arXiv:2508.04730) is not
implemented.

## Surprise-Based Detection

```julia
result = surprise_cp(A; max_iter=100)
```

Uses statistical surprise to find partitions that are unlikely under a random null model (Jeude et al., 2019). Surprise is higher when:

- Core-core and core-periphery edges are more than expected
- Periphery-periphery edges are fewer than expected

**Algorithm**: Greedy single-start node flipping from a mean-degree initial split,
accepting objective-improving moves. The statistic is the paper's exact joint-tail
probability, but this optimiser is a package variant; the paper's search procedure is
not reproduced. The tail sum is evaluated from the mode of the hypergeometric terms
outwards with constant-time ratio updates and truncated once the geometrically bounded
remainder is below double-precision resolution (checked against a `BigInt` oracle for
small graphs), so one evaluation costs roughly the product of the widths of the two
hypergeometric peaks rather than the full ``O(m^2)`` double sum. A greedy sweep is
therefore about ``n`` such evaluations plus the neighbour counts of each node
(``O(n)`` dense, ``O(\text{degree})`` sparse); a 400-node Erdős-Rényi graph at density
0.05 now fits in tens of milliseconds instead of seconds.

**Parameters**:

- `max_iter`: Maximum iterations (default: 100)

**When to use**: When you want a statistically-motivated definition of core-periphery structure.

## Label Switching Algorithm

```julia
result = label_switching_cp(A; max_iter=100, n_runs=1, directed=false, cp_block=:ones)
```

A seedable greedy algorithm that optimizes the discrete Borgatti-Everett
Pearson correlation while excluding degenerate all-core/all-periphery states. It
optimises the same objective as `borgatti_everett_discrete` (including the `directed`
and `cp_block` variants) and is the stronger optimiser of the two because the random
sweep order and restarts escape more local optima.

**Algorithm**:

1. Initialize based on node degrees
2. Process nodes in random order
3. Flip labels that improve the objective
4. Repeat from multiple seeded starts when requested

**Parameters**:

- `max_iter`: Maximum iterations (default: 100)
- `n_runs`: Number of starts (default: 1)
- `rng`: Explicit random-number generator
- `directed`: Accept asymmetric input and correlate over all ordered off-diagonal dyads
- `cp_block`: `:ones` (Borgatti-Everett equation 2) or `:missing` (equation 4)

**When to use**: Fast discrete classification, especially for large networks. Similar to Lip's algorithm but with different optimization strategy.

## Disconnected input and isolates

Core-periphery objectives are defined for a single connected network. On a union of
components the optimiser may assemble a core from several components or label a whole
component as periphery, so the package makes two guarantees and one check:

- **Isolates are never core.** `borgatti_everett_continuous`,
  `borgatti_everett_discrete`, `label_switching_cp`, `lip_discrete`,
  `rombach_continuous`, `rombach_core_score`, `minres_symmetric`, `minres_svd`,
  `minres_svd_directed`, and `multiple_cp_pairs` give every node with zero in- and
  out-strength coreness `0` and keep it out of `core_nodes` (`multiple_cp_pairs`
  reports pair label `0`). Rombach and the multi-pair detector drop isolates before
  fitting, so the fit equals that of the non-isolated subnetwork; the Borgatti-Everett
  correlations keep isolates as peripheral dyads of the sociogram. An edgeless network
  returns all-zero coreness and an empty core.
- **Disconnected input warns.** The same detectors emit a warning, once per detector
  per Julia session (`maxlog=1`), when the input has more than one connected component
  with at least two nodes. Isolates alone do not trigger it. The result is still
  returned; analyse each component separately when the components are substantively
  separate networks. The check is skipped for the null samples drawn inside
  `cp_significance`, which are random graphs and often disconnected; the observed
  network is still checked.
- **Rossa requires connectivity.** `random_walker_profiling` and
  `rossa_stationary_distribution` raise an `ArgumentError` unless undirected input is
  connected or directed input is strongly connected, because the persistence profile is
  undefined otherwise.

`surprise_cp` and `spectral_method` apply neither the isolate rule nor the warning:
Surprise scores the partition of the whole simple graph, and LowRank-Core's Find-Cut
is defined on the full score order. Inspect `coreness` for isolates before reporting
their partitions on fragmented input.

## Comparison Summary

The table lists one row per public estimator; `minres_svd` is the `CPResult`
projection of `minres_svd_directed`, and `multiple_cp_pairs_config` is the
configuration-null specialization of `multiple_cp_pairs`.

| Algorithm | Type | Speed | Directed | Key Feature |
|-----------|------|-------|----------|-------------|
| `borgatti_everett_continuous` | Continuous | Medium | With `directed=true` | Classic baseline |
| `borgatti_everett_discrete` | Discrete | Medium (exhaustive: exponential) | With `directed=true` | Discrete correlation; missing-block variant; exact option |
| `lip_discrete` | Discrete | Fast | No | Exact degree-prefix optimum |
| `rombach_continuous` | Continuous | Slow | No | Tunable template; `beta` sets the core size |
| `rombach_core_score` | Continuous | Slow (grid of fits) | No | Parameter-averaged aggregate core score |
| `spectral_method` | Continuous | Fast | No | No greedy iteration; sparse eigensolve |
| `random_walker_profiling` | Continuous | Medium | **Yes** | Stationary-flow persistence |
| `minres_symmetric` | Continuous | Medium | No | Symmetric rank-one residual |
| `minres_svd` | Continuous | Medium | **Yes** | Handles asymmetric networks |
| `multiple_cp_pairs` | Multi-pair | Medium | No | Multiple structures |
| `surprise_cp` | Discrete | Medium | No | Exact surprise statistic |
| `label_switching_cp` | Discrete | Fast | With `directed=true` | Efficient label updates |

## [References](@id algorithm-references)

1. Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures. *Social Networks*, 21(4), 375-395. [doi:10.1016/S0378-8733(99)00019-2](https://doi.org/10.1016/S0378-8733(99)00019-2)
2. Lip, S.Z.W. (2011). A fast algorithm for the discrete core/periphery bipartitioning problem. [arXiv:1102.5511](https://arxiv.org/abs/1102.5511). [doi:10.48550/arXiv.1102.5511](https://doi.org/10.48550/arXiv.1102.5511)
3. Rombach, M.P., Porter, M.A., Fowler, J.H., Mucha, P.J. (2017). Core-periphery structure in networks (revisited). *SIAM Review*, 59(3), 619-646. [doi:10.1137/17M1130046](https://doi.org/10.1137/17M1130046)
4. Cucuringu, M., Rombach, P., Lee, S.H., Porter, M.A. (2016). Detection of core-periphery structure in networks using spectral methods and geodesic paths. *European Journal of Applied Mathematics*, 27(6), 846-887. [doi:10.1017/S095679251600022X](https://doi.org/10.1017/S095679251600022X)
5. Della Rossa, F., Dercole, F., Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. *Scientific Reports*, 3, 1467. [doi:10.1038/srep01467](https://doi.org/10.1038/srep01467)
6. Boyd, J.P., Fitzgerald, W.J., Mahutga, M.C., Smith, D.A. (2010). Computing continuous core/periphery structures for social relations data with MINRES/SVD. *Social Networks*, 32(2), 125-137. [doi:10.1016/j.socnet.2009.09.003](https://doi.org/10.1016/j.socnet.2009.09.003)
7. Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks. *Physical Review E*, 96(5), 052313. [doi:10.1103/PhysRevE.96.052313](https://doi.org/10.1103/PhysRevE.96.052313)
8. Kojaku, S., Masuda, N. (2018). Core-periphery structure requires something else in the network. *New Journal of Physics*, 20, 043012. [doi:10.1088/1367-2630/aab547](https://doi.org/10.1088/1367-2630/aab547)
9. van Lidth de Jeude, J., Caldarelli, G., Squartini, T. (2019). Detecting core-periphery structures by surprise. *EPL (Europhysics Letters)*, 125(6), 68001. [doi:10.1209/0295-5075/125/68001](https://doi.org/10.1209/0295-5075/125/68001)
10. Yanchenko, E., Sengupta, S. (2025). A label-switching algorithm for fast core-periphery identification. [arXiv:2506.02069](https://arxiv.org/abs/2506.02069). [doi:10.48550/arXiv.2506.02069](https://doi.org/10.48550/arXiv.2506.02069)
11. Yanchenko, E., Sengupta, S., Mukherjee, D. (2025). Statistical inference for core-periphery structures. [arXiv:2508.04730](https://arxiv.org/abs/2508.04730). Related work; not implemented.
12. Milo, R., Kashtan, N., Itzkovitz, S., Newman, M.E.J., Alon, U. (2003). On the uniform generation of random graphs with prescribed degree sequences. [arXiv:cond-mat/0312028](https://arxiv.org/abs/cond-mat/0312028). Null-model reference for the switching chains.
