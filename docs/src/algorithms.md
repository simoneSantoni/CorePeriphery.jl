# Algorithms

CorePeriphery.jl implements 10 algorithms for detecting core-periphery structure. This page provides detailed descriptions of each algorithm, their theoretical foundations, and guidance on when to use them.

## Algorithm Categories

The algorithms can be grouped into three categories:

1. **Continuous methods**: Assign each node a coreness score in [0, 1]
2. **Discrete methods**: Classify nodes as either core (1) or periphery (0)
3. **Multi-pair methods**: Detect multiple non-overlapping core-periphery pairs

## Borgatti-Everett Continuous Model

```julia
result = borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, step=0.05, init=nothing)
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

**When to use**: When you need continuous coreness scores and want a well-established baseline method.

## Borgatti-Everett Discrete Model

```julia
result = borgatti_everett_discrete(A; max_iter=1000, init=nothing)
```

The discrete version of the Borgatti-Everett model. Finds a binary partition maximizing correlation with the ideal discrete pattern:

```math
\Delta_{ij} = \max(c_i, c_j)
```

where ``c_i \in \{0, 1\}``.

**Algorithm**: Greedy optimization that iteratively flips node labels to improve quality.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)
- `init`: Initial binary assignment (default: degree-based)

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

**When to use**: When you need fast discrete classification for large networks. Often faster than Borgatti-Everett discrete while achieving similar results.

## Rombach's Generalized Model

```julia
result = rombach_continuous(A; alpha=0.5, beta=0.5, n_runs=10, rng=Random.default_rng())
```

A generalized continuous model by Rombach et al. (2017) with parameters controlling the core-periphery structure:

- ``\alpha``: Controls the sharpness of the core-periphery boundary
- ``\beta``: Controls the relative size of the core

The parameters define a fixed rank template. The algorithm assigns a
permutation of that template to nodes.

**Algorithm**: Seedable random restarts with objective-improving pair swaps,
maximizing ``\sum_{i<j} A_{ij}c_ic_j``.

**Parameters**:

- `alpha`: Boundary sharpness, 0 to 1 (default: 0.5)
- `beta`: Core size parameter, 0 to 1 (default: 0.5)
- `n_runs`: Number of random restarts (default: 10)

**When to use**: When you want to tune the expected core-periphery structure or explore how results change with different boundary definitions.

## Spectral Method

```julia
result = spectral_method(A)
```

Implements Cucuringu et al.'s LowRank-Core algorithm: reconstruct the adjacency
matrix from its two eigenpairs of largest magnitude, threshold the rank-two
matrix, rank nodes by its row sums, and apply Find-Cut.

**Algorithm**: Rank-two symmetric eigendecomposition plus the published
core/core-periphery/periphery density cut objective. Input is binary,
undirected, and loop-free.

**Parameters**: `beta` sets the minimum fraction allowed in either side of the cut.

**When to use**: Fast method for continuous coreness scores. Works well when the core-periphery structure aligns with the network's spectral properties. No iterative optimization needed.

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

**When to use**: For directed (asymmetric) outer-product structure. The current
`CPResult` reports the average of in- and out-coreness; directed quality uses
all ordered off-diagonal dyads.

For a symmetric rank-one model, use:

```julia
result = minres_symmetric(A; max_iter=10_000, tol=1e-8)
```

This fits ``A_{ij}\approx w_iw_j`` and is a different estimand from the directed
``u_iv_j`` model.

## Multiple Core-Periphery Pairs

```julia
result = multiple_cp_pairs(A; max_pairs=10, min_pair_size=2, max_iter=100)
```

Detects multiple non-overlapping core-periphery pairs using the ``Q^{cp}`` quality function (Kojaku & Masuda, 2017).

**Algorithm**:

1. Create balanced seeded initial pair assignments
2. Jointly propose each node as core or periphery in its own or a neighboring pair
3. Accept only objective-improving moves across multiple restarts

Returns a `CPMultiResult` with pair labels and within-pair coreness.

**Parameters**:

- `max_pairs`: Maximum number of pairs (default: 10)
- `min_pair_size`: Minimum nodes per pair (default: 2)
- `max_iter`: Optimization iterations (default: 100)
- `n_runs`: Seeded restarts (default: 10)
- `null_model`: `:er` or `:configuration`
- `rng`: Explicit random-number generator
- `pair_selection`: `:objective` for compatibility or `:penalized` for complexity control
- `pair_penalty`: Multiplier for ``(k-1)\log(n)/n`` (default: 0.5)

**When to use**: When you expect the network to contain multiple distinct core-periphery structures rather than a single global one.

## Surprise-Based Detection

```julia
result = surprise_cp(A; max_iter=100)
```

Uses statistical surprise to find partitions that are unlikely under a random null model (Jeude et al., 2019). Surprise is higher when:

- Core-core and core-periphery edges are more than expected
- Periphery-periphery edges are fewer than expected

**Algorithm**: Greedy optimization of the surprise score.

**Parameters**:

- `max_iter`: Maximum iterations (default: 100)

**When to use**: When you want a statistically-motivated definition of core-periphery structure.

## Label Switching Algorithm

```julia
result = label_switching_cp(A; max_iter=100)
```

A seedable greedy algorithm that optimizes the discrete Borgatti-Everett
Pearson correlation while excluding degenerate all-core/all-periphery states.

**Algorithm**:

1. Initialize based on node degrees
2. Process nodes in random order
3. Flip labels that improve the objective
4. Repeat from multiple seeded starts when requested

**Parameters**:

- `max_iter`: Maximum iterations (default: 100)
- `n_runs`: Number of starts (default: 1)
- `rng`: Explicit random-number generator

**When to use**: Fast discrete classification, especially for large networks. Similar to Lip's algorithm but with different optimization strategy.

## Comparison Summary

| Algorithm | Type | Speed | Directed | Key Feature |
|-----------|------|-------|----------|-------------|
| `borgatti_everett_continuous` | Continuous | Medium | No | Classic baseline |
| `borgatti_everett_discrete` | Discrete | Medium | No | Optimal discrete correlation |
| `lip_discrete` | Discrete | Fast | No | Exact degree-prefix optimum |
| `rombach_continuous` | Continuous | Slow | No | Tunable parameters |
| `spectral_method` | Continuous | Fast | No | No iteration needed |
| `random_walker_profiling` | Continuous | Medium | **Yes** | Stationary-flow persistence |
| `minres_symmetric` | Continuous | Medium | No | Symmetric rank-one residual |
| `minres_svd` | Continuous | Medium | **Yes** | Handles asymmetric networks |
| `multiple_cp_pairs` | Multi-pair | Medium | No | Multiple structures |
| `surprise_cp` | Discrete | Medium | No | Statistical foundation |
| `label_switching_cp` | Discrete | Fast | No | Efficient label updates |
