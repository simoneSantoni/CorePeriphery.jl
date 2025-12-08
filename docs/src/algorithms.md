# Algorithms

CorePeriphery.jl implements 10 algorithms for detecting core-periphery structure. This page provides detailed descriptions of each algorithm, their theoretical foundations, and guidance on when to use them.

## Algorithm Categories

The algorithms can be grouped into three categories:

1. **Continuous methods**: Assign each node a coreness score in [0, 1]
2. **Discrete methods**: Classify nodes as either core (1) or periphery (0)
3. **Multi-pair methods**: Detect multiple non-overlapping core-periphery pairs

## Borgatti-Everett Continuous Model

```julia
result = borgatti_everett_continuous(A; max_iter=1000, tol=1e-6, init=nothing)
```

The classic continuous core-periphery model introduced by Borgatti and Everett (2000). It finds a coreness vector ``c`` that maximizes the correlation between the adjacency matrix ``A`` and an ideal core-periphery pattern:

```math
\Delta_{ij} = c_i \times c_j
```

**Algorithm**: Coordinate ascent optimization. Each node's coreness is updated iteratively while holding others fixed.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)
- `tol`: Convergence tolerance (default: 1e-6)
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
result = lip_discrete(A; max_iter=1000)
```

An efficient algorithm for discrete core-periphery detection introduced by Lip (2011). Optimizes a simpler objective: maximizing the number of edges involving at least one core node.

**Algorithm**: Swap-based optimization with O(n) updates per iteration. Maintains running statistics for efficient delta computations.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)

**When to use**: When you need fast discrete classification for large networks. Often faster than Borgatti-Everett discrete while achieving similar results.

## Rombach's Generalized Model

```julia
result = rombach_continuous(A; alpha=0.5, beta=0.5, max_iter=1000, tol=1e-6, n_runs=10)
```

A generalized continuous model by Rombach et al. (2017) with parameters controlling the core-periphery structure:

- ``\alpha``: Controls the sharpness of the core-periphery boundary
- ``\beta``: Controls the relative size of the core

The ideal pattern uses a transition function ``g(c_i; \alpha, \beta)`` to compute ``\Delta_{ij} = g(c_i) \times g(c_j)``.

**Algorithm**: Grid search optimization with random restarts.

**Parameters**:

- `alpha`: Boundary sharpness, 0 to 1 (default: 0.5)
- `beta`: Core size parameter, 0 to 1 (default: 0.5)
- `n_runs`: Number of random restarts (default: 10)

**When to use**: When you want to tune the expected core-periphery structure or explore how results change with different boundary definitions.

## Spectral Method

```julia
result = spectral_method(A)
```

Uses the eigenvector corresponding to the largest eigenvalue of the adjacency matrix to determine coreness (Cucuringu et al., 2016).

**Algorithm**: Eigendecomposition of the adjacency matrix. The leading eigenvector components (absolute values, normalized) give coreness scores.

**Parameters**: None

**When to use**: Fast method for continuous coreness scores. Works well when the core-periphery structure aligns with the network's spectral properties. No iterative optimization needed.

## Random Walker Profiling

```julia
result = random_walker_profiling(A; n_walks=1000, walk_length=10)
```

Based on the insight that random walkers spend more time at core nodes (Della Rossa et al., 2013). Coreness is determined by visit frequency.

**Algorithm**: Simulate random walks starting from random nodes. Count visits to each node and normalize to get coreness scores.

**Parameters**:

- `n_walks`: Number of random walks (default: 1000)
- `walk_length`: Steps per walk (default: 10)

**When to use**: When you want a dynamics-based definition of coreness. Useful for understanding network flow and accessibility.

## MINRES/SVD Method

```julia
result = minres_svd(A; max_iter=1000, tol=1e-6)
```

Minimizes the residual between the adjacency matrix and an outer product model (Boyd et al., 2010):

```math
\min_{u,v} \sum_{i \neq j} (A_{ij} - u_i v_j)^2
```

This yields separate in-coreness (``v``) and out-coreness (``u``) vectors, which are averaged for the final coreness scores.

**Algorithm**: Alternating least squares optimization.

**Parameters**:

- `max_iter`: Maximum iterations (default: 1000)
- `tol`: Convergence tolerance (default: 1e-6)

**When to use**: **The only algorithm that handles directed (asymmetric) networks**. For undirected networks, in-coreness equals out-coreness.

## Multiple Core-Periphery Pairs

```julia
result = multiple_cp_pairs(A; max_pairs=10, min_pair_size=2, max_iter=100)
```

Detects multiple non-overlapping core-periphery pairs using the ``Q^{cp}`` quality function (Kojaku & Masuda, 2017).

**Algorithm**:

1. Initialize each node as its own pair
2. Greedily merge pairs to maximize ``Q^{cp}``
3. Optimize core/periphery assignment within each pair

Returns a `CPMultiResult` with pair labels and within-pair coreness.

**Parameters**:

- `max_pairs`: Maximum number of pairs (default: 10)
- `min_pair_size`: Minimum nodes per pair (default: 2)
- `max_iter`: Optimization iterations (default: 100)

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

A fast greedy algorithm with efficient O(n) updates per iteration (Yanchenko & Sengupta, 2025). Optimizes the sum ``M = \sum_{ij} A_{ij} \Delta_{ij}`` with the discrete ideal pattern.

**Algorithm**:

1. Initialize based on node degrees
2. Process nodes in random order
3. Flip labels that improve the objective
4. Maintain running neighbor statistics for O(1) delta computations

**Parameters**:

- `max_iter`: Maximum iterations (default: 100)

**When to use**: Fast discrete classification, especially for large networks. Similar to Lip's algorithm but with different optimization strategy.

## Comparison Summary

| Algorithm | Type | Speed | Directed | Key Feature |
|-----------|------|-------|----------|-------------|
| `borgatti_everett_continuous` | Continuous | Medium | No | Classic baseline |
| `borgatti_everett_discrete` | Discrete | Medium | No | Optimal discrete correlation |
| `lip_discrete` | Discrete | Fast | No | Efficient swap optimization |
| `rombach_continuous` | Continuous | Slow | No | Tunable parameters |
| `spectral_method` | Continuous | Fast | No | No iteration needed |
| `random_walker_profiling` | Continuous | Medium | No | Dynamics-based |
| `minres_svd` | Continuous | Medium | **Yes** | Handles asymmetric networks |
| `multiple_cp_pairs` | Multi-pair | Medium | No | Multiple structures |
| `surprise_cp` | Discrete | Medium | No | Statistical foundation |
| `label_switching_cp` | Discrete | Fast | No | Efficient label updates |
