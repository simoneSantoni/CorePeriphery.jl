# CorePeriphery.jl

A Julia package implementing various algorithms for detecting core-periphery structure in networks.

## Overview

Core-periphery structure is a fundamental mesoscale pattern in complex networks where nodes are partitioned into:
- **Core**: A densely interconnected group of central nodes
- **Periphery**: Sparsely connected nodes that attach to the core but not to each other

This package provides multiple algorithms for detecting and quantifying core-periphery structure.

## Installation

```julia
# Add the src directory to your load path
push!(LOAD_PATH, "path/to/cp/src")
using CorePeriphery
```

## Algorithms Implemented

### 1. Borgatti-Everett Continuous Model
```julia
result = borgatti_everett_continuous(A)
```
The classic model that finds coreness values maximizing correlation with ideal pattern Δ[i,j] = c[i] × c[j].

### 2. Borgatti-Everett Discrete Model
```julia
result = borgatti_everett_discrete(A)
```
Binary classification optimizing correlation with ideal discrete pattern Δ[i,j] = max(c[i], c[j]).

### 3. Lip's Fast Discrete Algorithm
```julia
result = lip_discrete(A)
```
Efficient swap-based optimization for discrete core-periphery bipartitioning.

### 4. Rombach's Generalized Model
```julia
result = rombach_continuous(A; alpha=0.5, beta=0.5)
```
Continuous model with parameters controlling boundary sharpness (α) and core size (β).

### 5. Spectral Method
```julia
result = spectral_method(A)
```
Uses the leading eigenvector of the adjacency matrix to determine coreness.

### 6. Random Walker Profiling
```julia
result = random_walker_profiling(A; n_walks=1000, walk_length=10)
```
Nodes visited more frequently by random walks are classified as more core-like.

## Usage

### Basic Example

```julia
using CorePeriphery

# Create adjacency matrix from edge list
edges = [
    (1, 2), (1, 3), (1, 4), (1, 5),
    (2, 3), (2, 4), (2, 5),
    (3, 4), (3, 5),
    (6, 1), (7, 2)  # Periphery nodes
]
A = adjacency_to_matrix(edges, 7)

# Run detection algorithm
result = borgatti_everett_continuous(A)

# Access results
println("Coreness scores: ", result.coreness)
println("Core nodes: ", result.core_nodes)
println("Periphery nodes: ", result.periphery_nodes)
println("Quality: ", result.quality)
```

### Result Structure

All algorithms return a `CPResult` with:
- `coreness::Vector{Float64}`: Coreness score for each node (0 to 1)
- `core_nodes::Vector{Int}`: Indices of nodes classified as core
- `periphery_nodes::Vector{Int}`: Indices of nodes classified as periphery
- `quality::Float64`: Quality score (correlation with ideal pattern)
- `algorithm::String`: Name of the algorithm used

### Comparing Algorithms

```julia
# Run multiple algorithms on the same network
results = [
    borgatti_everett_continuous(A),
    borgatti_everett_discrete(A),
    lip_discrete(A),
    spectral_method(A),
    random_walker_profiling(A)
]

for r in results
    println("$(r.algorithm): quality = $(round(r.quality, digits=3))")
end
```

### Weighted Networks

The algorithms support weighted adjacency matrices:

```julia
weighted_edges = [
    (1, 2, 5.0),  # Strong connection
    (1, 3, 2.0),  # Weaker connection
    (2, 3, 3.0)
]
A = adjacency_to_matrix(weighted_edges, 3)
result = borgatti_everett_continuous(A)
```

## Utility Functions

```julia
# Convert edge list to adjacency matrix
A = adjacency_to_matrix(edges, n)

# Compute quality of a given coreness assignment
q = core_quality(A, coreness_vector)

# Generate ideal core-periphery matrix
Δ = ideal_cp_matrix(coreness_vector)

# Get coreness scores from result
scores = coreness_scores(result)
```

## Running Tests

```julia
include("test/test_coreperiphery.jl")
```

## Running Examples

```julia
include("examples/basic_usage.jl")
```

## References

1. Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures. *Social Networks*, 21(4), 375-395.

2. Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem. *arXiv:1102.5511*.

3. Rombach, M.P., Porter, M.A., Fowler, J.H., Mucha, P.J. (2017). Core-Periphery Structure in Networks (Revisited). *SIAM Review*, 59(3), 619-646.

4. Cucuringu, M., Rombach, P., Lee, S.H., Porter, M.A. (2016). Detection of core-periphery structure in networks using spectral methods and geodesic paths. *European Journal of Applied Mathematics*, 27(6), 846-887.

5. Della Rossa, F., Dercole, F., Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. *Scientific Reports*, 3, 1467.

6. Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks. *Physical Review E*, 96(5), 052313.

7. Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise. *EPL*, 125(6), 68001.

## License

MIT License - see LICENSE file for details.
