# CorePeriphery.jl

A Julia package for detecting core-periphery structure in networks.

## Overview

Core-periphery structure is a fundamental mesoscale pattern in complex networks where nodes are partitioned into two groups:

- **Core**: A densely interconnected group of central nodes
- **Periphery**: Sparsely connected nodes that attach to the core but not to each other

This pattern appears in many real-world networks including social networks, economic systems, biological networks, and transportation infrastructure. CorePeriphery.jl provides multiple algorithms for detecting and quantifying this structure.

## Features

- **10 detection algorithms** covering discrete, continuous, and multi-pair methods
- **Weighted network support** for analyzing networks with edge weights
- **Directed network support** via the MINRES/SVD algorithm
- **Pure Julia implementation** with no dependencies beyond the standard library
- **Consistent API** with unified result structures across all algorithms

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/simoneSantoni/CorePeriphery.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/CorePeriphery.jl")
```

## Quick Start

```julia
using CorePeriphery

# Create a network from an edge list
edges = [
    (1, 2), (1, 3), (1, 4), (1, 5),  # Node 1 connected to many
    (2, 3), (2, 4), (2, 5),           # Node 2 also well-connected
    (3, 4), (3, 5),
    (6, 1), (7, 2)                    # Peripheral nodes
]
A = adjacency_to_matrix(edges, 7)

# Run a detection algorithm
result = borgatti_everett_continuous(A)

# Access results
println("Coreness scores: ", result.coreness)
println("Core nodes: ", result.core_nodes)
println("Periphery nodes: ", result.periphery_nodes)
println("Quality: ", result.quality)
```

## Choosing an Algorithm

| Use Case | Recommended Algorithm |
|----------|----------------------|
| Fast binary classification | [`lip_discrete`](@ref) or [`label_switching_cp`](@ref) |
| Continuous coreness scores | [`borgatti_everett_continuous`](@ref) or [`spectral_method`](@ref) |
| Directed networks | [`minres_svd`](@ref) |
| Multiple CP structures | [`multiple_cp_pairs`](@ref) |
| Tunable core boundary | [`rombach_continuous`](@ref) |
| Statistical significance | [`surprise_cp`](@ref) |

## Documentation

```@contents
Pages = ["tutorial.md", "algorithms.md", "api.md"]
Depth = 2
```

## References

1. Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures. *Social Networks*, 21(4), 375-395.
2. Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem. *arXiv:1102.5511*.
3. Rombach, M.P., Porter, M.A., Fowler, J.H., Mucha, P.J. (2017). Core-Periphery Structure in Networks (Revisited). *SIAM Review*, 59(3), 619-646.
4. Cucuringu, M., et al. (2016). Detection of core-periphery structure in networks using spectral methods and geodesic paths. *European Journal of Applied Mathematics*, 27(6), 846-887.
5. Della Rossa, F., Dercole, F., Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. *Scientific Reports*, 3, 1467.
6. Boyd, J.P., et al. (2010). Computing continuous core/periphery structures for social relations data with MINRES/SVD. *Social Networks*, 32(2), 125-137.
7. Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks. *Physical Review E*, 96(5), 052313.
8. Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise. *EPL*, 125(6), 68001.
9. Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection in networks. *arXiv preprint*.
