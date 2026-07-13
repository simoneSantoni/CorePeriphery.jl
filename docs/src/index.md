# CorePeriphery.jl

A Julia package for detecting core-periphery structure in networks.

## Overview

Core-periphery structure is a fundamental mesoscale pattern in complex networks where nodes are partitioned into two groups:

- **Core**: A densely interconnected group of central nodes
- **Periphery**: Sparsely connected nodes that attach to the core but not to each other

This pattern appears in many real-world networks including social networks, economic systems, biological networks, and transportation infrastructure. CorePeriphery.jl provides multiple algorithms for detecting and quantifying this structure.

## Features

- **Published discrete, continuous, spectral, directed, and multi-pair methods**
- **Weighted network support** for analyzing networks with edge weights
- **Directed in/out coreness** via `CPDirectedResult`
- **Optional Graphs.jl integration** and sparse-matrix input
- **Monte Carlo significance tests** with undirected, directed, and weighted nulls
  plus explicit swap-completion diagnostics
- **Matrix-free sparse kernels** for spectral, directed, and greedy methods
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
| Continuous coreness scores | [`borgatti_everett_continuous`](@ref) or [`rombach_continuous`](@ref) |
| Published spectral partition | [`spectral_method`](@ref) |
| Directed networks | [`random_walker_profiling`](@ref) or [`minres_svd_directed`](@ref) |
| Symmetric MINRES | [`minres_symmetric`](@ref) |
| Multiple CP structures | [`multiple_cp_pairs`](@ref) |
| Tunable core boundary | [`rombach_continuous`](@ref) |
| Statistical significance | [`cp_significance`](@ref) |
| Configuration null convenience | [`multiple_cp_pairs_config`](@ref) |

## Documentation

```@contents
Pages = ["tutorial.md", "algorithms.md", "contracts.md", "performance.md", "development.md", "migration.md", "api.md"]
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

## Citation

```biblatex
@misc{SantoniCorePeripheryJL,
  author = {Santoni, Simone},
  title = {CorePeriphery.jl: Core-Periphery Detection in Julia},
  year = {2026},
  url = {https://github.com/simoneSantoni/CorePeriphery.jl},
  note = {Homepage: https://www.bayes.citystgeorges.ac.uk/faculties-and-research/experts/simone-santoni; GitHub: https://github.com/simoneSantoni}
}
```
