# CorePeriphery.jl

A Julia package for detecting core-periphery structure in networks.

## Overview

Core-periphery structure is a fundamental mesoscale pattern in complex networks where nodes are partitioned into two groups:

- **Core**: A densely interconnected group of central nodes
- **Periphery**: Sparsely connected nodes that attach to the core but not to each other

This pattern appears in many real-world networks including social networks, economic systems, biological networks, and transportation infrastructure. CorePeriphery.jl provides multiple algorithms for detecting and quantifying this structure.

## Features

- **Published discrete, continuous, spectral, directed, and multi-pair methods**,
  including the Borgatti-Everett missing-block model, an exact discrete optimum for
  small sociograms, and Rombach's aggregate core score
- **Weighted network support** for analyzing networks with edge weights
- **Directed in/out coreness** via `CPDirectedResult`, and `directed=true` in the
  Borgatti-Everett family
- **Labelled edge lists** (`adjacency_to_matrix(edges)` returns the matrix and node
  labels) and `block_densities` for reporting a fit
- **Optional Graphs.jl integration** and sparse-matrix input
- **Monte Carlo significance tests** with undirected, directed, and weighted nulls,
  effect sizes (`zscore`), sparse-aware sampling, and explicit swap-completion
  diagnostics
- **Matrix-free sparse kernels** for spectral, directed, and greedy methods
- **Consistent API** with unified result structures across all algorithms, explicit
  `partition` rules for continuous scores, and an isolate guarantee (isolates are
  never core)

## Installation

CorePeriphery.jl is not yet in the General registry. Install it from GitHub:

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

```jldoctest index-quickstart
julia> using CorePeriphery

julia> edges = [
           (1, 2), (1, 3), (1, 4), (1, 5),  # Node 1 connected to many
           (2, 3), (2, 4), (2, 5),           # Node 2 also well-connected
           (3, 4), (3, 5),
           (6, 1), (7, 2)                    # Peripheral nodes
       ];

julia> A = adjacency_to_matrix(edges, 7);

julia> result = borgatti_everett_continuous(A)
CPResult: Borgatti-Everett Continuous
  nodes = 7, core = 3, periphery = 4
  quality = 0.7818

julia> println("Coreness scores: ", result.coreness)
Coreness scores: [0.7, 0.7, 0.65, 0.5, 0.5, 0.2, 0.2]

julia> println("Core nodes: ", result.core_nodes)
Core nodes: [1, 2, 3]

julia> println("Periphery nodes: ", result.periphery_nodes)
Periphery nodes: [4, 5, 6, 7]

julia> println("Quality: ", result.quality)
Quality: 0.78178482001133

julia> block_densities(A, result)   # report the block densities with the fit
(core_core = 1.0, core_periphery = 0.6666666666666666, periphery_core = 0.6666666666666666, periphery_periphery = 0.0)
```

Node identifiers need not be integers: `A, labels = adjacency_to_matrix(edges)` accepts
any identifiers and returns the label of every row, so `labels[result.core_nodes]`
reports the core in the original names (see the [Tutorial](@ref)).

## Choosing an Algorithm

| Use Case | Recommended Algorithm |
|----------|----------------------|
| Fast binary classification | [`lip_discrete`](@ref) or [`label_switching_cp`](@ref) |
| Exact discrete optimum on a small sociogram | [`borgatti_everett_discrete`](@ref) with `method=:exhaustive` |
| Continuous coreness scores | [`borgatti_everett_continuous`](@ref) or [`rombach_continuous`](@ref) |
| Parameter-averaged Rombach core score | [`rombach_core_score`](@ref) |
| Published spectral partition | [`spectral_method`](@ref) |
| Directed networks | [`random_walker_profiling`](@ref), [`minres_svd_directed`](@ref), or the Borgatti-Everett family with `directed=true` |
| Symmetric MINRES | [`minres_symmetric`](@ref) |
| Multiple CP structures | [`multiple_cp_pairs`](@ref) |
| Tunable core size | [`rombach_continuous`](@ref) (`beta` is the periphery fraction) |
| Statistical significance and effect sizes | [`cp_significance`](@ref) |
| Configuration null convenience | [`multiple_cp_pairs_config`](@ref) |
| Block-density reporting | [`block_densities`](@ref) |

Raw `quality` values are algorithm-specific; do not use them to rank algorithms against each other.

## Performance snapshots

### CorePeriphery.jl compared with cpnet

Across three identical planted fixtures and all nine comparable algorithm
configurations, CorePeriphery.jl was faster in 25 of 27 warmed fits. The table reports
the median cpnet/Julia runtime ratio across the three fixtures; a value above 1 means
that the Julia call completed sooner. Fit budgets are not matched: cpnet runs five
stochastic starts for BE, MINRES, and Surprise where the Julia call is a single
deterministic fit, so the largest ratios compare configured estimators rather than
equal optimizer work.

| Algorithm | Median cpnet/Julia ratio | Julia faster |
|:----------|--------------------------:|-------------:|
| Borgatti--Everett discrete (BE, [`borgatti_everett_discrete`](@ref)) | 141.8× | 3/3 |
| Lip ([`lip_discrete`](@ref)) | 10.3× | 3/3 |
| LowRank-Core ([`spectral_method`](@ref)) | 20.0× | 3/3 |
| Rombach ([`rombach_continuous`](@ref)) | 0.7× | 1/3 |
| Rossa ([`random_walker_profiling`](@ref)) | 222.6× | 3/3 |
| symmetric MINRES ([`minres_symmetric`](@ref)) | 196.8× | 3/3 |
| Surprise ([`surprise_cp`](@ref)) | 152.9× | 3/3 |
| KM-ER ([`multiple_cp_pairs`](@ref), `null_model=:er`) | 23.6× | 3/3 |
| KM-configuration ([`multiple_cp_pairs_config`](@ref)) | 21.9× | 3/3 |

These are configured-estimator timings on small `n = 20, 40, 48` fixtures, not
equal-work or large-network speedups. Several cpnet configurations use more stochastic
starts than their deterministic Julia counterparts. See [Performance](@ref) for every
timing pair, fit budgets, fixture densities, methodology, and interpretation limits.

### Julia dense and sparse scaling

CorePeriphery.jl provides matrix-free sparse paths for its principal spectral,
directed, and greedy algorithms. On the package's synthetic `n = 512` scaling fixture,
the following warmed median runtimes were measured with Julia 1.12, one BLAS thread,
and seven repetitions:

| Algorithm | Dense | Sparse | Sparse speedup |
|:----------|------:|-------:|---------------:|
| [`lip_discrete`](@ref) | 0.906 ms | 0.132 ms | 6.8× |
| [`minres_symmetric`](@ref) (`max_iter = 50`) | 1.565 ms | 0.337 ms | 4.6× |
| [`spectral_method`](@ref) | 39.754 ms | 4.215 ms | 9.4× |
| [`random_walker_profiling`](@ref) | 18.480 ms | 0.832 ms | 22.2× |
| [`minres_svd_directed`](@ref) (`max_iter = 50`) | 7.096 ms | 0.675 ms | 10.5× |

Absolute timings are host-specific. The scaling results show the effect of sparse
storage on the benchmark fixture and are not intended to predict runtime for every
network.

## Documentation

```@contents
Pages = ["tutorial.md", "algorithms.md", "contracts.md", "performance.md", "development.md", "migration.md", "api.md"]
Depth = 2
```

## References

Full citations with DOIs for every implemented method are listed on the
[Algorithms](@ref) page.

1. Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures. *Social Networks*, 21(4), 375-395. [doi:10.1016/S0378-8733(99)00019-2](https://doi.org/10.1016/S0378-8733(99)00019-2)
2. Lip, S.Z.W. (2011). A fast algorithm for the discrete core/periphery bipartitioning problem. [arXiv:1102.5511](https://arxiv.org/abs/1102.5511).
3. Rombach, M.P., Porter, M.A., Fowler, J.H., Mucha, P.J. (2017). Core-periphery structure in networks (revisited). *SIAM Review*, 59(3), 619-646. [doi:10.1137/17M1130046](https://doi.org/10.1137/17M1130046)
4. Cucuringu, M., Rombach, P., Lee, S.H., Porter, M.A. (2016). Detection of core-periphery structure in networks using spectral methods and geodesic paths. *European Journal of Applied Mathematics*, 27(6), 846-887. [doi:10.1017/S095679251600022X](https://doi.org/10.1017/S095679251600022X)
5. Della Rossa, F., Dercole, F., Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. *Scientific Reports*, 3, 1467. [doi:10.1038/srep01467](https://doi.org/10.1038/srep01467)
6. Boyd, J.P., Fitzgerald, W.J., Mahutga, M.C., Smith, D.A. (2010). Computing continuous core/periphery structures for social relations data with MINRES/SVD. *Social Networks*, 32(2), 125-137. [doi:10.1016/j.socnet.2009.09.003](https://doi.org/10.1016/j.socnet.2009.09.003)
7. Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks. *Physical Review E*, 96(5), 052313. [doi:10.1103/PhysRevE.96.052313](https://doi.org/10.1103/PhysRevE.96.052313)
8. Kojaku, S., Masuda, N. (2018). Core-periphery structure requires something else in the network. *New Journal of Physics*, 20, 043012. [doi:10.1088/1367-2630/aab547](https://doi.org/10.1088/1367-2630/aab547)
9. van Lidth de Jeude, J., Caldarelli, G., Squartini, T. (2019). Detecting core-periphery structures by surprise. *EPL (Europhysics Letters)*, 125(6), 68001. [doi:10.1209/0295-5075/125/68001](https://doi.org/10.1209/0295-5075/125/68001)
10. Yanchenko, E., Sengupta, S. (2025). A label-switching algorithm for fast core-periphery identification. [arXiv:2506.02069](https://arxiv.org/abs/2506.02069).

## Citation

```biblatex
@software{SantoniCorePeripheryJL,
  author = {Santoni, Simone},
  title = {CorePeriphery.jl: Core-Periphery Detection in Julia},
  version = {0.2.0},
  year = {2026},
  url = {https://github.com/simoneSantoni/CorePeriphery.jl},
  note = {Julia package, version 0.2.0. Source code and documentation at the URL above.}
}
```

A machine-readable `CITATION.cff` is included in the repository.
