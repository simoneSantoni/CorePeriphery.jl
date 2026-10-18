# CorePeriphery.jl

[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/simoneSantoni/CorePeriphery.jl)
[![Build Status](https://github.com/simoneSantoni/CorePeriphery.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/simoneSantoni/CorePeriphery.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/simoneSantoni/CorePeriphery.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/simoneSantoni/CorePeriphery.jl)
<!-- Re-enable the stable-docs badge after the first tagged release:
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://simoneSantoni.github.io/CorePeriphery.jl/stable/)
-->
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://simoneSantoni.github.io/CorePeriphery.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.10+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia package implementing published algorithms for detecting core-periphery structure in networks.

## Overview

Core-periphery structure is a fundamental mesoscale pattern in complex networks where nodes are partitioned into:
- **Core**: A densely interconnected group of central nodes
- **Periphery**: Sparsely connected nodes that attach to the core but not to each other

CorePeriphery.jl implements ten published core-periphery detection methods (the MINRES
method has distinct symmetric and directed estimators) plus Monte Carlo significance
testing. Every detector has a documented contract (input requirements, meaning of its
scores, meaning of its quality value) and an oracle test suite.

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

Julia 1.10 or newer is required.

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

# Output:
# Coreness scores: [0.7, 0.7, 0.65, 0.5, 0.5, 0.2, 0.2]
# Core nodes: [1, 2, 3]
# Periphery nodes: [4, 5, 6, 7]
# Quality: 0.78178482001133

# Report the UCINET-style block densities alongside the fit
block_densities(A, result)
# (core_core = 1.0, core_periphery = 0.6666666666666666, periphery_core = 0.6666666666666666, periphery_periphery = 0.0)
```

Node identifiers need not be integers: `A, labels = adjacency_to_matrix(edges)` accepts
any identifiers (strings, symbols, ...) and returns the label of every row so that
`labels[result.core_nodes]` reports the core in the original names.

## Choosing an Algorithm

| Use Case | Recommended Algorithm |
|----------|----------------------|
| Fast binary classification | `lip_discrete` or `label_switching_cp` |
| Exact discrete optimum on a small sociogram | `borgatti_everett_discrete(A; method=:exhaustive)` |
| Continuous coreness scores | `borgatti_everett_continuous` or `rombach_continuous` |
| Parameter-averaged Rombach core score | `rombach_core_score` |
| Published spectral partition | `spectral_method` |
| Directed networks | `random_walker_profiling`, `minres_svd_directed`, or the Borgatti-Everett family with `directed=true` |
| Symmetric MINRES | `minres_symmetric` |
| Multiple CP structures | `multiple_cp_pairs` |
| Tunable core size | `rombach_continuous` (`beta` is the periphery fraction) |
| Statistical significance and effect sizes | `cp_significance` |
| Configuration null convenience | `multiple_cp_pairs_config` |
| Block-density reporting | `block_densities` |

Raw `quality` values are algorithm-specific; do not use them to rank algorithms against each other.

## Algorithms Implemented

### 1. Borgatti-Everett Continuous Model
```julia
result = borgatti_everett_continuous(A)
```
The classic model that finds coreness values maximizing correlation with ideal pattern Δ[i,j] = c[i] × c[j].

### 2. Borgatti-Everett Discrete Model
```julia
result = borgatti_everett_discrete(A)
result = borgatti_everett_discrete(A; cp_block=:missing, method=:exhaustive)
```
Binary classification optimizing correlation with the ideal discrete pattern
Δ[i,j] = max(c[i], c[j]) (`cp_block=:ones`, BE 2000 eq. 2) or with the missing-block
model that Borgatti and Everett recommend (`cp_block=:missing`, eq. 4, UCINET's
default). `n_runs` adds random restarts, `directed=true` accepts asymmetric input, and
`method=:exhaustive` (n ≤ 24) returns a global optimum.

### 3. Lip's Fast Discrete Algorithm
```julia
result = lip_discrete(A)
```
Lip's exact degree-prefix minimizer for simple, undirected binary graphs. The
returned `quality` is the discrete Pearson correlation of the returned partition;
the Lip loss is minimized internally but not returned.

### 4. Rombach's Generalized Model
```julia
result = rombach_continuous(A; alpha=0.5, beta=0.5)
score = rombach_core_score(A)
```
Continuous model whose eq. (4.9) template is controlled by `alpha` (boundary
sharpness) and `beta`, the fraction of nodes in the periphery (paper semantics): the
core holds the `N - ⌊βN⌋` top template ranks, which the default `partition=:beta`
reports as `core_nodes`. `coreness` sums to one as in the paper and `quality` is the
eq. (4.10) core quality of one `(alpha, beta)` point, so it must not be used to select
the parameters. `rombach_core_score` returns the paper's aggregate core score (eq.
4.11) averaged over an `(alpha, beta)` grid.

### 5. Spectral Method (LowRank-Core)
```julia
result = spectral_method(A)
```
Cucuringu et al.'s LowRank-Core rank-two reconstruction and Find-Cut
partitioning algorithm. The row-score coreness values come from the thresholded
rank-two reconstruction, while the partition comes from the best Find-Cut prefix.
`beta` is the paper's lower bound on the size of either side of the cut, not a core
fraction.

### 6. Random Walker Profiling
```julia
result = random_walker_profiling(A)
```
Della Rossa persistence profiling for connected undirected or strongly connected
directed, nonnegative networks. The reproducible default resolves final ties by node
index; use `tie_break=:paper` with an explicit `rng` for the paper's stochastic rule.
Directed inputs use the stationary distribution of their row-normalized random walk.
Use `rossa_profile_ensemble` for uncertainty across paper-mode tie realizations.

### 7. MINRES/SVD Method
```julia
result = minres_svd(A)
```
Minimizes residual to find in-coreness and out-coreness vectors. Works with asymmetric (directed) networks.
Use `minres_svd_directed` to retain both vectors and the residual; `minres_svd`
returns the backward-compatible combined `CPResult`.
Use `minres_symmetric` for the distinct symmetric rank-one model.

### 8. Multiple Core-Periphery Pairs
```julia
result = multiple_cp_pairs(A; max_pairs=10, min_pair_size=2)
```
Jointly switches pair and core/periphery labels using the Q^cp objective.
Erdős-Rényi (`null_model=:er`, Kojaku-Masuda 2017 eq. 3) and undirected degree
configuration (`:configuration`, Kojaku-Masuda 2018 eq. 13, implemented exactly)
nulls are available; these are the only two null models `multiple_cp_pairs`
accepts, and the papers recommend the configuration null. The number of pairs is
chosen with the penalised rule (`pair_selection=:penalized`, the default; the raw
objective over-partitions). Use `multiple_cp_pairs_config(A)` as a convenience
wrapper for the configuration null. The pair-count selection rule is a package
variant, not the procedure of the Kojaku-Masuda papers. Isolates get pair label 0.

### 9. Surprise-Based Detection
```julia
result = surprise_cp(A)
```
Optimizes the multivariate-hypergeometric joint-tail surprise in log space for
simple, undirected binary graphs. The statistic is exact (the tail is summed from its
mode outwards with rigorous truncation); the greedy single-start optimizer is a
package variant.

### 10. Label-Switching Algorithm
```julia
result = label_switching_cp(A)
```
Seedable, multi-start greedy optimization of the discrete Pearson objective (the same
objective as `borgatti_everett_discrete`, with the same `directed` and `cp_block`
options).

### 11. Statistical Significance
```julia
result = cp_significance(A, lip_discrete; null_model=:configuration, n_samples=999)
```
Compares one detector's observed quality against a Monte Carlo null
distribution (999 samples by default) and reports the plus-one p-value with the
effect sizes `null_mean`, `null_std`, and `zscore`. Supported nulls are `:er`,
`:configuration`, `:directed_configuration`, and `:weight_permutation`. The wrapper
can run threaded, forwards a child RNG to stochastic detectors by default
(`pass_rng=true`), samples sparse input in sparse form, and records swap
completion diagnostics for configuration-style nulls. The null model decides what
"significant" means: `:er` tests "more core-periphery-like than density alone",
`:configuration` tests "beyond degree heterogeneity"; report both.

## Usage

### Result Structure

The single-structure detectors return a `CPResult` with:
- `coreness::Vector{Float64}`: Coreness score for each node (0 to 1; isolates score 0)
- `core_nodes::Vector{Int}`: Indices of nodes classified as core (never isolates)
- `periphery_nodes::Vector{Int}`: Indices of nodes classified as periphery
- `quality::Float64`: Algorithm-specific quality score. Do not compare raw
  values across methods with different objectives.
- `algorithm::String`: Name of the algorithm used

For continuous methods, `core_nodes` is the split chosen by the `partition` keyword:
`:median` (default; nodes strictly above the median coreness, so ties at the median
stay peripheral), `:none`, `:threshold` with `threshold=`, or `:core_size` with
`core_size=`. Rombach defaults to `:beta` (the core size implied by `beta`), and
LowRank-Core uses its own Find-Cut partition. `minres_svd_directed` returns a
`CPDirectedResult`, `multiple_cp_pairs` a `CPMultiResult`, and `cp_significance` a
`CPSignificanceResult`; all result types define `==`, `isequal`, and `hash`. See the
contracts page in the documentation for each method's rule.

Every detector except `surprise_cp` and `spectral_method` warns once per session when
the input has more than one connected component with at least two nodes; core-periphery
objectives assume a single connected network. Solvers that fail to converge (ARPACK,
Rossa stationary solvers, incomplete configuration-null switching) raise
`ConvergenceError`.

### Comparing Algorithms

```julia
# Run multiple algorithms on the same network
results = [
    borgatti_everett_continuous(A),
    borgatti_everett_discrete(A),
    lip_discrete(A),
    spectral_method(A),
    random_walker_profiling(A),
    minres_svd(A),
    surprise_cp(A),
    label_switching_cp(A)
]

for r in results
    println("$(r.algorithm): quality = $(round(r.quality, digits=3))")
end

# For multiple CP pairs detection
result_multi = multiple_cp_pairs(A)
println("$(result_multi.algorithm): $(result_multi.n_pairs) pairs detected")
```

### Weighted Networks

Weighted nonnegative adjacency matrices are supported by the continuous
Borgatti-Everett, Rombach, random-walker, MINRES, and multi-pair methods.
Lip, LowRank-Core, and Surprise require binary input. Significance supports binary
ER/configuration nulls and a fixed-topology weighted permutation null. Check each
method's documentation for symmetry and connectivity requirements.

```julia
weighted_edges = [
    (1, 2, 5.0),  # Strong connection
    (1, 3, 2.0),  # Weaker connection
    (2, 3, 3.0)
]
A_weighted = adjacency_to_matrix(weighted_edges, 3)
result = borgatti_everett_continuous(A_weighted)
```

### Graphs.jl and sparse matrices

Loading Graphs.jl activates an optional extension, allowing graph objects to be
passed directly to every detector. `adjacency_to_matrix(graph)` returns a sparse
`Float64` matrix. Sparse inputs use nonzero-based validation and quality kernels,
matrix-free ARPACK eigenpairs for LowRank-Core, sparse matrix-vector products for
directed MINRES, and adjacency-list multi-pair updates.

### Significance testing

`cp_significance(A, detector)` compares observed quality with undirected ER,
undirected or directed degree-preserving, or weight-permutation null samples and returns a
`CPSignificanceResult` with the null distribution, the plus-one Monte Carlo
p-value, and the effect sizes `null_mean`, `null_std`, and `zscore`. Start Julia with
multiple threads and pass `threaded=true` to distribute samples reproducibly across
threads: `thread_schedule=:auto` (and `:greedy`) uses a work queue with one task per
thread, `:static` assigns samples round-robin, and both reproduce the serial result bit
for bit on every supported Julia version. For configuration-style nulls,
`swap_shortfall=:warn` or `:accept` makes incomplete proposal chains explicit, while
the default is to require completion.

## Utility Functions

```julia
# Convert edge list to adjacency matrix (conflicting duplicate dyads raise unless
# combine=:sum|:max|:min|:last; directed=true keeps (i, j) and (j, i) distinct)
A = adjacency_to_matrix(edges, n)

# Arbitrary node identifiers: returns the matrix and the label of every row
A, labels = adjacency_to_matrix([("a", "b"), ("b", "c")])

# Compute quality of a given coreness assignment
q = core_quality(A, coreness_vector)

# Generate ideal core-periphery matrix
Δ = ideal_cp_matrix(coreness_vector)

# Densities of the core-core, core-periphery, and periphery-periphery blocks
densities = block_densities(A, result)

# Get coreness scores from result
scores = coreness_scores(result)
```

## Running Tests

```julia
using Pkg
Pkg.test()
```

## Running Examples

```julia
include("examples/basic_usage.jl")
```

## Cross-Package Benchmark

See the website's [performance and benchmark section](docs/src/performance.md) for the
complete nine-configuration CorePeriphery.jl/`cpnet` timing summary, fit-budget caveats,
and Julia scaling results. [benchmark/COMPARISON.md](benchmark/COMPARISON.md) retains
the reproducible per-fixture recovery, concordance, and timing details.

The [development record and roadmap](docs/src/development.md) documents the scientific
audit, implementation and performance refactor, Rossa paper review, completed
validation, known limitations, and prioritized follow-up work.
Users upgrading from 0.1 should read the [migration guide](docs/src/migration.md) and
[changelog](CHANGELOG.md) before comparing stored scores.

## References

1. Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures. *Social Networks*, 21(4), 375-395. https://doi.org/10.1016/S0378-8733(99)00019-2
2. Lip, S.Z.W. (2011). A fast algorithm for the discrete core/periphery bipartitioning problem. arXiv:1102.5511. https://doi.org/10.48550/arXiv.1102.5511
3. Rombach, M.P., Porter, M.A., Fowler, J.H., Mucha, P.J. (2017). Core-periphery structure in networks (revisited). *SIAM Review*, 59(3), 619-646. https://doi.org/10.1137/17M1130046
4. Cucuringu, M., Rombach, P., Lee, S.H., Porter, M.A. (2016). Detection of core-periphery structure in networks using spectral methods and geodesic paths. *European Journal of Applied Mathematics*, 27(6), 846-887. https://doi.org/10.1017/S095679251600022X
5. Della Rossa, F., Dercole, F., Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. *Scientific Reports*, 3, 1467. https://doi.org/10.1038/srep01467
6. Boyd, J.P., Fitzgerald, W.J., Mahutga, M.C., Smith, D.A. (2010). Computing continuous core/periphery structures for social relations data with MINRES/SVD. *Social Networks*, 32(2), 125-137. https://doi.org/10.1016/j.socnet.2009.09.003
7. Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks. *Physical Review E*, 96(5), 052313. https://doi.org/10.1103/PhysRevE.96.052313
8. Kojaku, S., Masuda, N. (2018). Core-periphery structure requires something else in the network. *New Journal of Physics*, 20, 043012. https://doi.org/10.1088/1367-2630/aab547
9. van Lidth de Jeude, J., Caldarelli, G., Squartini, T. (2019). Detecting core-periphery structures by surprise. *EPL (Europhysics Letters)*, 125(6), 68001. https://doi.org/10.1209/0295-5075/125/68001
10. Yanchenko, E., Sengupta, S. (2025). A label-switching algorithm for fast core-periphery identification. arXiv:2506.02069. https://doi.org/10.48550/arXiv.2506.02069

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

A machine-readable [`CITATION.cff`](CITATION.cff) is included in the repository.

## License

MIT License - see LICENSE file for details.
