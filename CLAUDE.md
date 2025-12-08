# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CorePeriphery.jl is a Julia package implementing algorithms for detecting core-periphery structure in networks. Core-periphery structure is a mesoscale pattern where nodes partition into a densely interconnected core and a sparsely connected periphery that attaches to the core but not to each other.

## Common Commands

### Running Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Running Examples

```julia
include("examples/basic_usage.jl")
```

### Building Documentation

```bash
julia --project=docs docs/make.jl
```

### Interactive Development

```julia
using Pkg
Pkg.activate(".")
using CorePeriphery
```

## Architecture

The package is implemented as a single module in [src/CorePeriphery.jl](src/CorePeriphery.jl) (~1200 lines) with no external dependencies beyond Julia's standard library (LinearAlgebra, Statistics, Random).

### Result Types

**`CPResult`** - Standard result returned by most algorithms:

- `coreness::Vector{Float64}` - Coreness score for each node (normalized to [0,1])
- `core_nodes::Vector{Int}` - Indices of nodes classified as core
- `periphery_nodes::Vector{Int}` - Indices of nodes classified as periphery
- `quality::Float64` - Quality score (Pearson correlation with ideal pattern)
- `algorithm::String` - Name of the algorithm used

**`CPMultiResult`** - For multiple core-periphery pairs detection:

- `pair_labels::Vector{Int}` - Pair assignment for each node (1, 2, ..., K)
- `coreness::Vector{Float64}` - Binary coreness within each pair
- `n_pairs::Int` - Number of detected core-periphery pairs
- `quality::Float64` - Q^cp quality score
- `algorithm::String` - Name of the algorithm used

### Algorithms (10 total)

All algorithms accept an adjacency matrix `A::Matrix{Float64}`:

| Algorithm | Function | Type | Key Parameters |
|-----------|----------|------|----------------|
| Borgatti-Everett Continuous | `borgatti_everett_continuous(A)` | Continuous | `max_iter`, `tol`, `init` |
| Borgatti-Everett Discrete | `borgatti_everett_discrete(A)` | Discrete | `max_iter`, `init` |
| Lip's Fast Discrete | `lip_discrete(A)` | Discrete | `max_iter` |
| Rombach Generalized | `rombach_continuous(A)` | Continuous | `alpha`, `beta`, `n_runs` |
| Spectral Method | `spectral_method(A)` | Continuous | None |
| Random Walker Profiling | `random_walker_profiling(A)` | Continuous | `n_walks`, `walk_length` |
| MINRES/SVD | `minres_svd(A)` | Continuous | `max_iter`, `tol` |
| Multiple CP Pairs | `multiple_cp_pairs(A)` | Multi-pair | `max_pairs`, `min_pair_size` |
| Surprise-Based | `surprise_cp(A)` | Discrete | `max_iter` |
| Label Switching | `label_switching_cp(A)` | Discrete | `max_iter` |

**Algorithm Selection Guide:**

- **Fast discrete classification**: `lip_discrete` or `label_switching_cp`
- **Continuous coreness scores**: `borgatti_everett_continuous` or `spectral_method`
- **Directed/asymmetric networks**: `minres_svd`
- **Multiple CP structures**: `multiple_cp_pairs`
- **Tunable core boundary**: `rombach_continuous` (use `alpha` for sharpness, `beta` for core size)

### Utility Functions

- `adjacency_to_matrix(edges, n)` - Convert edge list (tuples) to symmetric adjacency matrix; supports weighted edges as `(i, j, weight)` tuples
- `ideal_cp_matrix(c; discrete=false)` - Generate ideal CP pattern: continuous uses `c[i]*c[j]`, discrete uses `max(c[i], c[j])`
- `core_quality(A, c; discrete=false)` - Compute Pearson correlation between adjacency and ideal pattern
- `coreness_scores(result)` - Extract coreness vector from CPResult

## Key Implementation Details

- Adjacency matrices are always `Matrix{Float64}`, symmetric for undirected networks
- Coreness scores are normalized to [0,1] range after optimization
- Core/periphery classification uses median coreness as threshold
- Weighted networks are supported - edge weights are preserved in the adjacency matrix
- Algorithms use coordinate ascent, greedy swaps, or spectral decomposition depending on the method

## Testing

Tests in [test/runtests.jl](test/runtests.jl) cover:

- Utility functions (matrix conversion, ideal patterns)
- All 10 algorithms with synthetic CP networks
- Edge cases (empty networks, single nodes, complete graphs)
- Consistency checks across algorithm runs
