# Tutorial

This tutorial walks through common use cases for CorePeriphery.jl.

## Creating Networks

### From Edge Lists

The most common way to create a network is from an edge list:

```julia
using CorePeriphery

# Unweighted edges as (source, target) tuples
edges = [
    (1, 2), (1, 3), (1, 4),
    (2, 3), (2, 4),
    (3, 4),
    (5, 1), (6, 2)  # Peripheral nodes
]
A = adjacency_to_matrix(edges, 6)
```

### Weighted Networks

For weighted networks, include weights as the third element:

```julia
weighted_edges = [
    (1, 2, 5.0),  # Strong connection
    (1, 3, 3.0),
    (2, 3, 4.0),
    (4, 1, 1.0),  # Weak peripheral connection
]
A = adjacency_to_matrix(weighted_edges, 4)
```

### From Existing Matrices

If you already have an adjacency matrix, you can use it directly:

```julia
using LinearAlgebra

# Create a 6x6 adjacency matrix
A = zeros(6, 6)
A[1, 2] = A[2, 1] = 1.0
A[1, 3] = A[3, 1] = 1.0
# ... etc
```

## Basic Core-Periphery Detection

### Running an Algorithm

All algorithms follow the same pattern:

```julia
result = borgatti_everett_continuous(A)
```

### Understanding Results

The result is a `CPResult` struct with several fields:

```julia
# Coreness scores (0 to 1, higher = more core-like)
println(result.coreness)

# Nodes classified as core
println(result.core_nodes)

# Nodes classified as periphery
println(result.periphery_nodes)

# Quality score (correlation with ideal pattern)
println(result.quality)

# Algorithm name
println(result.algorithm)
```

### Visualizing Results

```julia
# Sort nodes by coreness
sorted_idx = sortperm(result.coreness, rev=true)

println("Node ranking by coreness:")
for (rank, idx) in enumerate(sorted_idx)
    score = round(result.coreness[idx], digits=3)
    label = idx in result.core_nodes ? "CORE" : "PERIPHERY"
    println("  $rank. Node $idx: $score ($label)")
end
```

## Comparing Algorithms

Different algorithms may identify different core-periphery structures:

```julia
# Run multiple algorithms
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

# Compare quality scores
println("Algorithm Comparison:")
println("-" ^ 50)
for r in results
    q = round(r.quality, digits=4)
    n_core = length(r.core_nodes)
    println("$(rpad(r.algorithm, 30)) Quality: $q  Core: $n_core")
end
```

### Measuring Agreement

To compare which nodes different algorithms classify as core:

```julia
using Statistics

function jaccard_similarity(set1, set2)
    intersection = length(intersect(Set(set1), Set(set2)))
    union = length(union(Set(set1), Set(set2)))
    return union > 0 ? intersection / union : 1.0
end

# Compare core node sets
r1 = borgatti_everett_continuous(A)
r2 = spectral_method(A)

similarity = jaccard_similarity(r1.core_nodes, r2.core_nodes)
println("Jaccard similarity of core sets: $similarity")
```

## Working with Specific Algorithms

### Tuning Rombach's Model

The Rombach model has parameters to control the core-periphery structure:

```julia
# Sharp boundary, small core
result1 = rombach_continuous(A; alpha=0.9, beta=0.3)

# Soft boundary, large core
result2 = rombach_continuous(A; alpha=0.1, beta=0.7)

# Grid search over parameters
best_quality = -Inf
best_params = nothing

for alpha in 0.1:0.2:0.9
    for beta in 0.1:0.2:0.9
        r = rombach_continuous(A; alpha=alpha, beta=beta, n_runs=3)
        if r.quality > best_quality
            best_quality = r.quality
            best_params = (alpha=alpha, beta=beta)
        end
    end
end

println("Best parameters: α=$(best_params.alpha), β=$(best_params.beta)")
```

### Random Walker Parameters

Adjust walk parameters based on network size:

```julia
n = size(A, 1)

# For small networks (< 100 nodes)
result_small = random_walker_profiling(A; n_walks=1000, walk_length=10)

# For larger networks
result_large = random_walker_profiling(A; n_walks=n*100, walk_length=20)
```

### Directed Networks with MINRES/SVD

For directed networks, use `minres_svd`:

```julia
# Create a directed network (asymmetric matrix)
A_directed = zeros(5, 5)
A_directed[1, 2] = 1.0  # 1 → 2
A_directed[1, 3] = 1.0  # 1 → 3
A_directed[2, 3] = 1.0  # 2 → 3
A_directed[4, 1] = 1.0  # 4 → 1
A_directed[5, 2] = 1.0  # 5 → 2

result = minres_svd(A_directed)
```

Note: The result averages in-coreness and out-coreness. For separate values, you would need to access the algorithm internals.

## Multiple Core-Periphery Pairs

Some networks contain multiple distinct CP structures:

```julia
# Create a network with two communities, each with CP structure
A = zeros(12, 12)

# Community 1 (nodes 1-6): nodes 1-2 are core
for i in 1:2, j in 1:6
    if i != j
        A[i, j] = A[j, i] = 1.0
    end
end

# Community 2 (nodes 7-12): nodes 7-8 are core
for i in 7:8, j in 7:12
    if i != j
        A[i, j] = A[j, i] = 1.0
    end
end

# Detect multiple pairs
result = multiple_cp_pairs(A)

println("Number of CP pairs detected: $(result.n_pairs)")
println("Pair labels: $(result.pair_labels)")
println("Within-pair coreness: $(result.coreness)")
```

## Quality Assessment

### Computing Quality Manually

```julia
# Get a coreness vector (from algorithm or custom)
c = result.coreness

# Compute quality with continuous ideal pattern
q_cont = core_quality(A, c; discrete=false)

# Compute quality with discrete ideal pattern
c_binary = Float64.(c .>= median(c))
q_disc = core_quality(A, c_binary; discrete=true)

println("Continuous quality: $q_cont")
println("Discrete quality: $q_disc")
```

### Ideal Pattern Visualization

```julia
# Generate the ideal CP matrix
c = [1.0, 0.9, 0.8, 0.3, 0.2, 0.1]  # Example coreness
ideal_cont = ideal_cp_matrix(c; discrete=false)
ideal_disc = ideal_cp_matrix(c; discrete=true)

# Compare with actual adjacency
println("Continuous ideal pattern:")
display(round.(ideal_cont, digits=2))

println("\nDiscrete ideal pattern:")
display(round.(ideal_disc, digits=2))
```

## Synthetic Network Generation

For testing and validation, create networks with known CP structure:

```julia
using Random

function generate_cp_network(n_core, n_periphery;
                             p_cc=0.8, p_cp=0.4, p_pp=0.05, seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    n = n_core + n_periphery
    A = zeros(Float64, n, n)

    for i in 1:n, j in (i+1):n
        # Determine probability based on node types
        if i <= n_core && j <= n_core
            p = p_cc  # Core-core
        elseif i <= n_core || j <= n_core
            p = p_cp  # Core-periphery
        else
            p = p_pp  # Periphery-periphery
        end

        if rand() < p
            A[i, j] = A[j, i] = 1.0
        end
    end

    return A
end

# Generate and test
A = generate_cp_network(5, 15; seed=42)
result = borgatti_everett_continuous(A)

# True labels (first 5 are core)
true_core = Set(1:5)
detected_core = Set(result.core_nodes)

precision = length(intersect(true_core, detected_core)) / length(detected_core)
recall = length(intersect(true_core, detected_core)) / length(true_core)

println("Precision: $(round(precision, digits=3))")
println("Recall: $(round(recall, digits=3))")
```
