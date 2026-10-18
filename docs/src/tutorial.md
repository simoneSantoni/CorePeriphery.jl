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

A dyad listed twice with the same weight is stored once; a dyad listed with
*conflicting* weights raises an `ArgumentError` unless you say how to resolve it
(`combine=:sum`, `:max`, `:min`, or `:last`). Pass `directed=true` when `(i, j)` and
`(j, i)` are distinct arcs rather than the same undirected dyad.

### Labelled Node Identifiers

Relational data rarely arrives with nodes numbered `1:n`. Omit `n` and the edge list
may use any identifiers; the call returns the matrix together with the label of every
row/column, in order of first appearance (or sorted with `sort_labels=true`):

```jldoctest labelled
julia> using CorePeriphery

julia> friendships = [("alice", "bob"), ("alice", "carol"), ("bob", "carol"),
                      ("alice", "dave"), ("bob", "erin"), ("carol", "frank")];

julia> A, labels = adjacency_to_matrix(friendships);

julia> labels
6-element Vector{String}:
 "alice"
 "bob"
 "carol"
 "dave"
 "erin"
 "frank"

julia> result = lip_discrete(A);

julia> labels[result.core_nodes]
3-element Vector{String}:
 "alice"
 "bob"
 "carol"
```

Index `labels` with `core_nodes`, `periphery_nodes`, or `sortperm(result.coreness)` to
report results in the original identifiers.

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

# Algorithm-specific quality score
println(result.quality)

# Algorithm name
println(result.algorithm)
```

For the continuous detectors `core_nodes` is not part of the model: it is the split
chosen by the `partition` keyword. The default `:median` puts the nodes *strictly*
above the median coreness in the core (ties at the median stay peripheral, so a star's
leaves never join its hub). Make the rule an explicit analytic choice when you report
it:

```julia
borgatti_everett_continuous(A; partition=:none)                      # scores only
borgatti_everett_continuous(A; partition=:threshold, threshold=0.5)  # strictly above 0.5
borgatti_everett_continuous(A; partition=:core_size, core_size=3)    # three highest
```

Isolates (nodes with no ties) always receive coreness `0` and are never core, and a
network with more than one non-trivial connected component triggers a one-time warning:
core-periphery objectives assume one connected network.

### Reporting a Fit

The standard companion of a Borgatti-Everett fit is the density table of the four
blocks induced by the partition (UCINET's density table). `block_densities` computes it
from any result with `core_nodes` (or from an explicit core set):

```jldoctest quickstart
julia> using CorePeriphery

julia> A = adjacency_to_matrix([(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5),
                                (3, 4), (3, 5), (6, 1), (7, 2)], 7);

julia> result = borgatti_everett_continuous(A);

julia> result.core_nodes
3-element Vector{Int64}:
 1
 2
 3

julia> block_densities(A, result)
(core_core = 1.0, core_periphery = 0.6666666666666666, periphery_core = 0.6666666666666666, periphery_periphery = 0.0)
```

Report the fit (`quality`), the core size, the block densities, and a significance test
(below) together; a correlation of 0.6 with a dense periphery block means something
different from the same correlation with an empty one.

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

# Display quality scores. Their scales differ by method, so do not rank
# algorithms by these raw values.
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

```jldoctest jaccard
julia> function jaccard_similarity(set1, set2)
           intersection_size = length(intersect(Set(set1), Set(set2)))
           union_size = length(union(Set(set1), Set(set2)))
           return union_size > 0 ? intersection_size / union_size : 1.0
       end;

julia> jaccard_similarity([1, 2, 3], [2, 3, 4])
0.5
```

Apply it to the core sets returned by two detectors:

```julia
r1 = borgatti_everett_continuous(A)
r2 = spectral_method(A)

similarity = jaccard_similarity(r1.core_nodes, r2.core_nodes)
println("Jaccard similarity of core sets: $similarity")
```

## Working with Specific Algorithms

### Tuning Rombach's Model

The Rombach model has two template parameters: `alpha` sets the sharpness of the
core-periphery boundary and `beta` is the fraction of nodes in the *periphery* (paper
semantics): the first `⌊βN⌋` ranks of the template are peripheral and the remaining
`N − ⌊βN⌋` are core, so a smaller `beta` yields a larger core. The default
`partition=:beta` reports exactly those `N − ⌊βN⌋` nodes as `core_nodes`, so `beta` is
the core size you ask for, and `coreness` is the optimised template normalised to sum
to one:

```julia
using Random

# Sharp boundary, large core (30% of nodes peripheral)
result1 = rombach_continuous(A; alpha=0.9, beta=0.3)

# Soft boundary, small core (70% of nodes peripheral)
result2 = rombach_continuous(A; alpha=0.1, beta=0.7)
```

Do not select `(alpha, beta)` by maximising `result.quality`. The quality is the
equation-(4.10) objective of one template; it depends on the template as well as on the
data, and the paper never selects a single `(alpha, beta)` by it. Treat the parameters
as analyst choices and report a sensitivity sweep of the resulting core sets:

```julia
for beta in (0.2, 0.5, 0.8)
    r = rombach_continuous(A; alpha=0.5, beta=beta, n_runs=5, rng=MersenneTwister(1))
    println("beta = $beta: core_nodes = $(r.core_nodes)")
end
```

Rombach et al.'s own deliverable is the *aggregate core score* of their equation
(4.11): the template-optimal vectors averaged over a grid of `(alpha, beta)` values and
weighted by their objective, scaled so the largest score is one. `rombach_core_score`
computes it (a `10 × 10` grid by default; the paper uses `100 × 100`):

```julia
score = rombach_core_score(A; alphas=0.1:0.1:1.0, betas=0.1:0.1:1.0, n_runs=3,
                           rng=MersenneTwister(1))
score.coreness   # aggregate core score, maximum 1
score.core_nodes # strict-median split of the score (configurable via `partition`)
```

Compare the core sets across settings (for example with `jaccard_similarity` above)
or validate them against an external criterion.

### Random-Walker Requirements

`random_walker_profiling` accepts connected undirected or strongly connected directed,
loop-free networks with nonnegative weights. It uses lower node index for a
reproducible profile by default. To reproduce the original paper's stochastic
selection hierarchy, supply an explicit RNG:

```julia
using Random
paper_result = random_walker_profiling(
    A; tie_break=:paper, rng=MersenneTwister(42),
)
```

Directed inputs use total in-plus-out strength for ties and a lazy stationary
iteration for the row-normalized random walk. Its convergence controls are
`stationary_tol` and `stationary_max_iter`. The legacy `n_walks` and `walk_length`
keywords remain accepted for compatibility but do not affect the result.

For slowly mixing directed graphs, inspect the solver directly or select the robust
dense solve:

```julia
D = [0.0 2 0 1; 0 0 3 0; 4 0 0 1; 1 2 0 0]
stationary = rossa_stationary_distribution(D; method=:linear)
result = random_walker_profiling(D; stationary_method=:linear)
```

When paper-mode ties make node ranks unstable, summarize an ensemble:

```julia
ensemble = rossa_profile_ensemble(
    A; n_runs=200, rng=MersenneTwister(42), threaded=true,
)
ensemble.rank_stability
```

### Directed Networks with MINRES/SVD

For separate sender and receiver roles rather than a stationary-flow profile, use
`minres_svd_directed`:

```julia
# Create a directed network (asymmetric matrix)
A_directed = zeros(5, 5)
A_directed[1, 2] = 1.0  # 1 → 2
A_directed[1, 3] = 1.0  # 1 → 3
A_directed[2, 3] = 1.0  # 2 → 3
A_directed[4, 1] = 1.0  # 4 → 1
A_directed[5, 2] = 1.0  # 5 → 2

result = minres_svd_directed(A_directed)
```

`result.out_coreness` and `result.in_coreness` retain the two directed roles;
`result.coreness` is their combined ranking.

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

The two communities above are disconnected from each other, so this call also
emits the one-time disconnected-input warning; that is expected for a planted
two-pair fixture.

Prefer `multiple_cp_pairs_config(A)`, or
`multiple_cp_pairs(A; null_model=:configuration)`, for the degree-preserving
configuration null (Kojaku & Masuda 2018 equation 13): under the default ER null
any network with heterogeneous degrees exhibits core-periphery pairs, because
high-degree nodes look "core" relative to a uniform density.

The number of pairs is chosen with `pair_selection=:penalized` by default, which
discourages unsupported extra pairs; the raw-objective rule (`:objective`) tends to
over-partition a single planted pair. Inspect `candidate_pair_counts`,
`candidate_qualities`, and `selection_scores` in the result rather than treating
`n_pairs` as uncertainty-free. Isolates are reported with pair label `0`.

## Statistical Significance

```julia
using Random

test_result = cp_significance(
    A,
    lip_discrete;
    null_model=:configuration,
    n_samples=999,            # the default; 199 is enough for a quick look
    rng=MersenneTwister(42),
)
println(test_result.pvalue)      # plus-one Monte Carlo p-value
println(test_result.zscore)      # effect size: (observed - null mean) / null sd
println(test_result.null_mean, " ± ", test_result.null_std)
println(test_result.null_diagnostics.preserved)  # constraint preserved by the null
println(test_result.null_diagnostics.complete)   # did every switching chain complete?
```

Report the p-value together with the effect size and the null model: a tiny p-value
with a z-score of 1.5 on a large network says something different from a z-score of
12. The ER, configuration, and directed-configuration nulls preserve properties of a
simple graph and therefore require binary input; a weighted network raises an error
pointing to `:weight_permutation`.

Additional nulls support directed binary and weighted networks:

```julia
weighted_A = adjacency_to_matrix(
    [(1, 2, 0.5), (1, 3, 1.5), (2, 3, 2.5), (2, 4, 3.5)], 4,
)

directed_test = cp_significance(
    A_directed, minres_svd_directed;
    null_model=:directed_configuration,
    n_samples=199,
)

weighted_test = cp_significance(
    weighted_A, borgatti_everett_continuous;
    null_model=:weight_permutation,
    n_samples=199,
)
```

Configuration switching must complete by default. The returned
`null_diagnostics` records the preserved constraint, requested and accepted
swaps, attempts, and the attempt budget; accepting a shortfall requires an
explicit `swap_shortfall=:warn` or `:accept`.

### Which null model answers which question

The null model is not an implementation detail; it decides what "significant" means.

- `:er` keeps only the number of edges. A small p-value says the network is more
  core-periphery-like than a random graph of the same density.
- `:configuration` keeps every node's degree. A small p-value says the core-periphery
  signal goes beyond degree heterogeneity. A few high-degree nodes attached to many
  low-degree nodes already produce a high core-periphery quality under `:er`; the
  configuration null removes that explanation.
- `:directed_configuration` is the directed analogue (in- and out-degrees preserved,
  reciprocity not preserved).
- `:weight_permutation` keeps the topology and shuffles the weights; it asks whether
  the weights, not the edges, carry the structure.

A network with a planted core can therefore be significant under `:er` and not under
`:configuration`:

```julia
A_planted = generate_cp_network(12, 48; p_cc=0.8, p_cp=0.3, p_pp=0.05, seed=1)
er_test = cp_significance(A_planted, lip_discrete; null_model=:er,
                          n_samples=199, rng=MersenneTwister(42))
cm_test = cp_significance(A_planted, lip_discrete; null_model=:configuration,
                          n_samples=199, rng=MersenneTwister(42))
println("ER p = $(er_test.pvalue); configuration p = $(cm_test.pvalue)")
```

(`generate_cp_network` is defined in the synthetic-network section below.) A
disagreement between the two is the expected outcome when the planted structure is
mostly a degree effect, not a contradiction. Report both p-values and state which null
model each one uses.

Stochastic detectors (`label_switching_cp`, `rombach_continuous`,
`borgatti_everett_discrete` with restarts, ...) receive a fresh child RNG derived from
`rng` on every call because `pass_rng=true` is the default, so the whole test is
reproducible from `rng` alone. When you run many samples, `threaded=true` distributes
them across threads from a precomputed seed stream; serial and threaded runs give the
same result bit for bit. Sparse input is sampled directly into sparse form, which is
the fast path for large networks.

## Quality Assessment

### Computing Quality Manually

```julia
# Get a coreness vector (from algorithm or custom)
c = result.coreness

# Compute quality with continuous ideal pattern
q_cont = core_quality(A, c; discrete=false)

# Compute quality with discrete ideal pattern (strict median split, as in the
# package's default partition rule for continuous methods)
using Statistics
c_binary = Float64.(c .> median(c))
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
                             p_cc=0.8, p_cp=0.4, p_pp=0.05, seed=1)
    rng = MersenneTwister(seed)   # explicit RNG: no global seeding
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

        if rand(rng) < p
            A[i, j] = A[j, i] = 1.0
        end
    end

    return A
end

# Generate and test
A = generate_cp_network(5, 15; seed=42)
result = borgatti_everett_continuous(A)

# True labels (first 5 are core); `core_nodes` is the strict-median split
true_core = Set(1:5)
detected_core = Set(result.core_nodes)

precision = length(intersect(true_core, detected_core)) / length(detected_core)
recall = length(intersect(true_core, detected_core)) / length(true_core)

println("Precision: $(round(precision, digits=3))")
println("Recall: $(round(recall, digits=3))")
```

A sparse random periphery can leave a peripheral node without ties, or split off a
small component; the detectors then give that node coreness `0` and, for a second
non-trivial component, warn once. Check `sum(A, dims=2)` before interpreting such a
fixture.
