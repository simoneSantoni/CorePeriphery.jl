"""
Basic usage examples for CorePeriphery module.

This script demonstrates how to use different core-periphery detection algorithms.
"""

# Add source directory to path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using CorePeriphery
using LinearAlgebra
using Random

# Set seed for reproducibility
Random.seed!(42)

println("=" ^ 60)
println("Core-Periphery Detection Examples")
println("=" ^ 60)

# Example 1: Create a synthetic core-periphery network
println("\nExample 1: Synthetic Core-Periphery Network")
println("-" ^ 40)

function generate_cp_network(n_core::Int, n_periphery::Int;
                             p_cc::Float64=0.8,  # Core-core connection prob
                             p_cp::Float64=0.4,  # Core-periphery connection prob
                             p_pp::Float64=0.05) # Periphery-periphery connection prob
    n = n_core + n_periphery
    A = zeros(Float64, n, n)

    for i in 1:n
        for j in (i+1):n
            # Determine connection probability based on node types
            if i <= n_core && j <= n_core
                p = p_cc
            elseif i <= n_core || j <= n_core
                p = p_cp
            else
                p = p_pp
            end

            if rand() < p
                A[i, j] = 1.0
                A[j, i] = 1.0
            end
        end
    end

    return A
end

# Generate network with 10 core nodes and 20 periphery nodes
n_core, n_periphery = 10, 20
A = generate_cp_network(n_core, n_periphery)
n = n_core + n_periphery

println("Network: $n_core core nodes, $n_periphery periphery nodes")
println("Total nodes: $n")
println("Total edges: $(Int(sum(A)/2))")

# True labels (1 = core, 0 = periphery)
true_labels = vcat(ones(n_core), zeros(n_periphery))

# Helper function to compute accuracy
function compute_accuracy(result, true_labels)
    predicted = result.coreness .>= median(result.coreness)
    # Check both direct and inverted match
    acc1 = mean(predicted .== (true_labels .== 1.0))
    acc2 = mean(predicted .== (true_labels .== 0.0))
    return max(acc1, acc2)
end

# Run different algorithms
println("\nRunning algorithms...")
println()

# Algorithm 1: Borgatti-Everett Continuous
result_be_cont = borgatti_everett_continuous(A)
acc = compute_accuracy(result_be_cont, true_labels)
println("1. $(result_be_cont.algorithm)")
println("   Quality: $(round(result_be_cont.quality, digits=3))")
println("   Core nodes: $(length(result_be_cont.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 2: Borgatti-Everett Discrete
result_be_disc = borgatti_everett_discrete(A)
acc = compute_accuracy(result_be_disc, true_labels)
println("2. $(result_be_disc.algorithm)")
println("   Quality: $(round(result_be_disc.quality, digits=3))")
println("   Core nodes: $(length(result_be_disc.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 3: Lip's Fast Discrete
result_lip = lip_discrete(A)
acc = compute_accuracy(result_lip, true_labels)
println("3. $(result_lip.algorithm)")
println("   Quality: $(round(result_lip.quality, digits=3))")
println("   Core nodes: $(length(result_lip.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 4: Spectral Method
result_spectral = spectral_method(A)
acc = compute_accuracy(result_spectral, true_labels)
println("4. $(result_spectral.algorithm)")
println("   Quality: $(round(result_spectral.quality, digits=3))")
println("   Core nodes: $(length(result_spectral.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 5: Random Walker Profiling
result_rw = random_walker_profiling(A)
acc = compute_accuracy(result_rw, true_labels)
println("5. $(result_rw.algorithm)")
println("   Quality: $(round(result_rw.quality, digits=3))")
println("   Core nodes: $(length(result_rw.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 6: Rombach Continuous (slower, so use fewer runs)
println("6. Rombach Continuous (running...)")
result_rombach = rombach_continuous(A; n_runs=3)
acc = compute_accuracy(result_rombach, true_labels)
println("   Quality: $(round(result_rombach.quality, digits=3))")
println("   Core nodes: $(length(result_rombach.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 7: MINRES/SVD
result_minres = minres_svd(A)
acc = compute_accuracy(result_minres, true_labels)
println("7. $(result_minres.algorithm)")
println("   Quality: $(round(result_minres.quality, digits=3))")
println("   Core nodes: $(length(result_minres.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 8: Surprise CP
result_surprise = surprise_cp(A)
acc = compute_accuracy(result_surprise, true_labels)
println("8. $(result_surprise.algorithm)")
println("   Quality: $(round(result_surprise.quality, digits=3))")
println("   Core nodes: $(length(result_surprise.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 9: Label Switching CP
result_ls = label_switching_cp(A)
acc = compute_accuracy(result_ls, true_labels)
println("9. $(result_ls.algorithm)")
println("   Quality: $(round(result_ls.quality, digits=3))")
println("   Core nodes: $(length(result_ls.core_nodes))")
println("   Accuracy: $(round(acc * 100, digits=1))%")
println()

# Algorithm 10: Multiple CP Pairs
result_multi = multiple_cp_pairs(A)
println("10. $(result_multi.algorithm)")
println("   Quality: $(round(result_multi.quality, digits=3))")
println("   Number of pairs: $(result_multi.n_pairs)")
println()

# Example 2: Using edge list input
println("\nExample 2: Creating Network from Edge List")
println("-" ^ 40)

edges = [
    (1, 2), (1, 3), (1, 4), (1, 5),  # Node 1 connected to many (core-like)
    (2, 3), (2, 4), (2, 5),           # Node 2 also well-connected
    (3, 4), (3, 5),                   # Node 3
    (4, 5),                           # Node 4
    (6, 1), (7, 2), (8, 3)            # Periphery nodes with single connections
]

A_small = adjacency_to_matrix(edges, 8)
result = borgatti_everett_continuous(A_small)

println("Network with 8 nodes")
println("Coreness scores:")
for i in 1:8
    score = round(result.coreness[i], digits=3)
    label = i in result.core_nodes ? "CORE" : "PERIPHERY"
    println("  Node $i: $score ($label)")
end

# Example 3: Weighted network
println("\nExample 3: Weighted Network")
println("-" ^ 40)

weighted_edges = [
    (1, 2, 5.0), (1, 3, 4.0), (1, 4, 3.0),
    (2, 3, 4.0), (2, 4, 2.0),
    (3, 4, 3.0),
    (5, 1, 1.0), (6, 2, 1.0)
]

A_weighted = adjacency_to_matrix(weighted_edges, 6)
result_weighted = spectral_method(A_weighted)

println("Weighted network analysis:")
for i in 1:6
    score = round(result_weighted.coreness[i], digits=3)
    label = i in result_weighted.core_nodes ? "CORE" : "PERIPHERY"
    println("  Node $i: $score ($label)")
end

println("\n" * "=" ^ 60)
println("Examples completed!")
println("=" ^ 60)
