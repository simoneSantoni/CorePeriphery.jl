module CorePeripheryGraphsExt

using CorePeriphery
using Graphs
using SparseArrays

import CorePeriphery:
    adjacency_to_matrix,
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    cp_significance,
    label_switching_cp,
    lip_discrete,
    minres_svd,
    minres_symmetric,
    minres_svd_directed,
    multiple_cp_pairs,
    random_walker_profiling,
    rossa_profile_ensemble,
    rombach_continuous,
    spectral_method,
    surprise_cp

"""
    adjacency_to_matrix(graph::Graphs.AbstractGraph)

Convert a Graphs.jl graph to a sparse Float64 adjacency matrix. Vertex indices
follow Graphs.jl's `vertices(graph)` ordering.
"""
function adjacency_to_matrix(graph::Graphs.AbstractGraph)
    return SparseMatrixCSC{Float64,Int}(Graphs.adjacency_matrix(graph, Float64))
end

for detector in (
    :borgatti_everett_continuous,
    :borgatti_everett_discrete,
    :label_switching_cp,
    :lip_discrete,
    :minres_svd,
    :minres_symmetric,
    :minres_svd_directed,
    :multiple_cp_pairs,
    :random_walker_profiling,
    :rossa_profile_ensemble,
    :rombach_continuous,
    :spectral_method,
    :surprise_cp,
)
    @eval function $detector(graph::Graphs.AbstractGraph; kwargs...)
        return $detector(adjacency_to_matrix(graph); kwargs...)
    end
end

function cp_significance(graph::Graphs.AbstractGraph, detector; kwargs...)
    return cp_significance(adjacency_to_matrix(graph), detector; kwargs...)
end

end
