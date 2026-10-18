"""
    CorePeriphery

A Julia module implementing various core-periphery detection algorithms for network analysis.

Core-periphery structure is a mesoscale pattern in networks where nodes are divided into
a densely interconnected core and a sparsely connected periphery.

# Algorithms Implemented
- Borgatti-Everett continuous model (2000)
- Borgatti-Everett discrete model (2000)
- Lip's fast discrete algorithm (2011)
- Rombach's generalized model (2017)
- Spectral method (Cucuringu et al., 2016)
- Random walker profiling (Rossa et al., 2013)
- MINRES/SVD for asymmetric networks (Boyd et al., 2010)
- Multiple core-periphery pairs (Kojaku & Masuda, 2017)
- Surprise-based detection (Jeude et al., 2019)
- Label-switching algorithm (Yanchenko & Sengupta, 2025)

# References
- Borgatti, S.P., Everett, M.G. (2000). Models of core/periphery structures.
- Lip, S.Z.W. (2011). A Fast Algorithm for the Discrete Core/Periphery Bipartitioning Problem.
- Rombach, M.P., et al. (2017). Core-Periphery Structure in Networks (Revisited).
- Cucuringu, M., et al. (2016). Detection of core-periphery structure using spectral methods.
- Rossa, F.D., et al. (2013). Profiling core-periphery network structure by random walkers.
- Boyd, J.P., et al. (2010). Computing core/periphery structures and permutation tests.
- Kojaku, S., Masuda, N. (2017). Finding multiple core-periphery pairs in networks.
- Jeude, J., et al. (2019). Detecting Core-Periphery Structures by Surprise.
- Yanchenko, K., Sengupta, S. (2025). A fast label-switching algorithm for core-periphery detection.
"""
module CorePeriphery

using LinearAlgebra
using Statistics
using Random
using SparseArrays
using Arpack

export
    # Data structures
    ConvergenceError,
    CPResult,
    CPMultiResult,
    CPDirectedResult,
    CPSignificanceResult,
    StationaryDistributionResult,
    RossaEnsembleResult,
    # Algorithms
    borgatti_everett_continuous,
    borgatti_everett_discrete,
    lip_discrete,
    rombach_continuous,
    rombach_core_score,
    spectral_method,
    random_walker_profiling,
    rossa_stationary_distribution,
    rossa_profile_ensemble,
    minres_svd,
    minres_symmetric,
    minres_svd_directed,
    multiple_cp_pairs,
    multiple_cp_pairs_config,
    surprise_cp,
    label_switching_cp,
    cp_significance,
    # Utilities
    coreness_scores,
    core_quality,
    ideal_cp_matrix,
    adjacency_to_matrix,
    block_densities


include("errors.jl")
include("results.jl")
include("validation.jl")
include("quality.jl")
include("borgatti_everett.jl")
include("lip.jl")
include("rombach.jl")
include("rossa.jl")
include("surprise.jl")
include("label_switching.jl")
include("spectral.jl")
include("directed.jl")
include("multipair.jl")
include("significance.jl")

end # module
