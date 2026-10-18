# API Reference

## Module

```@docs
CorePeriphery
```

## Result Types

```@docs
CPResult
CPMultiResult
CPDirectedResult
CPSignificanceResult
CorePeriphery.NullDiagnostics
StationaryDistributionResult
RossaEnsembleResult
```

Every result type defines structural `==`, `isequal`, and `hash`, so results
can be compared, deduplicated, and used as `Dict`/`Set` keys. `==` follows the
field semantics (a `NaN` quality or z-score is not `==` to itself); use
`isequal` to compare results that may carry `NaN` fields.

## Exceptions

```@docs
ConvergenceError
```

## Detection Algorithms

### Continuous Methods

```@docs
borgatti_everett_continuous
rombach_continuous
rombach_core_score
spectral_method
random_walker_profiling
rossa_stationary_distribution
rossa_profile_ensemble
minres_svd
minres_symmetric
minres_svd_directed
```

### Discrete Methods

```@docs
borgatti_everett_discrete
lip_discrete
surprise_cp
label_switching_cp
```

### Multi-Pair Methods

```@docs
multiple_cp_pairs
multiple_cp_pairs_config
```

## Utility Functions

```@docs
adjacency_to_matrix
block_densities
ideal_cp_matrix
core_quality
coreness_scores
cp_significance
```

## Index

```@index
```
