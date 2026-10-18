# Causes of lower CorePeriphery.jl/cpnet concordance

This audit uses cpnet commit
`6aad458a6d434a3617d33e74f7163d514a27fecb`, CorePeriphery.jl 0.2.0,
the three planted benchmark networks, and the raw fitted vectors retained by the
comparison harness. The diagnostic script reproduces cpnet's stochastic paths with
the same seeds and then changes one mechanism at a time.

"Concordance" means similarity of exposed node scores or partitions; it is not an
assumption that cpnet is a reference implementation. Several methods expose different
mathematical quantities under the same algorithm name.

## Root-cause summary

| Algorithm | Observed rank concordance | Dominant causes | Causal evidence |
|---|---:|---|---|
| MINRES | 0.320 | cpnet selects the largest residual; its scorer is inconsistent with its gradient; loose stopping; symmetric `ww'` versus Julia `uv'` | Correct selection raises mean rank concordance to 0.595; correct selection plus tight convergence raises it to 0.918 |
| KM-configuration | 0.652 | Missing node-degree factor in cpnet's move gain; different initialization and pair constraints; diagonal-null convention | Cross-evaluation exposes poor local states and pair counts ranging from one to seven |
| KM-ER | 0.658 | Different initialization, feasible pair counts, and local basins; cpnet uses inconsistent densities in optimization and scoring | Both objective evaluators often prefer the same state, while the returned pair counts differ sharply |
| Rombach | 0.837 in the main comparison | cpnet's loop exits after its first improving sweep; fixed initial ordering; shifted template indexing | Completing best-swap search recovers Julia's ordering exactly on the noisy and two-pair graphs |
| LowRank-Core | 0.743 | Julia exposes reconstruction scores; cpnet exposes a binary Find-Cut result; unstable ordering of tied integer scores | Latent reconstruction rankings agree exactly on all three graphs |
| Rossa | 0.778 | Julia's benchmark used its deterministic tie mode; cpnet minimizes a surrogate, omits the paper's secondary strength filter, and reports a different centralization | Julia reproduces equations (5) and (8), and now also exposes the paper's stochastic tie rule and directed stationary formulation; cpnet violates the paper's selection rule repeatedly on the noisy graphs |

## MINRES

cpnet and CorePeriphery.jl do not initially solve the same numerical problem:

- cpnet fits a symmetric rank-one model `A[i,j] ≈ w[i]w[j]` using stochastic ADAM.
- Julia fits a directed rank-one model `A[i,j] ≈ u[i]v[j]` using deterministic
  alternating least squares, then averages normalized in- and out-coreness.
- cpnet stops when relative iterate change falls below `1e-2`; Julia defaults to
  `1e-6`.

There are also two defects in cpnet's implementation. Its gradient corresponds to
the off-diagonal residual

`||A||² - 2w'Aw + ||w||⁴ - sum(w[i]⁴)`,

but its scorer substitutes `||w||²` for `sum(w[i]⁴)`. It then selects
`max(residual)` across restarts even though MINRES must minimize the residual.

| Dataset | Reported ρ | Correct run selection ρ | Selection + `1e-6` ρ | Reported AUC | Corrected AUC |
|---|---:|---:|---:|---:|---:|
| ideal single-pair | 0.246 | 0.734 | 0.753 | 0.653 | 1.000 |
| noisy single-pair | 0.581 | 0.915 | 1.000 | 0.823 | 0.993 |
| two pairs | 0.134 | 0.136 | 1.000 | 0.697 | 0.752 |
| **Mean** | **0.320** | **0.595** | **0.918** | **0.724** | **0.915** |

After corrected selection and convergence, the noisy and two-pair score rankings are
identical to Julia's. The ideal graph retains lower rank correlation because its
structurally equivalent nodes admit tied or nearly tied scores; planted-core recovery
is nevertheless exact. Top-k Jaccard is `1.0` on all three graphs after the corrected
`1e-6` convergence run, compared with `0.429`, `0.429`, and `0.143` for cpnet's
reported MINRES fits.

## KM-ER and KM-configuration

Both implementations optimize a residual-density objective over a pair label and a
binary role for each node, but their state spaces and searches differ:

- cpnet begins with every node in its own pair and every node core-like. It permits
  singleton pairs and has no `max_pairs` or `min_pair_size` constraint.
- Julia uses balanced random assignments over feasible pair counts, initializes roles
  from degree, enforces `max_pairs` and `min_pair_size`, and compares those pair counts
  across restarts.
- Consequently the methods enter different local basins even when the underlying
  objective ordering agrees.

For KM-ER, cpnet optimizes with density `m / choose(n,2)` but scores restarts with
`2m/n²`. Its aggregate score also includes a self-null convention absent from Julia's
strict off-diagonal simple-graph objective.

For KM-configuration, cpnet's move gain subtracts

`(Dcore + x*Dperiphery)/(2m)`

without multiplying by the moving node's degree. The required degree factor appears
in an older commented formula in the same file, but not in the executed code. cpnet
therefore searches using a different expression from the one used to score its final
state. Julia's move gain uses `degree[node] * degree_sum / (2m)` and excludes diagonal
dyads consistently.

The following table cross-evaluates the Julia (J) and cpnet (P) fitted states. Values
are only comparable within one objective column.

| Dataset | Model | Julia objective J / P | cpnet score J / P | Pair count J / P |
|---|---|---:|---:|---:|
| ideal | KM-ER | 0.295 / 0.553 | 0.142 / 0.281 | 4 / 1 |
| noisy | KM-ER | 0.429 / 0.435 | 0.215 / 0.217 | 2 / 7 |
| two pairs | KM-ER | 0.641 / 0.610 | 0.320 / 0.303 | 3 / 4 |
| ideal | KM-config | 0.204 / 0.257 | 0.139 / 0.195 | 2 / 1 |
| noisy | KM-config | 0.254 / 0.087 | 0.227 / 0.058 | 3 / 1 |
| two pairs | KM-config | 0.533 / 0.542 | 0.507 / 0.517 | 2 / 2 |

This shows that KM divergence is mainly a search/state-space problem, amplified in
KM-configuration by an incorrect cpnet move gain. There is no single consistent
"better cpnet state": on the noisy configuration graph, Julia's state is substantially
better under both evaluators, whereas cpnet's state is preferred on the ideal graph.

## Rombach

The broad objective is shared: assign a fixed coreness template to nodes to maximize
`sum(A[i,j]c[i]c[j])`. The implementations differ in three details:

- cpnet uses ranks `0:(n-1)` with `floor(beta*n)`; Julia evaluates ranks `1:n`.
- cpnet always starts from the same ordering and randomizes only the node visitation
  order; Julia randomizes the full template across restarts.
- cpnet's loop condition is inverted: after any improving swap sets `isupdated=true`,
  the `while isupdated is False` condition terminates. Thus an improving run receives
  at most one sweep instead of continuing to a local optimum.

To isolate search from template scaling, both rankings were placed on Julia's template
and cpnet's ordering was completed with the same best-swap rule used by Julia.

| Dataset | Julia objective | cpnet ordering | Completed search | Extra swaps | ρ before / after |
|---|---:|---:|---:|---:|---:|
| ideal | 13.090 | 13.090 | 13.090 | 0 | 0.591 / 0.591 |
| noisy | 34.157 | 33.457 | 34.157 | 24 | 0.935 / 1.000 |
| two pairs | 27.848 | 27.269 | 27.848 | 25 | 0.987 / 1.000 |

The ideal graph has many objective-equivalent permutations, so lower rank concordance
there does not represent a worse optimum. On the other two graphs, early termination
accounts for the measurable difference.

## LowRank-Core

Both packages compute the same two-eigenpair reconstruction, threshold it at `0.5`,
sum the thresholded rows, and apply Find-Cut. The apparent score divergence is largely
an API difference:

- Julia returns normalized reconstruction row sums as `coreness` and stores the cut in
  `core_nodes`.
- cpnet replaces `x_` with the binary Find-Cut vector.

| Dataset | Reconstruction-score ρ | Exposed-score ρ | Find-Cut Jaccard |
|---|---:|---:|---:|
| ideal | 1.000 | 1.000 | 1.000 |
| noisy | 1.000 | 0.556 | 1.000 |
| two pairs | 1.000 | 0.673 | 0.625 |

The two-pair cut difference occurs despite identical reconstruction rankings. Scores
are integer-valued with many ties; Julia uses a stable node-order tie convention while
cpnet uses NumPy's default unstable `argsort`. The packages then inspect different
prefixes. cpnet additionally considers a reversed periphery orientation and excludes
one upper boundary cut that Julia includes. Its final quality scorer computes the
periphery-periphery term from the core vector again; that defect changes reported
quality but not the already selected cut.

## Rossa

The [original Scientific Reports paper](https://www.nature.com/articles/srep01467)
and its [supplementary algorithm and robustness analysis](https://static-content.springer.com/esm/art%3A10.1038%2Fsrep01467/MediaObjects/41598_2013_BFsrep01467_MOESM1_ESM.pdf)
define a three-stage greedy choice in equation (5):

1. start randomly from a minimum-strength node;
2. at each step, minimize the persistence probability of the enlarged set;
3. if the minimum is not unique, retain its minimum-strength candidates and choose
   randomly among them.

CorePeriphery.jl implements the exact persistence ratio, secondary minimum-strength
criterion, and equation (8) cp-centralization. Its default `tie_break=:deterministic`
replaces random choices with lower node index for reproducibility; `tie_break=:paper`
restores both stochastic choices using an explicit RNG. Directed inputs implement the
paper's total-strength tie rule and stationary-distribution persistence equation.

cpnet differs more substantially. For an undirected current peripheral set with volume
`a`, internal weight `b`, candidate strength `d`, and links `l`, equation (5) minimizes
`(b + 2l)/(a + d)`. cpnet instead minimizes `2a*l - b*d` divided by one global constant.
The missing candidate-dependent denominator means this is not an order-preserving
transformation. cpnet then samples across all surrogate minimizers without applying
the paper's secondary minimum-strength filter. Its claimed directed support uses the
same undirected volume calculation rather than the stationary-distribution equation.

| Dataset | cpnet ρ range over 200 seeds | Median ρ | Unique profiles | Mean paper-rule violations per run | Mean surrogate-argmin mismatches |
|---|---:|---:|---:|---:|---:|
| ideal | 0.948–0.997 | 0.974 | 99 | 0.000 | 0.000 |
| noisy | 0.435–0.998 | 0.663 | 200 | 11.385 | 1.120 |
| two pairs | 0.264–0.939 | 0.628 | 200 | 16.950 | 0.960 |

Most cpnet rule violations arise from omitting the minimum-strength tie filter; roughly
one step per noisy run additionally has a different surrogate minimizer from equation
(5). A separate deterministic implementation of equation (5), including the secondary
strength rule, reproduces Julia's complete score vectors exactly on all three graphs.

cpnet reports `1 - mean(alpha)` as graph quality, whereas the paper's equation (8) is
`1 - 2*sum(alpha[1:n-1])/(n-2)`.

| Dataset | Julia equation-(8) C | Equation-(8) C from cpnet profile | cpnet reported quality |
|---|---:|---:|---:|
| ideal | 0.726 | 0.726 | 0.827 |
| noisy | 0.587 | 0.559 | 0.766 |
| two pairs | 0.646 | 0.537 | 0.757 |

Thus Rossa's divergence is not merely random tie variation. Julia follows the paper's
recurrence and centralization, with deterministic and paper-stochastic tie modes and a
directed stationary implementation; cpnet changes the node-selection criterion, tie
hierarchy, directed formulation, and quality scaling.

## Overall conclusion

The low-concordance cases do not share one cause:

1. **MINRES:** mostly correctable cpnet run-selection, scoring, and stopping defects.
2. **KM:** different search spaces and initialization, plus a substantive
   configuration-gain defect and null-model conventions.
3. **Rombach:** premature cpnet termination and initialization strategy.
4. **LowRank-Core:** primarily different exposed outputs and tied-cut ordering.
5. **Rossa:** deterministic specialization in Julia versus a different selection
   surrogate, incomplete tie hierarchy, and non-paper centralization in cpnet.

None of these mechanisms is caused by Julia's sparse or high-performance kernels.
