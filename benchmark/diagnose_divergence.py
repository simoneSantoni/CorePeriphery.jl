#!/usr/bin/env python3
"""Reproduce targeted CorePeriphery.jl/cpnet divergence diagnostics."""

from pathlib import Path
import csv
import sys

import numpy as np
from scipy import sparse
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parent
UPSTREAM = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/tmp/core-periphery-detection-upstream"
)
sys.path.insert(0, str(UPSTREAM))

from cpnet.KM_ER import _score_ as km_er_score  # noqa: E402
from cpnet.KM_config import _score_ as km_config_score  # noqa: E402
from cpnet.Rossa import Rossa  # noqa: E402
from cpnet.adam import ADAM  # noqa: E402


DATASETS = ("ideal_single", "noisy_single", "two_pairs")


def read_rows(name):
    with (ROOT / name).open() as handle:
        return {
            (row["dataset"], row["algorithm"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def vector(value, dtype=float):
    return np.asarray([dtype(float(item)) for item in value.split(";") if item])


def truth(dataset):
    with (ROOT / "data" / f"{dataset}_truth.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return np.asarray([int(row["core"]) for row in rows])


def auc(scores, labels):
    positive = labels == 1
    n_positive = positive.sum()
    n_negative = len(labels) - n_positive
    ranks = rankdata(scores)
    return (
        ranks[positive].sum() - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)


def rho(left, right):
    value = spearmanr(left, right).statistic
    return float(value) if np.isfinite(value) else 0.0


def jaccard(left, right):
    left, right = set(left), set(right)
    return len(left & right) / len(left | right)


def top_k_jaccard(left, right, k):
    left_nodes = set(np.argsort(-left, kind="stable")[:k])
    right_nodes = set(np.argsort(-right, kind="stable")[:k])
    return len(left_nodes & right_nodes) / len(left_nodes | right_nodes)


def minres_runs(A, seed, tolerance, count=5):
    np.random.seed(seed)
    runs = []
    for _ in range(count):
        w = np.random.rand(A.shape[0])
        adam = ADAM()
        for iteration in range(10_000):
            norm = np.linalg.norm(w)
            gradient = A @ w - (norm**2 - w**2) * w
            candidate = adam.update(w, gradient, 0, False)
            difference = np.linalg.norm(candidate - w) / norm
            w = candidate.copy()
            if difference < tolerance:
                break

        norm = np.linalg.norm(w)
        cross = float(w @ A @ w)
        adjacency_squared = float(np.sum(A.data**2))
        cpnet_score = adjacency_squared - 2 * cross + norm**2 * (norm**2 - 1)
        true_residual = (
            adjacency_squared - 2 * cross + norm**4 - np.sum(w**4)
        )
        runs.append({
            "scores": w,
            "cpnet_score": cpnet_score,
            "true_residual": true_residual,
            "iterations": iteration + 1,
        })
    return runs


def lowrank_reconstruction_scores(A, seed):
    np.random.seed(seed)
    values, vectors = eigs(A, k=2, which="LM")
    reconstructed = np.dot(vectors * diags(values), vectors.T)
    return np.asarray((reconstructed > 0.5).astype(int).sum(axis=0)).reshape(-1)


def deterministic_rossa(A):
    n = len(A)
    degrees = A.sum(axis=1)
    in_periphery = np.zeros(n, dtype=bool)
    scores = np.zeros(n)
    first = min(range(n), key=lambda node: (degrees[node], node))
    in_periphery[first] = True
    volume = degrees[first]
    internal = 0.0
    links = A[:, first].copy()
    for _ in range(1, n):
        candidates = np.flatnonzero(~in_periphery)
        node = min(
            candidates,
            key=lambda candidate: (
                (internal + 2 * links[candidate]) /
                (volume + degrees[candidate]),
                degrees[candidate],
                candidate,
            ),
        )
        internal += 2 * links[node]
        volume += degrees[node]
        scores[node] = internal / volume
        in_periphery[node] = True
        links += A[:, node]
    return scores


def cpnet_rossa_trace(A, seed):
    """Replay cpnet and count deviations from Della Rossa et al. equation (5)."""
    np.random.seed(seed)
    n = len(A)
    degrees = A.sum(axis=0)
    selected = np.zeros(n, dtype=bool)
    node = np.random.choice(np.flatnonzero(degrees == degrees.min()))
    selected[node] = True
    volume = degrees[node]
    internal = 0.0
    profile = np.zeros(n)
    paper_rule_violations = 0
    surrogate_mismatches = 0

    for _ in range(1, n):
        links = A[:, selected].sum(axis=1)
        candidate_alpha = (internal + 2 * links) / (volume + degrees)
        candidate_alpha[selected] = np.inf
        minimum_alpha = candidate_alpha.min()
        alpha_minimizers = np.flatnonzero(candidate_alpha == minimum_alpha)
        minimum_strength = degrees[alpha_minimizers].min()
        paper_eligible = set(
            alpha_minimizers[degrees[alpha_minimizers] == minimum_strength]
        )

        # cpnet minimizes this numerator divided by one global constant rather
        # than equation (5), whose denominator depends on the candidate node.
        surrogate = 2 * volume * links - internal * degrees
        surrogate[selected] = np.inf
        surrogate_minimizers = np.flatnonzero(surrogate == surrogate.min())
        surrogate_mismatches += int(
            set(surrogate_minimizers) != set(alpha_minimizers)
        )
        node = np.random.choice(surrogate_minimizers)
        paper_rule_violations += int(node not in paper_eligible)

        selected[node] = True
        volume += degrees[node]
        internal = A[np.ix_(selected, selected)].sum()
        profile[node] = internal / max(1, volume)

    return profile, paper_rule_violations, surrogate_mismatches


def paper_centralization(profile):
    n = len(profile)
    return 1 - 2 * (profile.sum() - profile.max()) / (n - 2)


def rombach_objective(A, scores):
    return float(scores @ A @ scores / 2)


def greedy_best_swaps(A, scores):
    scores = scores.copy()
    steps = 0
    while True:
        neighbor_scores = A @ scores
        best_delta = 1e-12
        best_pair = None
        for left in range(len(scores) - 1):
            a = scores[left]
            for right in range(left + 1, len(scores)):
                b = scores[right]
                edge = A[left, right]
                delta = (b - a) * (
                    (neighbor_scores[left] - edge * b) -
                    (neighbor_scores[right] - edge * a)
                )
                if delta > best_delta:
                    best_delta = delta
                    best_pair = (left, right)
        if best_pair is None:
            return scores, steps
        left, right = best_pair
        scores[left], scores[right] = scores[right], scores[left]
        steps += 1


def julia_km_objective(A, pairs, roles, model):
    degrees = np.asarray(A.sum(axis=1)).reshape(-1)
    total_weight = degrees.sum() / 2
    density = total_weight / (len(roles) * (len(roles) - 1) / 2)
    quality = 0.0
    for left in range(len(roles) - 1):
        for right in range(left + 1, len(roles)):
            if pairs[left] != pairs[right]:
                continue
            observed = A[left, right]
            expected = density if model == "KM_ER" else (
                degrees[left] * degrees[right] / (2 * total_weight)
            )
            quality += (observed - expected) * max(roles[left], roles[right])
    return quality / total_weight


def main():
    current = read_rows("ablation_current_dense.tsv")
    current_full = read_rows("coreperiphery_results.tsv")
    cpnet = read_rows("cpnet_results.tsv")
    records = []

    def record(algorithm, dataset, metric, value):
        records.append({
            "algorithm": algorithm,
            "dataset": dataset,
            "metric": metric,
            "value": float(value),
        })

    for dataset_index, dataset in enumerate(DATASETS):
        matrix = np.loadtxt(ROOT / "data" / f"{dataset}.csv", delimiter=",")
        A = sparse.csr_matrix(matrix)
        labels = truth(dataset)

        # MINRES: reproduce cpnet, then correct run selection and convergence.
        julia_scores = vector(current[(dataset, "MINRES")]["coreness"])
        seed = 1000 + 100 * dataset_index + 5
        loose_runs = minres_runs(A, seed, 1e-2)
        converged_runs = minres_runs(A, seed, 1e-6)
        reported = max(loose_runs, key=lambda run: run["cpnet_score"])
        corrected = min(loose_runs, key=lambda run: run["true_residual"])
        converged = min(converged_runs, key=lambda run: run["true_residual"])
        raw_cpnet = vector(cpnet[(dataset, "MINRES")]["coreness"])
        if not np.allclose(reported["scores"], raw_cpnet):
            raise RuntimeError(f"failed to reproduce cpnet MINRES for {dataset}")
        for label, run in (
            ("reported", reported),
            ("corrected_selection", corrected),
            ("corrected_convergence", converged),
        ):
            record("MINRES", dataset, f"rho_{label}", rho(julia_scores, run["scores"]))
            record(
                "MINRES", dataset, f"topk_{label}",
                top_k_jaccard(julia_scores, run["scores"], int(labels.sum())),
            )
            record("MINRES", dataset, f"auc_{label}", auc(run["scores"], labels))
            record("MINRES", dataset, f"residual_{label}", run["true_residual"])
            record("MINRES", dataset, f"iterations_{label}", run["iterations"])

        # LowRank-Core: compare latent reconstruction scores and exposed cuts.
        julia_scores = vector(current[(dataset, "LowRankCore")]["coreness"])
        reconstruction = lowrank_reconstruction_scores(
            A, 1000 + 100 * dataset_index + 2,
        )
        cpnet_cut = vector(cpnet[(dataset, "LowRankCore")]["coreness"], int)
        julia_core = vector(current[(dataset, "LowRankCore")]["core_nodes"], int)
        cpnet_core = np.flatnonzero(cpnet_cut == 1) + 1
        if rho(julia_scores, reconstruction) < 1 - 1e-12:
            raise RuntimeError(f"LowRank reconstruction mismatch for {dataset}")
        record("LowRankCore", dataset, "rho_reconstruction", rho(julia_scores, reconstruction))
        record("LowRankCore", dataset, "rho_exposed", rho(julia_scores, cpnet_cut))
        record("LowRankCore", dataset, "cut_jaccard", jaccard(julia_core, cpnet_core))

        # Rossa: quantify stochastic tie sensitivity and deterministic equivalence.
        julia_scores = vector(current[(dataset, "Rossa")]["coreness"])
        correlations = []
        profiles = set()
        paper_violations = []
        surrogate_mismatches = []
        for seed in range(200):
            np.random.seed(seed)
            model = Rossa()
            model.detect(A)
            correlations.append(rho(julia_scores, model.x_))
            profiles.add(tuple(model.x_))
            replay, violations, mismatches = cpnet_rossa_trace(matrix, seed)
            if not np.allclose(model.x_, replay, atol=1e-14):
                raise RuntimeError(f"failed to reproduce cpnet Rossa for {dataset}")
            paper_violations.append(violations)
            surrogate_mismatches.append(mismatches)
        deterministic = deterministic_rossa(matrix)
        if not np.allclose(julia_scores, deterministic, atol=1e-14):
            raise RuntimeError(f"deterministic Rossa mismatch for {dataset}")
        record("Rossa", dataset, "rho_seed_min", min(correlations))
        record("Rossa", dataset, "rho_seed_median", np.median(correlations))
        record("Rossa", dataset, "rho_seed_max", max(correlations))
        record("Rossa", dataset, "unique_profiles_200_seeds", len(profiles))
        record("Rossa", dataset, "paper_rule_violations_mean", np.mean(paper_violations))
        record("Rossa", dataset, "paper_rule_violations_min", min(paper_violations))
        record("Rossa", dataset, "paper_rule_violations_max", max(paper_violations))
        record("Rossa", dataset, "surrogate_argmin_mismatches_mean", np.mean(
            surrogate_mismatches
        ))
        record("Rossa", dataset, "deterministic_max_delta", np.max(np.abs(
            julia_scores - deterministic
        )))
        cpnet_profile = vector(cpnet[(dataset, "Rossa")]["coreness"])
        record("Rossa", dataset, "paper_centralization_julia", paper_centralization(
            julia_scores
        ))
        record("Rossa", dataset, "paper_centralization_cpnet", paper_centralization(
            cpnet_profile
        ))
        record("Rossa", dataset, "cpnet_reported_quality", float(
            cpnet[(dataset, "Rossa")]["quality"]
        ))

        # Rombach: place both orderings on one template, then finish cpnet's search.
        julia_scores = vector(current_full[(dataset, "Rombach")]["coreness"])
        cpnet_scores = vector(cpnet[(dataset, "Rombach")]["coreness"])
        template = np.sort(julia_scores)
        cpnet_order = np.argsort(cpnet_scores, kind="stable")
        common_cpnet = np.empty(len(cpnet_scores))
        common_cpnet[cpnet_order] = template
        reoptimized, steps = greedy_best_swaps(matrix, common_cpnet)
        if not np.isclose(
            rombach_objective(matrix, julia_scores),
            rombach_objective(matrix, reoptimized),
            atol=1e-10,
        ):
            raise RuntimeError(f"Rombach completed-search mismatch for {dataset}")
        record("Rombach", dataset, "rho_common_template", rho(julia_scores, common_cpnet))
        record("Rombach", dataset, "rho_after_completed_search", rho(julia_scores, reoptimized))
        record("Rombach", dataset, "objective_julia", rombach_objective(matrix, julia_scores))
        record("Rombach", dataset, "objective_cpnet_order", rombach_objective(matrix, common_cpnet))
        record("Rombach", dataset, "objective_completed_search", rombach_objective(matrix, reoptimized))
        record("Rombach", dataset, "additional_best_swaps", steps)

        # KM: cross-evaluate both fitted states under both objective conventions.
        for algorithm, scorer in (
            ("KM_ER", km_er_score),
            ("KM_config", km_config_score),
        ):
            julia_pairs = vector(current[(dataset, algorithm)]["pairs"], int)
            julia_roles = vector(current[(dataset, algorithm)]["coreness"])
            cpnet_pairs = vector(cpnet[(dataset, algorithm)]["pairs"], int)
            cpnet_roles = vector(cpnet[(dataset, algorithm)]["coreness"])
            _, julia_pairs = np.unique(julia_pairs, return_inverse=True)
            _, cpnet_pairs = np.unique(cpnet_pairs, return_inverse=True)
            record(algorithm, dataset, "julia_objective_julia_state", julia_km_objective(
                matrix, julia_pairs, julia_roles, algorithm,
            ))
            record(algorithm, dataset, "julia_objective_cpnet_state", julia_km_objective(
                matrix, cpnet_pairs, cpnet_roles, algorithm,
            ))
            record(algorithm, dataset, "cpnet_score_julia_state", scorer(
                A.indptr, A.indices, A.data, julia_pairs, julia_roles,
                len(julia_roles), 0.5,
            )[0])
            record(algorithm, dataset, "cpnet_score_cpnet_state", scorer(
                A.indptr, A.indices, A.data, cpnet_pairs, cpnet_roles,
                len(cpnet_roles), 0.5,
            )[0])
            record(algorithm, dataset, "julia_pair_count", len(np.unique(julia_pairs)))
            record(algorithm, dataset, "cpnet_pair_count", len(np.unique(cpnet_pairs)))

    with (ROOT / "divergence_diagnostics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("algorithm", "dataset", "metric", "value"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
