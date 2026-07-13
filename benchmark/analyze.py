#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import numpy as np
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parent


def decode(value):
    return np.array([float(item) for item in value.split(";") if item != ""])


def read_results(path):
    with path.open() as handle:
        return {
            (row["dataset"], row["algorithm"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def choose2(value):
    return value * (value - 1) / 2


def adjusted_rand(left, right):
    left_values = {value: i for i, value in enumerate(np.unique(left))}
    right_values = {value: i for i, value in enumerate(np.unique(right))}
    table = np.zeros((len(left_values), len(right_values)), dtype=int)
    for a, b in zip(left, right):
        table[left_values[a], right_values[b]] += 1
    index = sum(choose2(value) for value in table.ravel())
    rows = sum(choose2(value) for value in table.sum(axis=1))
    cols = sum(choose2(value) for value in table.sum(axis=0))
    total = choose2(len(left))
    expected = rows * cols / total if total else 0
    maximum = (rows + cols) / 2
    return 1.0 if maximum == expected else (index - expected) / (maximum - expected)


def auc(scores, truth):
    positives = truth == 1
    n_pos = positives.sum()
    n_neg = len(truth) - n_pos
    ranks = rankdata(scores)
    return (ranks[positives].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def top_k_jaccard(left, right, k):
    left_set = set(np.argsort(-left, kind="stable")[:k])
    right_set = set(np.argsort(-right, kind="stable")[:k])
    return len(left_set & right_set) / len(left_set | right_set)


def fmt(value):
    return "NA" if not math.isfinite(value) else f"{value:.3f}"


def main():
    julia = read_results(ROOT / "coreperiphery_results.tsv")
    python = read_results(ROOT / "cpnet_results.tsv")
    rows = []
    for key in sorted(set(julia) & set(python)):
        dataset, algorithm = key
        left, right = julia[key], python[key]
        if not left["coreness"] or not right["coreness"]:
            continue
        jcore, pcore = decode(left["coreness"]), decode(right["coreness"])
        jpairs, ppairs = decode(left["pairs"]), decode(right["pairs"])
        truth_rows = list(csv.DictReader(
            (ROOT / "data" / f"{dataset}_truth.tsv").open(), delimiter="\t"))
        truth_core = np.array([int(row["core"]) for row in truth_rows])
        truth_pairs = np.array([int(row["pair"]) for row in truth_rows])
        correlation = spearmanr(jcore, pcore).statistic
        rows.append({
            "dataset": dataset,
            "algorithm": algorithm,
            "spearman": float(correlation) if np.isfinite(correlation) else 0.0,
            "topk_jaccard": top_k_jaccard(jcore, pcore, int(truth_core.sum())),
            "julia_auc": auc(jcore, truth_core),
            "cpnet_auc": auc(pcore, truth_core),
            "pair_ari": adjusted_rand(jpairs, ppairs),
            "julia_truth_ari": adjusted_rand(jpairs, truth_pairs),
            "cpnet_truth_ari": adjusted_rand(ppairs, truth_pairs),
            "julia_ms": float(left["runtime_ms"]),
            "cpnet_ms": float(right["runtime_ms"]),
            "julia_quality": float(left["quality"]),
            "cpnet_quality": float(right["quality"]),
        })

    with (ROOT / "comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with (ROOT / "COMPARISON.md").open("w") as handle:
        handle.write("# CorePeriphery.jl vs cpnet empirical comparison\n\n")
        julia_environment = (ROOT / "coreperiphery_environment.txt").read_text().strip()
        cpnet_environment = (ROOT / "cpnet_environment.txt").read_text().strip()
        handle.write("## Environment\n\n```text\n")
        handle.write(julia_environment + "\n" + cpnet_environment + "\n")
        handle.write("```\n\n")
        handle.write(
            "Quality values are reported by each package but are not directly "
            "comparable when objective definitions or scaling differ. Runtime "
            "is the median of five warmed fits; cpnet's Numba compilation is "
            "excluded. NumPy, Python, and Numba RNG state is reset before every "
            "cpnet fit. Both packages receive sparse matrices and BLAS is pinned "
            "to one thread by the benchmark launcher.\n\n"
        )
        handle.write(
            "| Dataset | Algorithm | Spearman | Top-k Jaccard | Julia AUC | "
            "cpnet AUC | Pair ARI | Julia truth ARI | cpnet truth ARI | "
            "Julia ms | cpnet ms |\n"
        )
        handle.write(
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        )
        for row in rows:
            handle.write(
                f"| {row['dataset']} | {row['algorithm']} | "
                f"{fmt(row['spearman'])} | {fmt(row['topk_jaccard'])} | "
                f"{fmt(row['julia_auc'])} | {fmt(row['cpnet_auc'])} | "
                f"{fmt(row['pair_ari'])} | {fmt(row['julia_truth_ari'])} | "
                f"{fmt(row['cpnet_truth_ari'])} | {fmt(row['julia_ms'])} | "
                f"{fmt(row['cpnet_ms'])} |\n"
            )
        speedups = np.array([
            row["cpnet_ms"] / row["julia_ms"] for row in rows
            if row["julia_ms"] > 0 and row["cpnet_ms"] > 0
        ])
        handle.write(
            "\n## Aggregate description\n\n"
            f"CorePeriphery.jl was faster in {(speedups > 1).sum()} of "
            f"{len(speedups)} warmed rows; the median cpnet/Julia runtime ratio "
            f"was {np.median(speedups):.1f}×. This ratio includes each "
            "implementation's configured search procedure: cpnet uses five "
            "random starts for BE, MINRES, Surprise, Rombach, and both KM "
            "methods, whereas CorePeriphery.jl uses five starts only for "
            "Rombach and KM and deterministic searches for BE and Surprise. "
            "The runtime ratio therefore describes these fitted estimators, not "
            "equal low-level work per restart.\n"
        )
        handle.write(
            "\n## Interpretation\n\n"
            "- Lip agrees exactly on the noisy and two-pair core rankings; its "
            "ideal-case difference is a tied boundary convention in cpnet.\n"
            "- LowRank-Core agrees exactly on the ideal graph. On noisy graphs, "
            "CorePeriphery.jl retains rank-two reconstruction scores whereas "
            "cpnet exposes the binary Find-Cut result as coreness.\n"
            "- Rombach, Rossa, and discrete BE show strong agreement on planted "
            "single-pair structure. The headline MINRES comparison now fits the "
            "same symmetric off-diagonal `w*w'` model in both packages. Its "
            "remaining divergence is optimizer behavior: cpnet selects the "
            "largest residual across restarts, uses inconsistent update/scoring "
            "equations, and stops loosely.\n"
            "- On the planted two-pair graph, configuration-null pair recovery "
            "agrees exactly in this seeded fit; ER-null recovery differs "
            "more. Both packages can over-partition single-pair graphs without "
            "a significance-filtering stage.\n"
            "- Surprise agrees on single-pair graphs but uses different "
            "statistics and search procedures, so raw quality values—including "
            "their signs—must not be compared.\n"
            "- CorePeriphery.jl is faster in 25 of 27 rows. Only Rombach is "
            "runtime-competitive here: cpnet is faster on the noisy and "
            "two-pair Rombach fits. The largest "
            "ratios occur for deterministic Julia fits compared with cpnet's "
            "five-start stochastic fits. These tiny-network timings are "
            "descriptive microbenchmarks, not hardware-independent claims.\n"
        )


if __name__ == "__main__":
    main()
