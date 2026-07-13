#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
RNG = np.random.default_rng(20260712)


def connect_components(A):
    n = len(A)
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in np.flatnonzero(A[i]):
            if int(j) not in seen:
                seen.add(int(j))
                stack.append(int(j))
    while len(seen) < n:
        representative = next(i for i in range(n) if i not in seen)
        A[0, representative] = A[representative, 0] = 1
        seen.add(representative)
        stack = [representative]
        while stack:
            i = stack.pop()
            for j in np.flatnonzero(A[i]):
                if int(j) not in seen:
                    seen.add(int(j))
                    stack.append(int(j))


def sample(name, core, pairs, probabilities):
    n = len(core)
    A = np.zeros((n, n), dtype=int)
    for i in range(n - 1):
        for j in range(i + 1, n):
            p = probabilities(i, j)
            if RNG.random() < p:
                A[i, j] = A[j, i] = 1
    connect_components(A)
    np.savetxt(DATA / f"{name}.csv", A, delimiter=",", fmt="%d")
    with (DATA / f"{name}_truth.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["node", "core", "pair"])
        for i, (is_core, pair) in enumerate(zip(core, pairs), start=1):
            writer.writerow([i, int(is_core), int(pair)])


def main():
    n = 20
    core = np.arange(n) < 5
    pairs = np.ones(n, dtype=int)
    A = np.zeros((n, n), dtype=int)
    for i in range(n - 1):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = int(core[i] or core[j])
    np.savetxt(DATA / "ideal_single.csv", A, delimiter=",", fmt="%d")
    with (DATA / "ideal_single_truth.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["node", "core", "pair"])
        for i in range(n):
            writer.writerow([i + 1, int(core[i]), 1])

    n = 40
    core = np.arange(n) < 10
    pairs = np.ones(n, dtype=int)
    sample(
        "noisy_single",
        core,
        pairs,
        lambda i, j: 0.80 if core[i] and core[j]
        else 0.35 if core[i] or core[j]
        else 0.05,
    )

    n = 48
    pairs = np.repeat([1, 2], n // 2)
    core = np.array([(i % 24) < 6 for i in range(n)])
    sample(
        "two_pairs",
        core,
        pairs,
        lambda i, j: 0.01 if pairs[i] != pairs[j]
        else 0.80 if core[i] and core[j]
        else 0.40 if core[i] or core[j]
        else 0.04,
    )


if __name__ == "__main__":
    main()
