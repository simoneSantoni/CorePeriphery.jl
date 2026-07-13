#!/usr/bin/env python3
from pathlib import Path
import csv
import random
import subprocess
import sys
import time
import numpy as np
import scipy
import numba
from numba import njit
from scipy import sparse

ROOT = Path(__file__).resolve().parent
UPSTREAM = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/tmp/core-periphery-detection-upstream"
)
sys.path.insert(0, str(UPSTREAM))
import cpnet  # noqa: E402

DATASETS = ("ideal_single", "noisy_single", "two_pairs")


@njit
def seed_numba(seed):
    np.random.seed(seed)


def factories():
    return {
        "BE": lambda: cpnet.BE(num_runs=5),
        "Lip": cpnet.Lip,
        "LowRankCore": lambda: cpnet.LowRankCore(beta=0.1),
        "Rombach": lambda: cpnet.Rombach(num_runs=5, alpha=0.5, beta=0.8),
        "Rossa": cpnet.Rossa,
        "MINRES": lambda: cpnet.MINRES(num_runs=5),
        "Surprise": lambda: cpnet.Surprise(num_runs=5),
        "KM_ER": lambda: cpnet.KM_ER(num_runs=5),
        "KM_config": lambda: cpnet.KM_config(num_runs=5),
    }


def scalar_quality(value):
    array = np.asarray(value, dtype=float).reshape(-1)
    return float(array.sum())


def fit(factory, A, seed):
    np.random.seed(seed)
    random.seed(seed)
    seed_numba(seed)
    warm = factory()
    warm.detect(A)
    timings = []
    model = None
    for _ in range(5):
        np.random.seed(seed)
        random.seed(seed)
        seed_numba(seed)
        model = factory()
        start = time.perf_counter()
        model.detect(A)
        timings.append((time.perf_counter() - start) * 1000)
    return model, float(np.median(timings))


def encode(values):
    return ";".join(str(float(value)) for value in np.asarray(values).reshape(-1))


def main():
    output = ROOT / "cpnet_results.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["package", "dataset", "algorithm", "quality", "runtime_ms",
             "coreness", "pairs"]
        )
        for dataset_index, dataset in enumerate(DATASETS):
            matrix = np.loadtxt(ROOT / "data" / f"{dataset}.csv", delimiter=",")
            A = sparse.csr_matrix(matrix)
            for algorithm_index, (name, factory) in enumerate(factories().items()):
                seed = 1000 + 100 * dataset_index + algorithm_index
                try:
                    model, elapsed = fit(factory, A, seed)
                    writer.writerow(
                        ["cpnet", dataset, name, scalar_quality(model.Q_), elapsed,
                         encode(model.x_), encode(model.c_)]
                    )
                except Exception as error:
                    writer.writerow(
                        ["cpnet", dataset, name, "nan", "nan",
                         "", f"ERROR:{type(error).__name__}:{error}"]
                    )

    commit = subprocess.check_output(
        ["git", "-C", str(UPSTREAM), "rev-parse", "HEAD"], text=True
    ).strip()
    describe = subprocess.check_output(
        ["git", "-C", str(UPSTREAM), "describe", "--tags", "--always"],
        text=True,
    ).strip()
    with (ROOT / "cpnet_environment.txt").open("w") as handle:
        handle.write(f"cpnet {describe}\n")
        handle.write(f"upstream commit {commit}\n")
        handle.write(f"Python {sys.version.split()[0]}\n")
        handle.write(f"NumPy {np.__version__}\n")
        handle.write(f"SciPy {scipy.__version__}\n")
        handle.write(f"Numba {numba.__version__}\n")


if __name__ == "__main__":
    main()
