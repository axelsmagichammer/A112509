"""
Run n=104 to n=111 simultaneously, 2 cores each (16 cores total).
Same parameters as the current MH_algorithm_bounded.py main() config.
"""
import signal
import time
import sys
import os

# Force unbuffered stdout so child-process prints appear immediately
os.environ["PYTHONUNBUFFERED"] = "1"

from multiprocessing import Process, cpu_count
from typing import List, Optional

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.algorithms.MH_algorithm import RunLengthOptimizer

CORES_PER_N = 2

# ── Shared parameters (same as current main() in MH_algorithm_bounded.py) ────

NUM_RESTARTS        = 15_000
ITERATIONS_PER_RUN  = 100_000

MIN_ONES: Optional[List[int]]   = [25, 22, 20, 1]
MAX_ONES: List[Optional[int]]   = [30, 27, 24, 3]
MIN_ZEROS: Optional[List[int]]  = [1, 2, 1, 1]
MAX_ZEROS: List[Optional[int]]  = [2, 3, 2, 2]
MIN_ONES_BLOCKS: Optional[int] = 15
MAX_ONES_BLOCKS: Optional[int] = 22
MIN_DENSITY                     = 0.82
MAX_DENSITY                     = 0.90
TARGET_DENSITY                  = 0.86

RUN_LEVEL_SET_WALK    = True
LEVEL_SET_CHAINS      = CORES_PER_N * 4
LEVEL_SET_STEPS       = 200_000
LEVEL_SET_EPSILON     = 0
LEVEL_SET_EPSILON_MAX = 3
LEVEL_SET_MAX_ROUNDS  = 4
LEVEL_SET_BITSWAP     = True


def run_one(n: int):
    _base = os.path.join(os.path.dirname(__file__), "..", "..")
    optimizer = RunLengthOptimizer(
        n=n,
        target_density=TARGET_DENSITY,
        results_dir=os.path.join(_base, "results", "mh_bounded"),
    )

    best_value, solutions = optimizer.multi_restart_search(
        num_restarts=NUM_RESTARTS,
        iterations_per_run=ITERATIONS_PER_RUN,
        min_ones=MIN_ONES,
        max_ones=MAX_ONES,
        min_zeros=MIN_ZEROS,
        max_zeros=MAX_ZEROS,
        min_density=MIN_DENSITY,
        max_density=MAX_DENSITY,
        num_processes=CORES_PER_N,
        min_ones_blocks=MIN_ONES_BLOCKS,
        max_ones_blocks=MAX_ONES_BLOCKS,
    )

    print(f"\n[n={n}] a({n}) = {best_value}  |  {len(solutions)} solution(s) found")

    if RUN_LEVEL_SET_WALK:
        print(f"\n[n={n}] Starting level-set walk...")
        optimizer.find_all_optimal(
            num_chains=LEVEL_SET_CHAINS,
            steps_per_chain=LEVEL_SET_STEPS,
            epsilon=LEVEL_SET_EPSILON,
            epsilon_max=LEVEL_SET_EPSILON_MAX,
            max_rounds=LEVEL_SET_MAX_ROUNDS,
            bitswap_expand=LEVEL_SET_BITSWAP,
        )


if __name__ == "__main__":
    N_VALUES = list(range(193, 201))  # 193-200 inclusive

    total_cores = len(N_VALUES) * CORES_PER_N
    print(f"Launching n={N_VALUES[0]} to n={N_VALUES[-1]} simultaneously")
    print(f"{CORES_PER_N} cores each × {len(N_VALUES)} values = {total_cores} cores total (of {cpu_count()} available)\n")

    processes = [Process(target=run_one, args=(n,), name=f"n={n}") for n in N_VALUES]

    for p in processes:
        p.start()

    try:
        for p in processes:
            p.join()
    except KeyboardInterrupt:
        print("\nInterrupted — terminating all workers...")
        for p in processes:
            p.terminate()
        for p in processes:
            p.join()

    print("\nAll runs finished.")
