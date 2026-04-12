"""
Benchmark timing script for OEIS A112509.

For each n from 1 to 80, runs the Metropolis-Hastings algorithm in repeated
restarts until the known value a(n) is matched, recording how long it takes.

Output: benchmark_timing_results.json
"""

import sys
import os
import time
import json
from datetime import datetime
from typing import Optional

# Allow importing from src/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.algorithms.MH_algorithm import RunLengthOptimizer


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_known_values(path: str) -> dict:
    """Parse known_values.txt into {n: a(n)} dict."""
    known = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                known[int(parts[0])] = int(parts[1])
    return known


def iterations_for_n(n: int) -> int:
    """Scale the per-restart iteration budget with n."""
    if n <= 10:
        return 5_000
    elif n <= 20:
        return 20_000
    elif n <= 40:
        return 50_000
    elif n <= 60:
        return 100_000
    else:
        return 200_000


def search_until_target(n: int, target: int, max_wall_seconds: float = 3600.0) -> dict:
    """
    Run MH restarts for n-bit strings until `target` distinct values are found
    (or `max_wall_seconds` elapses).

    Returns a result dict for this n.
    """
    iterations = iterations_for_n(n)
    optimizer = RunLengthOptimizer(n=n, target_density=None, results_dir="results")

    best_overall = 0
    restart_count = 0
    t_start = time.perf_counter()
    t_matched: Optional[float] = None

    print(f"  n={n:3d}  target={target:5d}  iter/restart={iterations:,}", flush=True)

    while True:
        elapsed = time.perf_counter() - t_start

        # Bail out if we've exceeded the wall-clock budget
        if elapsed > max_wall_seconds:
            print(f"    → timeout after {elapsed:.1f}s  (best={best_overall})", flush=True)
            break

        # Vary temperature/cooling to encourage diversity across restarts
        t_init  = 10.0 + (restart_count % 4) * 2.5
        cooling = 0.9995 - (restart_count % 3) * 0.0001

        count, _runs, _bits = optimizer.metropolis_hastings(
            max_iterations=iterations,
            T_initial=t_init,
            cooling_rate=cooling,
            use_adaptive_temp=True,
        )

        restart_count += 1

        if count > best_overall:
            best_overall = count

        if best_overall >= target and t_matched is None:
            t_matched = time.perf_counter() - t_start
            print(
                f"    → matched! best={best_overall}  restarts={restart_count}"
                f"  time={t_matched:.3f}s",
                flush=True,
            )
            break

    total_elapsed = time.perf_counter() - t_start

    return {
        "n": n,
        "best_value_calculated": best_overall,
        "known_value": target,
        "matches_known_value": best_overall >= target,
        "time_to_match_seconds": round(t_matched, 6) if t_matched is not None else None,
        "total_elapsed_seconds": round(total_elapsed, 6),
        "restarts_used": restart_count,
        "timed_out": t_matched is None,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    _base = os.path.join(os.path.dirname(__file__), "..", "..")
    known_path = os.path.join(_base, "data", "reference", "known_values.txt")
    output_path = os.path.join(_base, "output", "benchmarks", "benchmark_timing_results.json")

    known_values = load_known_values(known_path)

    # -----------------------------------------------------------------------
    # Configuration
    # -----------------------------------------------------------------------
    N_START = 1
    N_END   = 80

    # Max seconds to spend per n before giving up and moving on.
    # Increase this if you have more time / want guaranteed matches for large n.
    MAX_SECONDS_PER_N = 3600.0   # 1 hour per n (adjust as needed)
    # -----------------------------------------------------------------------

    results = []
    overall_start = time.perf_counter()

    print(f"A112509 Benchmark Timing  |  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Running n={N_START} to n={N_END}\n")

    for n in range(N_START, N_END + 1):
        if n not in known_values:
            print(f"  n={n}: no known value, skipping")
            continue

        target = known_values[n]
        result = search_until_target(n, target, max_wall_seconds=MAX_SECONDS_PER_N)
        results.append(result)

        # Save after every n so progress is not lost on interruption
        with open(output_path, "w") as f:
            json.dump(
                {
                    "description": "A112509 MH timing benchmark – time to first match a(n)",
                    "generated_at": datetime.now().isoformat(),
                    "results": results,
                },
                f,
                indent=2,
            )

    total_time = time.perf_counter() - overall_start
    print(f"\nAll done. Total wall time: {total_time:.1f}s  ({total_time/60:.1f} min)")
    print(f"Results saved to: {output_path}")

    # Summary table
    matched   = sum(1 for r in results if r["matches_known_value"])
    timed_out = sum(1 for r in results if r["timed_out"])
    print(f"\nMatched {matched}/{len(results)} known values  |  {timed_out} timeouts")


if __name__ == "__main__":
    main()
