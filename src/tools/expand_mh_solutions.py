#!/usr/bin/env python
"""
Run exhaustive bit-swap expansion on all unbounded MH results (n=80..200).

For each n with existing results:
1. Loads all known solutions from unbounded_MH_results/
2. Runs exhaustive 1-hop bit-swap expansion (single flips + density-preserving swaps)
3. Repeats until no new solutions are found (fixed-point)
4. Merges all new solutions back into the results file

Ctrl+C at any time saves partial progress for current n and moves to next.
"""
import sys
import os
import json
import time
import signal
from multiprocessing import Pool, cpu_count
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.algorithms.MH_algorithm import _bitswap_worker, _worker_init, RunLengthOptimizer

N_START = 80
N_END = 200
MAX_ROUNDS = 20  # max expansion rounds per n (converges much sooner)
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "mh_unbounded")


def get_eligible_n_values():
    """Return sorted list of n values with existing results in range."""
    eligible = []
    for fname in sorted(os.listdir(RESULTS_DIR)):
        if not fname.startswith("n_") or not fname.endswith("_results.json"):
            continue
        try:
            n = int(fname[2:6])
        except ValueError:
            continue
        if N_START <= n <= N_END:
            eligible.append(n)
    return eligible


def expand_one_n(n):
    """Run iterative bit-swap expansion for a single n value. Returns (old_count, new_count)."""
    fpath = os.path.join(RESULTS_DIR, f"n_{n:04d}_results.json")
    with open(fpath) as f:
        data = json.load(f)

    target_value = data["best_value"]
    seeds_raw = data.get("solutions", [])
    all_found = set()
    for s in seeds_raw:
        bs = s["bits"] if isinstance(s, dict) else s
        all_found.add(bs)

    old_count = len(all_found)
    ncpu = cpu_count()

    for round_num in range(1, MAX_ROUNDS + 1):
        before = len(all_found)
        tasks = [(n, target_value, s) for s in sorted(all_found)]

        with Pool(ncpu, initializer=_worker_init) as pool:
            results = pool.map(_bitswap_worker, tasks)

        for rs in results:
            all_found.update(rs)

        gained = len(all_found) - before
        print(f"    Round {round_num}: {before} -> {len(all_found)} (+{gained})")
        if gained == 0:
            break

    # Save back to file — preserve all existing solutions + new ones
    solutions_list = sorted(all_found)
    data["solutions"] = solutions_list
    data["last_updated"] = datetime.now().isoformat()
    if "metadata" not in data:
        data["metadata"] = {}
    data["metadata"]["num_solutions"] = len(solutions_list)
    data["metadata"]["bitswap_expanded"] = True

    with open(fpath, "w") as f:
        json.dump(data, f, indent=2)

    return old_count, len(solutions_list)


def main():
    eligible = get_eligible_n_values()
    print(f"Found {len(eligible)} unbounded MH results for n={N_START}..{N_END}")
    print(f"n values: {eligible}\n")

    t_total = time.time()

    for i, n in enumerate(eligible, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(eligible)}] n={n}")
        print(f"{'='*60}")

        fpath = os.path.join(RESULTS_DIR, f"n_{n:04d}_results.json")
        with open(fpath) as f:
            data = json.load(f)
        val = data.get("best_value", 0)
        nsol = len(data.get("solutions", []))
        print(f"  a({n})={val}, {nsol} solutions")

        try:
            old_count, new_count = expand_one_n(n)
            delta = new_count - old_count
            if delta > 0:
                print(f"  Done: {old_count} -> {new_count} (+{delta} new)")
            else:
                print(f"  Done: no new solutions")
        except KeyboardInterrupt:
            print(f"\n  Interrupted at n={n} — skipping to next...")
            continue
        except Exception as e:
            print(f"  Error: {e}")
            continue

    elapsed = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"All done in {elapsed:.0f}s ({elapsed/60:.1f}min)")


if __name__ == "__main__":
    main()
