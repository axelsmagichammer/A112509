#!/usr/bin/env python
"""
Compute pairwise Hamming distance statistics for optimal solutions of A112509.

For each n in cached_results.json, computes the min, max, and mean Hamming
distance across all pairs of optimal strings.  Single-solution cases are
recorded as min=max=mean=0.

Writes results to data/hamming_distances.json.
"""
import json
import os
from itertools import combinations

CACHED = os.path.join(os.path.dirname(__file__), "..", "..", "data", "cached_results.json")
OUTPUT = os.path.join(os.path.dirname(__file__), "..", "..", "data", "hamming_distances.json")


def hamming(a: str, b: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def main():
    with open(CACHED) as f:
        cached = json.load(f)

    results = {}
    for n_str in sorted(cached.keys(), key=int):
        entry = cached[n_str]
        sols = entry.get("optimal_strings", [])
        k = len(sols)

        if k <= 1:
            results[n_str] = {"num_solutions": k, "min": 0, "max": 0, "mean": 0.0}
            continue

        dists = [hamming(a, b) for a, b in combinations(sols, 2)]
        d_min = min(dists)
        d_max = max(dists)
        d_mean = sum(dists) / len(dists)

        results[n_str] = {
            "num_solutions": k,
            "min": d_min,
            "max": d_max,
            "mean": round(d_mean, 2),
        }
        print(f"n={n_str:>3}: {k:>4} sols, Hamming min={d_min} max={d_max} mean={d_mean:.2f}")

    with open(OUTPUT, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nWrote {len(results)} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
