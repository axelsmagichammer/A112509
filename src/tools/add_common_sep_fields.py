"""
Add K_common and common_seps to every entry in data/cached_results.json.

K_common  : length of the longest separator (0-block) prefix that is
            identical across 100% of optimal solutions for that n.
common_seps: the actual separator values at those positions, e.g. [1, 2, 1].

A 'separator' is a run of 0s that sits between two runs of 1s.
The structure of a solution is:
    [1-block, sep, 1-block, sep, 1-block, ...]
Separators are the 0-blocks in that alternating sequence.
"""

import json
import os

CACHE_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "data", "cached_results.json")


def find_runs(s):
    if not s:
        return []
    runs = []
    cur = s[0]
    cnt = 1
    for c in s[1:]:
        if c == cur:
            cnt += 1
        else:
            runs.append((cur, cnt))
            cur = c
            cnt = 1
    runs.append((cur, cnt))
    return runs


def extract_seps(s):
    """Return the list of separator (0-block) lengths between 1-blocks."""
    runs = find_runs(s)
    seps = []
    i = 0
    while i < len(runs):
        ch, ln = runs[i]
        if ch == '1':
            i += 1
            if i < len(runs) and runs[i][0] == '0':
                seps.append(runs[i][1])
                i += 1
        else:
            break
    return seps


def compute_common_sep_prefix(optimal_strings):
    """Return (K_common, common_seps) for a list of bit-strings."""
    if not optimal_strings:
        return 0, []

    all_seps = [extract_seps(s) for s in optimal_strings]

    if not all_seps or not all_seps[0]:
        return 0, []

    max_possible = max(len(s) for s in all_seps)
    if max_possible == 0:
        return 0, []

    K_common = 0
    for k in range(1, max_possible + 1):
        prefix = all_seps[0][:k]
        if all(len(s) >= k and s[:k] == prefix for s in all_seps):
            K_common = k
        else:
            break

    common_seps = list(all_seps[0][:K_common])
    return K_common, common_seps


def main():
    with open(CACHE_PATH, "r") as f:
        cache = json.load(f)

    updated = 0
    skipped = 0

    for key in sorted(cache.keys(), key=int):
        entry = cache[key]
        strings = entry.get("optimal_strings", [])

        if not strings:
            skipped += 1
            continue

        k, seps = compute_common_sep_prefix(strings)
        entry["K_common"] = k
        entry["common_seps"] = seps
        updated += 1

    with open(CACHE_PATH, "w") as f:
        json.dump(cache, f, indent=2)

    print(f"Updated {updated} entries, skipped {skipped} (no optimal_strings).")
    print(f"Saved to {CACHE_PATH}")

    # Show a sample
    print("\nSample (n=80 to n=90):")
    for key in sorted(cache.keys(), key=int):
        n = int(key)
        if n < 80 or n > 90:
            continue
        entry = cache[key]
        k = entry.get("K_common", "?")
        seps = entry.get("common_seps", "?")
        n_sols = entry.get("num_optimal", "?")
        print(f"  n={n:3d}: K_common={k}, common_seps={seps}  ({n_sols} solutions)")


if __name__ == "__main__":
    main()
