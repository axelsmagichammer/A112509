#!/usr/bin/env python3
"""
De Bruijn embedding analysis for optimal A112509 strings.
=========================================================

For each n with a known optimal solution in cached_results.json,
determines the minimum de Bruijn order k such that the optimal string
is a contiguous substring of a linear de Bruijn sequence B(2,k).

A linear de Bruijn sequence of order k has length 2^k + k - 1 and contains
every k-bit binary string exactly once as a contiguous substring.

The key insight: a string s can be embedded in B(2,k) if and only if
  1. All k-length substrings of s are distinct (no repeated k-grams).
  2. The residual de Bruijn graph (edges not used by s) is weakly connected.

Condition 2 is trivially satisfied when 2^k >> |s| (almost always the case).
Distinctness is monotone: if all k-substrings are distinct, so are all
(k+1)-substrings. So we binary-search for the minimum k.

Usage:
    python tools/debruijn_analysis.py          # analyze all n in cached_results
    python tools/debruijn_analysis.py 80       # analyze just n=80
    python tools/debruijn_analysis.py 50 100   # analyze n=50..100
"""

import json
import math
import sys
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Core analysis functions
# ---------------------------------------------------------------------------

def longest_run(s, ch):
    """Length of the longest consecutive run of character ch in s."""
    best = cur = 0
    for c in s:
        if c == ch:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def all_k_substrings_distinct(s, k):
    """Check if all k-length substrings of s are distinct."""
    n = len(s)
    if k > n:
        return True
    seen = set()
    for i in range(n - k + 1):
        sub = s[i:i+k]
        if sub in seen:
            return False
        seen.add(sub)
    return True


def first_repeated_k_substring(s, k):
    """Return the first repeated k-substring and its positions, or None."""
    n = len(s)
    seen = {}
    for i in range(n - k + 1):
        sub = s[i:i+k]
        if sub in seen:
            return sub, seen[sub], i
        seen[sub] = i
    return None


def min_k_all_distinct(s):
    """Find the minimum k such that all k-substrings of s are distinct.

    Uses the monotonicity property: if all k-substrings are distinct,
    then all (k+1)-substrings are too. This allows binary search.
    """
    n = len(s)

    # Lower bounds
    # 1. Combinatorial: need n - k + 1 <= 2^k
    k_lo = 1
    while n - k_lo + 1 > (1 << k_lo):
        k_lo += 1

    # 2. Longest run: a run of L identical chars means "c^k" repeats for k <= L-1
    #    So we need k >= L at minimum (binary search verifies the exact min).
    max_run = max(longest_run(s, '0'), longest_run(s, '1'))
    k_lo = max(k_lo, max_run)

    # Upper bound: k = n always works (only one n-substring)
    k_hi = n

    # Binary search for minimum k (distinctness is monotone)
    while k_lo < k_hi:
        mid = (k_lo + k_hi) // 2
        if all_k_substrings_distinct(s, mid):
            k_hi = mid
        else:
            k_lo = mid + 1

    return k_lo


def check_residual_connected(s, k):
    """Check if the residual de Bruijn graph is weakly connected after
    embedding s.  Assumes all k-substrings of s are distinct.
    """
    n = len(s)
    total_edges = 1 << k
    num_used = n - k + 1

    if num_used >= total_edges:
        return True  # s IS a de Bruijn sequence

    # For large graphs, removing a tiny fraction of edges can't disconnect
    if total_edges > 64 * num_used:
        return True

    # Explicit check for small graphs
    used = frozenset(s[i:i+k] for i in range(num_used))

    adj = defaultdict(set)
    nodes = set()
    for i in range(total_edges):
        edge = format(i, '0%db' % k)
        if edge not in used:
            u, v = edge[:-1], edge[1:]
            adj[u].add(v)
            adj[v].add(u)
            nodes.add(u)
            nodes.add(v)

    if not nodes:
        return True

    start_node = s[:k-1]
    end_node = s[-(k-1):]

    # All residual nodes must be reachable; if start != end, both must be in component
    required = set(nodes)
    if start_node != end_node:
        required.add(start_node)
        required.add(end_node)

    seed = end_node if end_node in nodes else (
        start_node if start_node in nodes else next(iter(nodes)))

    visited = {seed}
    stack = [seed]
    while stack:
        node = stack.pop()
        for nb in adj.get(node, ()):
            if nb not in visited:
                visited.add(nb)
                stack.append(nb)

    return required <= visited


def analyze_string(s, verbose=False):
    """Find the minimum de Bruijn order k for embedding string s.

    Returns dict with:
        n, min_k, debruijn_length, longest_ones_run, longest_zeros_run,
        num_distinct_at_k, embeddable
    """
    n = len(s)
    max_ones = longest_run(s, '1')
    max_zeros = longest_run(s, '0')

    k = min_k_all_distinct(s)
    num_distinct = n - k + 1  # all are distinct at this k

    embeddable = check_residual_connected(s, k)
    if not embeddable:
        # Search upward (shouldn't happen in practice)
        for k2 in range(k + 1, n + 1):
            if all_k_substrings_distinct(s, k2) and check_residual_connected(s, k2):
                k = k2
                num_distinct = n - k + 1
                embeddable = True
                break

    db_len = (1 << k) + k - 1

    result = {
        'n': n,
        'min_k': k,
        'debruijn_length': db_len,
        'longest_ones_run': max_ones,
        'longest_zeros_run': max_zeros,
        'num_k_substrings': num_distinct,
        'total_possible_k_substrings': 1 << k,
        'embeddable': embeddable,
    }

    if verbose:
        print("  longest 1-run=%d, longest 0-run=%d" % (max_ones, max_zeros))
        print("  min k for all-distinct k-substrings: %d" % k)
        print("  k-substrings used: %d / %d (%.1f%%)" % (
            num_distinct, 1 << k, 100.0 * num_distinct / (1 << k)))
        print("  de Bruijn length: 2^%d + %d = %s" % (k, k - 1, format(db_len, ',')))
        print("  ratio dB_len / n: %.1f" % (db_len / n))

        # Show the first repeated substring at k-1 (for insight)
        if k > 1:
            rep = first_repeated_k_substring(s, k - 1)
            if rep:
                sub, pos1, pos2 = rep
                print("  first repeat at k-1=%d: '%s' at positions %d and %d" % (
                    k - 1, sub, pos1, pos2))

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    data_path = Path(__file__).resolve().parent.parent.parent / "data" / "cached_results.json"
    if not data_path.exists():
        print("Error: %s not found" % data_path)
        sys.exit(1)

    cached = json.loads(data_path.read_text())

    # Parse command-line range
    n_min, n_max = 4, 999
    if len(sys.argv) == 2:
        n_min = n_max = int(sys.argv[1])
    elif len(sys.argv) >= 3:
        n_min, n_max = int(sys.argv[1]), int(sys.argv[2])

    results = []

    for key in sorted(cached.keys(), key=int):
        n = int(key)
        if n < n_min or n > n_max:
            continue

        entry = cached[key]
        strings = entry.get('optimal_strings', [])
        if not strings:
            continue

        valid = [s for s in strings if len(s) == n]
        if not valid:
            continue

        score = entry.get('a(n)', 0)
        num_opt = len(valid)

        # Analyze ALL optimal strings
        all_analyses = [analyze_string(s) for s in valid]
        all_ks = [a['min_k'] for a in all_analyses]
        best_k = min(all_ks)
        worst_k = max(all_ks)
        best_db = (1 << best_k) + best_k - 1
        worst_db = (1 << worst_k) + worst_k - 1

        # Use the string with the best (smallest) k for detailed verbose output
        best_idx = all_ks.index(best_k)
        best_analysis = all_analyses[best_idx]

        print("\nn=%d: a(n)=%d, %d optimal string%s" % (
            n, score, num_opt, 's' if num_opt > 1 else ''))
        print("  best  min_k=%d  (dB_length=%s)" % (best_k, format(best_db, ',')))
        if best_k != worst_k:
            print("  worst min_k=%d  (dB_length=%s)" % (worst_k, format(worst_db, ',')))
            # Count how many have each k
            from collections import Counter
            k_counts = Counter(all_ks)
            print("  k distribution: %s" % ', '.join(
                'k=%d: %d' % (k, c) for k, c in sorted(k_counts.items())))

        r = {
            'n': n,
            'a_n': score,
            'num_optimal': num_opt,
            'best_k': best_k,
            'worst_k': worst_k,
            'best_debruijn_length': best_db,
            'worst_debruijn_length': worst_db,
            'longest_ones_run_at_best': best_analysis['longest_ones_run'],
            'longest_zeros_run_at_best': best_analysis['longest_zeros_run'],
            'all_embeddable': all(a['embeddable'] for a in all_analyses),
            'ratio_best_dblen_n': best_db / n,
        }
        results.append(r)

    # Summary table
    if results:
        print("\n" + "=" * 105)
        print("%5s %8s %5s %6s %6s %5s %17s %10s" % (
            'n', 'a(n)', '#opt', 'bestk', 'wrstk', 'range', 'best_dB_len', 'dB_len/n'))
        print("-" * 105)
        for r in results:
            k_range = r['worst_k'] - r['best_k']
            print("%5d %8d %5d %6d %6d %5d %17s %10.1f" % (
                r['n'], r['a_n'], r['num_optimal'],
                r['best_k'], r['worst_k'], k_range,
                format(r['best_debruijn_length'], ','),
                r['ratio_best_dblen_n']))

        # Save results
        out_path = Path(__file__).resolve().parent.parent.parent / "results" / "debruijn_analysis.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(results, indent=2))
        print("\nResults saved to %s" % out_path)


if __name__ == "__main__":
    main()
