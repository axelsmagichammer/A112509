"""
evaluate.py - compute a(n) (distinct binary substring values) for a given bit string.

Uses a suffix-array approach: O(n log^2 n) rather than the naive O(n^2).

Key insight: two substrings have the same integer value iff they are identical
after stripping leading zeros. Therefore:
  distinct values = (1 if any 0-bit exists) + (distinct substrings starting with 1)

Usage:
    python evaluate.py <bitstring>
    python evaluate.py                  # prompts for input
    echo 11011010 | python evaluate.py  # pipe from stdin
"""

import sys


def compute_distinct_values(bits: list[int]) -> int:
    n = len(bits)
    if n == 0:
        return 0

    # Build suffix array via prefix doubling O(n log^2 n)
    sa = list(range(n))
    rank = bits[:]
    k = 1
    while k < n:
        key = [(rank[i], rank[i + k] if i + k < n else -1) for i in range(n)]
        sa.sort(key=lambda i: key[i])
        new_rank = [0] * n
        for i in range(1, n):
            new_rank[sa[i]] = new_rank[sa[i - 1]]
            if key[sa[i]] != key[sa[i - 1]]:
                new_rank[sa[i]] += 1
        rank = new_rank
        if rank[sa[-1]] == n - 1:
            break
        k *= 2

    # Build LCP array via Kasai's algorithm O(n)
    pos = [0] * n
    for i, v in enumerate(sa):
        pos[v] = i
    lcp = [0] * n
    h = 0
    for i in range(n):
        if pos[i] > 0:
            j = sa[pos[i] - 1]
            while i + h < n and j + h < n and bits[i + h] == bits[j + h]:
                h += 1
            lcp[pos[i]] = h
            if h > 0:
                h -= 1

    # Distinct substrings starting with 1
    count = sum((n - sa[i]) - lcp[i] for i in range(n) if bits[sa[i]] == 1)

    # Add value 0 if any zero bit present
    if 0 in bits:
        count += 1

    return count


def main():
    if len(sys.argv) > 1:
        s = sys.argv[1].strip()
    elif not sys.stdin.isatty():
        s = sys.stdin.read().strip()
    else:
        s = input("Enter bit string: ").strip()

    if not s:
        print("Error: empty string", file=sys.stderr)
        sys.exit(1)
    invalid = [c for c in s if c not in "01"]
    if invalid:
        print(f"Error: non-binary characters found: {set(invalid)}", file=sys.stderr)
        sys.exit(1)

    bits = [int(c) for c in s]
    n = len(bits)
    result = compute_distinct_values(bits)

    print(f"n = {n}")
    print(f"a(n) = {result}")


if __name__ == "__main__":
    main()
