"""
Brute force computation of OEIS A112509.

For each n, finds:
  a(n) = maximum number of distinct integers representable by
         contiguous substrings of an n-bit binary number.

Also reports all n-bit strings achieving this maximum.
"""

from itertools import product


def distinct_substring_values(s: str) -> set[int]:
    """Return the set of distinct integer values from all contiguous substrings of s."""
    n = len(s)
    values = set()
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub = s[i:j]
            values.add(int(sub, 2))
    return values


def a112509(n: int) -> tuple[int, list[str]]:
    """
    Compute a(n) and all n-bit strings achieving the maximum.

    An n-bit number has a leading 1, so we iterate over all (n-1)-bit suffixes.

    Returns:
        (max_count, list_of_optimal_strings)
    """
    if n == 1:
        return 1, ["0", "1"]

    best_count = 0
    best_strings = []

    for suffix_bits in product("01", repeat=n - 1):
        s = "1" + "".join(suffix_bits)
        count = len(distinct_substring_values(s))

        if count > best_count:
            best_count = count
            best_strings = [s]
        elif count == best_count:
            best_strings.append(s)

    return best_count, best_strings


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from data.reference.known_values import KNOWN_VALUES

    print(f"{'n':>3}  {'a(n)':>6}  {'check':>6}  {'#opt':>4}  optimal strings")
    print("-" * 80)

    for n in range(1, 21):
        count, strings = a112509(n)
        known = KNOWN_VALUES[n - 1] if n <= len(KNOWN_VALUES) else "?"
        check = "OK" if count == known else "FAIL"

        # Show strings (truncate if too many)
        if len(strings) <= 8:
            str_display = ", ".join(strings)
        else:
            str_display = ", ".join(strings[:4]) + f", ... ({len(strings)} total)"

        print(f"{n:>3}  {count:>6}  {check:>6}  {len(strings):>4}  {str_display}")