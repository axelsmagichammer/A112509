"""
Unit tests for SA construction, LCP (Kasai), and score computation.

Tests the core computational primitives in surgical_nudge.py:
  - build_sa_lcp:       suffix array + LCP array construction
  - score_from_sa_lcp:  A112509 score from pre-built SA/LCP
  - compute_score:      full pipeline (builds SA+LCP internally)

Correctness is verified three ways:
  1. Structural properties of SA and LCP arrays (sorted suffixes, correct
     common prefix lengths).
  2. Agreement with brute-force distinct-substring counting from brute_force.py
     for all bit-strings of length ≤ 8.
  3. Hand-computed ground-truth values for selected strings.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "algorithms"))

from surgical_nudge import build_sa_lcp, compute_score, score_from_sa_lcp
from brute_force import distinct_substring_values


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bits(s: str) -> np.ndarray:
    """Convert a binary string to an int32 numpy array."""
    return np.array([int(c) for c in s], dtype=np.int32)


def _brute_score(s: str) -> int:
    """Ground-truth score: number of distinct integer values of all substrings."""
    return len(distinct_substring_values(s))


def _suffix(bits: np.ndarray, pos: int) -> tuple:
    """Return suffix starting at pos as a tuple (for comparison)."""
    return tuple(bits[pos:].tolist())


def _common_prefix_len(bits: np.ndarray, i: int, j: int) -> int:
    """Compute the actual longest common prefix of suffixes starting at i and j."""
    n = len(bits)
    h = 0
    while i + h < n and j + h < n and bits[i + h] == bits[j + h]:
        h += 1
    return h


# ---------------------------------------------------------------------------
# Suffix array structural properties
# ---------------------------------------------------------------------------

class TestSuffixArrayProperties:
    """Verify that build_sa_lcp returns a valid suffix array + LCP array."""

    def _check_sa(self, bits: np.ndarray):
        n = len(bits)
        sa, lcp = build_sa_lcp(bits)
        assert len(sa) == n, "SA length must equal n"
        assert len(lcp) == n, "LCP length must equal n"
        # SA must be a permutation of [0, n)
        assert sorted(sa.tolist()) == list(range(n)), "SA must be a permutation"
        # Suffixes in SA order must be lexicographically sorted
        for i in range(1, n):
            assert _suffix(bits, sa[i - 1]) <= _suffix(bits, sa[i]), (
                f"SA not sorted at position {i}: "
                f"suffix[{sa[i-1]}]={_suffix(bits, sa[i-1])} > "
                f"suffix[{sa[i]}]={_suffix(bits, sa[i])}"
            )

    def _check_lcp(self, bits: np.ndarray):
        n = len(bits)
        sa, lcp = build_sa_lcp(bits)
        # lcp[0] is conventionally 0
        assert lcp[0] == 0, "lcp[0] must be 0"
        # Each lcp[i] must equal the actual common prefix length
        for i in range(1, n):
            expected = _common_prefix_len(bits, int(sa[i - 1]), int(sa[i]))
            assert lcp[i] == expected, (
                f"lcp[{i}] = {lcp[i]}, expected {expected} "
                f"for suffixes at {sa[i-1]} and {sa[i]}"
            )

    @pytest.mark.parametrize("s", [
        "0",
        "1",
        "00",
        "01",
        "10",
        "11",
        "010",
        "101",
        "110",
        "1010",
        "1100",
        "10110",
        "11010",
        "1011010",
        "10110100",
        "11111111",
        "00000000",
        "10101010",
    ])
    def test_sa_sorted(self, s):
        self._check_sa(_bits(s))

    @pytest.mark.parametrize("s", [
        "0",
        "1",
        "00",
        "01",
        "10",
        "11",
        "010",
        "101",
        "110",
        "1010",
        "1100",
        "10110",
        "11010",
        "1011010",
        "10110100",
        "11111111",
        "00000000",
        "10101010",
    ])
    def test_lcp_correct(self, s):
        self._check_lcp(_bits(s))

    def test_sa_all_ones(self):
        # "1111...1" — all suffixes are equal prefixes; SA must be n-1, n-2, ..., 0
        bits = _bits("1" * 10)
        self._check_sa(bits)
        self._check_lcp(bits)

    def test_sa_all_zeros(self):
        bits = _bits("0" * 10)
        self._check_sa(bits)
        self._check_lcp(bits)

    def test_sa_empty(self):
        bits = np.array([], dtype=np.int32)
        sa, lcp = build_sa_lcp(bits)
        assert len(sa) == 0
        assert len(lcp) == 0

    def test_sa_single(self):
        for c in ("0", "1"):
            bits = _bits(c)
            sa, lcp = build_sa_lcp(bits)
            assert sa.tolist() == [0]
            assert lcp.tolist() == [0]

    def test_lcp_monotone_for_repeated(self):
        # "aaaa" style: lcp[i] = n - i for all-equal string
        bits = _bits("0000000")
        sa, lcp = build_sa_lcp(bits)
        # For all-zeros, suffixes sorted: [6,5,4,3,2,1,0] (shortest first)
        # lcp between consecutive should be: 0,1,2,3,4,5,6
        self._check_lcp(bits)


# ---------------------------------------------------------------------------
# Score hand-computed ground truth
# ---------------------------------------------------------------------------

class TestScoreHandComputed:
    """Verify compute_score against manually derived values."""

    @pytest.mark.parametrize("s, expected", [
        # "0"  → substrings: {"0"} → {0} → 1 distinct value
        ("0", 1),
        # "1"  → substrings: {"1"} → {1} → 1 distinct value
        ("1", 1),
        # "10" → "1","10","0" → {1,2,0} → 3
        ("10", 3),
        # "11" → "1","11","1" → {1,3} → 2
        ("11", 2),
        # "00" → "0","00","0" → {0} → 1
        ("00", 1),
        # "110" → "1","11","110","1","10","0" → {1,3,6,2,0} → 5
        ("110", 5),
        # "101" → "1","10","101","0","01","1" → {1,2,5,0} → 4
        # Wait: "01" = 1 already counted. Values: 1,2,5,0,1 → {0,1,2,5} → 4
        ("101", 4),
        # "010" → "0","01","010","1","10","0" → int vals: 0,1,2,1,2,0 → {0,1,2} → 3
        ("010", 3),
    ])
    def test_hand_computed(self, s, expected):
        assert compute_score(_bits(s)) == expected

    def test_all_zeros_score_is_1(self):
        # All-zeros string: every distinct substring evaluates to 0, so only 1 distinct value
        for n in range(1, 10):
            assert compute_score(np.zeros(n, dtype=np.int32)) == 1

    def test_all_ones_score(self):
        # "1"*n → substrings are "1","11","111",... → values 1,3,7,15,...
        # These are all distinct (2^k - 1 for k=1..n), so score = n
        for n in range(1, 12):
            bits = np.ones(n, dtype=np.int32)
            assert compute_score(bits) == n


# ---------------------------------------------------------------------------
# Score agrees with brute-force for all short strings
# ---------------------------------------------------------------------------

class TestScoreVsBruteForce:
    """For all bit-strings of length ≤ 8, compute_score must match brute force."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7, 8])
    def test_all_strings_length_n(self, n):
        for mask in range(2 ** n):
            s = format(mask, f"0{n}b")
            bits = _bits(s)
            expected = _brute_score(s)
            got = compute_score(bits)
            assert got == expected, (
                f"n={n}, s={s!r}: compute_score={got}, brute_force={expected}"
            )


# ---------------------------------------------------------------------------
# API consistency: compute_score vs score_from_sa_lcp(build_sa_lcp(...))
# ---------------------------------------------------------------------------

class TestAPIConsistency:
    """compute_score and score_from_sa_lcp(build_sa_lcp(...)) must agree."""

    @pytest.mark.parametrize("s", [
        "10110100",
        "11010011",
        "10101010",
        "11001100",
        "1" * 20,
        "10" * 10,
        "110" * 7,
        "10110" * 5,
    ])
    def test_two_api_paths_agree(self, s):
        bits = _bits(s)
        sa, lcp = build_sa_lcp(bits)
        via_parts = score_from_sa_lcp(bits, sa, lcp)
        via_full = compute_score(bits)
        assert via_parts == via_full

    def test_score_from_sa_lcp_returns_int(self):
        bits = _bits("10110")
        sa, lcp = build_sa_lcp(bits)
        result = score_from_sa_lcp(bits, sa, lcp)
        assert isinstance(result, int)

    def test_compute_score_returns_int(self):
        assert isinstance(compute_score(_bits("1010")), int)

    def test_sa_dtype(self):
        bits = _bits("10110")
        sa, lcp = build_sa_lcp(bits)
        assert sa.dtype == np.int32
        assert lcp.dtype == np.int32


# ---------------------------------------------------------------------------
# Score monotonicity and bounds
# ---------------------------------------------------------------------------

class TestScoreBounds:
    """Basic sanity checks on score values."""

    @pytest.mark.parametrize("n", range(1, 16))
    def test_score_positive(self, n):
        # Every bit-string has at least 1 distinct substring value
        for mask in range(min(2 ** n, 64)):  # sample first 64
            bits = _bits(format(mask, f"0{n}b"))
            assert compute_score(bits) >= 1

    @pytest.mark.parametrize("n", range(1, 12))
    def test_score_at_most_n_choose_2_plus_1(self, n):
        # Maximum possible distinct substrings = n*(n+1)/2 + 1 (including 0)
        # but score counts integer values, not substrings; still upper-bounded by
        # the number of distinct non-negative integers up to 2^n - 1 plus zero = 2^n
        upper = 2 ** n
        for mask in range(2 ** n):
            bits = _bits(format(mask, f"0{n}b"))
            assert compute_score(bits) <= upper

    def test_score_nonnegative(self):
        # Zero-length string
        assert compute_score(np.array([], dtype=np.int32)) == 0


# ---------------------------------------------------------------------------
# Known A112509 values
# ---------------------------------------------------------------------------

class TestKnownA112509Values:
    """Spot-check that optimal strings achieve the known a(n) values from OEIS."""

    # Known: a(n) = max over all n-bit strings with leading 1 of compute_score
    # Source: OEIS A112509, confirmed by brute_force.a112509
    KNOWN = {
        1: 1,
        2: 3,
        3: 5,
        4: 7,
        5: 10,
        6: 13,
        7: 17,
        8: 22,
        9: 27,
        10: 33,
    }

    @pytest.mark.parametrize("n, expected_max", KNOWN.items())
    def test_max_over_all_strings(self, n, expected_max):
        # A112509 is the max over n-bit strings with a leading 1 (n-bit binary numbers)
        best = 0
        suffix_bits = 0 if n == 1 else (n - 1)
        for mask in range(2 ** suffix_bits):
            s = "1" + format(mask, f"0{suffix_bits}b") if suffix_bits > 0 else "1"
            v = compute_score(_bits(s))
            if v > best:
                best = v
            # also check "0"*n for n=1 (a(1) considers all 1-bit strings)
        if n == 1:
            best = max(best, compute_score(_bits("0")))
        assert best == expected_max, f"a({n}): got {best}, expected {expected_max}"
