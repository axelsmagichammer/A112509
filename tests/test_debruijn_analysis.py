"""Unit tests for src/tools/debruijn_analysis.py"""

import sys
from pathlib import Path

# Ensure src/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from tools.debruijn_analysis import (
    longest_run,
    all_k_substrings_distinct,
    first_repeated_k_substring,
    min_k_all_distinct,
    check_residual_connected,
    analyze_string,
)


# ---------------------------------------------------------------------------
# longest_run
# ---------------------------------------------------------------------------

class TestLongestRun:
    def test_empty(self):
        assert longest_run('', '1') == 0

    def test_all_ones(self):
        assert longest_run('11111', '1') == 5

    def test_all_zeros(self):
        assert longest_run('0000', '0') == 4

    def test_mixed(self):
        assert longest_run('10011101', '1') == 3
        assert longest_run('10011101', '0') == 2

    def test_single_char(self):
        assert longest_run('0', '0') == 1
        assert longest_run('0', '1') == 0
        assert longest_run('1', '1') == 1


# ---------------------------------------------------------------------------
# all_k_substrings_distinct
# ---------------------------------------------------------------------------

class TestAllKSubstringsDistinct:
    def test_k_larger_than_string(self):
        assert all_k_substrings_distinct('010', 5) is True

    def test_trivially_distinct(self):
        # "01" has 1-substrings '0','1' — distinct
        assert all_k_substrings_distinct('01', 1) is True

    def test_repeated(self):
        # "00" has 1-substrings '0','0' — repeated
        assert all_k_substrings_distinct('00', 1) is False

    def test_debruijn_b2_3(self):
        # B(2,3) = "0001011100" (length 10) — all 3-substrings distinct
        db = '0001011100'
        assert all_k_substrings_distinct(db, 3) is True
        # But 2-substrings have repeats (only 4 possible, string has 9)
        assert all_k_substrings_distinct(db, 2) is False

    def test_known_optimal_n5(self):
        # n=5 optimal: "11100" — 3-substrings: '111','110','100' distinct
        s = '11100'
        assert all_k_substrings_distinct(s, 3) is True
        # 2-substrings: '11','11','10','00' — '11' repeats
        assert all_k_substrings_distinct(s, 2) is False
        # 1-substrings: repeated
        assert all_k_substrings_distinct(s, 1) is False


# ---------------------------------------------------------------------------
# first_repeated_k_substring
# ---------------------------------------------------------------------------

class TestFirstRepeated:
    def test_no_repeat(self):
        assert first_repeated_k_substring('0123', 1) is None

    def test_finds_repeat(self):
        result = first_repeated_k_substring('aba', 1)
        assert result is not None
        sub, pos1, pos2 = result
        assert sub == 'a'
        assert pos1 == 0 and pos2 == 2

    def test_binary_repeat(self):
        result = first_repeated_k_substring('00', 1)
        assert result is not None
        assert result[0] == '0'


# ---------------------------------------------------------------------------
# min_k_all_distinct — the core function
# ---------------------------------------------------------------------------

class TestMinKAllDistinct:
    def test_n1_zero(self):
        # "0" — only 1-substring is '0', distinct at k=1
        assert min_k_all_distinct('0') == 1

    def test_n1_one(self):
        assert min_k_all_distinct('1') == 1

    def test_n2_debruijn(self):
        # "01" IS B(2,1) — all 1-substrings distinct
        assert min_k_all_distinct('01') == 1

    def test_n3_optimal(self):
        # "110" — 2-substrings: '11','10' distinct, but 1-substrings: '1','1','0' repeat
        # Actually: 1-subs '1','1','0' — NOT distinct. Need k=2.
        assert min_k_all_distinct('110') == 2

    def test_n5_optimal(self):
        # n=5 optimal strings: '11100' and '11101', both have max run=3
        assert min_k_all_distinct('11100') == 3
        assert min_k_all_distinct('11101') == 3

    def test_short_runs_lower_k(self):
        # '11001' has max run=2, so min_k=2
        assert min_k_all_distinct('11001') == 2

    def test_run_determines_k(self):
        # A string with a long run of 5 ones needs k >= 5
        s = '111110'
        k = min_k_all_distinct(s)
        assert k >= 5
        # Verify: at k, all substrings are distinct
        assert all_k_substrings_distinct(s, k) is True
        # And at k-1, they are NOT
        if k > 1:
            assert all_k_substrings_distinct(s, k - 1) is False

    def test_monotonicity(self):
        # For any string, if k works then k+1 works
        s = '1100101110'
        k = min_k_all_distinct(s)
        for k2 in range(k, len(s) + 1):
            assert all_k_substrings_distinct(s, k2) is True
        for k2 in range(1, k):
            assert all_k_substrings_distinct(s, k2) is False

    def test_known_values_from_cached(self):
        """Spot-check against known best_k from the analysis."""
        # These are (string, expected_min_k) from known optimal strings
        cases = [
            ('0', 1),        # n=1
            ('01', 1),       # n=2
            ('110', 2),      # n=3
            ('11100', 3),    # n=5 optimal
        ]
        for s, expected_k in cases:
            assert min_k_all_distinct(s) == expected_k, (
                "min_k(%r) = %d, expected %d" % (s, min_k_all_distinct(s), expected_k))


# ---------------------------------------------------------------------------
# check_residual_connected
# ---------------------------------------------------------------------------

class TestResidualConnected:
    def test_full_debruijn(self):
        # B(2,2) = "00110" (length 5) uses all 4 edges — trivially connected
        assert check_residual_connected('00110', 2) is True

    def test_small_string_large_k(self):
        # "01" with k=10: uses 0 of 1024 edges, graph still connected
        assert check_residual_connected('01', 10) is True

    def test_single_edge_used(self):
        # "01" with k=2: uses 1 of 4 edges ("01"), residual has "00","10","11"
        assert check_residual_connected('01', 2) is True


# ---------------------------------------------------------------------------
# analyze_string — integration tests
# ---------------------------------------------------------------------------

class TestAnalyzeString:
    def test_returns_required_keys(self):
        r = analyze_string('11001')
        required = {'n', 'min_k', 'debruijn_length', 'longest_ones_run',
                     'longest_zeros_run', 'num_k_substrings',
                     'total_possible_k_substrings', 'embeddable'}
        assert required <= set(r.keys())

    def test_n5(self):
        r = analyze_string('11100')
        assert r['n'] == 5
        assert r['min_k'] == 3
        assert r['debruijn_length'] == (1 << 3) + 3 - 1  # 10
        assert r['longest_ones_run'] == 3
        assert r['longest_zeros_run'] == 2
        assert r['embeddable'] is True

    def test_n3(self):
        r = analyze_string('110')
        assert r['min_k'] == 2
        assert r['debruijn_length'] == 5  # 2^2 + 2 - 1
        assert r['embeddable'] is True

    def test_debruijn_length_formula(self):
        # dB_length = 2^k + k - 1 for every result
        for s in ['0', '01', '110', '11001', '1100101110']:
            r = analyze_string(s)
            k = r['min_k']
            assert r['debruijn_length'] == (1 << k) + k - 1

    def test_all_embeddable_small(self):
        # Every binary string of length <= 7 should be embeddable in some dB
        import itertools
        for n in range(2, 8):
            for bits in itertools.product('01', repeat=n):
                s = ''.join(bits)
                r = analyze_string(s)
                assert r['embeddable'] is True, (
                    "String %r not embeddable?" % s)


# ---------------------------------------------------------------------------
# Cross-validation with cached_results.json
# ---------------------------------------------------------------------------

class TestAgainstCachedResults:
    """Verify analysis matches previously computed results."""

    @staticmethod
    def _load_cached():
        path = Path(__file__).resolve().parent.parent / "data" / "cached_results.json"
        if not path.exists():
            return None
        import json
        return json.loads(path.read_text())

    def test_best_k_matches_longest_run(self):
        """best_k = max longest run across 0 and 1 for the best string."""
        cached = self._load_cached()
        if cached is None:
            return  # skip if no data
        # Check a sample of small n values
        for n_str in ['3', '5', '8', '11', '14', '18', '23']:
            if n_str not in cached:
                continue
            entry = cached[n_str]
            strings = entry.get('optimal_strings', [])
            n = int(n_str)
            valid = [s for s in strings if len(s) == n]
            if not valid:
                continue
            ks = [min_k_all_distinct(s) for s in valid]
            best_k = min(ks)
            # Find the string that achieves best_k
            best_s = valid[ks.index(best_k)]
            max_run = max(longest_run(best_s, '0'), longest_run(best_s, '1'))
            assert best_k == max_run, (
                "n=%s: best_k=%d but max_run=%d" % (n_str, best_k, max_run))

    def test_spot_check_best_k_values(self):
        """Verify best_k for specific n against known results."""
        cached = self._load_cached()
        if cached is None:
            return
        # (n, expected_best_k) from our analysis
        expected = [
            (3, 2), (5, 3), (8, 4), (11, 5), (14, 6), (18, 7),
            (23, 8), (27, 9), (33, 10), (37, 11), (44, 12),
        ]
        for n, exp_k in expected:
            n_str = str(n)
            if n_str not in cached:
                continue
            entry = cached[n_str]
            strings = [s for s in entry.get('optimal_strings', []) if len(s) == n]
            if not strings:
                continue
            best_k = min(min_k_all_distinct(s) for s in strings)
            assert best_k == exp_k, (
                "n=%d: best_k=%d, expected %d" % (n, best_k, exp_k))

    def test_range2_cases(self):
        """n=51, 70, 116 are known to have k-range of 2."""
        cached = self._load_cached()
        if cached is None:
            return
        for n in [51, 70, 116]:
            n_str = str(n)
            if n_str not in cached:
                continue
            entry = cached[n_str]
            strings = [s for s in entry.get('optimal_strings', []) if len(s) == n]
            if not strings:
                continue
            ks = [min_k_all_distinct(s) for s in strings]
            k_range = max(ks) - min(ks)
            assert k_range == 2, (
                "n=%d: expected range 2, got %d" % (n, k_range))
