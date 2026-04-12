"""
OEIS A112509 - Template-Based Greedy Search
For large n (e.g., n=1000), exploits the discovered prefix structure pattern.
PREFIX STRUCTURE (discovered from n=150,200,300):
- First 8 zero blocks ALWAYS: [1, 2, 1, 1, 1, 1, 1, 1]
- First 8 ones blocks scale linearly with n
- Tail is more flexible/random
Strategy:
1. Fix the 8 zero blocks in prefix
2. Estimate the 8 ones blocks from linear scaling
3. Greedily optimize tail structure
4. Local search for refinement
"""
import random
import math
import json
import time
import numpy as np
from datetime import datetime
from typing import List, Tuple, Optional
from pathlib import Path
# ---------------------------------------------------------------------------
# Numba-JIT Kasai LCP — compiled once, reused every evaluation call.
# Falls back to a pure-Python version if numba is not installed.
# ---------------------------------------------------------------------------
try:
    from numba import njit as _njit

    @_njit(cache=True)
    def _kasai_jit(bits_arr, sa_arr, pos_arr, n):
        """Kasai LCP construction, JIT-compiled with numba."""
        lcp = np.zeros(n, dtype=np.int32)
        h = np.int32(0)
        for i in range(n):
            if pos_arr[i] > 0:
                j = sa_arr[pos_arr[i] - 1]
                while i + h < n and j + h < n and bits_arr[i + h] == bits_arr[j + h]:
                    h += np.int32(1)
                lcp[pos_arr[i]] = h
                if h > 0:
                    h -= np.int32(1)
        return lcp
    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

    def _kasai_jit(bits_arr, sa_arr, pos_arr, n):  # type: ignore[misc]
        """Pure-Python fallback Kasai when numba is absent."""
        lcp = np.zeros(n, dtype=np.int32)
        h = 0
        for i in range(n):
            if pos_arr[i] > 0:
                j = int(sa_arr[pos_arr[i] - 1])
                while i + h < n and j + h < n and bits_arr[i + h] == bits_arr[j + h]:
                    h += 1
                lcp[pos_arr[i]] = h
                if h > 0:
                    h -= 1
        return lcp

class TemplateGreedyOptimizer:
    """Template-based greedy search for large n"""

    def __init__(self, n: int, results_dir: str = "results", seed_solution: Optional[str] = None):
        self.n = n
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        # seeds/ stores starting bitstrings; results/ stores verified bests only
        self.seeds_dir = Path("seeds")
        self.seeds_dir.mkdir(exist_ok=True)
        self.max_score_cache_size = 50000
        self._score_cache = {}
        # Fixed prefix zeros structure (discovered pattern)
        self.prefix_zeros = [1, 2, 1, 1, 1, 1, 1, 1]
        # Predict prefix ones blocks from linear scaling
        self.prefix_ones_base = self._predict_prefix_ones(n)
        # Optional seed (full bitstring) from theory or prior runs
        self.seed_solution_bits = self._parse_seed_solution(seed_solution)
        if self.seed_solution_bits is None:
            self.seed_solution_bits = self._load_existing_best_solution()
        self.seed_prefix_ones = None
        self.seed_tail = None
        if self.seed_solution_bits is not None:
            self.seed_prefix_ones, self.seed_tail = self._extract_seed_template(self.seed_solution_bits)
            if self.seed_prefix_ones is not None:
                print("Using provided seed to initialize prefix/tail template")
            else:
                print("Seed accepted as baseline, but does not match fixed prefix template")

    def _parse_seed_solution(self, seed_solution: Optional[str]) -> Optional[List[int]]:
        """Parse and validate optional seed solution bitstring"""
        if seed_solution is None:
            return None
        cleaned = ''.join(seed_solution.split())
        if not cleaned:
            return None
        if any(ch not in '01' for ch in cleaned):
            print("Warning: seed_solution contains non-binary characters; ignoring seed")
            return None
        if len(cleaned) != self.n:
            print(f"Warning: seed_solution length {len(cleaned)} != n={self.n}; ignoring seed")
            return None
        return [int(ch) for ch in cleaned]

    def _extract_seed_template(self, seed_bits: List[int]) -> Tuple[Optional[List[int]], Optional[List[int]]]:
        """Extract prefix ones and tail runs from a full seed if template-compatible"""
        runs = self.bits_to_runs(seed_bits)
        if len(runs) < 16:
            return None, None
        seed_prefix_zeros = [runs[i] for i in range(1, 16, 2)]
        if seed_prefix_zeros != self.prefix_zeros:
            return None, None
        seed_prefix_ones = [runs[i] for i in range(0, 16, 2)]
        seed_tail = runs[16:]
        return seed_prefix_ones, seed_tail

    def _load_existing_best_solution(self) -> Optional[List[int]]:
        """Load starting bitstring from disk.
        Priority:
          1. seeds/n_XXXX_seed.txt  — plain bitstring, best available starting point
          2. results/n_XXXX_results.json solutions[0] — verified optimal (fallback)
        """
        # --- 1. seeds/ directory (preferred) ---
        seed_file = self.seeds_dir / f"n_{self.n:04d}_seed.txt"
        if seed_file.exists():
            candidate = seed_file.read_text().strip()
            if len(candidate) == self.n and all(ch in '01' for ch in candidate):
                print(f"Loaded seed from {seed_file.as_posix()}")
                return [int(ch) for ch in candidate]
        # --- 2. results/ fallback ---
        result_file = self.results_dir / f"n_{self.n:04d}_results.json"
        if not result_file.exists():
            return None
        try:
            with open(result_file, 'r') as f:
                data = json.load(f)
        except Exception:
            return None
        solutions = data.get("solutions", [])
        if not solutions:
            return None
        candidate = ''.join(str(solutions[0]).split())
        if len(candidate) != self.n or any(ch not in '01' for ch in candidate):
            return None
        print(f"Loaded seed from {result_file.as_posix()} (verified result)")
        return [int(ch) for ch in candidate]

    def _predict_prefix_ones(self, n: int) -> List[int]:
        """Predict the 8 ones blocks in the prefix based on linear scaling"""
        # Scaling coefficients from linear regression on n=70..300 trusted optimal solutions
        # Maintains the correct b0 > b1 > b2 > b4 > b6 ordering at all n
        # Format: (slope, intercept) for ones[0,2,4,6,8,10,12,14]
        scaling = [
            (0.090425, 10.2216),  # Position 0
            (0.089229,  7.4797),  # Position 2
            (0.087364,  4.6787),  # Position 4
            (0.0,       1.0),     # Position 6 (FIXED)
            (0.087088,  0.5342),  # Position 8
            (0.0,       2.0),     # Position 10 (FIXED)
            (0.081596, -2.3581),  # Position 12
            (0.0,       3.0),     # Position 14 (FIXED)
        ]
        ones = []
        for slope, intercept in scaling:
            value = round(slope * n + intercept)
            ones.append(max(1, value))  # Ensure at least 1
        return ones

    def runs_to_bits(self, runs: List[int]) -> List[int]:
        """Convert run-length encoding to bit array"""
        bits = []
        is_one = True
        for run_length in runs:
            bits.extend([1 if is_one else 0] * run_length)
            is_one = not is_one
        return bits

    def bits_to_runs(self, bits: List[int]) -> List[int]:
        """Convert bit array to run-length encoding"""
        if not bits:
            return []
        runs = []
        current_bit = bits[0]
        current_length = 1
        for bit in bits[1:]:
            if bit == current_bit:
                current_length += 1
            else:
                runs.append(current_length)
                current_bit = bit
                current_length = 1
        runs.append(current_length)
        return runs

    def compute_distinct_values(self, bits: List[int]) -> int:
        """Count distinct substring values using suffix array - O(n log^2 n).
        Key insight: two substrings have the same integer value iff they are
        identical after stripping leading zeros.  Therefore:
          distinct values = (1 if any 0-bit exists)
                          + (number of distinct substrings starting with 1)
        The second term is computed from the suffix array (SA) and LCP array:
          sum over i where bits[SA[i]]==1 of  (n - SA[i] - LCP[i])
        which is the standard SA distinct-substring formula restricted to
        suffixes whose first character is 1.
        Implementation uses numpy for SA construction (C-level argsort) and
        vectorised counting; only the O(n) Kasai LCP loop stays in Python.
        """
        n = len(bits)
        if n == 0:
            return 0
        # --- Build suffix array via prefix doubling using numpy argsort ---
        bits_np = np.asarray(bits, dtype=np.int32)
        idx = np.arange(n, dtype=np.int32)
        rank = bits_np.copy()
        k = 1
        while k < n:
            next_idx = np.minimum(idx + k, n - 1)
            r2 = np.where(idx + k < n, rank[next_idx], -1)
            max_r = int(rank.max()) + 2
            key = rank.astype(np.int64) * max_r + r2.astype(np.int64)
            sa = np.argsort(key, kind='stable')
            changed = np.empty(n, dtype=np.int32)
            changed[0] = 0
            changed[1:] = (key[sa[1:]] != key[sa[:-1]]).cumsum()
            rank = np.empty(n, dtype=np.int32)
            rank[sa] = changed
            if rank[sa[-1]] == n - 1:   # all ranks distinct → done
                break
            k *= 2
        # --- Build LCP array via Kasai (numba JIT or Python fallback) ---
        pos = np.empty(n, dtype=np.int32)
        pos[sa] = np.arange(n, dtype=np.int32)
        lcp = _kasai_jit(bits_np, sa, pos, n)
        # --- Count distinct substrings starting with 1 (vectorised) ---
        mask = bits_np[sa] == 1
        count = int(np.sum((n - sa[mask]) - lcp[mask]))
        # --- Add value 0 if any zero bit present ---
        if 0 in bits:
            count += 1
        return count

    def _evaluate_runs(self, runs: List[int]) -> int:
        """Evaluate objective for a run sequence with caching"""
        key = tuple(runs)
        cached = self._score_cache.get(key)
        if cached is not None:
            return cached
        value = self.compute_distinct_values(self.runs_to_bits(runs))
        if len(self._score_cache) >= self.max_score_cache_size:
            oldest_key = next(iter(self._score_cache))
            self._score_cache.pop(oldest_key)
        self._score_cache[key] = value
        return value

    def construct_from_template(self, prefix_ones: List[int], tail_runs: List[int]) -> List[int]:
        """Construct full run-length encoding from prefix template + tail"""
        runs = []
        # Interleave prefix ones and zeros
        for i in range(8):
            runs.append(prefix_ones[i])
            runs.append(self.prefix_zeros[i])
        # Add tail runs
        runs.extend(tail_runs)
        return runs

    def greedy_tail_construction(self, prefix_ones: List[int], target_total: int,
                                 max_tail_runs: int = 30,
                                 seed_tail: Optional[List[int]] = None,
                                 strategy_mode: str = "full",
                                 randomize_tail_ratio: float = 0.2,
                                 randomize_tail_variants: int = 4) -> List[int]:
        """Greedily construct the tail to maximize distinct values
        Args:
            prefix_ones: The 8 ones blocks for the prefix
            target_total: Total number of bits (n)
            max_tail_runs: Maximum number of runs to use in tail
        Returns:
            List of tail run lengths
        """
        prefix_bits = sum(prefix_ones) + sum(self.prefix_zeros)
        remaining_bits = target_total - prefix_bits
        print(f"Prefix uses {prefix_bits} bits, {remaining_bits} bits remaining for tail")
        if remaining_bits <= 0:
            return []
        # Start with a simple tail: try different structures
        best_tail = None
        best_value = 0
        # Try several different tail initialization strategies
        strategies = []
        valid_seed_tail = seed_tail is not None and sum(seed_tail) == remaining_bits and len(seed_tail) > 0
        if strategy_mode == "seed_only":
            if valid_seed_tail:
                assert seed_tail is not None
                print("  Using seed tail only (strategy_mode=seed_only)")
                return seed_tail.copy()
            print("  seed_only requested but no valid seed tail; falling back to decreasing")
            return self._tail_strategy_decreasing(remaining_bits)
        if strategy_mode == "decreasing_only":
            if valid_seed_tail:
                assert seed_tail is not None
                strategies.append(seed_tail.copy())
            strategies.append(self._tail_strategy_decreasing(remaining_bits))
        elif strategy_mode == "best_or_randomized_tail":
            if valid_seed_tail:
                assert seed_tail is not None
                strategies.append(seed_tail.copy())
                # Diversify around current best by randomizing last X% of full bits
                base_runs = self.construct_from_template(prefix_ones, seed_tail)
                base_bits = self.runs_to_bits(base_runs)
                for _ in range(max(0, randomize_tail_variants)):
                    candidate_bits = self._randomize_last_bits(base_bits, randomize_tail_ratio)
                    candidate_runs = self.bits_to_runs(candidate_bits)
                    if len(candidate_runs) >= 16:
                        candidate_tail = candidate_runs[16:]
                        if len(candidate_tail) > 0 and sum(candidate_tail) == remaining_bits:
                            strategies.append(candidate_tail)
                # Keep one deterministic fallback
                strategies.append(self._tail_strategy_decreasing(remaining_bits))
            else:
                # No valid seed tail available; fall back to standard exploration
                strategies.extend([
                    self._tail_strategy_decreasing(remaining_bits),
                    self._tail_strategy_balanced(remaining_bits),
                    self._tail_strategy_fibonacci(remaining_bits),
                    self._tail_strategy_random(remaining_bits, max_tail_runs),
                ])
        else:
            if valid_seed_tail:
                assert seed_tail is not None
                strategies.append(seed_tail.copy())
            strategies.extend([
                self._tail_strategy_decreasing(remaining_bits),
                self._tail_strategy_balanced(remaining_bits),
                self._tail_strategy_fibonacci(remaining_bits),
                self._tail_strategy_random(remaining_bits, max_tail_runs),
            ])
        for i, tail in enumerate(strategies):
            runs = self.construct_from_template(prefix_ones, tail)
            if sum(runs) == target_total:
                value = self._evaluate_runs(runs)
                print(f"  Strategy {i+1}: {value} distinct values, {len(tail)} tail runs")
                if value > best_value:
                    best_value = value
                    best_tail = tail
        return best_tail if best_tail is not None else strategies[0]

    def _randomize_last_bits(self, bits: List[int], ratio: float) -> List[int]:
        """Return a copy of bits with last ratio fraction randomized"""
        if not bits:
            return []
        clamped_ratio = min(max(ratio, 0.0), 1.0)
        random_count = max(1, int(len(bits) * clamped_ratio))
        start = len(bits) - random_count
        randomized = bits.copy()
        for i in range(start, len(randomized)):
            randomized[i] = random.randint(0, 1)
        return randomized

    def _tail_strategy_decreasing(self, remaining_bits: int) -> List[int]:
        """Tail with decreasing run lengths"""
        tail = []
        used = 0
        length = min(20, remaining_bits // 2)
        while used < remaining_bits and length > 0:
            if used + length <= remaining_bits:
                tail.append(length)
                used += length
                length -= 1
            else:
                length -= 1
        # Add final bits
        if used < remaining_bits:
            tail.append(remaining_bits - used)
        return tail

    def _tail_strategy_balanced(self, remaining_bits: int) -> List[int]:
        """Tail with roughly equal-sized runs"""
        num_runs = min(20, remaining_bits // 10)
        if num_runs == 0:
            return [remaining_bits]
        base_length = remaining_bits // num_runs
        extra = remaining_bits % num_runs
        tail = [base_length] * num_runs
        for i in range(extra):
            tail[i] += 1
        return tail

    def _tail_strategy_fibonacci(self, remaining_bits: int) -> List[int]:
        """Tail with Fibonacci-like spacing"""
        tail = []
        used = 0
        a, b = 1, 2
        while used < remaining_bits:
            if used + a <= remaining_bits:
                tail.append(a)
                used += a
                a, b = b, a + b
            else:
                tail.append(remaining_bits - used)
                break
        return tail

    def _tail_strategy_random(self, remaining_bits: int, max_runs: int) -> List[int]:
        """Random tail strategy"""
        num_runs = random.randint(max(1, max_runs // 2), max_runs)
        tail = []
        used = 0
        for i in range(num_runs - 1):
            max_len = (remaining_bits - used) - (num_runs - i - 1)
            if max_len > 0:
                length = random.randint(1, max_len)
                tail.append(length)
                used += length
        # Last run takes remaining
        tail.append(remaining_bits - used)
        return tail

    def local_search_tail(self, prefix_ones: List[int], initial_tail: List[int],
                          max_iterations: int = 10000,
                          max_seconds: Optional[float] = None) -> Tuple[List[int], int]:
        """Local search to refine the tail structure using iterated local search.
        Operations:
        - Adjust individual run lengths (±1, ±2, ±5)
        - Split a run into two
        - Merge adjacent runs
        - Swap lengths between runs
        When stagnant for stagnation_kick_threshold iterations, applies a random
        perturbation kick (ILS) and resumes rather than stopping early.
        """
        current_tail = initial_tail.copy()
        current_runs = self.construct_from_template(prefix_ones, current_tail)
        current_value = self._evaluate_runs(current_runs)
        print(f"\nStarting local search from {current_value} distinct values")
        print(f"Initial tail: {len(current_tail)} runs, {sum(current_tail)} bits")
        best_value = current_value
        best_tail = current_tail.copy()
        stagnant_iterations = 0
        kick_count = 0
        # Scale kick threshold with n so small n kicks more frequently
        stagnation_kick_threshold = max(2000, self.n * 2)
        kick_strength = max(10, len(initial_tail) // 3)
        start_time = time.time()
        try:
          for iteration in range(max_iterations):
            if max_seconds is not None and (time.time() - start_time) >= max_seconds:
                print(f"  Time cap reached after {time.time() - start_time:.1f}s")
                break
            # ILS kick: when stuck, perturb from best and restart
            if stagnant_iterations >= stagnation_kick_threshold:
                kick_count += 1
                print(f"  Iter {iteration}: kicking (kick #{kick_count}, stagnant={stagnant_iterations})"
                      f" -> restarting from best={best_value}")
                current_tail = best_tail.copy()
                # Apply kick_strength random moves unconditionally
                for _ in range(kick_strength):
                    move_type = random.choice(['adjust', 'split', 'merge', 'swap'])
                    candidate = self._apply_tail_move(current_tail, move_type)
                    if candidate is not None and sum(candidate) == sum(initial_tail):
                        current_tail = candidate
                current_runs = self.construct_from_template(prefix_ones, current_tail)
                current_value = self._evaluate_runs(current_runs)
                stagnant_iterations = 0
                continue
            # Try a random move
            move_type = random.choice(['adjust', 'split', 'merge', 'swap'])
            new_tail = self._apply_tail_move(current_tail, move_type)
            if new_tail is None or sum(new_tail) != sum(initial_tail):
                continue
            if new_tail == current_tail:
                continue
            # Evaluate
            new_runs = self.construct_from_template(prefix_ones, new_tail)
            if sum(new_runs) != self.n:
                continue
            new_value = self._evaluate_runs(new_runs)
            # Accept if better
            if new_value > current_value:
                current_tail = new_tail
                current_value = new_value
                stagnant_iterations = 0
                if new_value > best_value:
                    best_value = new_value
                    best_tail = new_tail.copy()
                    print(f"  Iter {iteration}: NEW BEST = {best_value}")
            else:
                stagnant_iterations += 1
            # Progress update
            if iteration % 1000 == 0 and iteration > 0:
                print(f"  Iter {iteration}: best={best_value}, current={current_value}, "
                      f"stagnant={stagnant_iterations}, kicks={kick_count}")
        except KeyboardInterrupt:
            print(f"\n  [Interrupted at iter {iteration}] Saving best so far: {best_value}")
            return best_tail, best_value
        print(f"\nLocal search complete: {best_value} distinct values ({kick_count} kicks)")
        return best_tail, best_value

    def _apply_tail_move(self, tail: List[int], move_type: str) -> Optional[List[int]]:
        """Apply a random move to the tail"""
        if len(tail) == 0:
            return None
        new_tail = tail.copy()
        changed = False
        if move_type == 'adjust':
            if len(new_tail) < 2:
                return None
            # Adjust one run length
            i = random.randint(0, len(new_tail) - 1)
            delta = random.choice([-5, -2, -1, 1, 2, 5])
            # Find another run to balance
            j = random.randint(0, len(new_tail) - 1)
            while j == i:
                j = random.randint(0, len(new_tail) - 1)
            if new_tail[i] + delta > 0 and new_tail[j] - delta > 0:
                new_tail[i] += delta
                new_tail[j] -= delta
                changed = True
        elif move_type == 'split':
            # Split one run into two
            i = random.randint(0, len(new_tail) - 1)
            if new_tail[i] >= 2:
                split_point = random.randint(1, new_tail[i] - 1)
                new_tail[i:i+1] = [split_point, new_tail[i] - split_point]
                changed = True
        elif move_type == 'merge':
            # Merge two adjacent runs
            if len(new_tail) >= 2:
                i = random.randint(0, len(new_tail) - 2)
                merged = new_tail[i] + new_tail[i+1]
                new_tail[i:i+2] = [merged]
                changed = True
        elif move_type == 'swap':
            # Swap lengths between two runs
            if len(new_tail) >= 2:
                i, j = random.sample(range(len(new_tail)), 2)
                new_tail[i], new_tail[j] = new_tail[j], new_tail[i]
                changed = True
        if not changed or new_tail == tail:
            return None
        return new_tail

    def optimize_prefix_ones(self, initial_ones: List[int], tail: List[int],
                             tolerance: int = 5) -> Tuple[List[int], List[int], int]:
        """Fine-tune the prefix ones blocks (within ±tolerance of initial values)"""
        print(f"\nOptimizing prefix ones blocks (+/-{tolerance} from predictions)...")
        best_ones = initial_ones.copy()
        best_tail = tail.copy()
        runs = self.construct_from_template(best_ones, best_tail)
        best_value = self._evaluate_runs(runs)
        print(f"  Initial: {best_value} distinct values")
        # Try small adjustments to each ones block
        for pos in range(8):
            # Skip fixed positions
            if initial_ones[pos] in [1, 2, 3]:
                continue
            for delta in range(-tolerance, tolerance + 1):
                if delta == 0:
                    continue
                test_ones = best_ones.copy()
                test_ones[pos] += delta
                if test_ones[pos] < 1:
                    continue
                # Adjust tail to maintain total bits
                test_tail = tail.copy()
                if len(test_tail) > 0:
                    test_tail[0] -= delta
                    if test_tail[0] < 1:
                        continue
                test_runs = self.construct_from_template(test_ones, test_tail)
                if sum(test_runs) != self.n:
                    continue
                test_value = self._evaluate_runs(test_runs)
                if test_value > best_value:
                    best_value = test_value
                    best_ones = test_ones.copy()
                    best_tail = test_tail.copy()
                    print(f"    Position {pos*2}: delta={delta:+d} -> {test_value} (NEW BEST)")
        return best_ones, best_tail, best_value

    def run(self, num_restarts: int = 10, local_search_iterations: int = 10000,
            max_seconds_per_restart: Optional[float] = None,
            tail_strategy_mode: str = "full",
            tail_randomize_ratio: float = 0.2,
            tail_randomize_variants: int = 4):
        """Run the template-based greedy optimization"""
        print("=" * 80)
        print(f"Template-Based Greedy Search for n={self.n}")
        print("=" * 80)
        print()
        print("PREFIX TEMPLATE:")
        print(f"  Zeros blocks (FIXED): {self.prefix_zeros}")
        print(f"  Ones blocks (predicted): {self.prefix_ones_base}")
        print(f"  Prefix total: {sum(self.prefix_ones_base) + sum(self.prefix_zeros)} bits")
        print()
        overall_best_value = 0
        overall_best_solution = None
        overall_best_runs = None
        carried_prefix_ones = self.seed_prefix_ones.copy() if self.seed_prefix_ones is not None else None
        carried_tail = self.seed_tail.copy() if self.seed_tail is not None else None
        if self.seed_solution_bits is not None:
            seed_runs = self.bits_to_runs(self.seed_solution_bits)
            seed_value = self._evaluate_runs(seed_runs)
            overall_best_value = seed_value
            overall_best_solution = self.seed_solution_bits.copy()
            overall_best_runs = self.bits_to_runs(overall_best_solution)
            print(f"Seed baseline value: {seed_value}")
            # Persist the seed as the starting point; only update results/ if
            # seed_value actually beats the stored best (handled inside save_result).
            self.save_result(overall_best_solution, overall_best_value)
        try:
          for restart in range(num_restarts):
            restart_start = time.time()
            print(f"\n{'='*80}")
            print(f"RESTART {restart + 1}/{num_restarts}")
            print(f"{'='*80}")
            if max_seconds_per_restart is not None:
                print(f"Time cap for this restart: {max_seconds_per_restart:.1f}s")
            # Prefix selection: anchor on incumbent best when available
            base_prefix = carried_prefix_ones if carried_prefix_ones is not None else self.prefix_ones_base
            # Keep prefix stable once a seeded/incumbent prefix exists, so tail bit-budget stays compatible
            variation = 0 if carried_prefix_ones is not None else 3
            prefix_ones = [
                max(1, x + random.randint(-variation, variation)) if x > 3 else x
                for x in base_prefix
            ]
            print(f"Prefix ones (with variation): {prefix_ones}")
            # Greedy tail construction
            initial_tail = self.greedy_tail_construction(
                prefix_ones,
                self.n,
                seed_tail=carried_tail,
                strategy_mode=tail_strategy_mode,
                randomize_tail_ratio=tail_randomize_ratio,
                randomize_tail_variants=tail_randomize_variants
            )
            elapsed = time.time() - restart_start
            remaining_time = None
            if max_seconds_per_restart is not None:
                remaining_time = max(0.0, max_seconds_per_restart - elapsed)
            # Reserve time for prefix optimisation at the end (runs ~80 evaluations).
            # Local search gets the budget minus this reserve; prefix opt always runs.
            PREFIX_OPT_RESERVE = 120  # seconds
            local_search_budget = None
            if remaining_time is not None:
                local_search_budget = max(0.0, remaining_time - PREFIX_OPT_RESERVE)
            # Local search on tail
            if local_search_budget is not None and local_search_budget <= 0:
                print("Restart time budget exhausted before local search; using initial tail")
                best_tail = initial_tail.copy()
                runs = self.construct_from_template(prefix_ones, best_tail)
                tail_value = self._evaluate_runs(runs)
            else:
                best_tail, tail_value = self.local_search_tail(
                    prefix_ones,
                    initial_tail,
                    max_iterations=local_search_iterations,
                    max_seconds=local_search_budget
                )
            elapsed = time.time() - restart_start
            remaining_time = None
            if max_seconds_per_restart is not None:
                remaining_time = max(0.0, max_seconds_per_restart - elapsed)
            # Fine-tune prefix ones (always runs — local search budget reserved time for this)
            final_ones, final_tail, final_value = self.optimize_prefix_ones(
                prefix_ones, best_tail, tolerance=5
            )
            # Check if this is the best overall
            if final_value > overall_best_value:
                overall_best_value = final_value
                overall_best_runs = self.construct_from_template(final_ones, final_tail)
                overall_best_solution = self.runs_to_bits(overall_best_runs)
                carried_prefix_ones = final_ones.copy()
                carried_tail = final_tail.copy()
                print(f"\n*** NEW OVERALL BEST: {final_value} distinct values ***")
                # Save immediately
                self.save_result(overall_best_solution, overall_best_value)
            else:
                print("No overall improvement this restart; keeping incumbent best as next restart seed")
            # Checkpoint after every restart so long runs are resumable
            if overall_best_solution is not None:
                self.save_result(overall_best_solution, overall_best_value)
            if max_seconds_per_restart is not None:
                print(f"Restart elapsed: {time.time() - restart_start:.1f}s")
        except KeyboardInterrupt:
            print("\n\nInterrupted! Saving current best before exit...")
            if overall_best_solution is not None:
                self.save_result(overall_best_solution, overall_best_value)
                print(f"Saved best={overall_best_value}")
            return overall_best_solution, overall_best_value
        print("\n" + "=" * 80)
        print("OPTIMIZATION COMPLETE")
        print("=" * 80)
        print(f"Best value found: {overall_best_value}")
        if overall_best_solution is not None:
            print(f"Best solution: {''.join(map(str, overall_best_solution))}")
        return overall_best_solution, overall_best_value

    def _save_seed(self, solution: List[int]):
        """Write the best available starting bitstring to seeds/.
        Makes no claim about optimality — just a good starting point."""
        seed_file = self.seeds_dir / f"n_{self.n:04d}_seed.txt"
        seed_file.write_text(''.join(map(str, solution)))

    def save_result(self, solution: List[int], value: int):
        """Save a verified best result.
        Rules:
          1. Always write the bitstring to seeds/ as the best available seed.
          2. Compute the actual score of *this* bitstring (to catch any mismatch).
          3. Only update results/ when the verified score strictly beats the
             current stored best — so results/best_value always matches
             results/solutions[0] and is never downgraded.
        """
        # --- 1. Always persist as seed ---
        self._save_seed(solution)
        # --- 2. Verify value integrity ---
        actual_value = self.compute_distinct_values(solution)
        if actual_value != value:
            print(f"  [integrity] passed value {value} != computed {actual_value}; "
                  f"using computed value")
            value = actual_value
        # --- 3. Only update results if this is a genuine new best ---
        result_file = self.results_dir / f"n_{self.n:04d}_results.json"
        existing_best = 0
        if result_file.exists():
            try:
                existing_best = json.load(open(result_file)).get("best_value", 0)
            except Exception:
                pass
        if value <= existing_best:
            return  # don't overwrite a better verified record
        solution_str = ''.join(map(str, solution))
        data = {
            "n": self.n,
            "best_value": value,
            "solutions": [solution_str],
            "last_updated": datetime.now().isoformat(),
            "metadata": {
                "algorithm": "template_greedy",
                "prefix_zeros": self.prefix_zeros,
                "prefix_ones_base": self.prefix_ones_base
            }
        }
        with open(result_file, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"\nNew verified best saved to {result_file.as_posix()} (value={value})")

def main():
    """Main entry point"""
    # CONFIGURATION
    n = 20000
    num_restarts = 1
    local_search_iterations = 5000000
    max_seconds_per_restart = 18000  # 5 hours per restart
    tail_strategy_mode = "seed_only"  # options: full, decreasing_only, seed_only, best_or_randomized_tail
    tail_randomize_ratio = 0.2
    tail_randomize_variants = 4
    # Optional: explicit seed (set None to auto-load latest best from results/n_3000_results.json)
    seed_solution = None
    # Run optimization
    optimizer = TemplateGreedyOptimizer(n, seed_solution=seed_solution)
    solution, value = optimizer.run(
        num_restarts=num_restarts,
        local_search_iterations=local_search_iterations,
        max_seconds_per_restart=max_seconds_per_restart,
        tail_strategy_mode=tail_strategy_mode,
        tail_randomize_ratio=tail_randomize_ratio,
        tail_randomize_variants=tail_randomize_variants
    )
    print(f"\nFinal result: {value} distinct values for n={n}")
if __name__ == "__main__":
    main()
