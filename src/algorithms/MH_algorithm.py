#!/Users/benclare/.pyenv/versions/3.12.0/bin/python
"""
Metropolis-Hastings Algorithm for OEIS A112509

This module implements a sophisticated Metropolis-Hastings search algorithm to find
n-bit binary strings that maximize the number of distinct integers representable by
their substrings (sequence A112509).

Algorithm Overview:
==================

The algorithm uses run-length encoding to efficiently represent and manipulate binary
strings, employing various optimization techniques:

1. Run-Length Encoding:
   - Binary strings represented as alternating runs of 1s and 0s
   - Example: "1110011" → [3, 2, 2] (3 ones, 2 zeros, 2 ones)
   - Reduces search space and enables efficient structural moves

2. Metropolis-Hastings MCMC:
   - Proposes random modifications to current solution
   - Always accepts improvements (higher distinct count)
   - Sometimes accepts worse solutions with probability exp(Δ/T)
   - Temperature T decreases over time (simulated annealing)

3. Move Types:
   - Transfer: Move bits between adjacent runs
   - Split: Split a run into two parts with separator
   - Merge: Merge two runs by removing separator
   - Swap: Exchange non-adjacent run lengths of same type
   - Multi-transfer: Move bits between non-adjacent runs

4. Enhanced Features:
   - Adaptive temperature with reheating to escape local optima
   - Phase-aware move selection (exploration → exploitation)
   - Periodic hill-climbing for greedy refinement
   - Constraint enforcement on block structure
   - Parallel multi-restart with parameter variation

5. Multi-Restart Strategy:
   - Runs many independent searches in parallel
   - Each restart uses slightly varied parameters
   - Explores different regions of search space
   - Utilizes all CPU cores for speed

Key Parameters:
==============
- n: Length of binary string
- target_density: Proportion of 1s (None = no bias, or 0.0 to 1.0 for biased search)
- T_initial: Initial temperature (controls exploration, default 10.0)
- cooling_rate: Temperature decay (typically 0.9995)
- num_restarts: Number of parallel searches
- iterations_per_run: Search length per restart (default 100,000)

Author: A112509 Research Project
Date: February 2026
"""

import random
import math
import json
import time
import signal
from datetime import datetime
from typing import List, Set, Tuple, Dict, Optional
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import os
from pathlib import Path

import numpy as np

try:
    from numba import njit
    _NUMBA_OK = True
except ImportError:
    _NUMBA_OK = False
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        def _wrap(fn):
            return fn
        return _wrap


# ─────────────────────────────────────────────────────────────────────
# Numba-accelerated suffix-automaton distinct-count kernel
# ─────────────────────────────────────────────────────────────────────

@njit(cache=True)
def _sam_distinct_count(bits):
    """Count distinct integer values representable by substrings of *bits*.

    Pure-integer suffix-automaton build + topological path count, compiled
    to native code by Numba.  The input must be a 1-D ``np.int8`` array.

    Returns the count as an ``int64``.
    """
    n = bits.shape[0]
    if n == 0:
        return 0

    max_states = 2 * n + 2
    sam_len  = np.zeros(max_states, dtype=np.int64)
    sam_link = np.full(max_states, -1, dtype=np.int64)
    sam_t0   = np.full(max_states, -1, dtype=np.int64)
    sam_t1   = np.full(max_states, -1, dtype=np.int64)
    size = 1
    last = 0
    has_zero = False

    for i in range(n):
        c = bits[i]
        if c == 0:
            has_zero = True
        cur = size
        sam_len[cur] = sam_len[last] + 1
        size += 1
        p = last
        while p != -1:
            if c == 0:
                if sam_t0[p] >= 0:
                    break
                sam_t0[p] = cur
            else:
                if sam_t1[p] >= 0:
                    break
                sam_t1[p] = cur
            p = sam_link[p]
        if p == -1:
            sam_link[cur] = 0
        else:
            q = sam_t0[p] if c == 0 else sam_t1[p]
            if sam_len[p] + 1 == sam_len[q]:
                sam_link[cur] = q
            else:
                clone = size
                sam_len[clone] = sam_len[p] + 1
                sam_link[clone] = sam_link[q]
                sam_t0[clone] = sam_t0[q]
                sam_t1[clone] = sam_t1[q]
                size += 1
                while p != -1:
                    if c == 0:
                        if sam_t0[p] == q:
                            sam_t0[p] = clone
                        else:
                            break
                    else:
                        if sam_t1[p] == q:
                            sam_t1[p] = clone
                        else:
                            break
                    p = sam_link[p]
                sam_link[q] = clone
                sam_link[cur] = clone
        last = cur

    # --- topological path count ---
    max_l = 0
    for s in range(size):
        if sam_len[s] > max_l:
            max_l = sam_len[s]

    bucket = np.zeros(max_l + 1, dtype=np.int64)
    for s in range(size):
        bucket[sam_len[s]] += 1
    total = 0
    for i in range(max_l, -1, -1):
        c = bucket[i]
        bucket[i] = total
        total += c
    order = np.zeros(size, dtype=np.int64)
    for s in range(size):
        order[bucket[sam_len[s]]] = s
        bucket[sam_len[s]] += 1

    cnt = np.zeros(size, dtype=np.int64)
    for oi in range(size):
        s = order[oi]
        val = np.int64(0)
        if sam_t0[s] >= 0:
            val += 1 + cnt[sam_t0[s]]
        if sam_t1[s] >= 0:
            val += 1 + cnt[sam_t1[s]]
        cnt[s] = val

    count = np.int64(0)
    if sam_t1[0] >= 0:
        count = 1 + cnt[sam_t1[0]]
    if has_zero:
        count += 1
    return count


@njit(cache=True)
def _bitswap_count_flips(bits, target_value):
    """Try all single-bit flips and density-preserving swaps (Numba kernel).

    Returns a 2-D int8 array where each row is a solution bitstring that
    achieves *target_value*.  If none are found, returns an empty (0, n)
    array.
    """
    n = bits.shape[0]
    max_results = n + n * n
    results = np.empty((max_results, n), dtype=np.int8)
    count = 0

    work = bits.copy()

    # --- 1. Single bit flips ---
    for i in range(n):
        original = work[i]
        work[i] = np.int8(1 - original)
        if work[0] == 1:
            v = _sam_distinct_count(work)
            if v == target_value:
                results[count] = work
                count += 1
        work[i] = original

    # --- 2. Density-preserving bit swaps ---
    n_zeros = 0
    n_ones = 0
    for i in range(n):
        if work[i] == 0:
            n_zeros += 1
        elif i > 0:
            n_ones += 1

    zeros = np.empty(n_zeros, dtype=np.int64)
    ones = np.empty(n_ones, dtype=np.int64)
    zi = 0
    oi = 0
    for i in range(n):
        if work[i] == 0:
            zeros[zi] = i
            zi += 1
        elif i > 0:
            ones[oi] = i
            oi += 1

    for zi in range(n_zeros):
        j0 = zeros[zi]
        for oi in range(n_ones):
            j1 = ones[oi]
            work[j0] = np.int8(1)
            work[j1] = np.int8(0)
            v = _sam_distinct_count(work)
            if v == target_value:
                results[count] = work
                count += 1
            work[j0] = np.int8(0)
            work[j1] = np.int8(1)

    return results[:count]


class RunLengthOptimizer:
    """
    Search for optimal n-bit strings using run-length encoding and Metropolis-Hastings.

    This class implements a sophisticated search algorithm that uses run-length encoding
    to represent binary strings efficiently. The encoding represents consecutive runs of
    1s and 0s, starting with 1s (since valid n-bit numbers must start with 1).

    For example, the binary string "1110011" is encoded as [3, 2, 2]:
    - 3 ones, then 2 zeros, then 2 ones

    Attributes:
        n (int): Length of binary strings to search
        target_density (float | None): Target proportion of 1s (None = no bias)
        target_ones (int): Target number of 1-bits (for initialization only)
        best_value (int): Best count of distinct integers found
        best_solutions (List[Tuple]): List of (runs, bits) tuples for optimal solutions
        results_dir (Path): Directory for saving results
        start_time (float): Timestamp when search started
    """

    def __init__(self, n: int, target_density: Optional[float] = None, results_dir: str = "results"):
        """
        Initialize the optimizer for n-bit strings.

        Args:
            n: Length of binary strings to optimize
            target_density: Target proportion of 1s (None = no bias, free exploration)
            results_dir: Directory to save results (default "results")
        """
        self.n = n
        self.target_density = target_density
        self.target_ones = int(n * target_density) if target_density is not None else n // 2
        self.best_value = 0
        self.best_solutions = []
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        self.start_time = None
        self._interrupted = False  # set True if search was stopped by Ctrl+C

    def runs_to_bits(self, runs: List[int]) -> np.ndarray:
        """
        Convert run-length encoding to bit array (numpy int8).

        The run-length encoding alternates between runs of 1s and 0s, starting with 1s.
        For example, [3, 2, 2] -> array([1, 1, 1, 0, 0, 1, 1], dtype=int8)

        Args:
            runs: List of run lengths, alternating 1s and 0s

        Returns:
            np.ndarray of 0s and 1s (dtype=int8) representing the binary string
        """
        vals = np.zeros(len(runs), dtype=np.int8)
        vals[::2] = 1  # even indices are 1-runs
        return np.repeat(vals, runs)

    def bits_to_runs(self, bits) -> List[int]:
        """
        Convert bit array to run-length encoding.

        Args:
            bits: List of 0s and 1s

        Returns:
            List of run lengths, starting with 1s (adds 0-length prefix if needed)
        """
        if len(bits) == 0:
            return []

        runs = []
        current_bit = int(bits[0])
        current_length = 1

        for bit in bits[1:]:
            if int(bit) == current_bit:
                current_length += 1
            else:
                runs.append(current_length)
                current_bit = int(bit)
                current_length = 1
        runs.append(current_length)

        # Ensure we start with 1s (add 0-length prefix if string starts with 0)
        if current_bit == 0 and len(runs) == 1:
            runs = [0] + runs
        elif int(bits[0]) == 0:
            runs = [0] + runs

        return runs

    def clean_runs(self, runs: List[int]) -> List[int]:
        """
        Remove zero-length runs and merge adjacent runs of the same type.

        This handles edge cases where operations might create invalid run encodings.

        Args:
            runs: Run-length encoding that may contain zeros

        Returns:
            Cleaned run-length encoding
        """
        if not runs:
            return []

        cleaned = []
        i = 0

        while i < len(runs):
            if runs[i] == 0:
                # Skip zero-length run
                # If there are runs before and after, they're now adjacent and same type
                if cleaned and i + 1 < len(runs):
                    # Merge the last cleaned run with the next run
                    cleaned[-1] += runs[i + 1]
                    i += 2  # Skip both the zero and the next run (already merged)
                else:
                    i += 1  # Just skip the zero
            else:
                cleaned.append(runs[i])
                i += 1

        return cleaned

    def compute_distinct_values(self, bits: List[int]) -> Set[int]:
        """
        Compute all distinct integer values representable by substrings.

        This is the objective function we're trying to maximize. Valid substrings are:
        - The single bit '0' (value 0)
        - Any substring starting with '1' (no leading zeros)

        Args:
            bits: Binary string as list of 0s and 1s

        Returns:
            Set of all distinct integer values from valid substrings

        Example:
            bits = [1, 1, 0] represents "110"
            Valid substrings: "1", "11", "110", "1", "10", "0"
            Distinct values: {0, 1, 2, 6} (0, 1, 2=10₂, 6=110₂)
            Returns: {0, 1, 2, 6} - count is 4
        """
        values = set()
        n = len(bits)
        add_value = values.add

        for i in range(n):
            if bits[i] == 0:
                # Substrings starting with 0 are only valid if single 0
                add_value(0)
                continue

            val = 0
            for j in range(i, n):
                val = val * 2 + bits[j]
                add_value(val)

        return values

    def compute_distinct_count(self, bits) -> int:
        """
        Count distinct integer values representable by substrings of `bits`.

        Delegates to a Numba-JIT compiled suffix automaton kernel for speed.
        Accepts either a Python list or a numpy int8 array.
        """
        if isinstance(bits, np.ndarray):
            return int(_sam_distinct_count(bits))
        return int(_sam_distinct_count(np.array(bits, dtype=np.int8)))

    def structured_initialization(self,
                                  min_ones: Optional[List[int]] = None,
                                  max_ones: Optional[List[Optional[int]]] = None,
                                  min_zeros: Optional[List[int]] = None,
                                  max_zeros: Optional[List[Optional[int]]] = None,
                                  min_ones_blocks: Optional[int] = None,
                                  max_ones_blocks: Optional[int] = None) -> List[int]:
        """
        Initialize a run-length encoding with explicit constraints on block structure.

        This creates an initial solution that satisfies specified constraints on the
        first 4 blocks of ones and first 4 blocks of zeros. The run structure alternates:
        [1s, 0s, 1s, 0s, 1s, 0s, 1s, 0s, ...]

        The function distributes bits to meet the specified constraints while trying to
        achieve the target density.

        Args:
            min_ones: [min_1st_ones, min_2nd_ones, min_3rd_ones, min_4th_ones]
            max_ones: [max_1st_ones, max_2nd_ones, max_3rd_ones, max_4th_ones]
            min_zeros: [min_1st_zeros, min_2nd_zeros, min_3rd_zeros, min_4th_zeros]
            max_zeros: [max_1st_zeros, max_2nd_zeros, max_3rd_zeros, max_4th_zeros]
            min_ones_blocks: Minimum total number of 1-runs (None = no constraint)
            max_ones_blocks: Maximum total number of 1-runs (None = no constraint);
                implies max 0-runs = max_ones_blocks (ends in 0) or max_ones_blocks-1 (ends in 1)

        Returns:
            Run-length encoding satisfying the constraints

        Raises:
            ValueError: If constraints are impossible to satisfy for the given n
        """
        # Default bounds if not specified
        if min_ones is None:
            min_ones = [0, 0, 0, 0]
        if max_ones is None:
            max_ones = [None, None, None, None]
        if min_zeros is None:
            min_zeros = [0, 0, 0, 0]
        if max_zeros is None:
            max_zeros = [None, None, None, None]

        # Build initial structure for first 4 blocks
        runs = []
        total_used = 0

        # Add initial ones and zeros blocks (up to 4 pairs, capped by max_ones_blocks)
        initial_blocks = min(4, max_ones_blocks) if max_ones_blocks is not None else 4
        for i in range(initial_blocks):
            # Add i-th block of ones
            should_add_ones = False
            if min_ones and min_ones[i] > 0:
                should_add_ones = True
            elif max_ones:
                max_val_ones = max_ones[i]
                if max_val_ones is not None and max_val_ones > 0:
                    should_add_ones = True
            
            if should_add_ones:
                ones_val = min_ones[i] if min_ones else 0
                if max_ones and max_ones[i] is not None:
                    max_val = max_ones[i]
                    if max_val is not None:
                        ones_val = min(ones_val, max_val)
                runs.append(ones_val)
                total_used += ones_val
            else:
                runs.append(0)  # Will be validated later

            # Add i-th block of zeros (if we're not done yet)
            should_add_zeros = False
            if i < 3:
                should_add_zeros = True
            elif min_zeros and min_zeros[i] > 0:
                should_add_zeros = True
            elif max_zeros:
                max_val_zeros = max_zeros[i]
                if max_val_zeros is not None and max_val_zeros > 0:
                    should_add_zeros = True
            
            if should_add_zeros:
                zeros_val = min_zeros[i] if min_zeros else 0
                if max_zeros and max_zeros[i] is not None:
                    max_val = max_zeros[i]
                    if max_val is not None:
                        zeros_val = min(zeros_val, max_val)
                runs.append(zeros_val)
                total_used += zeros_val

        # Calculate remaining bits to distribute
        remaining = self.n - total_used

        if remaining < 0:
            raise ValueError(f"Constraints exceed n={self.n}: used {total_used}, need {remaining} more")

        # Distribute remaining bits while respecting constraints
        while remaining > 0:
            # Decide whether to add to ones or zeros based on target density
            if self.target_density is not None:
                # Use target density to guide distribution
                current_ones = sum(runs[i] for i in range(0, len(runs), 2))
                current_zeros = sum(runs[i] for i in range(1, len(runs), 2))
                total_bits = current_ones + current_zeros

                # Avoid division by zero - if no bits yet, compare ones count directly
                if total_bits > 0:
                    current_density = current_ones / total_bits
                else:
                    current_density = 0.0

                add_to_ones = current_density < self.target_density
            else:
                # No target density - randomly choose or distribute evenly
                current_ones = sum(runs[i] for i in range(0, len(runs), 2))
                current_zeros = sum(runs[i] for i in range(1, len(runs), 2))
                # Distribute randomly (50/50 chance)
                add_to_ones = random.random() < 0.5

            if add_to_ones:
                # Add to ones - find a block that can accept more
                added = False
                for i in range(0, min(8, len(runs)), 2):  # First 4 ones blocks
                    block_idx = i // 2
                    if max_ones and block_idx < len(max_ones):
                        max_val = max_ones[block_idx]
                        if max_val is None or runs[i] < max_val:
                            runs[i] += 1
                            added = True
                            break
                    else:
                        runs[i] += 1
                        added = True
                        break

                if not added:
                    # All first 4 blocks at max, add to a new block
                    current_ones_blocks = (len(runs) + 1) // 2
                    if max_ones_blocks is None or current_ones_blocks < max_ones_blocks:
                        if len(runs) % 2 == 0:  # Last is zeros
                            runs.append(1)  # New ones block
                        else:
                            runs[-1] += 1  # Add to last ones block
                    else:
                        # At block count limit: grow the last ones block instead
                        last_ones_idx = len(runs) - 1 if len(runs) % 2 == 1 else len(runs) - 2
                        if last_ones_idx >= 0:
                            runs[last_ones_idx] += 1
            else:
                # Add to zeros - find a block that can accept more
                added = False
                for i in range(1, min(8, len(runs)), 2):  # First 4 zeros blocks
                    block_idx = i // 2
                    if max_zeros and block_idx < len(max_zeros):
                        max_val = max_zeros[block_idx]
                        if max_val is None or runs[i] < max_val:
                            runs[i] += 1
                            added = True
                            break
                    else:
                        runs[i] += 1
                        added = True
                        break

                if not added:
                    # All first 4 blocks at max, add to a new block
                    if len(runs) % 2 == 1:  # Last is ones
                        runs.append(1)  # New zeros block
                    else:
                        runs[-1] += 1  # Add to last zeros block

            remaining -= 1

        # Final validation - ensure first run is non-zero
        if len(runs) > 0 and runs[0] == 0:
            # Move bits from elsewhere
            if len(runs) > 2 and runs[2] > min_ones[1]:
                transfer = min(runs[2] - min_ones[1], 1)
                runs[0] += transfer
                runs[2] -= transfer
            else:
                runs[0] = 1
                if len(runs) > 1:
                    runs[1] = max(0, runs[1] - 1)

        # Enforce minimum ones-block count by splitting large blocks if needed
        if min_ones_blocks is not None:
            current_ones_blocks = (len(runs) + 1) // 2
            while current_ones_blocks < min_ones_blocks:
                # Find the largest ones block that can be split (needs > 2 bits: part + sep + part)
                best_idx = -1
                best_size = 2
                for k in range(0, len(runs), 2):
                    if runs[k] > best_size:
                        best_idx = k
                        best_size = runs[k]
                if best_idx == -1:
                    break  # Cannot split further
                half = runs[best_idx] // 2
                tail = runs[best_idx] - half - 1
                if tail < 1:
                    break
                runs[best_idx:best_idx + 1] = [half, 1, tail]
                current_ones_blocks = (len(runs) + 1) // 2

        return runs

    def propose_move(self, runs: List[int], iteration: int = 0, max_iterations: int = 10000,
                    min_ones: Optional[List[int]] = None, max_ones: Optional[List[Optional[int]]] = None,
                    min_zeros: Optional[List[int]] = None, max_zeros: Optional[List[Optional[int]]] = None,
                    min_density: Optional[float] = None, max_density: Optional[float] = None,
                    min_ones_blocks: Optional[int] = None, max_ones_blocks: Optional[int] = None) -> List[int]:
        """
        Propose a modification to the run-length encoding.

        This is the key function for the Metropolis-Hastings algorithm. It proposes
        random changes to the current solution while maintaining validity and constraints.

        Move Types:
        -----------
        1. Transfer: Move 1-3 bits from one run to an adjacent run
        2. Split: Split a run into two parts with a separator
        3. Merge: Merge two runs by removing their separator
        4. Swap: Exchange the lengths of two non-adjacent runs of the same type
        5. Multi-transfer: Move bits between non-adjacent runs

        The selection probability for each move type adapts based on search phase:
        - Early phase (0-30%): Aggressive structural changes (more splits/merges)
        - Mid phase (30-70%): Balanced exploration
        - Late phase (70-100%): Fine-tuning (more transfers, less structural change)

        Args:
            runs: Current run-length encoding
            iteration: Current iteration number (for phase detection)
            max_iterations: Total iterations (for phase detection)
            min_ones: Minimum sizes for first 4 ones blocks
            max_ones: Maximum sizes for first 4 ones blocks
            min_zeros: Minimum sizes for first 4 zeros blocks
            max_zeros: Maximum sizes for first 4 zeros blocks
            min_density: Minimum proportion of 1s
            max_density: Maximum proportion of 1s
            min_ones_blocks: Minimum total number of 1-runs (None = no constraint)
            max_ones_blocks: Maximum total number of 1-runs (None = no constraint)

        Returns:
            New run-length encoding (or original if move invalid)
        """
        # Fast path: no constraints active — skip all normalization and
        # the validation loop at the bottom (the common unconstrained case).
        _no_constraints = (
            min_ones is None and max_ones is None and
            min_zeros is None and max_zeros is None and
            min_density is None and max_density is None and
            min_ones_blocks is None and max_ones_blocks is None
        )

        if not _no_constraints:
            if min_ones is None:
                min_ones = [0, 0, 0, 0]
            if max_ones is None:
                max_ones = [None, None, None, None]
            if min_zeros is None:
                min_zeros = [0, 0, 0, 0]
            if max_zeros is None:
                max_zeros = [None, None, None, None]

        run_count = len(runs)
        if run_count == 0:
            return runs

        new_runs = runs.copy()
        changed = False

        # Weight moves differently based on search phase
        early_phase = iteration < max_iterations * 0.3
        mid_phase = max_iterations * 0.3 <= iteration < max_iterations * 0.7

        if early_phase:
            # transfer, split, merge, swap, multi_transfer cumulative thresholds
            thresholds = (0.25, 0.55, 0.75, 0.90)
        elif mid_phase:
            thresholds = (0.30, 0.55, 0.75, 0.90)
        else:
            thresholds = (0.40, 0.60, 0.75, 0.90)

        move_roll = random.random()

        if move_roll < thresholds[0] and run_count >= 2:
            # Transfer 1-3 units from one run to an adjacent run
            idx = random.randint(0, run_count - 2)
            transfer_amount = min(new_runs[idx] - 1, random.randint(1, 3))
            if transfer_amount > 0:
                new_runs[idx] -= transfer_amount
                new_runs[idx + 1] += transfer_amount
                changed = True

        elif move_roll < thresholds[1] and run_count >= 1:
            # Split a run into two with a separator
            idx = random.randint(0, run_count - 1)
            if new_runs[idx] >= 3:  # Need at least 3 to split meaningfully
                split_point = random.randint(1, new_runs[idx] - 1)
                new_runs[idx:idx+1] = [split_point, 1, new_runs[idx] - split_point - 1]
                changed = True

        elif move_roll < thresholds[2] and run_count >= 3:
            # Merge two runs by removing separator
            idx = random.randint(0, run_count - 3)
            if idx % 2 == 0:  # Merging two 1-runs
                new_runs[idx:idx+3] = [new_runs[idx] + new_runs[idx+1] + new_runs[idx+2]]
                changed = True

        elif move_roll < thresholds[3] and run_count >= 4:
            # Swap two non-adjacent runs (of same type)
            idx1 = random.randint(0, run_count - 1)
            idx2 = random.randint(0, run_count - 1)
            # Ensure same type (both even or both odd indices)
            if idx1 % 2 == idx2 % 2 and idx1 != idx2:
                new_runs[idx1], new_runs[idx2] = new_runs[idx2], new_runs[idx1]
                changed = True

        elif run_count >= 3:
            # Transfer between non-adjacent runs
            idx1 = random.randint(0, run_count - 1)
            idx2 = random.randint(0, run_count - 1)
            if idx1 != idx2 and new_runs[idx1] > 1:
                transfer = min(new_runs[idx1] - 1, random.randint(1, 2))
                new_runs[idx1] -= transfer
                new_runs[idx2] += transfer
                changed = True

        # ~5% chance: double-pair cross-swap.  Simultaneously move 1 bit from run A
        # to run B and 1 bit from run C to run D (all same parity), keeping total
        # length constant.  This covers the symmetric class of moves where two
        # independent adjacent-pair transfers are needed at the same time — a key
        # move for crossing connectivity gaps in the optimal-solution graph.
        if not changed and run_count >= 5 and random.random() < 0.05:
            same_parity = [i for i in range(run_count) if i % 2 == 0]
            if len(same_parity) >= 4:
                random.shuffle(same_parity)
                a, b, c, d = same_parity[:4]
                if new_runs[a] > 1 and new_runs[c] > 1:
                    new_runs[a] -= 1
                    new_runs[b] += 1
                    new_runs[c] -= 1
                    new_runs[d] += 1
                    changed = True

        if not changed:
            return runs  # Invalid move, return original

        # Ensure first run is 1s (non-zero)
        if len(new_runs) > 0 and new_runs[0] == 0:
            return runs

        # No constraints — return immediately without any validation loop.
        if _no_constraints:
            return new_runs

        # Enforce constraints on first 4 blocks
        # Must count dynamically because splits/merges change indices
        ones_block_count = 0
        zeros_block_count = 0
        for idx, run_length in enumerate(new_runs):
            if idx % 2 == 0:  # Ones block (even indices)
                if ones_block_count < 4:
                    # Check minimum constraint
                    if min_ones is not None and run_length < min_ones[ones_block_count]:
                        return runs
                    # Check maximum constraint
                    if max_ones is not None:
                        max_val = max_ones[ones_block_count]
                        if max_val is not None and run_length > max_val:
                            return runs
                ones_block_count += 1
            else:  # Zeros block (odd indices)
                if zeros_block_count < 4:
                    # Check minimum constraint
                    if min_zeros is not None and run_length < min_zeros[zeros_block_count]:
                        return runs
                    # Check maximum constraint
                    if max_zeros is not None:
                        max_val = max_zeros[zeros_block_count]
                        if max_val is not None and run_length > max_val:
                            return runs
                zeros_block_count += 1

            # Early stop if no density/block-count checks and first 4 block constraints done
            if (min_density is None and max_density is None and
                    min_ones_blocks is None and max_ones_blocks is None and
                    ones_block_count >= 4 and zeros_block_count >= 4):
                break

        # Enforce density constraints if specified
        if min_density is not None or max_density is not None:
            ones_count = sum(new_runs[i] for i in range(0, len(new_runs), 2))
            current_density = ones_count / self.n

            if min_density is not None and current_density < min_density:
                return runs
            if max_density is not None and current_density > max_density:
                return runs

        # Enforce total ones-block count constraints
        if min_ones_blocks is not None or max_ones_blocks is not None:
            num_ones_blocks = (len(new_runs) + 1) // 2
            if min_ones_blocks is not None and num_ones_blocks < min_ones_blocks:
                return runs
            if max_ones_blocks is not None and num_ones_blocks > max_ones_blocks:
                return runs

        return new_runs

    def metropolis_hastings(self, max_iterations: int = 10000,
                           T_initial: float = 10.0,
                           cooling_rate: float = 0.9995,
                           min_ones: Optional[List[int]] = None,
                           max_ones: Optional[List[Optional[int]]] = None,
                           min_zeros: Optional[List[int]] = None,
                           max_zeros: Optional[List[Optional[int]]] = None,
                           min_density: Optional[float] = None,
                           max_density: Optional[float] = None,
                           use_adaptive_temp: bool = True,
                           min_ones_blocks: Optional[int] = None,
                           max_ones_blocks: Optional[int] = None) -> Tuple[int, List[int], List[int]]:
        """
        Run Metropolis-Hastings search with adaptive temperature and hill-climbing.

        The Metropolis-Hastings Algorithm:
        ----------------------------------
        1. Start with an initial solution
        2. At each iteration:
           a. Propose a random modification
           b. Compute the change in objective (Δ = new_value - current_value)
           c. Accept the proposal if:
              - Δ > 0 (improvement), OR
              - random() < exp(Δ/T) (probabilistic acceptance of worse solution)
        3. Decrease temperature T over time (reduces acceptance of worse solutions)

        Enhancements:
        ------------
        - Adaptive Temperature: Reheat when stuck to escape local optima
        - Phase-Aware Moves: Different move probabilities in early/mid/late phases
        - Hill-Climbing: Periodic greedy search phases for exploitation
        - Constraint Enforcement: All moves respect structural constraints

        Args:
            max_iterations: Number of iterations to run
            T_initial: Initial temperature (higher = more exploration)
            cooling_rate: Temperature multiplier each iteration (typically 0.995-0.9995)
            min_ones: Minimum sizes for first 4 ones blocks
            max_ones: Maximum sizes for first 4 ones blocks
            min_zeros: Minimum sizes for first 4 zeros blocks
            max_zeros: Maximum sizes for first 4 zeros blocks
            min_density: Minimum proportion of 1s
            max_density: Maximum proportion of 1s
            use_adaptive_temp: Enable adaptive temperature with reheating
            min_ones_blocks: Minimum total number of 1-runs (None = no constraint)
            max_ones_blocks: Maximum total number of 1-runs (None = no constraint)

        Returns:
            Tuple of (max_distinct_count, best_runs, best_bits)
        """
        # Initialize
        runs = self.structured_initialization(min_ones, max_ones, min_zeros, max_zeros,
                                              min_ones_blocks, max_ones_blocks)
        bits = self.runs_to_bits(runs)
        current_count = self.compute_distinct_count(bits)

        best_count = current_count
        best_runs = runs.copy()
        best_bits = bits.copy()

        T = T_initial
        no_improvement_count = 0
        last_best = best_count

        for iteration in range(max_iterations):
            # Propose move with phase awareness and constraint enforcement
            new_runs = self.propose_move(runs, iteration, max_iterations,
                                        min_ones, max_ones, min_zeros, max_zeros,
                                        min_density, max_density,
                                        min_ones_blocks, max_ones_blocks)

            # Fast path: invalid/no-op proposal returns unchanged state
            if new_runs is runs:
                if use_adaptive_temp:
                    if best_count > last_best:
                        T *= cooling_rate
                        no_improvement_count = 0
                        last_best = best_count
                    else:
                        no_improvement_count += 1
                        T *= cooling_rate
                        if no_improvement_count > 0 and no_improvement_count % 2000 == 0:
                            T = min(T * 5.0, T_initial * 0.5)
                else:
                    T *= cooling_rate
                continue

            new_bits = self.runs_to_bits(new_runs)
            new_count = self.compute_distinct_count(new_bits)

            # Metropolis-Hastings acceptance criterion
            delta = new_count - current_count

            if delta > 0 or (T > 0 and random.random() < math.exp(delta / T)):
                runs = new_runs
                bits = new_bits
                current_count = new_count

                if current_count > best_count:
                    best_count = current_count
                    best_runs = runs.copy()
                    best_bits = bits.copy()
                    no_improvement_count = 0

            # Adaptive temperature with reheating
            if use_adaptive_temp:
                if best_count > last_best:
                    # Found improvement, cool normally
                    T *= cooling_rate
                    no_improvement_count = 0
                    last_best = best_count
                else:
                    no_improvement_count += 1
                    T *= cooling_rate

                    # Reheat if stuck for too long (escape local optimum)
                    if no_improvement_count > 0 and no_improvement_count % 2000 == 0:
                        T = min(T * 5.0, T_initial * 0.5)  # Reheat but not to initial
            else:
                T *= cooling_rate

            # Periodic hill-climbing phase (pure greedy search for exploitation)
            if iteration > 0 and iteration % 5000 == 0 and iteration < max_iterations - 1000:
                # Do 100 iterations of pure hill climbing
                for _ in range(100):
                    hc_runs = self.propose_move(runs, iteration, max_iterations,
                                               min_ones, max_ones, min_zeros, max_zeros,
                                               min_density, max_density,
                                               min_ones_blocks, max_ones_blocks)
                    if hc_runs is runs:
                        continue
                    hc_bits = self.runs_to_bits(hc_runs)
                    hc_count = self.compute_distinct_count(hc_bits)
                    if hc_count > current_count:
                        runs, bits, current_count = hc_runs, hc_bits, hc_count
                        if current_count > best_count:
                            best_count = current_count
                            best_runs = runs.copy()
                            best_bits = bits.copy()

        return best_count, best_runs, best_bits.tolist()

    # ------------------------------------------------------------------
    # Level-set sampler: find all / many optimal solutions
    # ------------------------------------------------------------------

    def _propose_long_range_swap(self, bits: np.ndarray) -> Optional[np.ndarray]:
        """Swap a random 0-bit with a random 1-bit anywhere in the string.

        Inspired by ``surgical_nudge.py`` Phase 1b (targeted relocation).
        Preserves total density (ones count).  Position 0 is excluded as a
        swap candidate for the 1-bit in order to maintain the leading-1
        invariant (the string must start with 1 to represent valid binary
        numbers without degenerate leading zeros).

        Returns the modified bits array, or ``None`` if no valid swap exists.
        """
        zeros = np.where(bits == 0)[0]
        ones  = np.where(bits[1:] == 1)[0] + 1  # exclude pos 0
        if len(zeros) == 0 or len(ones) == 0:
            return None
        j0 = zeros[random.randint(0, len(zeros) - 1)]
        j1 = ones[random.randint(0, len(ones) - 1)]
        new_bits = bits.copy()
        new_bits[j0], new_bits[j1] = new_bits[j1], new_bits[j0]
        return new_bits

    def level_set_walk(self,
                       target_value: int,
                       seed_bits: List[int],
                       steps: int = 500_000,
                       epsilon: int = 2,
                       ) -> Set[str]:
        """
        Random walk restricted to near ``target_value`` with basin-hopping.

        Starting from ``seed_bits`` (which must already achieve ``target_value``),
        propose moves and accept states that are within ``epsilon`` of the target.
        Every distinct state that achieves exactly ``target_value`` is recorded.

        The ``epsilon`` tolerance allows the walk to briefly cross connectivity
        barriers — low-value "valleys" that separate distinct connected components
        of the optimal-solution graph.  Without this, solutions in disconnected
        components can never be found from a single seed.

        Args:
            target_value:  The optimal value a(n) to maximise.
            seed_bits:     A starting solution achieving ``target_value``.
            steps:         Number of proposals to make.
            epsilon:       Maximum allowed drop below ``target_value``.
                           epsilon=0 recovers the original strict level-set walk.
                           epsilon=2 is typically sufficient for n up to ~200.

        Returns:
            Set of bit-strings (as ``str``) that achieve ``target_value``.
        """
        found: Set[str] = set()
        runs = self.bits_to_runs(seed_bits)
        bits = np.array(seed_bits, dtype=np.int8) if not isinstance(seed_bits, np.ndarray) else seed_bits.copy()
        current_count = target_value
        found.add(''.join(str(int(b)) for b in bits))
        floor = target_value - epsilon

        for _ in range(steps):
            # 25% chance: long-range bit-swap (surgical_nudge.py Phase 1b).
            # Covers connectivity gaps that RLE moves cannot bridge — e.g. when
            # the optimal solution differs from the current state at exactly two
            # positions (one 0↔1 swap) that are far apart in the run-length
            # representation.  25% (raised from 15%) improves coverage of
            # disconnected solution components seen at n≥130.
            if random.random() < 0.25:
                new_bits = self._propose_long_range_swap(bits)
                if new_bits is None:
                    continue
                new_runs = self.bits_to_runs(new_bits)
                if not new_runs or new_runs[0] == 0:
                    continue
            else:
                new_runs = self.propose_move(runs)
                if new_runs is runs:
                    continue
                new_bits = self.runs_to_bits(new_runs)
            new_count = self.compute_distinct_count(new_bits)
            if new_count >= floor:
                runs = new_runs
                bits = new_bits
                current_count = new_count
                if current_count == target_value:
                    found.add(''.join(str(int(b)) for b in bits))

        return found

    def find_all_optimal(self,
                         num_chains: Optional[int] = None,
                         steps_per_chain: int = 500_000,
                         epsilon: int = 0,
                         epsilon_max: int = 3,
                         max_rounds: int = 10,
                         bitswap_expand: bool = True,
                         ) -> List[str]:
        """
        Run ``num_chains`` independent level-set walks and aggregate results.

        Call this after multi_restart_search has run — it reads the best value
        and seed solutions directly from the saved results file, so no prior
        knowledge of the optimal value is needed.

        Uses an **iterative expansion** strategy: each round runs ``num_chains``
        parallel walks seeded from all known solutions found so far (spreading
        seeds evenly).  Any newly discovered solutions are added to the seed
        pool and a new round begins.  Rounds stop when no new solutions are
        found or ``max_rounds`` is reached.

        After each round of random walks, a deterministic **exhaustive bit-swap
        expansion** step is performed: for every currently-known solution, all
        density-preserving single bit-swaps are tried, and any neighbours that
        achieve the target value are added to the found set.  This guarantees
        that no solution which is exactly one density-preserving bit-swap away
        from any found solution can be missed.  The expansion is parallelised
        across CPU cores; it is O(|solutions| × n²) evaluations per round but
        is fast enough for n ≤ ~200.  For larger n, set ``bitswap_expand=False``
        to skip it.

        This ensures that when the optimal set is disconnected (i.e. some
        solutions can only be reached from others that weren't known at the
        start), they are still discovered as long as any single walk can cross
        the gap — and new seeds are distributed to all cores in subsequent rounds.

        Parallel via ``multiprocessing.Pool`` using all available cores.

        Args:
            num_chains:       Chains per round (default: ``cpu_count() × 4``).
            steps_per_chain:  Proposals per chain.
            epsilon:          Maximum allowed drop below the optimal value during
                              the walk (basin-hopping tolerance).  epsilon=0
                              (default) gives the strict level-set walk, which
                              works well because solutions are directly connected
                              by valid moves.  Increase only if you believe the
                              optimal set has deep barriers.
            max_rounds:       Maximum number of iterative expansion rounds.
            bitswap_expand:   If True (default), run exhaustive bit-swap
                              neighbourhood expansion after each walk round.
                              Guaranteed to find solutions one bit-swap away
                              from any found solution.  May be slow for n > 200.
            epsilon_max:      Upper bound for automatic epsilon escalation.
                              When a round finds no new solutions the walk epsilon
                              is incremented by 1 (up to ``epsilon_max``) before
                              trying again, enabling the walk to cross barriers
                              between disconnected solution components.  Set to 0
                              to disable escalation.

        Returns:
            Sorted list of distinct bit-strings achieving the best found value.
        """
        if num_chains is None:
            num_chains = cpu_count() * 4

        existing = self.load_existing_results()
        if existing is None or not existing.get("solutions"):
            raise ValueError(
                "No existing solutions found. Run multi_restart_search first."
            )

        # Read the best value from the file — no prior knowledge needed.
        target_value = existing["best_value"]

        seeds_raw = existing["solutions"]
        all_found: Set[str] = set()
        for s in seeds_raw:
            all_found.add(s["bits"] if isinstance(s, dict) else s)

        seed_list = [
            [int(b) for b in s] for s in sorted(all_found)
        ]

        print(f"\nLevel-set walk: a({self.n}) = {target_value}  (epsilon={epsilon}, epsilon_max={epsilon_max})")
        print(f"Seeds: {len(seed_list)} known solution(s)")
        print(f"Chains per round: {num_chains}  ×  {steps_per_chain:,} steps")
        print(f"Max rounds: {max_rounds}  |  Using {cpu_count()} CPU cores\n")

        t0 = time.time()
        current_epsilon = epsilon

        for round_num in range(1, max_rounds + 1):
            # Distribute chains evenly across all currently known seeds
            worker_args = [
                (self.n, target_value, seed_list[i % len(seed_list)], steps_per_chain, i, current_epsilon)
                for i in range(num_chains)
            ]

            round_found: Set[str] = set()
            lsw_pool = Pool(processes=cpu_count(), initializer=_worker_init)
            try:
                for result in lsw_pool.imap_unordered(_level_set_worker, worker_args):
                    round_found.update(result)
                lsw_pool.close()
                lsw_pool.join()
            except KeyboardInterrupt:
                print("\nKeyboard interrupt — saving results found so far...")
                lsw_pool.terminate()
                self._interrupted = True
                all_found.update(round_found)
                break

            new_this_round = round_found - all_found
            all_found.update(round_found)

            elapsed = time.time() - t0
            print(f"Round {round_num} (ε={current_epsilon}): {len(round_found)} found this round, "
                  f"{len(new_this_round)} new, "
                  f"total {len(all_found)}, elapsed {elapsed:.1f}s")

            # ---- Exhaustive bit-swap neighbourhood expansion ----
            # For every currently-known optimal solution, try every possible
            # density-preserving single bit-swap (swap a 0-bit with a 1-bit at
            # any positions, excluding position 0).  Add any neighbours that
            # achieve target_value.  This deterministically closes the gap for
            # solutions that are exactly 1 swap away from any known solution but
            # unreachable via the RLE random walk.
            if bitswap_expand:
                expand_args = [(self.n, target_value, s) for s in sorted(all_found)]
                expand_pool = Pool(processes=cpu_count(), initializer=_worker_init)
                bitswap_new: Set[str] = set()
                try:
                    for new_strs in expand_pool.imap_unordered(_bitswap_worker, expand_args):
                        bitswap_new.update(new_strs - all_found)
                    expand_pool.close()
                    expand_pool.join()
                except KeyboardInterrupt:
                    expand_pool.terminate()
                    self._interrupted = True
                if bitswap_new:
                    all_found.update(bitswap_new)
                    new_this_round.update(bitswap_new)
                    elapsed = time.time() - t0
                    print(f"  Bit-swap expansion: +{len(bitswap_new)} new, "
                          f"total {len(all_found)}, elapsed {elapsed:.1f}s")

            if not new_this_round:
                if current_epsilon < epsilon_max:
                    current_epsilon += 1
                    elapsed = time.time() - t0
                    print(f"  No new solutions — escalating epsilon to {current_epsilon} "
                          f"(elapsed {elapsed:.1f}s)")
                    # Continue to next round with higher epsilon; don't break
                else:
                    print("No new solutions found — converged.")
                    break
            else:
                # Progress made: step epsilon back toward base (avoid over-exploration)
                if current_epsilon > epsilon:
                    current_epsilon = max(epsilon, current_epsilon - 1)

            # Rebuild seed list from ALL found solutions so new ones are explored next round
            seed_list = [[int(b) for b in s] for s in sorted(all_found)]

        result_list = sorted(all_found)
        elapsed = time.time() - t0
        print(f"\nTotal: {len(result_list)} distinct optimal strings in {elapsed:.1f}s")

        # Merge new solutions into the saved results file.
        current_bits_set = {
            (s["bits"] if isinstance(s, dict) else s) for s in seeds_raw
        }
        # Populate self.best_solutions with ALL existing seeds first so that
        # save_results() writes them out alongside any newly found solutions.
        self.best_value = target_value
        for s in seeds_raw:
            bits_str = s["bits"] if isinstance(s, dict) else s
            bits = [int(b) for b in bits_str]
            runs = self.bits_to_runs(bits)
            self.best_solutions.append((runs, bits))
        new_count = 0
        for bits_str in result_list:
            if bits_str not in current_bits_set:
                bits = [int(b) for b in bits_str]
                runs = self.bits_to_runs(bits)
                self.best_solutions.append((runs, bits))
                current_bits_set.add(bits_str)
                new_count += 1

        if new_count:
            self.best_value = target_value
            self.save_results({"source": "level_set_walk",
                               "steps_per_chain": steps_per_chain,
                               "num_chains": num_chains,
                               "epsilon": epsilon,
                               "max_rounds": max_rounds,
                               "bitswap_expand": bitswap_expand})
            print(f"Added {new_count} new solution(s) to results file.")
        else:
            print("No new solutions beyond existing ones.")

        return result_list

    def bitswap_expand(self, max_rounds: int = 20) -> int:
        """Iterative bit-swap expansion until convergence.

        For every solution in ``self.best_solutions``, tries all single-bit
        flips and all density-preserving bit swaps.  Any new solution at
        ``self.best_value`` is added to the pool and the process repeats
        until no new solutions are found (fixed-point).

        Runs in parallel across all CPU cores.  Typically converges in
        2-3 rounds and costs O(rounds × |solutions| × n²) evaluations.

        Args:
            max_rounds: Safety cap on expansion rounds (default 20;
                        convergence usually happens in 2-3).

        Returns:
            Number of new solutions discovered.
        """
        if not self.best_solutions or self.best_value <= 0:
            return 0

        target_value = self.best_value
        all_found: set = set()
        for _, bits in self.best_solutions:
            all_found.add(''.join(map(str, bits)))
        initial_count = len(all_found)

        ncpu = cpu_count()
        print(f"\nBit-swap expansion: a({self.n})={target_value}, "
              f"{initial_count} seed(s), up to {max_rounds} rounds")

        for round_num in range(1, max_rounds + 1):
            before = len(all_found)
            tasks = [(self.n, target_value, s) for s in sorted(all_found)]
            pool = Pool(ncpu, initializer=_worker_init)
            try:
                for result_set in pool.imap_unordered(_bitswap_worker, tasks):
                    all_found.update(result_set)
                pool.close()
                pool.join()
            except Exception:
                pool.terminate()
            gained = len(all_found) - before
            print(f"  Round {round_num}: {before} -> {len(all_found)} (+{gained})")
            if gained == 0:
                break

        # Add any new solutions to self.best_solutions
        existing = {''.join(map(str, b)) for _, b in self.best_solutions}
        new_count = 0
        for bits_str in sorted(all_found - existing):
            bits_list = [int(c) for c in bits_str]
            self.best_solutions.append((self.bits_to_runs(bits_list), bits_list))
            new_count += 1

        total_new = len(all_found) - initial_count
        print(f"  Bit-swap expansion done: +{total_new} new solution(s)")
        return total_new

    def load_existing_results(self) -> Optional[Dict]:
        """
        Load existing results for this n value if they exist.

        Returns:
            Dictionary containing previous results or None if no file exists
        """
        results_file = self.results_dir / f"n_{self.n:04d}_results.json"
        if results_file.exists():
            with open(results_file, 'r') as f:
                return json.load(f)
        return None

    def save_results(self, metadata: Optional[Dict] = None):
        """
        Save current best results to JSON file.

        Args:
            metadata: Optional dictionary of metadata to include in the file
        """
        results_file = self.results_dir / f"n_{self.n:04d}_results.json"

        # Convert solutions to serializable format (just the bit strings).
        # Deduplicate while preserving order (a set() then sort is fine since
        # bit strings have a natural lexicographic order and callers don't
        # rely on insertion order).
        seen: set = set()
        solutions_data = []
        for _, bits in self.best_solutions:
            s = ''.join(map(str, bits))
            if s not in seen:
                seen.add(s)
                solutions_data.append(s)

        # Add num_solutions to metadata
        if metadata is None:
            metadata = {}
        metadata['num_solutions'] = len(solutions_data)

        data = {
            "n": self.n,
            "best_value": self.best_value,
            "solutions": solutions_data,
            "last_updated": datetime.now().isoformat(),
            "metadata": metadata
        }

        with open(results_file, 'w') as f:
            json.dump(data, f, indent=2)

    def merge_with_existing(self):
        """
        Load and merge with existing results from disk.

        This allows incremental improvement: new searches can build on previous results.
        - If new results are better, they replace the old ones
        - If they're equal, solutions are merged (avoiding duplicates)
        - If old results are better, they're kept
        """
        existing = self.load_existing_results()
        if existing is None:
            return

        existing_value = existing["best_value"]

        if existing_value > self.best_value:
            # Existing results are better, use them
            self.best_value = existing_value
            # Handle both old format (dict with runs/bits) and new format (just bit strings)
            self.best_solutions = []
            for sol in existing["solutions"]:
                if isinstance(sol, dict):
                    # Old format
                    bits = [int(b) for b in sol["bits"]]
                    runs = sol.get("runs", self.bits_to_runs(bits))
                else:
                    # New format (just bit string)
                    bits = [int(b) for b in sol]
                    runs = self.bits_to_runs(bits)
                self.best_solutions.append((runs, bits))
            print(f"Loaded existing results: a({self.n}) = {self.best_value} with {len(self.best_solutions)} solution(s)")
        elif existing_value == self.best_value:
            # Merge solutions, avoiding duplicates
            # Handle both formats
            existing_bits_set = set()
            for sol in existing["solutions"]:
                if isinstance(sol, dict):
                    existing_bits_set.add(sol["bits"])
                else:
                    existing_bits_set.add(sol)

            current_bits_set = {''.join(map(str, bits)) for _, bits in self.best_solutions}

            # Add existing solutions that we don't have
            for sol in existing["solutions"]:
                if isinstance(sol, dict):
                    bits_str = sol["bits"]
                    bits = [int(b) for b in bits_str]
                    runs = sol.get("runs", self.bits_to_runs(bits))
                else:
                    bits_str = sol
                    bits = [int(b) for b in bits_str]
                    runs = self.bits_to_runs(bits)

                if bits_str not in current_bits_set:
                    self.best_solutions.append((runs, bits))

            print(f"Merged with existing results: a({self.n}) = {self.best_value} with {len(self.best_solutions)} solution(s)")
        else:
            # New results are better, will replace existing
            print(f"New best value found! Previous: {existing_value}, New: {self.best_value}")

    def multi_restart_search(self, num_restarts: int = 20,
                            iterations_per_run: int = 50000,
                            min_ones: Optional[List[int]] = None,
                            max_ones: Optional[List[Optional[int]]] = None,
                            min_zeros: Optional[List[int]] = None,
                            max_zeros: Optional[List[Optional[int]]] = None,
                            min_density: Optional[float] = None,
                            max_density: Optional[float] = None,
                            num_processes: Optional[int] = None,
                            min_ones_blocks: Optional[int] = None,
                            max_ones_blocks: Optional[int] = None):
        """
        Run multiple independent searches in parallel.

        Multi-Restart Strategy:
        -----------------------
        Running multiple independent searches helps:
        1. Explore different regions of the search space
        2. Overcome the stochastic nature of the algorithm
        3. Utilize multiple CPU cores for speed

        Each restart uses slightly varied parameters (density, temperature, constraints)
        to encourage diversity in the search.

        Args:
            num_restarts: Number of independent searches to run
            iterations_per_run: Iterations per search
            min_ones: [min_1st, min_2nd, min_3rd, min_4th] ones blocks
            max_ones: [max_1st, max_2nd, max_3rd, max_4th] ones blocks
            min_zeros: [min_1st, min_2nd, min_3rd, min_4th] zeros blocks
            max_zeros: [max_1st, max_2nd, max_3rd, max_4th] zeros blocks
            min_density: Minimum density of ones (0.0 to 1.0, None = no constraint)
            max_density: Maximum density of ones (0.0 to 1.0, None = no constraint)
            num_processes: Number of parallel processes (None = use all cores)

        Returns:
            Tuple of (best_value, list of solutions)
        """
        if num_processes is None:
            num_processes = cpu_count()

        self.start_time = time.time()

        # Load existing best now so we can show it as a reference target during progress.
        # merge_with_existing() at the end will preserve this value even if this run
        # doesn't match it — so the saved result can never regress.
        existing = self.load_existing_results()
        existing_best = existing["best_value"] if existing else None

        print(f"\n{'='*70}")
        print(f"Multi-Restart Search for n={self.n}")
        if existing_best is not None:
            print(f"Previous best (from file): {existing_best}  <-- target to match/beat")
        else:
            print(f"No existing results found for n={self.n}")
        if self.target_density is not None:
            print(f"Target density: {self.target_density:.2f} ({self.target_ones} ones)")
        else:
            print(f"Target density: None (free exploration)")
        if min_density is not None or max_density is not None:
            density_str = f"[{min_density if min_density is not None else '0.0'}, {max_density if max_density is not None else '1.0'}]"
            print(f"Density constraints: {density_str}")
        print(f"Constraints on ones blocks: min={min_ones}, max={max_ones}")
        print(f"Constraints on zeros blocks: min={min_zeros}, max={max_zeros}")
        print(f"Ones-block count constraints: min={min_ones_blocks}, max={max_ones_blocks}")
        print(f"Using {num_processes} processes for {num_restarts} restarts")
        print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*70}\n")

        # Create worker arguments for each restart
        worker_args = [
            (self.n, self.target_density, iterations_per_run,
             min_ones, max_ones, min_zeros, max_zeros, min_density, max_density,
             min_ones_blocks, max_ones_blocks, restart_id)
            for restart_id in range(num_restarts)
        ]

        # Run searches in parallel with progress tracking
        all_results = []

        if num_processes > 1:
            # Create pool without 'with' statement so we control cleanup precisely.
            # Using 'with Pool(...) as pool:' causes a double terminate()+join() on
            # Ctrl+C (once in the except block, once in __exit__), which hangs on
            # Windows when workers are killed mid-startup.
            pool = Pool(processes=num_processes, initializer=_worker_init)
            async_results = []
            try:
                # Submit all jobs asynchronously
                print("Submitting jobs to worker processes...")
                async_results = [pool.apply_async(_run_single_search, args) for args in worker_args]
                print(f"All {num_restarts} jobs submitted. Waiting for results...\n")

                # Poll for completion and show progress
                completed = 0
                next_milestone = 5
                last_update_time = self.start_time

                while completed < num_restarts:
                    # Check how many are done
                    current_completed = sum(1 for ar in async_results if ar.ready())

                    # Update progress
                    elapsed = time.time() - self.start_time
                    progress = (current_completed / num_restarts) * 100
                    time_since_update = time.time() - last_update_time

                    # Show progress at 5% intervals or every 30 seconds
                    should_update = (progress >= next_milestone) or (time_since_update >= 30)

                    if should_update:
                        # Get results from completed tasks
                        finished_results = [ar.get() for ar in async_results if ar.ready()]
                        target_str = f" (target: {existing_best})" if existing_best is not None else ""
                        if finished_results:
                            current_best = max([r[0] for r in finished_results])
                            distinct_best = len({tuple(r[2]) for r in finished_results if r[0] == current_best})
                            restarts_at_best = sum(1 for r in finished_results if r[0] == current_best)
                            print(f"Progress: {current_completed}/{num_restarts} | "
                                  f"Best: {current_best}{target_str} | "
                                  f"Distinct solutions: {distinct_best} ({restarts_at_best} hits) | "
                                  f"Elapsed: {elapsed:.1f}s")
                        else:
                            # No results yet, just show we're working
                            print(f"Progress: {current_completed}/{num_restarts} | "
                                  f"Working...{target_str} | "
                                  f"Elapsed: {elapsed:.1f}s")

                        if progress >= next_milestone:
                            next_milestone += 5
                        last_update_time = time.time()

                    completed = current_completed

                    # Sleep briefly to avoid busy-waiting
                    if completed < num_restarts:
                        time.sleep(0.5)

                # Collect all results
                all_results = [ar.get() for ar in async_results]
                pool.close()
                pool.join()

            except KeyboardInterrupt:
                print("\n\nKeyboard interrupt received! Terminating workers...")
                pool.terminate()
                # Do NOT call pool.join() here — it hangs on Windows when workers
                # are killed mid-startup (before _worker_init has run). The OS will
                # clean up child processes when this process exits.
                print("Workers terminated. Saving partial results...\n")
                self._interrupted = True

                # Collect whatever results we have
                for ar in async_results:
                    if ar.ready():
                        try:
                            all_results.append(ar.get(timeout=0.1))
                        except Exception:
                            pass

                if not all_results:
                    print("No completed results to save. Exiting.")
                    return 0, []

                # Quick bit-swap expansion on partial results before saving.
                # This costs only O(|solutions| × n²) evaluations and catches
                # all Hamming-distance-1 solutions that the interrupted MH run
                # may have narrowly missed.  Runs in the main process serially
                # (no pool available after terminate) but is fast enough for
                # n ≤ ~200 even with hundreds of partial solutions.
                best_partial = max(r[0] for r in all_results)
                best_partial_bits = {
                    tuple(r[2]) for r in all_results if r[0] == best_partial
                }
                print(f"Running quick 1-hop bit-swap expansion on {len(best_partial_bits)} "
                      f"partial solutions (value={best_partial})...")
                expand_pool_int = Pool(processes=cpu_count(), initializer=_worker_init)
                expand_args_int = [
                    (self.n, best_partial, ''.join(map(str, b)))
                    for b in best_partial_bits
                ]
                extra_found: Set[str] = set()
                try:
                    for new_strs in expand_pool_int.imap_unordered(_bitswap_worker, expand_args_int):
                        extra_found.update(new_strs)
                    expand_pool_int.close()
                    expand_pool_int.join()
                except Exception:
                    expand_pool_int.terminate()
                if extra_found:
                    extra_filtered = extra_found - {''.join(map(str, b)) for b in best_partial_bits}
                    print(f"  +{len(extra_filtered)} extra solution(s) from bit-swap expansion")
                    for bits_str in extra_filtered:
                        bits_list = [int(c) for c in bits_str]
                        all_results.append((best_partial, self.bits_to_runs(bits_list), bits_list))
        else:
            # Single process mode for debugging
            completed = 0
            next_milestone = 5

            try:
                for args in worker_args:
                    result = _run_single_search(*args)
                    all_results.append(result)
                    completed += 1
                    progress = (completed / num_restarts) * 100

                    if progress >= next_milestone:
                        elapsed = time.time() - self.start_time
                        current_best = max([r[0] for r in all_results])
                        distinct_best = len({tuple(r[2]) for r in all_results if r[0] == current_best})
                        restarts_at_best = sum(1 for r in all_results if r[0] == current_best)
                        target_str = f" (target: {existing_best})" if existing_best is not None else ""
                        print(f"Progress: {completed}/{num_restarts} | "
                              f"Best: {current_best}{target_str} | "
                              f"Distinct solutions: {distinct_best} ({restarts_at_best} hits) | "
                              f"Elapsed: {elapsed:.1f}s")
                        next_milestone += 5
            except KeyboardInterrupt:
                print("\n\nKeyboard interrupt received! Saving partial results...\n")
                if not all_results:
                    print("No completed results to save. Exiting.")
                    return 0, []

        # Process results using set for faster duplicate checking
        seen_solutions = set()
        for restart_id, (count, runs, bits) in enumerate(all_results):
            if count > self.best_value:
                self.best_value = count
                self.best_solutions = [(runs, bits)]
                seen_solutions = {tuple(bits)}
            elif count == self.best_value:
                bits_tuple = tuple(bits)
                if bits_tuple not in seen_solutions:
                    self.best_solutions.append((runs, bits))
                    seen_solutions.add(bits_tuple)

        # Merge with existing results and save
        total_elapsed = time.time() - self.start_time
        self.merge_with_existing()

        # Iterative bit-swap expansion until convergence
        self.bitswap_expand()

        metadata = {
            "num_restarts": num_restarts,
            "iterations_per_run": iterations_per_run,
            "min_ones": min_ones,
            "max_ones": max_ones,
            "min_zeros": min_zeros,
            "max_zeros": max_zeros,
            "min_density": min_density,
            "max_density": max_density,
            "target_density": self.target_density,
            "total_time_seconds": total_elapsed,
            "min_ones_blocks": min_ones_blocks,
            "max_ones_blocks": max_ones_blocks
        }
        self.save_results(metadata)

        return self.best_value, self.best_solutions


def _worker_init():
    """
    Initializer for Pool workers: ignore SIGINT so only the main process
    handles Ctrl+C.  When the main process catches KeyboardInterrupt it
    calls pool.terminate(), which sends SIGTERM to workers — a clean
    shutdown path that avoids pool.join() hanging on Windows.

    Also redirects worker stderr to /dev/null to suppress BrokenPipeError
    noise that appears when the pool is shut down while workers are still
    trying to send results back through the now-closed pipes.
    Worker computation errors are caught in _run_single_search and printed
    to stdout, so nothing meaningful is lost.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    import sys
    try:
        sys.stderr = open(os.devnull, 'w')
    except Exception:
        pass


def _level_set_worker(args: Tuple) -> Set[str]:
    """
    Module-level worker for find_all_optimal parallel execution.

    Args:
        args: (n, target_value, seed_bits, steps, chain_id, epsilon)

    Returns:
        Set of bit-strings (str) found that achieve ``target_value``.
    """
    n, target_value, seed_bits, steps, chain_id, epsilon = args
    import struct
    base = struct.unpack('>I', os.urandom(4))[0]
    random.seed(base ^ (chain_id * 98317 + n * 13))

    opt = RunLengthOptimizer(n=n)
    return opt.level_set_walk(target_value, seed_bits, steps, epsilon=epsilon)


def _bitswap_worker(args: Tuple) -> Set[str]:
    """
    Module-level worker for exhaustive bit-neighbourhood expansion.

    Uses the Numba-JIT ``_bitswap_count_flips`` kernel to evaluate all
    single-bit flips and density-preserving swaps in compiled native code,
    then converts the resulting numpy rows back to Python strings.

    Args:
        args: (n, target_value, seed_str)

    Returns:
        Set of bit-strings (str) that are one bit-level step away from the seed
        and achieve ``target_value``.
    """
    n, target_value, seed_str = args
    bits_arr = np.array([int(c) for c in seed_str], dtype=np.int8)
    hits = _bitswap_count_flips(bits_arr, target_value)
    found: Set[str] = set()
    for row_idx in range(hits.shape[0]):
        found.add(''.join(str(int(hits[row_idx, j])) for j in range(n)))
    return found


def _run_single_search(n: int, target_density: Optional[float], max_iterations: int,
                      min_ones: Optional[List[int]], max_ones: Optional[List[Optional[int]]],
                      min_zeros: Optional[List[int]], max_zeros: Optional[List[Optional[int]]],
                      min_density: Optional[float], max_density: Optional[float],
                      min_ones_blocks: Optional[int], max_ones_blocks: Optional[int],
                      restart_id: int) -> Tuple[int, List[int], List[int]]:
    """
    Worker function for parallel execution.

    This function runs a single Metropolis-Hastings search and returns the best result.
    It must be at module level for multiprocessing to work properly.

    Parameter Diversification:
    -------------------------
    To encourage exploration of different regions, this function varies parameters
    across restarts:
    - Density: ±0.05 variation from target
    - Constraints: ±1 variation in min/max bounds
    - Temperature: Different initial temperatures and cooling rates

    Args:
        n: Length of binary strings
        target_density: Base target density (will be varied)
        max_iterations: Iterations per search
        min_ones, max_ones, min_zeros, max_zeros: Constraint arrays
        min_density, max_density: Density constraints
        restart_id: Index of this restart (used to seed variation)

    Returns:
        Tuple of (count, runs, bits) for the best solution found
    """
    try:
        # Seed the PRNG so workers are (a) different from each other within a run, and
        # (b) different across runs.  Using only restart_id gives a fixed deterministic
        # seed, so every run finds the same solution.  XOR-ing with OS entropy (4 random
        # bytes) ensures uniqueness across runs while restart_id keeps workers distinct.
        import struct
        base = struct.unpack('>I', os.urandom(4))[0]
        random.seed(base ^ (restart_id * 12347 + n * 7))

        # Diversify search parameters across restarts
        # Vary density slightly to explore different regions (only if target density specified)
        if target_density is not None:
            density_variation = (restart_id % 5) * 0.02 - 0.05  # -0.05 to +0.05
            varied_density = max(0.65, min(0.95, target_density + density_variation))
        else:
            # No target density - keep it None for free exploration
            varied_density = None

        # Vary constraints across restarts (±1 for diversity) if they exist
        if min_ones is not None:
            varied_min_ones = [max(0, min_ones[i] + ((restart_id >> i) % 3) - 1) for i in range(4)]
        else:
            varied_min_ones = None

        if max_ones is not None:
            varied_max_ones = []
            for i in range(4):
                max_val = max_ones[i]
                if max_val is None:
                    varied_max_ones.append(None)
                else:
                    varied_max_ones.append(max_val + ((restart_id >> (i+4)) % 3) - 1)
        else:
            varied_max_ones = None

        if min_zeros is not None:
            varied_min_zeros = [max(0, min_zeros[i] + ((restart_id >> (i+8)) % 3) - 1) for i in range(4)]
        else:
            varied_min_zeros = None

        if max_zeros is not None:
            varied_max_zeros = []
            for i in range(4):
                max_val = max_zeros[i]
                if max_val is None:
                    varied_max_zeros.append(None)
                else:
                    varied_max_zeros.append(max_val + ((restart_id >> (i+12)) % 3) - 1)
        else:
            varied_max_zeros = None

        # Vary ones-block count constraints (±1 for diversity)
        if min_ones_blocks is not None:
            varied_min_ones_blocks = max(1, min_ones_blocks + (restart_id % 3) - 1)
        else:
            varied_min_ones_blocks = None

        if max_ones_blocks is not None:
            varied_max_ones_blocks = max_ones_blocks + (restart_id % 3) - 1
            if varied_min_ones_blocks is not None:
                varied_max_ones_blocks = max(varied_min_ones_blocks, varied_max_ones_blocks)
        else:
            varied_max_ones_blocks = None

        # Create a new optimizer instance for this worker
        optimizer = RunLengthOptimizer(n=n, target_density=varied_density)

        # Original temperature and cooling parameters (proven to work well up to n~200)
        T_initial = 10.0 + (restart_id % 4) * 2.5  # 10.0, 12.5, 15.0, 17.5
        cooling = 0.9995 - (restart_id % 3) * 0.0001  # Slightly different cooling rates

        # Run the search (silently)
        count, runs, bits = optimizer.metropolis_hastings(
            max_iterations=max_iterations,
            T_initial=T_initial,
            cooling_rate=cooling,
            min_ones=varied_min_ones,
            max_ones=varied_max_ones,
            min_zeros=varied_min_zeros,
            max_zeros=varied_max_zeros,
            min_density=min_density,
            max_density=max_density,
            use_adaptive_temp=True,
            min_ones_blocks=varied_min_ones_blocks,
            max_ones_blocks=varied_max_ones_blocks
        )

        return count, runs, bits
        
    except Exception as e:
        print(f"Error in worker {restart_id}: {e}")
        import traceback
        traceback.print_exc()
        return 0, [], []

def main():
    """
    Example usage demonstrating the full search pipeline.

    This shows how to:
    1. Set up the optimizer
    2. Configure parallel search
    3. Apply constraints
    4. Save and report results
    """
    # ========================================================================
    # USER CONFIGURATION - CHANGE THESE VALUES AS NEEDED
    # ========================================================================
    
    # Basic parameters
    n = 194                    # Length of binary string to optimize
    target_density = None      # None = free exploration (no density bias)

    # Search parameters
    # num_restarts is set large so the script fills the full run time;
    # interrupt with Ctrl+C at any point — partial results are saved.
    num_restarts = 15_000      # Effectively unlimited; Ctrl+C saves progress
    iterations_per_run = 100_000  # Short restarts → more families explored per hour
    num_processes = None           # 8 of 16 cores (other 8 used by _mh_n129.py)

    # Structural constraints (None = no constraint)
    # Constraints on first 4 blocks of ones: [1st, 2nd, 3rd, 4th]
    # min_ones: Optional[List[int]] = [17, 15, 11, 1]            # Example: [18, 12, 10, 5]
    # max_ones: List[Optional[int]] = [22, 18, 15, 10]             # Example: [25, 20, 15, 12]
    min_ones: Optional[List[int]] = None            # Example: [18, 12, 10, 5]
    max_ones: List[Optional[int]] = [None, None, None, None]       
    
    # Constraints on first 4 blocks of zeros: [1st, 2nd, 3rd, 4th]
    # min_zeros: Optional[List[int]] = [1, 1, 1, 1]           # Example: [1, 1, 0, 0]
    # max_zeros: List[Optional[int]] = [3, 3, 3, 3]            # Example: [3, 5, 3, 3]
    min_zeros: Optional[List[int]] = None           # Example: [1, 1, 0, 0]
    max_zeros: List[Optional[int]] = [None, None, None, None]            # Example: [3, 5, 3, 3]
    
    # Density constraints (None = no constraint)
    min_density = None          # Minimum density, e.g., 0.75
    max_density = None          # Maximum density, e.g., 0.88 

    # # Structural constraints (None = no constraint)
    # # Constraints on first 4 blocks of ones: [1st, 2nd, 3rd, 4th]
    # min_ones: Optional[List[int]] = [26, 24, 20, 1]
    # max_ones: List[Optional[int]] = [30, 27, 24, 3]
    
    # # Constraints on first 4 blocks of zeros: [1st, 2nd, 3rd, 4th]
    # min_zeros: Optional[List[int]] = [1, 2, 1, 1]
    # max_zeros: List[Optional[int]] = [2, 3, 2, 2]
    
    # # Density constraints (None = no constraint)
    # min_density = 0.83        # Minimum density, e.g., 0.75
    # max_density = 0.88          # Maximum density, e.g., 0.88

    # Ones-block count constraints (None = no constraint)
    # Constrains the total number of distinct runs of 1s in the binary string.
    # Because the string starts with 1s and alternates:
    #   k ones-blocks  →  k-1 zeros-blocks (ends in 1s)  or  k zeros-blocks (ends in 0s)
    # e.g. "111 00 11 0 1" has 3 ones-blocks and 2 zeros-blocks.
    min_ones_blocks: Optional[int] = None   # e.g. 6  → at least 6 runs of 1s
    max_ones_blocks: Optional[int] = None   # e.g. 10 → at most  10 runs of 1s

    # Solution exploration (run AFTER the MH search above has completed)
    # Set to True to run additional passes that discover more distinct solutions
    # at whatever best value was found.  Does not require knowing the optimum.
    run_level_set_walk  = True  # Parallel random walks pinned to the best value
    level_set_chains    = None  # None = cpu_count() × 4
    level_set_steps     = 200_000  # Proposals per chain
    level_set_epsilon   = 0          # Starting epsilon (0 = strict level-set)
    level_set_epsilon_max = 3        # Auto-escalate up to this if stuck (0 = disable)
                                    # Allows walk to cross barriers between disconnected
                                    # solution components; escalates automatically per round
    level_set_max_rounds = 4        # Iterative expansion rounds (stops early if converged)
    level_set_bitswap_expand = True  # Exhaustive 1-hop bit-swap expansion after each round

    # ========================================================================
    # END USER CONFIGURATION
    # ========================================================================
    
    # Set up signal handler for clean exit
    original_sigint = signal.getsignal(signal.SIGINT)

    try:
        # Initialize optimizer
        optimizer = RunLengthOptimizer(n=n, target_density=target_density, results_dir="results/mh_unbounded")

        # Get number of available CPU cores
        num_cores = cpu_count()
        print(f"Detected {num_cores} CPU cores")
        
        # Set default num_restarts if not specified
        if num_restarts is None:
            num_restarts = num_cores * 64

        # Run the search
        best_value, solutions = optimizer.multi_restart_search(
            num_restarts=num_restarts,
            iterations_per_run=iterations_per_run,
            min_ones=min_ones,
            max_ones=max_ones,
            min_zeros=min_zeros,
            max_zeros=max_zeros,
            min_density=min_density,
            max_density=max_density,
            num_processes=num_processes,
            min_ones_blocks=min_ones_blocks,
            max_ones_blocks=max_ones_blocks
        )

        total_time = time.time() - (optimizer.start_time or time.time())

        print(f"\n{'='*70}")
        print(f"SEARCH COMPLETE")
        print(f"{'='*70}")
        print(f"n = {n}")
        print(f"Best value: a({n}) = {best_value}")
        print(f"Total optimal solutions found: {len(solutions)}")
        print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        print(f"Results saved to: {optimizer.results_dir}/n_{n:04d}_results.json")

        # Show first few solutions
        print(f"\nOptimal solutions:")
        for idx, (runs, bits) in enumerate(solutions[:3]):  # Show first 3
            cleaned_runs = optimizer.clean_runs(runs)
            print(f"\nSolution {idx + 1}:")
            print(f"  Runs: {cleaned_runs}")
            print(f"  Bits: {''.join(map(str, bits))}")
            print(f"  Density: {sum(bits)/len(bits):.3f} ({sum(bits)}/{n} ones)")

        if len(solutions) > 3:
            print(f"\n... and {len(solutions) - 3} more solution(s)")

        # ----------------------------------------------------------------
        # Post-search: find more distinct solutions at the best found value
        # (no knowledge of the optimum needed — value is read from file)
        # ----------------------------------------------------------------
        if run_level_set_walk:
            if optimizer._interrupted:
                print("\nNote: MH search was interrupted — running level-set walk on "
                      "partial results. The best value found may not be the true optimum.")
            print(f"\n{'='*70}")
            print("Level-set walk (parallel)")
            print(f"{'='*70}")
            optimizer.find_all_optimal(
                num_chains=level_set_chains,
                steps_per_chain=level_set_steps,
                epsilon=level_set_epsilon,
                epsilon_max=level_set_epsilon_max,
                max_rounds=level_set_max_rounds,
                bitswap_expand=level_set_bitswap_expand,
            )

        print(f"\n{'='*70}")

    except KeyboardInterrupt:
        print("\n\nProgram interrupted by user.")
    finally:
        # Restore original signal handler
        signal.signal(signal.SIGINT, original_sigint)

if __name__ == "__main__":
    # Required for Windows multiprocessing
    main()
