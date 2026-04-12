"""
OEIS A112509 - Surgical Nudge Optimizer
========================================
Starts from the best known seed for a given n and refines it by targeting the
highest-LCP collision in the suffix array — the worst source of "wasted"
substring space — and making the smallest possible bit-swap to break it.
Algorithm
---------
1.  Build the suffix array (SA) via numpy prefix-doubling   [O(n log² n) time,
    O(n) extra space; fast in practice due to C-level numpy argsort].
2.  Build the LCP array via Kasai's algorithm (numba JIT when available).
3.  Score  =  (distinct substrings starting with '1')  +  (1 if any '0' exists)
              i.e. the A112509 value.
4.  Surgical Nudge loop:
      a. Find the top-K LCP entries (worst collisions).
      b. For each collision (SA[i-1], SA[i]) with shared prefix length lcp[i]:
           - Collect candidate positions: every '0' bit within the shared
             region of either suffix (clamped to string bounds).
           - For each candidate position j, try swapping bits[j] (=0) with
             bits[j+d] (=1) for d in {-2, -1, +1, +2}.
           - Evaluate each swap; keep the best that improves the score.
      c. If no surgical move helps, fall back to random bit-swaps (ILS kick).
5.  Save improved solutions to seeds/ and results/ using the same conventions
    as template_greedy_search.py.
"""
import random
import json
import os
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from typing import List, Optional, Tuple
from pathlib import Path
# ---------------------------------------------------------------------------
# Optional C-speed suffix array via pydivsufsort (libdivsufsort wrapper).
# If available, SA construction is ~40x faster than the numpy prefix-doubling.
# Install with:  pip install pydivsufsort
# ---------------------------------------------------------------------------
try:
    from pydivsufsort import divsufsort as _divsufsort
    _DIVSUFSORT_AVAILABLE = True
except ImportError:
    _DIVSUFSORT_AVAILABLE = False
# ---------------------------------------------------------------------------
# Numba JIT Kasai — identical to template_greedy_search.py so the JIT cache
# is already warm when both scripts are used in the same session.
# ---------------------------------------------------------------------------
try:
    from numba import njit as _njit

    @_njit(cache=True)
    def _kasai_jit(bits_arr, sa_arr, pos_arr, n):
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

    @_njit(cache=True)
    def _score_jit(bits_arr, sa_arr, lcp_arr, n):  # type: ignore[misc]
        count = np.int64(0)
        has_zero = False
        for i in range(n):
            if bits_arr[sa_arr[i]] == 1:
                count += np.int64(n - sa_arr[i]) - np.int64(lcp_arr[i])
            if not has_zero and bits_arr[i] == 0:
                has_zero = True
        if has_zero:
            count += np.int64(1)
        return count

    @_njit(cache=True)
    def _hot_positions_jit(sa_arr, lcp_arr, top_k_arr, n, window):
        mask = np.zeros(n, dtype=np.bool_)
        for t in range(len(top_k_arr)):
            idx = top_k_arr[t]
            if idx == 0:
                continue
            L = lcp_arr[idx]
            if L == 0:
                continue
            for origin_sel in range(2):
                if origin_sel == 0:
                    origin = sa_arr[idx - 1]
                else:
                    origin = sa_arr[idx]
                offset = L - window
                if offset < 0:
                    offset = 0
                lo = origin + offset
                if lo < 0:
                    lo = 0
                hi = origin + L + 5
                if hi > n:
                    hi = n
                for pos in range(lo, hi):
                    mask[pos] = True
        return mask

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

    def _kasai_jit(bits_arr, sa_arr, pos_arr, n):  # type: ignore[misc]
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

    def _score_jit(bits_arr, sa_arr, lcp_arr, n):  # type: ignore[misc]
        count = 0
        has_zero = False
        for i in range(n):
            if bits_arr[sa_arr[i]] == 1:
                count += (n - int(sa_arr[i])) - int(lcp_arr[i])
            if not has_zero and bits_arr[i] == 0:
                has_zero = True
        if has_zero:
            count += 1
        return count

    _hot_positions_jit = None  # type: ignore[assignment]
# ---------------------------------------------------------------------------
# Core computational primitives
# ---------------------------------------------------------------------------

def build_sa_lcp(bits: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (suffix_array, lcp_array) for a binary numpy array.
    Uses pydivsufsort (C library, ~40x faster) when available, falling back
    to the pure-numpy O(n log²n) prefix-doubling implementation.
    """
    n = len(bits)
    if n == 0:
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)
    if _DIVSUFSORT_AVAILABLE:
        # divsufsort expects uint8; returns int64 SA → cast to int32
        sa = _divsufsort(bits.astype(np.uint8)).astype(np.int32)
    else:
        # --- Fallback: suffix array via prefix doubling (numpy argsort) ---
        idx = np.arange(n, dtype=np.int32)
        rank = bits.copy().astype(np.int32)
        k = 1
        while k < n:
            next_idx = np.minimum(idx + k, n - 1)
            r2 = np.where(idx + k < n, rank[next_idx], -1)
            max_r = int(rank.max()) + 2
            key = rank.astype(np.int64) * max_r + r2.astype(np.int64)
            sa = np.argsort(key, kind='stable').astype(np.int32)
            changed = np.empty(n, dtype=np.int32)
            changed[0] = 0
            changed[1:] = (key[sa[1:]] != key[sa[:-1]]).cumsum()
            rank = np.empty(n, dtype=np.int32)
            rank[sa] = changed
            if int(rank[sa[-1]]) == n - 1:
                break
            k *= 2
    # --- LCP array via Kasai (numba JIT when available) ---
    pos = np.empty(n, dtype=np.int32)
    pos[sa] = np.arange(n, dtype=np.int32)
    lcp = _kasai_jit(bits, sa, pos, n)
    return sa, lcp

def score_from_sa_lcp(bits: np.ndarray, sa: np.ndarray, lcp: np.ndarray) -> int:
    """Compute A112509 score from a pre-built SA/LCP pair.
    Uses numba JIT single-pass loop (no intermediate arrays) when available."""
    return int(_score_jit(bits, sa, lcp, len(bits)))

def compute_score(bits: np.ndarray) -> int:
    """Full score computation (builds SA+LCP internally)."""
    sa, lcp = build_sa_lcp(bits)
    return score_from_sa_lcp(bits, sa, lcp)

# ---------------------------------------------------------------------------
# Parallel batch evaluation
#
# For large n (>1M), each SA build takes 0.1–7s.  Instead of evaluating
# candidates sequentially, we generate a batch of swap specifications and
# evaluate them in parallel using a process pool.  Workers read the base
# bitstring from shared memory, apply their swap, build SA+LCP, score,
# and return just the score integer — avoiding pickling large arrays.
# ---------------------------------------------------------------------------

import multiprocessing as mp
from multiprocessing import shared_memory as _shm

_POOL: Optional[ProcessPoolExecutor] = None
_NUM_WORKERS: int = 1
_SHM_NAME: Optional[str] = None
_SHM_SIZE: int = 0
_SHM_OBJ: Optional[_shm.SharedMemory] = None  # keep reference alive on Windows


def _init_pool(num_workers: int, n: int = 0) -> None:
    """Initialize global process pool and shared memory for base bitstring."""
    global _POOL, _NUM_WORKERS, _SHM_NAME, _SHM_SIZE, _SHM_OBJ
    _cleanup_shm()  # clean up any existing shared memory
    if num_workers <= 1:
        _POOL = None
        _NUM_WORKERS = 1
        return
    _NUM_WORKERS = num_workers
    if n > 0:
        _SHM_SIZE = n * 4  # int32 = 4 bytes
        _SHM_OBJ = _shm.SharedMemory(create=True, size=_SHM_SIZE)
        _SHM_NAME = _SHM_OBJ.name
    _POOL = ProcessPoolExecutor(max_workers=num_workers)
    print(f"Process pool initialized: {num_workers} workers"
          + (f", shared memory: {_SHM_SIZE / 1e6:.0f}MB" if _SHM_NAME else ""))


def _update_shared_bits(bits: np.ndarray) -> None:
    """Copy current bitstring into shared memory."""
    if _SHM_OBJ is None:
        return
    arr = np.ndarray(len(bits), dtype=np.int32, buffer=_SHM_OBJ.buf)
    arr[:] = bits


def _cleanup_shm() -> None:
    """Release shared memory."""
    global _SHM_NAME, _SHM_OBJ
    if _SHM_OBJ is not None:
        try:
            _SHM_OBJ.close()
            _SHM_OBJ.unlink()
        except Exception:
            pass
        _SHM_OBJ = None
        _SHM_NAME = None


def _shutdown_pool() -> None:
    """Shutdown global process pool and shared memory."""
    global _POOL
    if _POOL is not None:
        _POOL.shutdown(wait=False, cancel_futures=True)
        _POOL = None
    _cleanup_shm()


def _worker_eval_swap(shm_name: str, n: int, swaps: List[Tuple[int, int]]) -> int:
    """Worker: read base bits from shared memory, apply swaps, score.
    swaps is a list of (position, new_value) pairs.
    Returns score.
    """
    shm = _shm.SharedMemory(name=shm_name)
    bits = np.ndarray(n, dtype=np.int32, buffer=shm.buf).copy()  # copy to own memory
    shm.close()
    for pos, val in swaps:
        bits[pos] = val
    sa, lcp = build_sa_lcp(bits)
    return score_from_sa_lcp(bits, sa, lcp)


def parallel_eval_best(
    bits: np.ndarray,
    swap_specs: List[List[Tuple[int, int]]],
    threshold: int = 0,
) -> Tuple[int, int, Optional[np.ndarray], Optional[np.ndarray]]:
    """Evaluate a batch of swap candidates, return the best.

    bits: current bitstring (used as base for all candidates).
    swap_specs: list of [(pos, new_val), ...] diffs from bits.
    threshold: minimum score to beat.

    Returns (best_score, best_index, best_sa, best_lcp).
    If no candidate beats threshold, returns (threshold, -1, None, None).
    """
    if not swap_specs:
        return threshold, -1, None, None

    n = len(bits)

    # Single-threaded: sequential with early exit
    if _POOL is None or _SHM_NAME is None or len(swap_specs) == 1:
        best_val = threshold
        best_idx = -1
        best_sa = None
        best_lcp = None
        for i, swaps in enumerate(swap_specs):
            # Apply swaps in-place, evaluate, revert
            old_vals = [(pos, int(bits[pos])) for pos, val in swaps]
            for pos, val in swaps:
                bits[pos] = val
            sa_c, lcp_c = build_sa_lcp(bits)
            val = score_from_sa_lcp(bits, sa_c, lcp_c)
            if val > best_val:
                best_val = val
                best_idx = i
                best_sa = sa_c
                best_lcp = lcp_c
                # Revert and break (first improvement)
                for pos, old_val in old_vals:
                    bits[pos] = old_val
                break
            # Revert
            for pos, old_val in old_vals:
                bits[pos] = old_val
        return best_val, best_idx, best_sa, best_lcp

    # Multi-process: update shared memory, submit swap specs
    _update_shared_bits(bits)
    futures = {}
    for i, swaps in enumerate(swap_specs):
        fut = _POOL.submit(_worker_eval_swap, _SHM_NAME, n, swaps)
        futures[fut] = i

    best_val = threshold
    best_idx = -1
    for fut in as_completed(futures):
        i = futures[fut]
        try:
            val = fut.result()
        except Exception:
            continue
        if val > best_val:
            best_val = val
            best_idx = i

    # Rebuild SA/LCP for the winner
    best_sa = None
    best_lcp = None
    if best_idx >= 0:
        winner = bits.copy()
        for pos, val in swap_specs[best_idx]:
            winner[pos] = val
        best_sa, best_lcp = build_sa_lcp(winner)

    return best_val, best_idx, best_sa, best_lcp

# ---------------------------------------------------------------------------
# Surgical nudge primitives
# ---------------------------------------------------------------------------

def find_top_collision_indices(lcp: np.ndarray, k: int = 20) -> np.ndarray:
    """Return the k SA-ranks with highest LCP values (worst collisions first)."""
    if k >= len(lcp):
        return np.argsort(lcp)[::-1]
    return np.argpartition(lcp, -k)[-k:][ np.argsort(lcp[np.argpartition(lcp, -k)[-k:]])[::-1] ]

def surgical_candidates(
    bits: np.ndarray,
    sa: np.ndarray,
    lcp: np.ndarray,
    lcp_idx: int,
    deltas: Tuple[int, ...] = (-2, -1, 1, 2),
    window_extra: int = 5,
    max_window: int = 500,
) -> List[Tuple[int, int]]:
    """Return (j, d) swap pairs that move a '0' bit near a collision.
    The collision at rank lcp_idx involves suffixes starting at p=SA[lcp_idx-1]
    and q=SA[lcp_idx].  They share L=lcp[lcp_idx] characters.  We look for '0'
    bits in the shared regions near the END of the overlap (where a small change
    most efficiently breaks the collision), capped to max_window characters so
    that very long collisions don't explode the candidate count.
    Specifically the search window is:
        [origin + max(0, L - max_window), origin + L + window_extra)
    """
    n = len(bits)
    L = int(lcp[lcp_idx])
    p = int(sa[lcp_idx - 1])
    q = int(sa[lcp_idx])
    candidates = []
    seen = set()
    for origin in (p, q):
        lo = max(0, origin + max(0, L - max_window))
        hi = min(n, origin + L + window_extra)
        for j in range(lo, hi):
            if bits[j] != 0:
                continue
            for d in deltas:
                jd = j + d
                if jd < 0 or jd >= n:
                    continue
                if bits[jd] != 1:
                    continue   # swap only flips if the two positions differ
                key = (min(j, jd), max(j, jd))
                if key not in seen:
                    seen.add(key)
                    candidates.append((j, d))
    return candidates

def apply_swap(bits: np.ndarray, j: int, d: int) -> np.ndarray:
    """Return a copy of bits with positions j and j+d swapped."""
    new_bits = bits.copy()
    new_bits[j], new_bits[j + d] = new_bits[j + d], new_bits[j]
    return new_bits

def _swap_inplace(bits: np.ndarray, j: int, jd: int) -> None:
    """Swap bits[j] and bits[jd] in-place (no allocation)."""
    bits[j], bits[jd] = bits[jd], bits[j]

def _hot_positions(
    sa: np.ndarray, lcp: np.ndarray, n: int, k: int = 30, window: int = 500
) -> np.ndarray:
    """Return a boolean mask of length n: True at positions that fall inside
    at least one of the top-k collision windows.  Zeros at these positions
    directly contribute to the worst collisions, so swapping them is more
    likely to improve the score than swapping a random zero elsewhere."""
    top_k = find_top_collision_indices(lcp, k)
    if _hot_positions_jit is not None:
        return _hot_positions_jit(sa, lcp, top_k, n, window)
    mask = np.zeros(n, dtype=bool)
    for idx in top_k:
        if idx == 0 or int(lcp[idx]) == 0:
            continue
        L = int(lcp[idx])
        for origin in (int(sa[idx - 1]), int(sa[idx])):
            lo = max(0, origin + max(0, L - window))
            hi = min(n, origin + L + 5)
            mask[lo:hi] = True
    return mask
# ---------------------------------------------------------------------------
# RLE-level structural moves
#
# These complement bit-level swaps by making changes that require multiple
# correlated bit flips — impossible to achieve one swap at a time.  Inspired
# by the MH algorithm's RLE representation (split, merge, transfer moves).
# ---------------------------------------------------------------------------

def _get_runs(bits: np.ndarray) -> List[Tuple[int, int, int]]:
    """Return list of (value, start_pos, length) for each run in bits."""
    n = len(bits)
    runs = []
    i = 0
    while i < n:
        v = int(bits[i])
        j = i
        while j < n and bits[j] == v:
            j += 1
        runs.append((v, i, j - i))
        i = j
    return runs


def structural_block_split(
    bits: np.ndarray, tries: int = 30
) -> Optional[np.ndarray]:
    """Split a large ones-block into [big_part, 0, small_block, 0, rest].

    Mimics the exact-solution pattern where small ones-blocks (1, 2, 3, ...)
    are interleaved between large ones-blocks.  Needs to steal 2 ones
    (to become the new separator zeros) plus the small-block ones from the
    big block, requiring multiple correlated bit changes.
    """
    runs = _get_runs(bits)
    ones_runs = [(v, pos, length, idx) for idx, (v, pos, length) in enumerate(runs) if v == 1]
    if len(ones_runs) < 2:
        return None
    # Sort by length descending — prefer splitting the biggest blocks
    big_runs = sorted(ones_runs, key=lambda x: -x[2])

    for _ in range(tries):
        # Pick a large ones-block to split
        pick_idx = min(int(np.random.exponential(2)), len(big_runs) - 1)
        _, pos, length, run_idx = big_runs[pick_idx]
        if length < 10:
            continue  # too small to split meaningfully

        # Choose small block size (1, 2, 3, ... with exponential preference for small)
        small_size = max(1, min(int(np.random.exponential(2)) + 1, length // 4))
        # Need: 2 zeros (new separators) + small_size ones from the block
        # Total ones consumed: small_size + 2 (turned to zeros)
        needed = small_size + 2
        if length - needed < 4:
            continue  # remainder too small

        # Pick split point within the block
        # Leave at least 3 ones on each side of the main parts
        margin = 3
        split_at = random.randint(margin, length - needed - margin)

        # Build new bits
        new_bits = bits.copy()
        # Original block at [pos, pos+length) all ones
        # New layout: [ones*split_at][0][ones*small_size][0][ones*(rest)]
        rest = length - split_at - needed
        cursor = pos + split_at
        new_bits[cursor] = 0          # first new separator
        cursor += 1
        # small block (already ones, leave them)
        cursor += small_size
        new_bits[cursor] = 0          # second new separator
        # Remaining ones after cursor+1 are already ones
        # But we consumed 2 ones (turned to 0), so we need to verify total
        # Actually: original block had `length` ones at [pos, pos+length).
        # We're setting bits[pos+split_at]=0 and bits[pos+split_at+1+small_size]=0
        # That gives: split_at ones, 0, small_size ones, 0, (length-split_at-small_size-2) ones
        return new_bits

    return None


def structural_block_merge(
    bits: np.ndarray, tries: int = 30
) -> Optional[np.ndarray]:
    """Merge a small ones-block with an adjacent block by removing the separator.

    Turns [big, 0, small, 0, big2] → [big+small+2, 0, big2] or similar,
    removing a separator and expanding an adjacent block.  Preserves total
    bits by turning 2 zeros into ones.
    """
    runs = _get_runs(bits)
    # Find small ones-blocks (size 1-5) that have zero-separators on both sides
    candidates = []
    for idx in range(2, len(runs) - 2):
        v, pos, length = runs[idx]
        if v != 1 or length > 5:
            continue
        # Check: left neighbor is zero-sep, left-left is ones
        if runs[idx - 1][0] == 0 and runs[idx - 2][0] == 1:
            # Check: right neighbor is zero-sep, right-right is ones
            if idx + 2 < len(runs) and runs[idx + 1][0] == 0 and runs[idx + 2][0] == 1:
                candidates.append(idx)
    if not candidates:
        return None

    for _ in range(tries):
        idx = random.choice(candidates)
        v, pos, length = runs[idx]
        left_sep_v, left_sep_pos, left_sep_len = runs[idx - 1]
        right_sep_v, right_sep_pos, right_sep_len = runs[idx + 1]

        new_bits = bits.copy()
        # Remove left separator: turn its zeros into ones → merges left block + small block
        if left_sep_len == 1:
            new_bits[left_sep_pos] = 1
            return new_bits
        # If left sep > 1, just shrink it by 1
        new_bits[left_sep_pos + left_sep_len - 1] = 1
        return new_bits

    return None


def structural_separator_relocate(
    bits: np.ndarray, tries: int = 30
) -> Optional[np.ndarray]:
    """Move a zero from a multi-zero separator to create a new single-zero
    separator elsewhere (splitting a large ones-block).

    This redistributes separators more evenly, turning e.g.
    [big1, 00, big2, big3] → [big1, 0, big2, 0, big3_left, big3_right].
    """
    runs = _get_runs(bits)
    # Find multi-zero separators (length >= 2) — donors
    donors = [(idx, v, pos, length) for idx, (v, pos, length) in enumerate(runs)
              if v == 0 and length >= 2]
    # Find large ones-blocks — recipients for a new separator
    large_ones = [(idx, v, pos, length) for idx, (v, pos, length) in enumerate(runs)
                  if v == 1 and length >= 8]
    if not donors or not large_ones:
        return None

    for _ in range(tries):
        # Pick a donor separator
        d_idx, _, d_pos, d_len = random.choice(donors)
        # Pick a large ones-block to split
        r_idx, _, r_pos, r_len = random.choice(large_ones)
        if d_idx == r_idx:
            continue

        new_bits = bits.copy()
        # Shrink donor: turn last zero of donor into a one
        new_bits[d_pos + d_len - 1] = 1
        # Insert new separator: turn a one in the middle of the recipient block into a zero
        split_point = r_pos + random.randint(3, r_len - 4)
        new_bits[split_point] = 0
        return new_bits

    return None


def structural_run_swap(
    bits: np.ndarray, tries: int = 30
) -> Optional[np.ndarray]:
    """Swap the lengths of two ones-blocks by physically rearranging bits.

    E.g. if block A has 100 ones and block B has 3 ones, swap so A gets 3
    and B gets 100.  This is the MH's non-adjacent RLE swap — it can move a
    small "punctuation" block from the tail into the early structure.
    """
    runs = _get_runs(bits)
    ones_runs = [(idx, pos, length) for idx, (v, pos, length) in enumerate(runs) if v == 1]
    if len(ones_runs) < 4:
        return None

    for _ in range(tries):
        # Pick two ones-blocks with different sizes
        i1, i2 = random.sample(range(len(ones_runs)), 2)
        idx_a, pos_a, len_a = ones_runs[i1]
        idx_b, pos_b, len_b = ones_runs[i2]
        if len_a == len_b:
            continue
        # Both must have same length of ones to swap — but they don't.
        # We need to redistribute: steal from the bigger, give to the smaller.
        # Just swap the run lengths by adjusting the boundary bits.
        # Simpler approach: fill both regions, then set correct lengths.
        new_bits = bits.copy()
        # Clear both regions to 1 first (they already are 1)
        # Set the shorter region: first len_b ones at pos_a (or len_a at pos_b)
        # But we can't change the region sizes — the separators are fixed positions.
        # Instead: we can only move ones between the two blocks.
        # Transfer: move (len_a - len_b) // 2 ones from bigger to smaller
        # by flipping border ones to zeros and border zeros to ones.
        # Actually this is complex — let's just do a direct content swap.
        if len_a > len_b:
            big_pos, big_len, small_pos, small_len = pos_a, len_a, pos_b, len_b
        else:
            big_pos, big_len, small_pos, small_len = pos_b, len_b, pos_a, len_a

        # Shrink the big block by turning some trailing ones into zeros
        # and expand the small block by turning adjacent zeros into ones
        # This only works if there's room. Try a partial transfer.
        transfer = min(big_len - small_len, max(2, big_len // 8))
        if transfer < 2:
            continue

        # Shrink big block: turn last `transfer` ones to zeros
        for k in range(transfer):
            new_bits[big_pos + big_len - 1 - k] = 0
        # Grow small block: turn `transfer` zeros after the block to ones
        # (or before, depending on what's available)
        grown = 0
        # Try growing rightward into the separator after small block
        pos_after = small_pos + small_len
        while grown < transfer and pos_after < len(new_bits) and new_bits[pos_after] == 0:
            new_bits[pos_after] = 1
            pos_after += 1
            grown += 1
        # Try growing leftward into the separator before small block
        pos_before = small_pos - 1
        while grown < transfer and pos_before >= 0 and new_bits[pos_before] == 0:
            new_bits[pos_before] = 1
            pos_before -= 1
            grown += 1
        # Any remaining transfer: just flip random zeros to ones elsewhere
        if grown < transfer:
            remaining_zeros = np.where(new_bits == 0)[0]
            if len(remaining_zeros) > 0:
                for k in range(min(transfer - grown, len(remaining_zeros))):
                    new_bits[remaining_zeros[k]] = 1
                    grown += 1
        return new_bits

    return None

# ---------------------------------------------------------------------------
# I/O helpers  (same conventions as template_greedy_search.py)
# ---------------------------------------------------------------------------

def generate_seed(n: int, results_dir: Path) -> np.ndarray:
    """Generate a high-quality seed by scaling the structure of the best
    available reference solution with n_ref <= n.

    Strategy (in order of preference):
    1. If a reference solution exists for some n_ref <= n, extract its full
       RLE structure (ones-blocks + separator lengths) and scale it up to n
       by proportionally enlarging block sizes while preserving separators.
    2. If no reference with n_ref <= n exists but a larger one does, scale
       it down similarly.
    3. If no reference at all, build from scratch using the known structural
       pattern from exact solutions: [big, big, big, 1, med, 2, med, 3, ...]
       with separator prefix [1, 2, 1, 1, 1, ...].
    """
    # ── Collect all available references ──────────────────────────────────
    refs = {}  # n_ref -> bitstring
    if results_dir.exists():
        for fname in sorted(results_dir.iterdir()):
            if not fname.name.endswith('_results.json'):
                continue
            try:
                data = json.loads(fname.read_text())
                ref_n = data.get('n', 0)
                sols = data.get('solutions', [])
                if sols and ref_n > 0 and ref_n != n:
                    s = sols[0] if isinstance(sols[0], str) else sols[0].get('bits', '')
                    if len(s) == ref_n and s.count('1') / len(s) > 0.90:
                        refs[ref_n] = s
            except Exception:
                continue

    # Also check cached_results.json for exact solutions (n <= 150)
    cached_path = results_dir.parent.parent / "data" / "cached_results.json"
    if cached_path.exists():
        try:
            cached = json.loads(cached_path.read_text())
            for key, entry in cached.items():
                ref_n = int(key)
                if ref_n == n or ref_n < 30:
                    continue
                strs = entry.get('optimal_strings', [])
                if strs and ref_n not in refs:
                    refs[ref_n] = strs[0]
        except Exception:
            pass

    # ── Pick reference: prefer largest n_ref <= n, else smallest n_ref > n ─
    leq = sorted([k for k in refs if k <= n], reverse=True)
    gt  = sorted([k for k in refs if k > n])
    ref_n = leq[0] if leq else (gt[0] if gt else None)
    best_ref = refs.get(ref_n) if ref_n is not None else None

    if best_ref is not None and ref_n is not None:
        return _scale_reference(best_ref, ref_n, n)
    else:
        return _build_from_pattern(n)


def _extract_rle(s: str) -> Tuple[List[int], List[int]]:
    """Extract (ones_blocks, sep_lengths) from a bitstring."""
    runs = []
    i = 0
    while i < len(s):
        c = s[i]; j = i
        while j < len(s) and s[j] == c:
            j += 1
        runs.append((c, j - i))
        i = j
    ones_blocks = [length for c, length in runs if c == '1']
    sep_lengths = [length for c, length in runs if c == '0']
    return ones_blocks, sep_lengths


def _scale_reference(ref_str: str, ref_n: int, target_n: int) -> np.ndarray:
    """Scale a reference solution's RLE structure to a different n.

    Enforces the known optimal structural pattern:
      [big0, big1, big2, 1, big3, 2, big4, 3, big5, 4, ...]
    Counter blocks at odd indices >= 3 are always set to (idx-3)//2 + 1,
    regardless of the reference's values there.  Only the separator pattern
    and the large-block relative sizes come from the reference.
    """
    ones_blocks, sep_lengths = _extract_rle(ref_str)
    ratio = target_n / ref_n
    print(f"Scaling reference n={ref_n} -> n={target_n} (ratio={ratio:.3f})")
    print(f"  Reference: {len(ones_blocks)} ones-blocks, {len(sep_lengths)} seps, "
          f"density={ref_str.count('1')/len(ref_str):.6f}")

    # Counter positions: odd indices >= 3 → values 1, 2, 3, 4, ...
    # These are structural and always enforced, never scaled.
    def is_counter_pos(idx: int) -> bool:
        return idx >= 3 and idx % 2 == 1

    def counter_value(idx: int) -> int:
        return (idx - 3) // 2 + 1

    # ── Determine separator list ──────────────────────────────────────────
    if abs(ratio - 1.0) < 0.01:
        new_seps = list(sep_lengths)
    elif ratio > 1.0:
        target_num_seps = max(len(sep_lengths), int(len(sep_lengths) * (ratio ** 0.5)))
        new_seps = list(sep_lengths)
        if target_num_seps > len(new_seps) and len(sep_lengths) > 3:
            tail_seps = sep_lengths[3:]
            while len(new_seps) < target_num_seps:
                new_seps.append(tail_seps[len(new_seps) % len(tail_seps)] if tail_seps else 1)
    else:
        target_num_seps = max(5, int(len(sep_lengths) * (ratio ** 0.5)))
        new_seps = sep_lengths[:target_num_seps]

    num_blocks = len(new_seps) + 1

    # ── Build large-block sizes as a decreasing arithmetic progression ────
    # Exact solutions show large blocks (even indices + 0,1,2) follow a
    # strictly decreasing pattern: [25, 21, 18, 14, 11, 7, 5, 6, 2].
    # We generate a decreasing arithmetic progression that sums to the
    # required budget, with all terms distinct and positive.
    num_large = sum(1 for i in range(num_blocks) if not is_counter_pos(i))
    counter_total = sum(counter_value(i) for i in range(num_blocks) if is_counter_pos(i))
    large_budget = target_n - sum(new_seps) - counter_total

    if num_large <= 1:
        large_sizes = [max(1, large_budget)]
    else:
        k = num_large
        # Arithmetic progression: a, a-d, a-2d, ..., a-(k-1)*d
        # Sum = k*a - d*k*(k-1)/2 = large_budget
        # Last term = a - (k-1)*d = last_val (must be >= 1)
        # Solve: last_val = max(1, large_budget // (k * 3))
        # a = last_val + (k-1)*d
        # d = 2*(large_budget - k*last_val) / (k*(k-1))
        last_val = max(1, large_budget // (k * 3))
        d = 2.0 * (large_budget - k * last_val) / (k * (k - 1))
        if d < 0:
            d = 0
        first_val = last_val + d * (k - 1)

        large_sizes = []
        for j in range(k):
            val = max(1, round(first_val - d * j))
            large_sizes.append(val)

        # Ensure all distinct: walk backwards and bump any duplicates down
        for j in range(len(large_sizes) - 2, -1, -1):
            if large_sizes[j] <= large_sizes[j + 1]:
                large_sizes[j] = large_sizes[j + 1] + 1

        # Adjust sum to exactly match large_budget — O(k) not O(diff)
        current = sum(large_sizes)
        diff = large_budget - current
        if diff > 0:
            per_block = diff // k
            leftover = diff % k
            for i in range(k):
                large_sizes[i] += per_block + (1 if i < leftover else 0)
        elif diff < 0:
            to_remove = abs(diff)
            per_block = to_remove // k
            leftover = to_remove % k
            for i in range(k):
                large_sizes[i] = max(1, large_sizes[i] - per_block - (1 if i < leftover else 0))

    # ── Assemble final ones-blocks: interleave large + counter ────────────
    new_ones = []
    large_iter = iter(large_sizes)
    for i in range(num_blocks):
        if is_counter_pos(i):
            new_ones.append(counter_value(i))
        else:
            new_ones.append(next(large_iter, 1))

    # ── Final length adjustment — O(num_blocks) not O(diff) ─────────────
    total = sum(new_ones) + sum(new_seps)
    diff = target_n - total
    if diff != 0:
        large_indices = [i for i in range(len(new_ones)) if not is_counter_pos(i)]
        if not large_indices:
            large_indices = list(range(len(new_ones)))
        m = len(large_indices)
        if diff > 0:
            per_block = diff // m
            leftover = diff % m
            for rank, bi in enumerate(large_indices):
                new_ones[bi] += per_block + (1 if rank < leftover else 0)
        else:
            to_remove = abs(diff)
            per_block = to_remove // m
            leftover = to_remove % m
            for rank, bi in enumerate(large_indices):
                new_ones[bi] = max(1, new_ones[bi] - per_block - (1 if rank < leftover else 0))

    # ── Build bitstring using numpy (fast for large n) ────────────────────
    result = _assemble_bitstring(new_ones, new_seps, target_n)
    ones = int(np.sum(result == 1))
    zeros = int(np.sum(result == 0))
    print(f"Generated seed: n={len(result)}, ones={ones}, zeros={zeros}, "
          f"density={ones/len(result):.6f}, blocks={len(new_ones)}+{len(new_seps)}")
    return result


def _assemble_bitstring(ones_blocks: List[int], seps: List[int], n: int) -> np.ndarray:
    """Build an n-bit numpy array from ones-block sizes and separator sizes.
    Uses numpy slicing instead of Python list.extend — ~100x faster at n=10M.
    """
    result = np.ones(n, dtype=np.int32)
    cursor = 0
    for i, bsize in enumerate(ones_blocks):
        cursor += bsize  # skip over ones (already 1)
        if i < len(seps) and cursor < n:
            end = min(cursor + seps[i], n)
            result[cursor:end] = 0
            cursor = end
    return result


def _build_from_pattern(n: int) -> np.ndarray:
    """Build a seed from scratch using the known structural pattern from
    exact solutions: [big, big, big, 1, med, 2, med, 3, med, 4, ...]
    with separator prefix [1, 2, 1, 1, 1, ...].

    The good-quality n=100M optimised seed shows all large blocks in a
    gently decreasing arithmetic sequence (8382, 8381, 8380, …) with
    counter blocks interleaved (1, 2, 3, …). Separators are almost all 1.
    This function replicates that pattern exactly.
    """
    # Target zeros ~ 2.3 * sqrt(n), with sep=1 dominant
    target_zeros = int(2.3 * (n ** 0.5))
    print(f"Building seed from structural pattern: ~{target_zeros} zeros for n={n}")

    # Separator pattern: common prefix then mostly 1s
    common_seps_prefix = [1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3]
    seps = list(common_seps_prefix)
    remaining_seps = target_zeros - len(seps)
    for _ in range(max(0, remaining_seps)):
        r = random.random()
        if r < 0.975:
            seps.append(1)
        elif r < 0.995:
            seps.append(2)
        else:
            seps.append(3)

    total_sep_bits = sum(seps)
    remaining_ones = n - total_sep_bits

    # Cap seps so counter blocks (1+2+...+k) don't exceed ~30% of remaining_ones.
    import math
    max_counter_budget = int(0.30 * remaining_ones)
    max_counters = int((-1 + math.sqrt(1 + 8 * max_counter_budget)) / 2)
    max_blocks = 3 + max_counters * 2
    if len(seps) + 1 > max_blocks:
        seps = seps[:max_blocks - 1]
        total_sep_bits = sum(seps)
        remaining_ones = n - total_sep_bits

    num_blocks = len(seps) + 1

    # Identify which block indices are counters vs large blocks.
    # Pattern: [large, large, large, counter=1, large, counter=2, large, ...]
    # Counter positions: indices 3, 5, 7, 9, ... (odd indices >= 3)
    def is_counter(i: int) -> bool:
        return i >= 3 and i % 2 == 1

    def counter_val(i: int) -> int:
        return (i - 3) // 2 + 1

    counter_total = sum(counter_val(i) for i in range(num_blocks) if is_counter(i))
    num_large = sum(1 for i in range(num_blocks) if not is_counter(i))
    large_budget = remaining_ones - counter_total

    # Build large blocks as a gently decreasing arithmetic sequence,
    # mirroring the optimised n=100M seed structure.
    # Sequence: first, first-1, first-2, ... (step of 1)
    # Sum = k*first - k*(k-1)/2 = large_budget
    # => first = (large_budget + k*(k-1)/2) / k
    k = num_large
    first_val = (large_budget + k * (k - 1) // 2) / k
    first_val = int(round(first_val))

    large_sizes = [max(1, first_val - j) for j in range(k)]

    # Adjust sum to exactly match budget — O(k)
    current = sum(large_sizes)
    diff = large_budget - current
    if diff > 0:
        per = diff // k
        leftover = diff % k
        for i in range(k):
            large_sizes[i] += per + (1 if i < leftover else 0)
    elif diff < 0:
        to_remove = abs(diff)
        per = to_remove // k
        leftover = to_remove % k
        for i in range(k):
            large_sizes[i] = max(1, large_sizes[i] - per - (1 if i < leftover else 0))

    # Assemble blocks: interleave large + counter
    blocks: List[int] = []
    large_iter = iter(large_sizes)
    for i in range(num_blocks):
        if is_counter(i):
            blocks.append(counter_val(i))
        else:
            blocks.append(next(large_iter, 1))

    # Final length adjustment — O(num_blocks)
    total = sum(blocks) + total_sep_bits
    diff = n - total
    if diff != 0:
        large_idx = [i for i in range(len(blocks)) if not is_counter(i)]
        if not large_idx:
            large_idx = list(range(len(blocks)))
        m = len(large_idx)
        if diff > 0:
            per = diff // m
            leftover = diff % m
            for rank, bi in enumerate(large_idx):
                blocks[bi] += per + (1 if rank < leftover else 0)
        else:
            to_remove = abs(diff)
            per = to_remove // m
            leftover = to_remove % m
            for rank, bi in enumerate(large_idx):
                blocks[bi] = max(1, blocks[bi] - per - (1 if rank < leftover else 0))

    # Build bitstring using numpy (fast for large n)
    result = _assemble_bitstring(blocks, seps, n)
    ones = int(np.sum(result == 1))
    zeros = int(np.sum(result == 0))
    print(f"Generated seed: n={len(result)}, ones={ones}, zeros={zeros}, "
          f"density={ones/len(result):.6f}, blocks={len(blocks)}+{len(seps)}")
    return result


def load_seed(n: int, seeds_dir: Path, results_dir: Path) -> Optional[np.ndarray]:
    """Load best available bitstring for this n.
    Tries numpy binary (.npy) first, then text seed, then JSON results.
    """
    # 1. numpy binary (fast for large n)
    npy_file = seeds_dir / f"n_{n}_seed.npy"
    if npy_file.exists():
        try:
            bits = np.load(str(npy_file))
            if len(bits) == n:
                print(f"Loaded seed from {npy_file} (numpy binary)")
                return bits.astype(np.int32)
        except Exception:
            pass
    # 2. text seed (legacy / small n)
    seed_file = seeds_dir / f"n_{n:04d}_seed.txt"
    if seed_file.exists():
        try:
            txt = seed_file.read_text().strip()
            if len(txt) == n:
                bits = np.frombuffer(txt.encode('ascii'), dtype=np.uint8) - ord('0')
                if np.all((bits == 0) | (bits == 1)):
                    print(f"Loaded seed from {seed_file}")
                    return bits.astype(np.int32)
        except Exception:
            pass
    # 3. JSON results
    result_file = results_dir / f"n_{n:04d}_results.json"
    if result_file.exists():
        try:
            data = json.loads(result_file.read_text())
            sols = data.get("solutions", [])
            if sols:
                txt = sols[0]
                if len(txt) == n:
                    bits = np.frombuffer(txt.encode('ascii'), dtype=np.uint8) - ord('0')
                    if np.all((bits == 0) | (bits == 1)):
                        print(f"Loaded seed from {result_file} (verified result)")
                        return bits.astype(np.int32)
        except Exception:
            pass
    return None

def save_result(
    bits: np.ndarray,
    value: int,
    n: int,
    seeds_dir: Path,
    results_dir: Path,
):
    """Persist improved bitstring.
    Always writes seeds/ (.npy binary for fast I/O at any n).
    Only updates results/ JSON when strictly better.
    Skips integrity recompute for large n (too expensive).
    """
    # 1. Always save seed as numpy binary
    npy_file = seeds_dir / f"n_{n}_seed.npy"
    np.save(str(npy_file), bits)
    # 2. Update results only if strictly better
    result_file = results_dir / f"n_{n:04d}_results.json"
    existing_best = 0
    if result_file.exists():
        try:
            existing_best = json.loads(result_file.read_text()).get("best_value", 0)
        except Exception:
            pass
    if value <= existing_best:
        return
    # For n > 1M, don't store the full bitstring in JSON (too large).
    # Store just the score and metadata; the seed .npy file has the bits.
    if n > 1_000_000:
        data = {
            "n": n,
            "best_value": value,
            "solutions": [],  # full bits in .npy file
            "seed_file": str(npy_file),
            "last_updated": datetime.now().isoformat(),
            "metadata": {"algorithm": "surgical_nudge"},
        }
    else:
        solution_str = (bits.astype(np.uint8) + ord('0')).tobytes().decode('ascii')
        data = {
            "n": n,
            "best_value": value,
            "solutions": [solution_str],
            "last_updated": datetime.now().isoformat(),
            "metadata": {"algorithm": "surgical_nudge"},
        }
    result_file.write_text(json.dumps(data, indent=2))
    print(f"  New verified best saved → {result_file} (value={value})")
# ---------------------------------------------------------------------------
# Main optimisation loop
# ---------------------------------------------------------------------------

def optimize(
    n: int,
    max_seconds: float = 3600.0,
    top_k_collisions: int = 30,
    random_fallback_swaps: int = 200,
    zero_migration_swaps: int = 200,
    zero_count_adjust_tries: int = 500,
    zero_count_interleave_every: int = 5,
    zero_count_interleave_tries: int = 15,
    targeted_relocation_tries: int = 30,
    random_swap_interleave_every: int = 3,
    random_swap_interleave_tries: int = 50,
    multi_swap_kick_every: int = 10,
    multi_swap_kick_size: int = 3,
    multi_swap_kick_tries: int = 30,
    structural_move_every: int = 5,
    structural_move_tries: int = 20,
    stagnation_threshold: int = 50,
    save_interval_seconds: float = 300.0,
    progress_interval_seconds: float = 30.0,
    num_workers: int = 1,
    seeds_dir: Path = Path("seeds"),
    results_dir: Path = Path("results"),
) -> Tuple[np.ndarray, int]:
    """Run the surgical-nudge hill climber.
    Parameters
    ----------
    n                          : string length
    max_seconds                : total wall-clock budget
    top_k_collisions           : how many worst collisions to inspect per step
    random_fallback_swaps      : random swaps to try when surgical nudge stalls
    zero_migration_swaps       : body<->tail zero relocations tried when surgical stalls
    zero_count_adjust_tries    : single bit-flips tried to change zero count when deep stall
    zero_count_interleave_every: run a speculative zero-count probe every N iterations
                                 (regardless of stagnation), 0 to disable
    zero_count_interleave_tries: max SA builds per interleaved probe (small = low overhead)
    targeted_relocation_tries  : per-iteration budget for targeted relocation (Phase 1b):
                                 zeros from collision window swapped with any 1 anywhere
    random_swap_interleave_every: run random long-range 0<->1 swaps every N iterations
                                 regardless of surgical progress (0 to disable);
                                 mirrors the MH's ~25%% random proposal rate
    random_swap_interleave_tries: SA builds budget per interleaved random swap probe
    multi_swap_kick_every      : every N iterations, try a multi-swap perturbation
                                 (flip multiple bits simultaneously for deeper exploration)
    multi_swap_kick_size       : number of simultaneous 0<->1 swaps per perturbation
    multi_swap_kick_tries      : how many multi-swap candidates to evaluate per kick
    structural_move_every      : every N iterations, try RLE-level structural moves
                                 (block split, merge, separator relocation, run swap)
    structural_move_tries      : candidates to try per structural move type
    stagnation_threshold       : stagnant iters before attempting zero-count adjustment
    save_interval_seconds      : checkpoint interval (seed always; results/ if new best)
    """
    seeds_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)
    # --- Initialize process pool for parallel candidate evaluation ---
    _init_pool(num_workers, n=n)
    # --- Load initial solution ---
    bits = load_seed(n, seeds_dir, results_dir)
    if bits is None:
        print(f"No seed found for n={n}; generating structured seed...")
        bits = generate_seed(n, results_dir)
    # --- Initial evaluation ---
    sa_method = "pydivsufsort (C)" if _DIVSUFSORT_AVAILABLE else "numpy prefix-doubling"
    print(f"Building initial SA+LCP for n={n}…  [{sa_method}]")
    sa, lcp = build_sa_lcp(bits)
    best_value = score_from_sa_lcp(bits, sa, lcp)
    best_bits = bits.copy()
    n_sq = float(n) * float(n)
    print(f"Starting score: {best_value}  (a/n²={best_value/n_sq:.5f})")
    print(f"Max LCP (worst collision): {int(lcp.max())}  at rank {int(lcp.argmax())}")
    print(f"Budget: {max_seconds:.0f}s\n")
    start_time = time.time()
    last_save_time = start_time
    last_progress_time = start_time
    iteration = 0
    surgical_improvements = 0
    targeted_improvements = 0
    migration_improvements = 0
    random_improvements = 0
    swap_interleave_improvements = 0
    multi_kick_improvements = 0
    structural_improvements = 0
    zcount_improvements = 0
    stagnant_steps = 0
    # Partition boundary for zero-migration moves (matches build_seed.py)
    tail_boundary = round(0.85 * n)
    try:
        while True:
            elapsed = time.time() - start_time
            if elapsed >= max_seconds:
                print(f"\nTime budget exhausted after {elapsed:.1f}s")
                break
            iteration += 1
            assert sa is not None and lcp is not None  # always set; helps Pylance
            # ------------------------------------------------------------------
            # Recovery: if SA-accepted structural moves have degraded `bits`
            # below best_value for too long, snap back to best_bits so we
            # don't waste time polishing a dead-end structure.
            # ------------------------------------------------------------------
            current_value = score_from_sa_lcp(bits, sa, lcp)
            if current_value < best_value and stagnant_steps > 0 and stagnant_steps % 20 == 0:
                bits = best_bits.copy()
                sa, lcp = build_sa_lcp(bits)
                print(f"  Iter {iteration}: Recovery — snapping back to best ({best_value})")
            # ------------------------------------------------------------------
            # Phase 0: Speculative zero-count probe (interleaved, low overhead).
            # Every zero_count_interleave_every iterations, try a handful of
            # single bit-flips (add or remove one zero) to see if a different
            # zero count is immediately better.  Uses first-improvement so it
            # rarely costs more than 1-3 SA builds.  Fires even while surgical
            # is making progress, so the zero count can drift toward optimal
            # without waiting for stagnation.
            # ------------------------------------------------------------------
            if (
                zero_count_interleave_every > 0
                and zero_count_interleave_tries > 0
                and iteration % zero_count_interleave_every == 0
            ):
                half = max(1, zero_count_interleave_tries // 2)
                # Generate batch of swap specs: half add-zero, half remove-zero
                zp_specs: List[List[Tuple[int, int]]] = []
                ones_idx = np.where(bits == 1)[0]
                zeros_idx = np.where(bits == 0)[0]
                for _ in range(half):
                    if len(ones_idx) == 0:
                        break
                    j = int(ones_idx[np.random.randint(len(ones_idx))])
                    zp_specs.append([(j, 0)])
                for _ in range(half):
                    if len(zeros_idx) == 0:
                        break
                    j = int(zeros_idx[np.random.randint(len(zeros_idx))])
                    zp_specs.append([(j, 1)])
                if zp_specs:
                    zp_val, zp_idx, zp_sa, zp_lcp = parallel_eval_best(
                        bits, zp_specs, threshold=best_value)
                    if zp_idx >= 0 and zp_sa is not None and zp_lcp is not None:
                        for pos, val in zp_specs[zp_idx]:
                            bits[pos] = val
                        sa, lcp = zp_sa, zp_lcp
                        best_value = zp_val
                        best_bits = bits.copy()
                        zcount_improvements += 1
                        stagnant_steps = 0
                        new_zcount = int(np.sum(bits == 0))
                        print(
                            f"  Iter {iteration}: NEW BEST = {best_value}"
                            f"  [zero-count probe #{zcount_improvements}, zeros={new_zcount}]"
                            f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                        )
                        now = time.time()
                        if now - last_save_time >= save_interval_seconds:
                            save_result(best_bits, best_value, n, seeds_dir, results_dir)
                            last_save_time = now
                        last_progress_time = now
                        continue
            # ------------------------------------------------------------------
            # Phase 0.5: Interleaved random long-range 0↔1 swaps.
            #
            # Mirrors the MH's ~25% random proposal rate.  At n=2M the top-K
            # collision windows cover <1% of positions, so most beneficial
            # swaps are invisible to surgical phases.  Running a small batch
            # of fully random swaps every few iterations maintains exploration
            # diversity even while surgical moves are succeeding.
            # ------------------------------------------------------------------
            if (
                random_swap_interleave_every > 0
                and random_swap_interleave_tries > 0
                and iteration % random_swap_interleave_every == 0
            ):
                rs_zeros = np.where(bits == 0)[0]
                rs_ones = np.where(bits == 1)[0]
                if len(rs_zeros) > 0 and len(rs_ones) > 0:
                    # Generate batch of swap specs for parallel eval
                    ri_specs: List[List[Tuple[int, int]]] = []
                    for _ in range(random_swap_interleave_tries):
                        j0 = int(rs_zeros[np.random.randint(len(rs_zeros))])
                        j1 = int(rs_ones[np.random.randint(len(rs_ones))])
                        ri_specs.append([(j0, int(bits[j1])), (j1, int(bits[j0]))])
                    ri_val, ri_idx, ri_sa, ri_lcp = parallel_eval_best(
                        bits, ri_specs, threshold=best_value)
                    if ri_idx >= 0 and ri_sa is not None and ri_lcp is not None:
                        for pos, val in ri_specs[ri_idx]:
                            bits[pos] = val
                        sa, lcp = ri_sa, ri_lcp
                        best_value = ri_val
                        best_bits = bits.copy()
                        swap_interleave_improvements += 1
                        stagnant_steps = 0
                        print(
                            f"  Iter {iteration}: NEW BEST = {best_value}"
                            f"  [random interleave #{swap_interleave_improvements}]"
                            f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                        )
                        now = time.time()
                        if now - last_save_time >= save_interval_seconds:
                            save_result(best_bits, best_value, n, seeds_dir, results_dir)
                            last_save_time = now
                        last_progress_time = now
            # ------------------------------------------------------------------
            # Phase 0.75: Multi-swap perturbation kick.
            #
            # Inspired by MH's double-pair cross-swap.  Flip multiple 0↔1
            # pairs simultaneously to reach neighbor states unreachable by
            # single swaps (combinatorial "tunneling").  Fires periodically
            # regardless of stagnation for deeper exploration.
            # ------------------------------------------------------------------
            if (
                multi_swap_kick_every > 0
                and multi_swap_kick_tries > 0
                and iteration % multi_swap_kick_every == 0
            ):
                mk_zeros = np.where(bits == 0)[0]
                mk_ones = np.where(bits == 1)[0]
                k_size = min(multi_swap_kick_size, len(mk_zeros), len(mk_ones))
                if k_size >= 2:
                    for _ in range(multi_swap_kick_tries):
                        # Pick k_size random 0↔1 pairs and swap them all
                        z_idx = np.random.choice(len(mk_zeros), size=k_size, replace=False)
                        o_idx = np.random.choice(len(mk_ones), size=k_size, replace=False)
                        swap_pairs = [(int(mk_zeros[z_idx[i]]), int(mk_ones[o_idx[i]]))
                                      for i in range(k_size)]
                        # Apply all swaps
                        for j0, j1 in swap_pairs:
                            _swap_inplace(bits, j0, j1)
                        csa, clcp = build_sa_lcp(bits)
                        cval = score_from_sa_lcp(bits, csa, clcp)
                        if cval > best_value:
                            sa, lcp = csa, clcp
                            best_value = cval
                            best_bits = bits.copy()
                            multi_kick_improvements += 1
                            stagnant_steps = 0
                            print(
                                f"  Iter {iteration}: NEW BEST = {best_value}"
                                f"  [multi-kick #{multi_kick_improvements}, size={k_size}]"
                                f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                            )
                            now = time.time()
                            if now - last_save_time >= save_interval_seconds:
                                save_result(best_bits, best_value, n, seeds_dir, results_dir)
                                last_save_time = now
                            last_progress_time = now
                            break
                        # Revert all swaps (reverse order)
                        for j0, j1 in reversed(swap_pairs):
                            _swap_inplace(bits, j0, j1)
            # ------------------------------------------------------------------
            # Phase 0.9: RLE-level structural moves with SA acceptance.
            #
            # These make multi-bit correlated changes that single swaps can
            # never achieve: splitting large ones-blocks to insert small
            # "punctuation" blocks (1, 2, 3, ...), merging tiny blocks,
            # relocating separators, and swapping run lengths between blocks.
            #
            # Unlike bit-level phases (strict hill climb), structural moves use
            # simulated annealing: worse solutions accepted with probability
            # exp(-delta/T).  This lets the optimizer pass through temporarily
            # worse states to reach the correct block structure.  The current
            # (possibly worse) solution is tracked in `bits`, while `best_bits`
            # always holds the best-ever.
            # ------------------------------------------------------------------
            if (
                structural_move_every > 0
                and structural_move_tries > 0
                and iteration % structural_move_every == 0
            ):
                # Temperature schedule: starts at ~1% of score, decays over time
                progress_frac = min(1.0, elapsed / max_seconds)
                struct_temp = best_value * 0.01 * (1.0 - progress_frac) + 1.0

                struct_candidates = []
                for move_fn, label in [
                    (structural_block_split, "block-split"),
                    (structural_block_merge, "block-merge"),
                    (structural_separator_relocate, "sep-relocate"),
                    (structural_run_swap, "run-swap"),
                ]:
                    candidate = move_fn(bits, tries=structural_move_tries)
                    if candidate is not None:
                        struct_candidates.append((candidate, label))

                if struct_candidates:
                    # Evaluate all, pick the best candidate
                    current_value = score_from_sa_lcp(bits, sa, lcp)
                    struct_best_val = -1
                    struct_best_bits = None
                    struct_best_sa = None
                    struct_best_lcp = None
                    struct_label = ""
                    for cand_bits, label in struct_candidates:
                        csa, clcp = build_sa_lcp(cand_bits)
                        cval = score_from_sa_lcp(cand_bits, csa, clcp)
                        if cval > struct_best_val:
                            struct_best_val = cval
                            struct_best_bits = cand_bits
                            struct_best_sa = csa
                            struct_best_lcp = clcp
                            struct_label = label

                    delta = struct_best_val - current_value
                    # Accept if better, or with SA probability if worse
                    accept = False
                    if delta > 0:
                        accept = True
                    elif struct_temp > 0:
                        prob = np.exp(delta / struct_temp)
                        accept = random.random() < prob

                    if accept and struct_best_bits is not None and struct_best_sa is not None and struct_best_lcp is not None:
                        bits = struct_best_bits
                        sa, lcp = struct_best_sa, struct_best_lcp
                        if struct_best_val > best_value:
                            best_value = struct_best_val
                            best_bits = bits.copy()
                            structural_improvements += 1
                            stagnant_steps = 0
                            print(
                                f"  Iter {iteration}: NEW BEST = {best_value}"
                                f"  [structural #{structural_improvements} ({struct_label})]"
                                f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                            )
                            now = time.time()
                            if now - last_save_time >= save_interval_seconds:
                                save_result(best_bits, best_value, n, seeds_dir, results_dir)
                                last_save_time = now
                            last_progress_time = now
                        else:
                            structural_improvements += 1
                            print(
                                f"  Iter {iteration}: SA-accept structural ({struct_label})"
                                f"  delta={delta:+d}  T={struct_temp:.0f}"
                                f"  current={struct_best_val}  best={best_value}"
                            )
                        continue
            # ------------------------------------------------------------------
            # Phase 1: Surgical nudge — target the top-K worst collisions
            # Phase 1a: local ±2 swap within collision window (fast)
            # Phase 1b: targeted relocation — zeros from collision window moved
            #           to any 1 anywhere in the string (no distance limit)
            # ------------------------------------------------------------------
            collision_indices = find_top_collision_indices(lcp, k=top_k_collisions)
            ones_all = np.where(bits == 1)[0]  # reused across 1b
            improved = False
            for lcp_idx in collision_indices:
                if lcp_idx == 0:
                    continue
                if int(lcp[lcp_idx]) == 0:
                    break
                # --- Phase 1a: \u00b12 local swap ---
                candidates = surgical_candidates(bits, sa, lcp, lcp_idx)
                if candidates:
                    random.shuffle(candidates)
                    for j, d in candidates:
                        jd = j + d
                        _swap_inplace(bits, j, jd)
                        csa, clcp = build_sa_lcp(bits)
                        cval = score_from_sa_lcp(bits, csa, clcp)
                        if cval > best_value:
                            sa, lcp = csa, clcp
                            best_value = cval
                            best_bits = bits.copy()
                            surgical_improvements += 1
                            stagnant_steps = 0
                            improved = True
                            print(
                                f"  Iter {iteration}: NEW BEST = {best_value}"
                                f"  [surgical #{surgical_improvements}]"
                                f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                            )
                            break
                        else:
                            _swap_inplace(bits, j, jd)
                if improved:
                    break
                # --- Phase 1b: targeted relocation (any-distance swap) ---
                # Collect zeros in the collision window (same region as 1a)
                if targeted_relocation_tries <= 0 or len(ones_all) == 0:
                    continue
                L = int(lcp[lcp_idx])
                p, q = int(sa[lcp_idx - 1]), int(sa[lcp_idx])
                window_zeros = []
                for origin in (p, q):
                    lo = max(0, origin + max(0, L - 500))
                    hi = min(n, origin + L + 5)
                    for jj in range(lo, hi):
                        if bits[jj] == 0:
                            window_zeros.append(jj)
                if not window_zeros:
                    continue
                random.shuffle(window_zeros)
                tries_left = targeted_relocation_tries
                for j0 in window_zeros:
                    if tries_left <= 0:
                        break
                    # sample a random 1 from anywhere
                    j1 = int(ones_all[np.random.randint(len(ones_all))])
                    if j1 == j0:
                        continue
                    tries_left -= 1
                    _swap_inplace(bits, j0, j1)
                    csa, clcp = build_sa_lcp(bits)
                    cval = score_from_sa_lcp(bits, csa, clcp)
                    if cval > best_value:
                        sa, lcp = csa, clcp
                        best_value = cval
                        best_bits = bits.copy()
                        targeted_improvements += 1
                        stagnant_steps = 0
                        improved = True
                        ones_all = np.where(bits == 1)[0]  # refresh after change
                        print(
                            f"  Iter {iteration}: NEW BEST = {best_value}"
                            f"  [targeted reloc #{targeted_improvements}]"
                            f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                        )
                        break
                    else:
                        _swap_inplace(bits, j0, j1)
                if improved:
                    break
            now = time.time()
            if improved:
                # Checkpoint if enough time has passed
                if now - last_save_time >= save_interval_seconds:
                    save_result(best_bits, best_value, n, seeds_dir, results_dir)
                    last_save_time = now
                last_progress_time = now  # reset progress timer on improvement
                continue
            # ------------------------------------------------------------------
            # Phase 1.5: Zero-migration — relocate zeros between body and tail.
            #
            # Small-n optimal solutions have ~37% of zeros in the last 15% of
            # the string (escalating run lengths toward the end).  Surgical
            # moves (delta ±5) can never bridge the body/tail gap, so we add
            # dedicated long-range swaps: body-zero ↔ tail-one and vice versa.
            # ------------------------------------------------------------------
            stagnant_steps += 1
            # Compute hot-position mask: zeros inside active collision windows
            # are preferentially sampled in migration and random phases.
            hot = _hot_positions(sa, lcp, n, k=top_k_collisions)
            body_zeros = np.where(bits[:tail_boundary] == 0)[0]
            tail_ones  = np.where(bits[tail_boundary:] == 1)[0] + tail_boundary
            tail_zeros = np.where(bits[tail_boundary:] == 0)[0] + tail_boundary
            body_ones  = np.where(bits[:tail_boundary] == 1)[0]
            # Prefer hot zeros when available
            body_hot_zeros = body_zeros[hot[body_zeros]] if len(body_zeros) else body_zeros
            tail_hot_zeros = tail_zeros[hot[tail_zeros]] if len(tail_zeros) else tail_zeros
            body_zeros_src = body_hot_zeros if len(body_hot_zeros) > 0 else body_zeros
            tail_zeros_src = tail_hot_zeros if len(tail_hot_zeros) > 0 else tail_zeros

            # --- Phase 1.5: Zero-migration (parallel batch) ---
            mig_specs: List[List[Tuple[int, int]]] = []
            if len(body_zeros) > 0 and len(tail_ones) > 0 and zero_migration_swaps > 0:
                n_each = zero_migration_swaps // 2
                for _ in range(n_each):
                    j0 = int(body_zeros_src[np.random.randint(len(body_zeros_src))])
                    j1 = int(tail_ones[np.random.randint(len(tail_ones))])
                    mig_specs.append([(j0, int(bits[j1])), (j1, int(bits[j0]))])
                if len(tail_zeros) > 0 and len(body_ones) > 0:
                    for _ in range(n_each):
                        j0 = int(tail_zeros_src[np.random.randint(len(tail_zeros_src))])
                        j1 = int(body_ones[np.random.randint(len(body_ones))])
                        mig_specs.append([(j0, int(bits[j1])), (j1, int(bits[j0]))])
            if mig_specs:
                mig_val, mig_idx, mig_sa, mig_lcp = parallel_eval_best(
                    bits, mig_specs, threshold=best_value)
                if mig_idx >= 0 and mig_sa is not None and mig_lcp is not None:
                    for pos, val in mig_specs[mig_idx]:
                        bits[pos] = val
                    sa, lcp = mig_sa, mig_lcp
                    best_value = mig_val
                    best_bits = bits.copy()
                    migration_improvements += 1
                    stagnant_steps = 0
                    print(
                        f"  Iter {iteration}: NEW BEST = {best_value}"
                        f"  [migration #{migration_improvements}]"
                        f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                    )
                    now = time.time()
                    if now - last_save_time >= save_interval_seconds:
                        save_result(best_bits, best_value, n, seeds_dir, results_dir)
                        last_save_time = now
                    last_progress_time = now
                    continue

            # --- Phase 2: random swap (parallel batch) ---
            ones_pos  = np.where(bits == 1)[0]
            zeros_pos = np.where(bits == 0)[0]
            if len(ones_pos) > 0 and len(zeros_pos) > 0:
                hot_zeros = zeros_pos[hot[zeros_pos]]
                zeros_src = hot_zeros if len(hot_zeros) > 0 else zeros_pos
                rnd_specs: List[List[Tuple[int, int]]] = []
                for _ in range(random_fallback_swaps):
                    j0 = int(zeros_src[np.random.randint(len(zeros_src))])
                    j1 = int(ones_pos[np.random.randint(len(ones_pos))])
                    rnd_specs.append([(j0, int(bits[j1])), (j1, int(bits[j0]))])
                rnd_val, rnd_idx, rnd_sa, rnd_lcp = parallel_eval_best(
                    bits, rnd_specs, threshold=best_value)
                if rnd_idx >= 0 and rnd_sa is not None and rnd_lcp is not None:
                    for pos, val in rnd_specs[rnd_idx]:
                        bits[pos] = val
                    sa, lcp = rnd_sa, rnd_lcp
                    best_value = rnd_val
                    best_bits = bits.copy()
                    random_improvements += 1
                    stagnant_steps = 0
                    print(
                        f"  Iter {iteration}: NEW BEST = {best_value}"
                        f"  [random kick #{random_improvements}]"
                        f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                    )
                    now = time.time()
                    if now - last_save_time >= save_interval_seconds:
                        save_result(best_bits, best_value, n, seeds_dir, results_dir)
                        last_save_time = now
                    last_progress_time = now

            # --- Phase 3: Zero-count adjustment (parallel batch) ---
            if stagnant_steps >= stagnation_threshold and zero_count_adjust_tries > 0:
                zca_specs: List[List[Tuple[int, int]]] = []
                ones_idx = np.where(bits == 1)[0]
                zeros_idx = np.where(bits == 0)[0]
                for _ in range(zero_count_adjust_tries // 2):
                    j = int(ones_idx[np.random.randint(len(ones_idx))])
                    zca_specs.append([(j, 0)])
                for _ in range(zero_count_adjust_tries // 2):
                    if len(zeros_idx) == 0:
                        break
                    j = int(zeros_idx[np.random.randint(len(zeros_idx))])
                    zca_specs.append([(j, 1)])
                if zca_specs:
                    zca_val, zca_idx, zca_sa, zca_lcp = parallel_eval_best(
                        bits, zca_specs, threshold=best_value)
                    if zca_idx >= 0 and zca_sa is not None and zca_lcp is not None:
                        for pos, val in zca_specs[zca_idx]:
                            bits[pos] = val
                        sa, lcp = zca_sa, zca_lcp
                        best_value = zca_val
                        best_bits = bits.copy()
                        zcount_improvements += 1
                        stagnant_steps = 0
                        new_zcount = int(np.sum(bits == 0))
                        print(
                            f"  Iter {iteration}: NEW BEST = {best_value}"
                            f"  [zero-count adjust #{zcount_improvements}, zeros now={new_zcount}]"
                            f"  max_lcp={int(lcp.max())}  a/n²={best_value/n_sq:.5f}"
                        )
                        now = time.time()
                        if now - last_save_time >= save_interval_seconds:
                            save_result(best_bits, best_value, n, seeds_dir, results_dir)
                            last_save_time = now
                        last_progress_time = now
                        continue
            # ------------------------------------------------------------------
            # Regular progress line (time-based, like the greedy search)
            # ------------------------------------------------------------------
            now = time.time()
            if now - last_progress_time >= progress_interval_seconds:
                print(
                    f"  Iter {iteration}: best={best_value}, stagnant={stagnant_steps},"
                    f" surgical={surgical_improvements}, random_kicks={random_improvements},"
                    f" max_lcp={int(lcp.max())}, a/n²={best_value/n_sq:.5f}, elapsed={now - start_time:.0f}s"
                )
                last_progress_time = now
    except KeyboardInterrupt:
        print(f"\n[Interrupted at iter {iteration}] Saving best so far: {best_value}")
    # --- Final save + cleanup ---
    save_result(best_bits, best_value, n, seeds_dir, results_dir)
    _shutdown_pool()
    total_elapsed = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"DONE  n={n}  score={best_value}")
    print(f"Surgical improvements (\u00b12):      {surgical_improvements}")
    print(f"Targeted reloc improvements:     {targeted_improvements}")
    print(f"Random interleave improvements:  {swap_interleave_improvements}")
    print(f"Multi-kick improvements:         {multi_kick_improvements}")
    print(f"Structural improvements:         {structural_improvements}")
    print(f"Migration improvements:          {migration_improvements}")
    print(f"Random kick improvements:        {random_improvements}")
    print(f"Zero-count improvements:         {zcount_improvements}")
    print(f"Final zero count:            {int(np.sum(best_bits == 0))}")
    print(f"Total iterations:            {iteration}")
    print(f"Elapsed:               {total_elapsed:.1f}s")
    print(f"{'='*70}")
    return best_bits, best_value
# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # ── CONFIGURATION ──────────────────────────────────────────────────────
    n                      = 100_000_000  # 100M — ~7s per SA build, ~1.2GB per worker
    max_seconds            = 36000      # wall-clock budget in seconds (10 hours)
    num_workers            = 14         # parallel SA evaluations (leave 2 of 16 logical cores for OS)
    top_k_collisions       = 30         # worst collisions to inspect per step
    zero_migration_swaps   = 50         # body<->tail zero relocation candidates per batch
    random_fallback_swaps  = 50         # random swap candidates per batch
    zero_count_adjust_tries = 100       # bit flip candidates when stuck
    zero_count_interleave_every = 1     # run speculative zero-count probe every N iters
    zero_count_interleave_tries = 20    # candidates per interleaved probe batch
    targeted_relocation_tries  = 10     # candidates for Phase 1b (any-distance swap)
    random_swap_interleave_every = 3    # random long-range swaps every N iters (MH-style)
    random_swap_interleave_tries = 20   # candidates per random interleave batch
    multi_swap_kick_every      = 10     # multi-swap perturbation every N iters
    multi_swap_kick_size       = 3      # simultaneous 0<->1 swaps per kick
    multi_swap_kick_tries      = 10     # candidates to evaluate per multi-kick
    structural_move_every  = 5          # RLE-level structural moves every N iters
    structural_move_tries  = 10         # candidates per structural move type
    stagnation_threshold   = 50         # stagnant iters before zero-count adjustment fires
    save_interval_seconds  = 120        # checkpoint every 2 minutes
    progress_interval_seconds = 60      # status line every minute
    # ───────────────────────────────────────────────────────────────────────
    seeds_dir   = Path("seeds")
    results_dir = Path(os.path.join(os.path.dirname(__file__), "..", "..", "results", "large_n"))
    optimize(
        n=n,
        max_seconds=max_seconds,
        top_k_collisions=top_k_collisions,
        zero_migration_swaps=zero_migration_swaps,
        random_fallback_swaps=random_fallback_swaps,
        zero_count_adjust_tries=zero_count_adjust_tries,
        zero_count_interleave_every=zero_count_interleave_every,
        zero_count_interleave_tries=zero_count_interleave_tries,
        targeted_relocation_tries=targeted_relocation_tries,
        random_swap_interleave_every=random_swap_interleave_every,
        random_swap_interleave_tries=random_swap_interleave_tries,
        multi_swap_kick_every=multi_swap_kick_every,
        multi_swap_kick_size=multi_swap_kick_size,
        multi_swap_kick_tries=multi_swap_kick_tries,
        structural_move_every=structural_move_every,
        structural_move_tries=structural_move_tries,
        stagnation_threshold=stagnation_threshold,
        save_interval_seconds=save_interval_seconds,
        progress_interval_seconds=progress_interval_seconds,
        num_workers=num_workers,
        seeds_dir=seeds_dir,
        results_dir=results_dir,
    )
if __name__ == "__main__":
    main()
