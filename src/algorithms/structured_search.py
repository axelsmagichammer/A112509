"""
Compute A112509(n) using pre-analyzed structural bounds.

Requires:
    - Pre-analyzed bounds in config/learned_bounds.json (use analyze_structure.py)
    - Cached results in data/cached_results.json

Algorithm:
    For each n, loads learned structural bounds (K common blocks + separators)
    and enumerates all permutations of block sizes, then all valid tail combinations.
"""

import itertools
import json
import math
import os
import sys
import threading
import time
import traceback
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import freeze_support, cpu_count

try:
    import numpy as np
    import numba
    _NUMBA_OK = True
except ImportError:
    _NUMBA_OK = False


# ---------------------------------------------------------------------------
# Core algorithm (optimized with bit operations)
# ---------------------------------------------------------------------------

def distinct_substring_count_int(num: int, n: int) -> int:
    """Count distinct integer values from all contiguous substrings of an n-bit integer.
    
    Args:
        num: The integer representation of the binary string
        n: Number of bits
    
    Returns:
        Count of distinct substring values
    """
    # Pre-compute masks to avoid repeated computation
    masks = [(1 << length) - 1 for length in range(n + 1)]
    values = set()
    
    # For each starting position
    for i in range(n):
        # For each length from that position
        max_length = n - i
        for length in range(1, max_length + 1):
            # Extract substring: shift right to align, then mask
            shift = max_length - length
            substring_val = (num >> shift) & masks[length]
            values.add(substring_val)
    
    return len(values)


def distinct_substring_count(s: str) -> int:
    """Return the count of distinct integer values from all contiguous substrings.
    Converts string to int and uses optimized bit operations.
    """
    n = len(s)
    num = int(s, 2)
    return distinct_substring_count_int(num, n)


# ---------------------------------------------------------------------------
# Fast evaluator: Numba JIT suffix array (O(n log n)), ~38-750x faster
# Falls back to naive if numba/numpy unavailable.
# ---------------------------------------------------------------------------

if _NUMBA_OK:
    @numba.njit(cache=True)
    def _counting_sort_sa(perm, keys, num_keys, n):
        """Stable counting sort: reorder perm[0..n) so that keys[perm[i]] is
        non-decreasing.  key values must be in [0, num_keys).
        Returns a new int32 array of length n.
        """
        counts = np.zeros(num_keys, dtype=np.int32)
        for i in range(n):
            counts[keys[perm[i]]] += 1
        # convert counts to start positions (prefix sum in-place)
        total = np.int32(0)
        for i in range(num_keys):
            c = counts[i]
            counts[i] = total
            total += c
        out = np.empty(n, dtype=np.int32)
        for i in range(n):
            kv = keys[perm[i]]
            out[counts[kv]] = perm[i]
            counts[kv] += np.int32(1)
        return out

    @numba.njit(cache=True)
    def _numba_eval_filtered(bits, max_zeros):
        """Count distinct substrings via suffix array + Kasai LCP.

        max_zeros: if >= 0, return -1 immediately if the longest zero run exceeds
        this value (early-exit filter at JIT speed, avoids the full SA computation).
        SA prefix-doubling uses O(n) 2-pass radix sort instead of O(n log n) argsort.
        """
        n = len(bits)
        if n == 0:
            return np.int32(0)
        # ── Zero-run early-exit filter ────────────────────────────────────
        if max_zeros >= 0:
            max_run = np.int32(0)
            run = np.int32(0)
            for i in range(n):
                if bits[i] == 0:
                    run += np.int32(1)
                    if run > max_run:
                        max_run = run
                        if max_run > max_zeros:
                            return np.int32(-1)
                else:
                    run = np.int32(0)
        # ── Suffix array (prefix-doubling with 2-pass radix sort) ─────────
        # ranks are always in [0, n-1]; with the +1 sentinel offset they fit
        # in [0, n], so num_keys = n+1 for both counting-sort passes.
        sa   = np.arange(n, dtype=np.int32)
        rank = bits.astype(np.int32)
        tmp  = np.empty(n, dtype=np.int32)
        sec  = np.empty(n, dtype=np.int32)  # secondary key per suffix
        pri  = np.empty(n, dtype=np.int32)  # primary   key per suffix
        num_keys = n + 1
        k = 1
        while k < n:
            for i in range(n):
                r2 = rank[i + k] if i + k < n else np.int32(-1)
                sec[i] = np.int32(r2 + 1)       # sentinel -1 → 0; ranks+1 ∈ [1,n]
                pri[i] = np.int32(rank[i] + 1)  # ranks+1 ∈ [1, n]
            # Pass 1: sort identity by secondary key (least-significant)
            init_perm = np.arange(n, dtype=np.int32)
            sa_tmp = _counting_sort_sa(init_perm, sec, num_keys, n)
            # Pass 2: stable-sort result by primary key (most-significant)
            sa = _counting_sort_sa(sa_tmp, pri, num_keys, n)
            # Update ranks based on (pri, sec) pair comparison
            tmp[sa[0]] = 0
            for i in range(1, n):
                tmp[sa[i]] = tmp[sa[i - 1]]
                if pri[sa[i]] != pri[sa[i - 1]] or sec[sa[i]] != sec[sa[i - 1]]:
                    tmp[sa[i]] += 1
            rank = tmp.copy()
            if rank[sa[n - 1]] == n - 1:
                break
            k *= 2
        # ── Kasai LCP ─────────────────────────────────────────────────────
        pos = np.empty(n, dtype=np.int32)
        for i in range(n):
            pos[sa[i]] = i
        lcp = np.zeros(n, dtype=np.int32)
        h = 0
        for i in range(n):
            if pos[i] > 0:
                j = int(sa[pos[i] - 1])
                while i + h < n and j + h < n and bits[i + h] == bits[j + h]:
                    h += 1
                lcp[pos[i]] = h
                if h > 0:
                    h -= 1
        count = 0
        for i in range(n):
            if bits[sa[i]] == 1:
                count += (n - int(sa[i])) - int(lcp[i])
        for i in range(n):
            if bits[i] == 0:
                count += 1
                break
        return count

    @numba.njit(cache=True)
    def _numba_eval_batch(bits_batch, max_zeros, template_len=0):
        """Evaluate a batch of bit strings via suffix automaton (SAM).

        Builds an online suffix automaton for each string in O(n) time,
        then counts distinct substrings starting with '1' via topological
        path counting on the SAM's transition DAG.

        When template_len > 0, builds the SAM once for the shared template
        prefix (first template_len bits, assumed identical across all rows),
        then clones and extends for each string's unique tail.  This saves
        ~template_len/n of the SAM build work per string.

        bits_batch: 2-D int8 array (batch, n). Returns int32 array of counts.
        """
        m = bits_batch.shape[0]
        n = bits_batch.shape[1]
        out = np.empty(m, dtype=np.int32)

        max_states = 2 * n + 2  # SAM has at most 2n-1 states

        # Pre-allocate SAM arrays once for the entire batch
        sam_len = np.empty(max_states, dtype=np.int32)
        sam_link = np.empty(max_states, dtype=np.int32)
        sam_t0 = np.empty(max_states, dtype=np.int32)     # transition on bit 0
        sam_t1 = np.empty(max_states, dtype=np.int32)     # transition on bit 1

        # Workspace for topological path counting
        order_buf = np.empty(max_states, dtype=np.int32)
        cnt_buf = np.empty(max_states, dtype=np.int32)
        bucket_buf = np.empty(n + 1, dtype=np.int32)

        # ── Build template SAM if template_len > 0 ───────────────────
        t_size = np.int32(0)
        t_last = np.int32(0)
        t_has_zero = False
        if template_len > 0 and m > 0:
            # Template SAM arrays (separate from working arrays)
            t_len = np.empty(max_states, dtype=np.int32)
            t_link = np.empty(max_states, dtype=np.int32)
            t_t0 = np.empty(max_states, dtype=np.int32)
            t_t1 = np.empty(max_states, dtype=np.int32)
            t_len[0] = np.int32(0)
            t_link[0] = np.int32(-1)
            t_t0[0] = np.int32(-1)
            t_t1[0] = np.int32(-1)
            t_size = np.int32(1)
            t_last = np.int32(0)
            bits0 = bits_batch[0]
            for i in range(template_len):
                c = bits0[i]
                if c == 0:
                    t_has_zero = True
                cur = t_size
                t_len[cur] = t_len[t_last] + np.int32(1)
                t_link[cur] = np.int32(-1)
                t_t0[cur] = np.int32(-1)
                t_t1[cur] = np.int32(-1)
                t_size += np.int32(1)
                p = t_last
                while p != -1:
                    if c == 0:
                        if t_t0[p] >= 0:
                            break
                        t_t0[p] = cur
                    else:
                        if t_t1[p] >= 0:
                            break
                        t_t1[p] = cur
                    p = t_link[p]
                if p == -1:
                    t_link[cur] = np.int32(0)
                else:
                    q = t_t0[p] if c == 0 else t_t1[p]
                    if t_len[p] + 1 == t_len[q]:
                        t_link[cur] = q
                    else:
                        cl = t_size
                        t_len[cl] = t_len[p] + np.int32(1)
                        t_link[cl] = t_link[q]
                        t_t0[cl] = t_t0[q]
                        t_t1[cl] = t_t1[q]
                        t_size += np.int32(1)
                        while p != -1:
                            if c == 0:
                                if t_t0[p] == q:
                                    t_t0[p] = cl
                                else:
                                    break
                            else:
                                if t_t1[p] == q:
                                    t_t1[p] = cl
                                else:
                                    break
                            p = t_link[p]
                        t_link[q] = cl
                        t_link[cur] = cl
                t_last = cur

        for idx in range(m):
            bits = bits_batch[idx]

            # ── Zero-run early-exit filter ────────────────────────────
            if max_zeros >= 0:
                max_run = np.int32(0)
                run = np.int32(0)
                skip = False
                for i in range(n):
                    if bits[i] == 0:
                        run += np.int32(1)
                        if run > max_run:
                            max_run = run
                            if max_run > max_zeros:
                                skip = True
                                break
                    else:
                        run = np.int32(0)
                if skip:
                    out[idx] = np.int32(-1)
                    continue

            # ── Build suffix automaton ────────────────────────────────
            if template_len > 0:
                # Clone template SAM, then extend with tail
                for s in range(t_size):
                    sam_len[s] = t_len[s]
                    sam_link[s] = t_link[s]
                    sam_t0[s] = t_t0[s]
                    sam_t1[s] = t_t1[s]
                size = t_size
                last_st = t_last
                has_zero = t_has_zero
                build_start = template_len
            else:
                sam_len[0] = np.int32(0)
                sam_link[0] = np.int32(-1)
                sam_t0[0] = np.int32(-1)
                sam_t1[0] = np.int32(-1)
                size = np.int32(1)
                last_st = np.int32(0)
                has_zero = False
                build_start = 0

            max_l = sam_len[last_st]  # track max len during build

            for i in range(build_start, n):
                c = bits[i]
                if c == 0:
                    has_zero = True
                cur = size
                new_len = sam_len[last_st] + np.int32(1)
                sam_len[cur] = new_len
                if new_len > max_l:
                    max_l = new_len
                sam_link[cur] = np.int32(-1)
                sam_t0[cur] = np.int32(-1)
                sam_t1[cur] = np.int32(-1)
                size += np.int32(1)

                p = last_st
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
                    sam_link[cur] = np.int32(0)
                else:
                    q = sam_t0[p] if c == 0 else sam_t1[p]

                    if sam_len[p] + 1 == sam_len[q]:
                        sam_link[cur] = q
                    else:
                        clone = size
                        sam_len[clone] = sam_len[p] + np.int32(1)
                        sam_link[clone] = sam_link[q]
                        sam_t0[clone] = sam_t0[q]
                        sam_t1[clone] = sam_t1[q]
                        size += np.int32(1)

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

                last_st = cur

            # ── Count substrings starting with '1' via path counting ──
            # Topological sort by sam_len descending (bucket sort)
            for i in range(max_l + 1):
                bucket_buf[i] = np.int32(0)
            for s in range(size):
                bucket_buf[sam_len[s]] += np.int32(1)

            total_tmp = np.int32(0)
            for i in range(max_l, -1, -1):
                c = bucket_buf[i]
                bucket_buf[i] = total_tmp
                total_tmp += c

            for s in range(size):
                order_buf[bucket_buf[sam_len[s]]] = np.int32(s)
                bucket_buf[sam_len[s]] += np.int32(1)

            # Propagate path counts (distinct non-empty paths from each state)
            for s in range(size):
                cnt_buf[s] = np.int32(0)

            for oi in range(size):
                s = order_buf[oi]
                val = np.int32(0)
                if sam_t0[s] >= 0:
                    val += np.int32(1) + cnt_buf[sam_t0[s]]
                if sam_t1[s] >= 0:
                    val += np.int32(1) + cnt_buf[sam_t1[s]]
                cnt_buf[s] = val

            # Distinct substrings starting with '1' = paths from root's '1' child
            if sam_t1[0] >= 0:
                count = np.int32(1) + cnt_buf[sam_t1[0]]
            else:
                count = np.int32(0)

            # Add 1 for integer 0 if string contains a '0' bit
            if has_zero:
                count += np.int32(1)

            out[idx] = count
        return out

    def _fast_eval(num: int, n: int, max_zeros: int = -1) -> int:
        """Evaluate using Numba JIT suffix array + Kasai LCP.

        max_zeros: if >= 0, returns -1 if longest zero run exceeds this (filter,
        no SA computed). Default -1 = no filtering.
        Uses int.to_bytes + np.unpackbits for fast bit extraction (~6x vs list comp).
        """
        nbytes = (n + 7) >> 3
        bits = np.unpackbits(
            np.frombuffer(num.to_bytes(nbytes, 'big'), dtype=np.uint8),
            bitorder='big'
        )[nbytes * 8 - n:].astype(np.int8)
        return int(_numba_eval_filtered(bits, max_zeros))  # type: ignore[arg-type]

else:
    def _fast_eval(num: int, n: int, max_zeros: int = -1) -> int:
        """Fallback to naive O(n²) evaluator (max_zeros filter ignored)."""
        return distinct_substring_count_int(num, n)


# ---------------------------------------------------------------------------
# Known-best propagation: shared counter accessible in every worker process
# ---------------------------------------------------------------------------
_global_best = None   # set to a multiprocessing.Value by _init_worker
_BATCH_SIZE  = 4096   # tails per Numba call (larger = less Python overhead)


def _init_worker(shared_best):
    """Pool initializer: bind the shared best Value into this worker process."""
    global _global_best
    _global_best = shared_best


# ---------------------------------------------------------------------------
# Combinatorial enumeration via Gosper's hack
#
# Instead of generating all 2^T tails and SWAR-filtering by popcount, we
# enumerate only the C(T,k) tails with exactly k ones for each valid k.
# This eliminates the popcount filter entirely and reduces SA evaluations
# by the fraction of tails that would have been rejected.
# ---------------------------------------------------------------------------

# Pre-computed C(n,k) table for 0 <= n,k <= 64.
# Entries that would overflow int64 are clamped to 2^62.
def _build_c_table():
    C = np.zeros((65, 65), dtype=np.int64)
    C[0, 0] = 1
    for i in range(1, 65):
        C[i, 0] = 1
        for j in range(1, i + 1):
            v = int(C[i - 1, j - 1]) + int(C[i - 1, j])
            C[i, j] = min(v, np.iinfo(np.int64).max // 2)
    return C

_C_TABLE = _build_c_table()

if _NUMBA_OK:
    @numba.njit(cache=True)
    def _comb_at_rank(T, k, rank, C):
        """Return the T-bit integer with exactly k ones at 0-based combinatorial rank.

        Uses the combinatorial number system: the rank is decomposed as
            rank = C(b_k, k) + C(b_{k-1}, k-1) + ... + C(b_1, 1)
        with b_k > b_{k-1} > ... > b_1 >= 0.  Bit b_i is set in the result.
        This corresponds to ascending integer order (= Gosper traversal order).
        """
        result = np.int64(0)
        rem = rank
        for i in range(k, 0, -1):
            b = i - 1
            while b + 1 <= T and C[b + 1, i] <= rem:
                b += 1
            rem -= C[b, i]
            result |= np.int64(1) << np.int64(b)
        return result

    @numba.njit(cache=True)
    def _fill_gosper_slice(T, k, rank_start, count, C, out):
        """Fill out[0:count] with T-bit integers having exactly k ones,
        beginning at combinatorial rank rank_start (ascending integer order).
        Uses Gosper's hack to iterate after finding the starting combination.
        """
        if count <= 0:
            return
        if k == 0:
            out[0] = np.int64(0)
            return
        x = _comb_at_rank(T, k, rank_start, C)
        for i in range(count):
            out[i] = x
            if i < count - 1:
                # Gosper's hack: next integer with same popcount
                c_val = x & (-x)
                r = x + c_val
                x = ((r ^ x) >> np.int64(2)) // c_val | r


def _evaluate_template_chunk(args):
    """Worker: evaluate a chunk of tail combinations for a template prefix.

    args = (n, template_str, k, rank_start, rank_end, tail_bits,
            max_zeros_block, min_one_blocks, max_one_blocks)

    k           : exact number of ones required in the tail (popcount)
    rank_start  : 0-based start index into the sorted C(tail_bits, k) combinations
    rank_end    : exclusive end index
    tail_bits   : n - len(template_str)

    Returns (best_count, best_strings).
    """
    n, template_str, k, rank_start, rank_end, tail_bits, \
        max_zeros_block, min_one_blocks, max_one_blocks = args

    template_len  = len(template_str)
    template_int  = int(template_str, 2) if template_str else 0
    max_zeros_arg = max_zeros_block if max_zeros_block is not None else -1
    chunk_count   = rank_end - rank_start

    best_count   = 0
    best_strings = []

    if _NUMBA_OK and tail_bits > 0 and chunk_count > 0:
        # ── Build template bit row ────────────────────────────────────────
        if template_len > 0:
            nbytes_t = (template_len + 7) >> 3
            tmpl_bits = np.unpackbits(
                np.frombuffer(template_int.to_bytes(nbytes_t, 'big'), dtype=np.uint8),
                bitorder='big'
            )[nbytes_t * 8 - template_len:].astype(np.int8)
        else:
            tmpl_bits = np.empty(0, dtype=np.int8)

        shifts = np.arange(tail_bits - 1, -1, -1, dtype=np.int64)

        # ── Enumerate valid tails by exact popcount via Gosper ────────────
        raw_tails = np.empty(chunk_count, dtype=np.int64)
        if k == 0:
            raw_tails[0] = np.int64(0)
        else:
            _fill_gosper_slice(tail_bits, k, np.int64(rank_start), chunk_count,
                               _C_TABLE, raw_tails)

        # ── Pre-allocate full-string buffer (template columns constant) ───
        bits_buf = np.empty((_BATCH_SIZE, n), dtype=np.int8)
        if template_len > 0:
            bits_buf[:, :template_len] = tmpl_bits[None, :]

        # Pre-compute 1-block count of the template prefix
        if (min_one_blocks is not None or max_one_blocks is not None) and template_len > 0:
            _td = np.diff(tmpl_bits)
            template_one_blocks = int((_td == 1).sum()) + int(tmpl_bits[0] == 1)
        else:
            template_one_blocks = None

        for b_start in range(0, chunk_count, _BATCH_SIZE):
            batch_np = raw_tails[b_start: b_start + _BATCH_SIZE]
            bsz      = len(batch_np)

            # Write tail bits into pre-allocated buffer
            bits_buf[:bsz, template_len:] = (
                (batch_np[:, None] >> shifts[None, :]) & np.int64(1)
            ).astype(np.int8)

            # ── Vectorised 1-block count filter ──────────────────────────
            if template_one_blocks is not None and (min_one_blocks is not None or max_one_blocks is not None):
                tail_part = bits_buf[:bsz, template_len:]
                if template_len > 0:
                    boundary = (tail_part[:, 0] - tmpl_bits[-1]).reshape(bsz, 1)
                else:
                    boundary = np.zeros((bsz, 1), dtype=np.int8)
                tail_diff = np.diff(tail_part, axis=1)
                joined = np.concatenate([boundary, tail_diff], axis=1)
                tail_new_blocks = (joined == 1).sum(axis=1).astype(np.int32)
                total_blk = tail_new_blocks + template_one_blocks
                blk_mask = np.ones(bsz, dtype=bool)
                if min_one_blocks is not None:
                    blk_mask &= total_blk >= min_one_blocks
                if max_one_blocks is not None:
                    blk_mask &= total_blk <= max_one_blocks
                if not np.any(blk_mask):
                    continue
                if blk_mask.all():
                    eval_batch = bits_buf[:bsz]
                    eval_tails = batch_np
                else:
                    eval_batch = np.ascontiguousarray(bits_buf[:bsz][blk_mask])
                    eval_tails = batch_np[blk_mask]
            else:
                eval_batch = bits_buf[:bsz]
                eval_tails = batch_np

            counts = _numba_eval_batch(eval_batch, np.int32(max_zeros_arg),
                                       template_len)

            valid_mask = counts >= 0
            if np.any(valid_mask):
                batch_max = int(np.max(counts[valid_mask]))
                if batch_max >= best_count:
                    top_tails = eval_tails[counts == batch_max]
                    if batch_max > best_count:
                        best_count   = batch_max
                        best_strings = [format((template_int << tail_bits) | int(t), f'0{n}b')
                                        for t in top_tails]
                        if _global_best is not None and batch_max > _global_best.value:
                            _global_best.value = batch_max
                    else:
                        best_strings.extend(
                            format((template_int << tail_bits) | int(t), f'0{n}b')
                            for t in top_tails
                        )

    elif tail_bits == 0:
        # No tail: the template is the full string
        count = _fast_eval(template_int, n, max_zeros=max_zeros_arg)
        if count > 0:
            best_count   = count
            best_strings = [template_str]

    else:
        # ── Scalar fallback (no numpy / numba) ───────────────────────────
        from itertools import combinations as _combs
        for combo in _combs(range(tail_bits), k):
            tail_int = sum(1 << b for b in combo)
            full_num = (template_int << tail_bits) | tail_int
            count    = distinct_substring_count_int(full_num, n)
            if count > best_count:
                best_count   = count
                best_strings = [format(full_num, f'0{n}b')]
            elif count == best_count:
                best_strings.append(format(full_num, f'0{n}b'))

    return best_count, best_strings


# ---------------------------------------------------------------------------
# MH seed loader — seed shared_best from existing results/n_XXXX_results.json
# ---------------------------------------------------------------------------

def _load_mh_seed(n: int) -> int:
    """Return best best_value from results/ and bounded_MH_results/ if they exist, else 0.

    The MH algorithm stores its best-found value there; seeding the exhaustive
    search with it causes aggressive pruning from the very first task.
    """
    base = os.path.join(os.path.dirname(__file__), '..', '..')
    best = 0
    for subdir in ('results/mh_unbounded', 'results/mh_bounded', 'results/large_n'):
        fname = os.path.join(base, subdir, f'n_{n:04d}_results.json')
        if not os.path.exists(fname):
            continue
        try:
            with open(fname, 'r') as f:
                data = json.load(f)
            val = int(data.get('best_value', 0))
            if val > best:
                best = val
        except Exception:
            pass
    return best


# ---------------------------------------------------------------------------
# Cache I/O
# ---------------------------------------------------------------------------

CACHE_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'cached_results.json')


def _extract_seps(s: str) -> list[int]:
    """Return separator (0-block) lengths between 1-blocks in a binary string."""
    seps = []
    i = 0
    n = len(s)
    while i < n:
        if s[i] == '1':
            while i < n and s[i] == '1':
                i += 1
            if i < n and s[i] == '0':
                run = 0
                while i < n and s[i] == '0':
                    run += 1
                    i += 1
                seps.append(run)
        else:
            break
    return seps


def _compute_common_sep_prefix(optimal_strings: list[str]) -> tuple[int, list[int]]:
    """Return (K_common, common_seps) for a list of optimal bit-strings.

    K_common is the length of the longest separator prefix that is identical
    across 100% of solutions.  common_seps is that prefix as a list of ints.
    """
    if not optimal_strings:
        return 0, []
    all_seps = [_extract_seps(s) for s in optimal_strings]
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
    return K_common, list(all_seps[0][:K_common])


def load_cache() -> dict:
    """Load cached results. Keys are string versions of n."""
    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, 'r') as f:
            return json.load(f)
    return {}


def save_cache(cache: dict):
    os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)
    with open(CACHE_PATH, 'w') as f:
        json.dump(cache, f, indent=2)


# ---------------------------------------------------------------------------
# Learned bounds loading
# ---------------------------------------------------------------------------

def _load_learned_bounds(n: int) -> dict | None:
    """Load pre-analyzed bounds for target n.

    Search order:
      1. config/learned_bounds.json  — active/current bounds (list under 'bounds' key)
      2. config/search_constraints.json — archive of per-n bounds (under '{n}_bounds' key)

    Returns dict with K_common, common_seps, block_ranges, and optional constraint
    fields, or None if no bounds found for this n.
    """
    cfg_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'config')

    def _parse_entry(entry):
        return {
            'common_seps':    entry['common_seps'],
            'K_common':       entry['K_common'],
            'block_ranges':   [tuple(r) for r in entry['block_ranges']],
            'min_total_1s':   entry.get('min_total_1s'),
            'max_total_1s':   entry.get('max_total_1s'),
            'max_zeros_block': entry.get('max_zeros_block'),
            'min_one_blocks': entry.get('min_one_blocks'),
            'max_one_blocks': entry.get('max_one_blocks'),
        }

    # 1. learned_bounds.json (most recent / active)
    learned_file = os.path.join(cfg_dir, 'learned_bounds.json')
    if os.path.exists(learned_file):
        try:
            with open(learned_file, 'r') as f:
                data = json.load(f)
            for entry in reversed(data.get('bounds', [])):
                if entry.get('target_n') == n:
                    return _parse_entry(entry)
        except Exception:
            pass

    # 2. search_constraints.json (archive)
    versions_file = os.path.join(cfg_dir, 'search_constraints.json')
    if os.path.exists(versions_file):
        try:
            with open(versions_file, 'r') as f:
                data = json.load(f)
            key = f'{n}_bounds'
            entries = data.get(key, [])
            if entries:
                # Take the last entry (most recent) for this n
                entry = entries[-1]
                if entry.get('target_n') == n:
                    return _parse_entry(entry)
        except Exception:
            pass

    return None


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def compute_a112509(n: int, num_workers: int, cache: dict) -> tuple[int, list[str]]:
    """Compute a(n) using learned bounds.
    
    Returns (count, list_of_optimal_strings)
    """
    if n == 1:
        return 1, ["0", "1"]
    
    learned = _load_learned_bounds(n)
    if learned is None:
        raise ValueError(f"No learned bounds found for n={n}. Add bounds to config/learned_bounds.json")
    
    # Trigger Numba JIT compilation in the main process before spawning workers.
    # This also pre-populates the on-disk cache so workers skip recompilation.
    if _NUMBA_OK:
        print("    [JIT] Warming up Numba (may take 30-60s on first run)...", flush=True)
        _t_jit = time.time()
        _fast_eval(0b10110101, 8)  # tiny dummy call to force compilation
        print(f"    [JIT] Warmup complete in {time.time()-_t_jit:.1f}s", flush=True)
    
    return _compute_with_learned_bounds(n, learned, num_workers, cache)


def _build_tasks(n: int, learned: dict) -> tuple[list, int, int]:
    """Generate evaluation tasks from learned bounds.

    Returns (tasks, total_valid_combs, num_templates).
    Each task is a tuple:
        (n, template_str, k, rank_start, rank_end, tail_bits,
         max_zeros_block, min_one_blocks, max_one_blocks)
    """
    K_common       = learned['K_common']
    common_seps    = learned['common_seps']
    block_ranges   = learned['block_ranges']
    min_total_1s   = learned.get('min_total_1s')
    max_total_1s   = learned.get('max_total_1s')
    max_zeros_block = learned.get('max_zeros_block')
    min_one_blocks  = learned.get('min_one_blocks')
    max_one_blocks  = learned.get('max_one_blocks')

    COMB_CHUNK = 2_000_000  # max combinations per task

    tasks = []
    total_valid_combs = 0
    num_templates  = 0
    num_skipped    = 0

    common_ranges = [range(lo, hi + 1) for lo, hi in block_ranges[:K_common]]

    for common_blocks in itertools.product(*common_ranges):
        common_len = sum(common_blocks) + sum(common_seps)
        if common_len > n:
            num_skipped += 1
            continue

        tail_bits = n - common_len
        tmpl_1s   = sum(common_blocks)

        template = ''
        for i, b in enumerate(common_blocks):
            template += '1' * b
            if i < len(common_seps):
                template += '0' * common_seps[i]

        num_templates += 1

        if tail_bits == 0:
            tasks.append((n, template, 0, 0, 1, 0,
                          max_zeros_block, min_one_blocks, max_one_blocks))
            total_valid_combs += 1
            continue

        k_min = max(0, min_total_1s - tmpl_1s) if min_total_1s is not None else 0
        k_max = min(tail_bits, max_total_1s - tmpl_1s) if max_total_1s is not None else tail_bits

        if k_min > k_max:
            num_skipped += 1
            num_templates -= 1
            continue

        for k in range(k_min, k_max + 1):
            N = int(_C_TABLE[tail_bits, k]) if tail_bits <= 64 else 0
            if N == 0:
                continue
            total_valid_combs += N
            for start in range(0, N, COMB_CHUNK):
                end = min(start + COMB_CHUNK, N)
                tasks.append((n, template, k, start, end, tail_bits,
                              max_zeros_block, min_one_blocks, max_one_blocks))

    print(f"    Skipped (len>n or k impossible): {num_skipped}")
    print(f"    Valid templates: {num_templates}")
    if min_total_1s is not None:
        print(f"    Tail popcount range: encoded per task (no SWAR filter)")
    if max_zeros_block is not None:
        print(f"    Max zeros block: {max_zeros_block}")
    if min_one_blocks is not None:
        print(f"    1-block count range: {min_one_blocks}–{max_one_blocks}")
    if total_valid_combs > 0:
        print(f"    Total valid tail combinations: {total_valid_combs:,}"
              f"  (~2^{math.log2(total_valid_combs):.1f})")
    print(f"    Parallel tasks created: {len(tasks)}", flush=True)

    return tasks, total_valid_combs, num_templates


def _run_parallel(tasks, total_valid_combs, n, num_workers):
    """Execute tasks across worker processes with progress monitoring."""
    start_time = time.time()
    print(f"    Spawning {num_workers} workers...", flush=True)
    mh_seed = _load_mh_seed(n)
    if mh_seed > 0:
        print(f"    Seeded best from MH results: {mh_seed}", flush=True)
    shared_best = multiprocessing.Value('i', mh_seed)

    with ProcessPoolExecutor(max_workers=num_workers,
                             initializer=_init_worker,
                             initargs=(shared_best,)) as executor:
        all_results  = []
        total_tasks  = len(tasks)
        running_best = mh_seed
        combs_done   = 0
        completed_i  = 0
        lock         = threading.Lock()
        stop_monitor = threading.Event()

        MONITOR_INTERVAL = 60

        def _monitor():
            while not stop_monitor.wait(MONITOR_INTERVAL):
                with lock:
                    cd, ci, best = combs_done, completed_i, running_best
                elapsed = time.time() - start_time
                frac = cd / total_valid_combs if total_valid_combs > 0 else 0
                eta  = (elapsed / frac - elapsed) if frac > 1e-9 else 0
                print(f"    [monitor] {frac*100:.1f}% combs ({cd:,}/{total_valid_combs:,}) "
                      f"| tasks {ci}/{total_tasks} "
                      f"| best={best} "
                      f"| {elapsed/3600:.2f}h elapsed "
                      f"| ETA {eta/3600:.2f}h", flush=True)

        mon_thread = threading.Thread(target=_monitor, daemon=True)
        mon_thread.start()

        MAX_INFLIGHT = num_workers * 64
        task_iter = iter(tasks)
        pending = {}

        for task in itertools.islice(task_iter, MAX_INFLIGHT):
            f = executor.submit(_evaluate_template_chunk, task)
            pending[f] = task

        while pending:
            done_set, _ = wait(list(pending.keys()), return_when=FIRST_COMPLETED)
            for future in done_set:
                task = pending.pop(future)
                try:
                    count, strings = future.result()
                except Exception as exc:
                    stop_monitor.set()
                    print(f"\n[ERROR] Worker raised: {type(exc).__name__}: {exc}", flush=True)
                    traceback.print_exc()
                    raise
                all_results.append((count, strings))
                with lock:
                    completed_i += 1
                    combs_done  += task[4] - task[3]
                    new_best = count > running_best
                    if new_best:
                        running_best = count
                        if count > shared_best.value:
                            shared_best.value = count
                    ci, cd, best = completed_i, combs_done, running_best
                if new_best:
                    elapsed = time.time() - start_time
                    frac = cd / total_valid_combs if total_valid_combs > 0 else 0
                    eta  = (elapsed / frac - elapsed) if frac > 1e-9 else 0
                    print(f"  *** NEW BEST: {best} "
                          f"| {frac*100:.1f}% combs ({cd:,}/{total_valid_combs:,}) "
                          f"| {ci}/{total_tasks} tasks "
                          f"| {elapsed:.0f}s elapsed "
                          f"| ETA {eta/3600:.2f}h ***", flush=True)
                try:
                    next_task = next(task_iter)
                    f = executor.submit(_evaluate_template_chunk, next_task)
                    pending[f] = next_task
                except StopIteration:
                    pass

        stop_monitor.set()
        mon_thread.join(timeout=2)

    return all_results


def _compute_with_learned_bounds(n: int, learned: dict, num_workers: int, cache: dict) -> tuple[int, list[str]]:
    """Compute using learned K-block structure + Gosper tail enumeration."""
    tasks, total_valid_combs, _ = _build_tasks(n, learned)

    if not tasks:
        return 0, []

    if total_valid_combs <= 100_000 or len(tasks) <= 2:
        all_results = [_evaluate_template_chunk(t) for t in tasks]
    else:
        all_results = _run_parallel(tasks, total_valid_combs, n, num_workers)

    best_count   = 0
    best_strings = []
    for count, strings in all_results:
        if count > best_count:
            best_count   = count
            best_strings = strings
        elif count == best_count:
            best_strings.extend(strings)

    best_strings = list(set(best_strings))
    return best_count, best_strings


def main():
    # ========================================================================
    # USER CONFIGURATION - CHANGE THESE VALUES AS NEEDED
    # ========================================================================

    n_start = 1                # first n to compute
    n_end   = 30               # last n to compute
    force_recompute = False    # True = ignore cache and recompute all
    num_workers = cpu_count()  # number of parallel worker processes

    # ========================================================================
    # END USER CONFIGURATION
    # ========================================================================

    # Also accept CLI overrides: [end] or [start end] [--force]
    args = sys.argv[1:]
    if '--force' in args or '-f' in args:
        force_recompute = True
        args = [a for a in args if a not in ('--force', '-f')]
    if len(args) == 1:
        n_start, n_end = 1, int(args[0])
    elif len(args) >= 2:
        n_start, n_end = int(args[0]), int(args[1])

    print(f"Computing A112509 for n = {n_start}..{n_end}")
    print(f"Using {num_workers} worker processes")
    print(f"Cache file: {os.path.abspath(CACHE_PATH)}")
    if force_recompute:
        print(f"Force recompute: YES (ignoring cache)")
    print()
    
    cache = load_cache()
    
    for n in range(n_start, n_end + 1):
        key = str(n)
        if key in cache and not force_recompute:
            print(f"n={n:>3}  a(n)={cache[key]['a(n)']:>5}  "
                  f"#opt={len(cache[key]['optimal_strings']):>4}  [cached]")
            continue
        
        learned = _load_learned_bounds(n)
        t0 = time.time()
        count, strings = compute_a112509(n, num_workers, cache)
        elapsed = time.time() - t0
        
        k_common, common_seps = _compute_common_sep_prefix(strings)

        cache[key] = {
            "a(n)": count,
            "num_optimal": len(strings),
            "optimal_strings": strings,
            "K_common": k_common,
            "common_seps": common_seps,
        }
        save_cache(cache)
        
        if learned:
            K = learned['K_common']
            seps = learned['common_seps']
            ranges_str = ' '.join(f"{chr(97+i)}=[{lo},{hi}]"
                                 for i, (lo, hi) in enumerate(learned['block_ranges'][:K]))
            print(f"n={n:>3}  a(n)={count:>5}  #opt={len(strings):>4}  "
                  f"time={elapsed:>7.2f}s  K={K} seps={seps}  {ranges_str}")
        else:
            print(f"n={n:>3}  a(n)={count:>5}  #opt={len(strings):>4}  "
                  f"time={elapsed:>7.2f}s  [brute force]")
    
    print("\nDone. Results saved to", os.path.abspath(CACHE_PATH))


if __name__ == '__main__':
    freeze_support()
    main()
