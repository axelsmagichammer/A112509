"""
Exact branch-and-bound search for OEIS A112509 values.

Overview
--------
For a given n, a(n) is the maximum number of distinct integers representable
as substrings of any n-bit binary string.  This module finds a(n) exactly and
enumerates ALL n-bit strings that achieve it.

The search is provably correct: the upper bound on every node is a rigorous
over-estimate of the best score achievable from that prefix, so no optimal
string is ever pruned.

------------------------------------------------------------------------
Algorithm: parallel iterative DFS branch-and-bound (primary mode)
------------------------------------------------------------------------

This is the mode used in production (collect_all=True).  It does NOT use
a priority queue.  Instead each worker runs an iterative depth-first search
with an explicit stack, bounded by a rigorous upper bound (UB) at every node.

Phase 1 -- Parallel enumeration of root tasks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BFS from the root prefix (1,) to depth split_depth (~n//3, clamped 10-25).
Every prefix whose UB >= incumbent is a candidate root for a parallel task.
"Hard" subtrees (UB <= incumbent + hard_split_gap) are expanded extra_levels
more BFS levels so that individual tasks are roughly equal in cost.
This produces O(thousands) of fine-grained independent tasks.

Phase 2 -- Parallel DFS branch-and-bound workers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each task is an independent DFS rooted at one BFS-frontier prefix.
Workers use an explicit stack of (prefix, precomputed_ub) pairs.
Key properties:
  - Memory: O(n^2) per worker -- stack depth is at most n+1 at any time,
    regardless of the number of optimal strings in the subtree.
    This is immune to the heap explosion that afflicts best-first search
    when collect_all=True (which must keep every UB==optimal node alive).
  - Each node's UB is computed ONCE at push time and stored in the stack
    entry.  On pop, the stored UB is compared against `best` (which can
    only increase) -- no redundant re-evaluation.
  - When a leaf is reached its stored UB equals the exact score (m_rem=0,
    lookahead is capped to L_use=0).  No re-scoring is needed.
  - If best increases, previously pushed nodes with stale UB < new best
    are correctly pruned on pop.

Phase 2b -- Tail re-split
~~~~~~~~~~~~~~~~~~~~~~~~~
When the external task queue drops to <=10% of its original size, the
remaining QUEUED (not yet submitted) tasks are each BFS-expanded 4 more
levels.  This converts the last few potentially huge tasks into many smaller
ones, keeping all workers busy through the expensive tail of the run where
a few near-optimal subtrees dominate wall time.

------------------------------------------------------------------------
Upper bound: L-step lookahead with incremental SAM + savings tightening
------------------------------------------------------------------------

For prefix p_k of length k:

  UB_L(p_k) = max over all 2^L bit-extensions e of length L:
                 score(p_k + e) + remaining(p_k + e, n)

where score(s) = number of distinct integer values represented by substrings
of s, and remaining is a provable upper bound on the additional values any
(n - k - L)-bit suffix can contribute.

Remaining bound (tight decomposition)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For extended prefix q = p_k + e of length k+L with m = n-k-L bits remaining:

  remaining(q, m) = a(m) + popcount(q) * m - SAM_savings

where:
  a(m)            -- tight bound on values contributed by pure-suffix substrings
                     (all m-bit substrings of the remaining suffix)
  popcount(q) * m -- bound on cross-boundary values: only substrings starting
                     at a '1'-bit in q can exceed all pure-suffix values; there
                     are popcount(q) such starting positions each contributing
                     at most m new values.  ('0'-leading cross-boundary substrings
                     equal a pure-suffix value and contribute nothing new.)
  SAM_savings     -- subtracts cross-boundary values already proven to be
                     in the SAM of q (i.e., both 1-step and 2-step continuations
                     of the SAM node for each '1'-bit suffix already exist),
                     so they cannot contribute NEW distinct integers.

SAM savings computation (Numba JIT, inside _sam_score_extensions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The prefix SAM is built once for p_k, then each of 2^L extensions is evaluated
by cloning the prefix SAM and extending it L more steps.  After extension the
suffix-link chain from s_last is traversed in O(k+L) to find, for each '1'-bit
at position i in q, the representative SAM node v_i.  A BFS from v_i up to
JMAX=2 steps checks whether all 2^j j-step continuations already exist in the
SAM; each fully-covered step subtracts 1 from the remaining cross-boundary
budget.  Total cost per UB evaluation: O(k + 2^L * (L + popcount(q) * JMAX)).

Correctness guarantee
~~~~~~~~~~~~~~~~~~~~~
For any n-bit string s, the UB at its length-k prefix is always >= score(s).
Proof: score(s) = score(p_k) + (cross-boundary values) + (suffix-only values).
  - score(p_k + e) >= score(p_k) (taking e = first L bits of s)
  - remaining(q, m) >= suffix-only values + cross-boundary values not yet in SAM
  Hence UB_L >= score(s) for the optimal extension e.
  Pruning only occurs when UB < incumbent <= a(n), so no optimal string is pruned.

------------------------------------------------------------------------
Sequential best-first branch-and-bound (collect_all=False)
------------------------------------------------------------------------
When collect_all=False, workers instead use a max-priority queue (sorted by UB)
rather than a DFS stack.  This visits nodes in best-first order, finds a(n)
quickly with few node expansions, but the heap can grow very large when many
near-optimal nodes exist -- causing OOM at large n with collect_all=True.
It is therefore used only for value-finding (not exhaustive enumeration).
"""

from __future__ import annotations

import heapq
import json
import multiprocessing
import os
import time
from collections import deque
from itertools import product as _product
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Global UB cap: per-length bound sum_b min(n-b+1, 2^(b-1)) + 1 (for value 0).
# ---------------------------------------------------------------------------

def _compute_global_max(n: int) -> int:
    """Return the tightest global UB on a(n), independent of any prefix.

    Computed as 1 + sum_{b=1}^{crossover} 2^(b-1) + sum_{b=crossover+1}^{n} (n-b+1)
    where crossover = max{b : 2^(b-1) <= n-b+1}.
    """
    crossover = max(b for b in range(1, n + 1) if (1 << (b - 1)) <= n - b + 1)
    return (1 + sum(1 << (b - 1) for b in range(1, crossover + 1))
            + sum(n - b + 1 for b in range(crossover + 1, n + 1)))


# ---------------------------------------------------------------------------
# Evaluator: reuse the SAM-based single-string evaluator from structured_search
# ---------------------------------------------------------------------------

try:
    from src.algorithms.structured_search import (
        _numba_eval_batch,
        _numba_eval_filtered,
        _NUMBA_OK as _SA_OK,
    )
    _eval_available = _SA_OK
except ImportError:
    _eval_available = False
    _numba_eval_batch = None

# ---------------------------------------------------------------------------
# Incremental SAM for lookahead UB: build the prefix SAM once, then clone
# and extend for each of 2^L extensions in a single numba-JIT function.
# This replaces 2^L full-string SAM builds with 1 prefix-SAM build +
# 2^L cheap L-step extensions per node expansion -- critical speedup.
# ---------------------------------------------------------------------------
try:
    import numba as _numba

    _BASE_SAM_JMAX = 2
    _REFINE_SAM_JMAX = 3

    @_numba.njit(cache=True)
    def _sam_score_extensions(
        prefix_bits,           # int8[k] -- the prefix
        extensions,            # int8[2^L, L] -- all extensions to evaluate
        full_len,              # n -- total string length
        known_a,               # int64[max_m+1] -- known_a[m] = a(m), 0 if unknown
        max_known,             # int -- len(known_a)-1
        sam_jmax,              # max depth of SAM savings BFS
    ):
        """Build prefix SAM once, then evaluate each extension by cloning it.

        Returns the maximum value of
            score(prefix+ext) + remaining(prefix+ext, n)
        over all rows of `extensions`.

        Uses the same SAM build logic as _numba_eval_batch in structured_search.
        """
        k = prefix_bits.shape[0]
        L = extensions.shape[0]
        ext_len = extensions.shape[1]
        full_ext = k + ext_len
        m_rem = full_len - full_ext

        max_states = 2 * full_len + 4

        # ── Arrays for the frozen prefix SAM ─────────────────────────────
        # Use flat arrays -- numba doesn't support list of lists for cloning.
        p_len_arr  = np.zeros(max_states, dtype=np.int32)
        p_link_arr = np.full(max_states, -1, dtype=np.int32)
        p_t0_arr   = np.full(max_states, -1, dtype=np.int32)
        p_t1_arr   = np.full(max_states, -1, dtype=np.int32)
        p_size     = np.int32(1)
        p_last     = np.int32(0)
        p_has_zero = False

        # Online SAM extend for prefix
        for i in range(k):
            c = prefix_bits[i]
            if c == 0:
                p_has_zero = True
            cur = p_size
            p_len_arr[cur] = p_len_arr[p_last] + np.int32(1)
            p_link_arr[cur] = np.int32(-1)
            p_t0_arr[cur]   = np.int32(-1)
            p_t1_arr[cur]   = np.int32(-1)
            p_size += np.int32(1)
            p = p_last
            while p != -1:
                if c == 0:
                    if p_t0_arr[p] >= 0:
                        break
                    p_t0_arr[p] = cur
                else:
                    if p_t1_arr[p] >= 0:
                        break
                    p_t1_arr[p] = cur
                p = p_link_arr[p]
            if p == -1:
                p_link_arr[cur] = np.int32(0)
            else:
                q = p_t0_arr[p] if c == 0 else p_t1_arr[p]
                if p_len_arr[p] + 1 == p_len_arr[q]:
                    p_link_arr[cur] = q
                else:
                    clone = p_size
                    p_len_arr[clone]  = p_len_arr[p] + np.int32(1)
                    p_link_arr[clone] = p_link_arr[q]
                    p_t0_arr[clone]   = p_t0_arr[q]
                    p_t1_arr[clone]   = p_t1_arr[q]
                    p_size += np.int32(1)
                    while p != -1:
                        if c == 0:
                            if p_t0_arr[p] == q:
                                p_t0_arr[p] = clone
                            else:
                                break
                        else:
                            if p_t1_arr[p] == q:
                                p_t1_arr[p] = clone
                            else:
                                break
                        p = p_link_arr[p]
                    p_link_arr[q]   = clone
                    p_link_arr[cur] = clone
            p_last = cur

        # ── Work arrays for cloned (per-extension) SAM ───────────────────
        s_len  = np.empty(max_states, dtype=np.int32)
        s_link = np.empty(max_states, dtype=np.int32)
        s_t0   = np.empty(max_states, dtype=np.int32)
        s_t1   = np.empty(max_states, dtype=np.int32)
        state_cov = np.empty(max_states, dtype=np.int32)
        order_buf  = np.empty(max_states, dtype=np.int32)
        cnt_buf    = np.empty(max_states, dtype=np.int32)
        bucket_buf = np.empty(full_len + 2, dtype=np.int32)

        best_result = np.float64(0.0)
        prefix_ones = np.int32(0)
        for i in range(k):
            if prefix_bits[i] == 1:
                prefix_ones += np.int32(1)

        for eidx in range(L):
            # Clone prefix SAM into working arrays
            sz = p_size
            for s in range(sz):
                s_len[s] = p_len_arr[s]
                s_link[s] = p_link_arr[s]
                s_t0[s] = p_t0_arr[s]
                s_t1[s] = p_t1_arr[s]
            s_size = sz
            s_last = p_last
            s_has_zero = p_has_zero

            ext_ones = prefix_ones
            # Extend with this extension's bits
            for j in range(ext_len):
                c = extensions[eidx, j]
                if c == 0:
                    s_has_zero = True
                else:
                    ext_ones += np.int32(1)
                cur = s_size
                s_len[cur]  = s_len[s_last] + np.int32(1)
                s_link[cur] = np.int32(-1)
                s_t0[cur]   = np.int32(-1)
                s_t1[cur]   = np.int32(-1)
                s_size += np.int32(1)
                p = s_last
                while p != -1:
                    if c == 0:
                        if s_t0[p] >= 0:
                            break
                        s_t0[p] = cur
                    else:
                        if s_t1[p] >= 0:
                            break
                        s_t1[p] = cur
                    p = s_link[p]
                if p == -1:
                    s_link[cur] = np.int32(0)
                else:
                    q = s_t0[p] if c == 0 else s_t1[p]
                    if s_len[p] + 1 == s_len[q]:
                        s_link[cur] = q
                    else:
                        clone = s_size
                        s_len[clone]  = s_len[p] + np.int32(1)
                        s_link[clone] = s_link[q]
                        s_t0[clone]   = s_t0[q]
                        s_t1[clone]   = s_t1[q]
                        s_size += np.int32(1)
                        while p != -1:
                            if c == 0:
                                if s_t0[p] == q:
                                    s_t0[p] = clone
                                else:
                                    break
                            else:
                                if s_t1[p] == q:
                                    s_t1[p] = clone
                                else:
                                    break
                            p = s_link[p]
                        s_link[q]   = clone
                        s_link[cur] = clone
                s_last = cur

            # Count distinct substrings via topological path count
            max_l = s_len[s_last]
            for i in range(max_l + 1):
                bucket_buf[i] = 0
            for s in range(s_size):
                bucket_buf[s_len[s]] += np.int32(1)
            total_tmp = np.int32(0)
            for i in range(max_l, -1, -1):
                c2 = bucket_buf[i]
                bucket_buf[i] = total_tmp
                total_tmp += c2
            for s in range(s_size):
                order_buf[bucket_buf[s_len[s]]] = np.int32(s)
                bucket_buf[s_len[s]] += np.int32(1)
            for s in range(s_size):
                cnt_buf[s] = 0
            for oi in range(s_size):
                s = order_buf[oi]
                val = np.int32(0)
                if s_t0[s] >= 0:
                    val += np.int32(1) + cnt_buf[s_t0[s]]
                if s_t1[s] >= 0:
                    val += np.int32(1) + cnt_buf[s_t1[s]]
                cnt_buf[s] = val
            if s_t1[0] >= 0:
                score = np.int32(1) + cnt_buf[s_t1[0]]
            else:
                score = np.int32(0)
            if s_has_zero:
                score += np.int32(1)

            # Remaining budget (tightened by subtracting cross-boundary pairs
            # whose values are guaranteed to already be in the SAM).
            if m_rem <= 0:
                rem = np.int64(0)
            elif m_rem <= max_known and known_a[m_rem] > 0:
                # ── Step 1: Build suffix-link chain ──────────────────────
                MAX_CHAIN = full_ext + np.int32(2)
                chain_nodes = np.empty(MAX_CHAIN, dtype=np.int32)
                chain_link_lens = np.empty(MAX_CHAIN, dtype=np.int32)
                chain_len = np.int32(0)
                v_ch = s_last
                while v_ch != np.int32(0):
                    chain_nodes[chain_len] = v_ch
                    lnk = s_link[v_ch]
                    chain_link_lens[chain_len] = s_len[lnk] if lnk >= 0 else np.int32(0)
                    chain_len += np.int32(1)
                    v_ch = lnk

                # ── Step 2: For each 1-bit, look up v_i and compute savings ──────
                JMAX = np.int32(sam_jmax)
                MAX_FRONT = np.int32(1 << sam_jmax)
                savings = np.int64(0)
                frontier_cur = np.full(MAX_FRONT, np.int32(-1), dtype=np.int32)
                frontier_nxt = np.full(MAX_FRONT, np.int32(-1), dtype=np.int32)
                for s in range(s_size):
                    state_cov[s] = np.int32(-1)

                chain_pos = np.int32(0)
                for i in range(full_ext):
                    if i < k:
                        bit_i = prefix_bits[i]
                    else:
                        bit_i = extensions[eidx, i - k]
                    if bit_i != np.int8(1):
                        continue

                    suf_len = np.int32(full_ext - i)
                    while chain_pos < chain_len and chain_link_lens[chain_pos] >= suf_len:
                        chain_pos += np.int32(1)
                    if chain_pos >= chain_len:
                        continue
                    v_i = chain_nodes[chain_pos]

                    jmax_use = np.int32(JMAX if JMAX < m_rem else m_rem)
                    cached_cov = state_cov[v_i]
                    if cached_cov >= np.int32(0):
                        covered_depth = cached_cov
                    else:
                        frontier_cur[0] = v_i
                        front_size = np.int32(1)
                        covered_depth = np.int32(0)
                        for j in range(np.int32(1), jmax_use + np.int32(1)):
                            all_covered = True
                            nxt_size = np.int32(0)
                            for fi in range(front_size):
                                v = frontier_cur[fi]
                                for c2 in range(np.int32(2)):
                                    nxt2 = s_t0[v] if c2 == 0 else s_t1[v]
                                    if nxt2 < np.int32(0):
                                        all_covered = False
                                        break
                                    found_nxt = False
                                    for ni in range(nxt_size):
                                        if frontier_nxt[ni] == nxt2:
                                            found_nxt = True
                                            break
                                    if not found_nxt:
                                        if nxt_size < MAX_FRONT:
                                            frontier_nxt[nxt_size] = nxt2
                                            nxt_size += np.int32(1)
                                        else:
                                            all_covered = False
                                            break
                                if not all_covered:
                                    break
                            if not all_covered:
                                break
                            covered_depth += np.int32(1)
                            for fi in range(nxt_size):
                                frontier_cur[fi] = frontier_nxt[fi]
                                frontier_nxt[fi] = np.int32(-1)
                            front_size = nxt_size
                        state_cov[v_i] = covered_depth
                    savings += np.int64(covered_depth)
                    if covered_depth <= np.int32(ext_len):
                        _sp = v_i
                        for _sj in range(covered_depth):
                            _sc = extensions[eidx, _sj]
                            _sp = s_t0[_sp] if _sc == np.int8(0) else s_t1[_sp]
                        _sp_jend = (np.int32(ext_len) if np.int32(ext_len) < jmax_use
                                    else jmax_use)
                        for _sj2 in range(covered_depth, _sp_jend):
                            _sc2 = extensions[eidx, _sj2]
                            _sn2 = s_t0[_sp] if _sc2 == np.int8(0) else s_t1[_sp]
                            if _sn2 < np.int32(0):
                                break
                            savings += np.int64(1)
                            _sp = _sn2

                # Cross-boundary UB: ext_ones * m_rem bounds cross values,
                # a(m_rem) bounds suffix-only values, SAM savings subtracted.
                cross_new = np.int64(ext_ones) * np.int64(m_rem) - savings
                if cross_new < np.int64(0):
                    cross_new = np.int64(0)
                rem = known_a[m_rem] + cross_new
            else:
                rem = (np.int64(full_len) * (full_len + 1) // 2
                       - np.int64(full_ext) * (full_ext + 1) // 2)

            ub_val = np.float64(score) + np.float64(rem)
            if eidx == 0 or ub_val > best_result:
                best_result = ub_val

        return best_result

    @_numba.njit(cache=True)
    def _sam_score_extensions_split_by_first_bit(
        prefix_bits,
        extensions,
        full_len,
        known_a,
        max_known,
        sam_jmax,
    ):
        """As _sam_score_extensions, but return the best UB by first extension bit."""
        k = prefix_bits.shape[0]
        L = extensions.shape[0]
        ext_len = extensions.shape[1]
        full_ext = k + ext_len
        m_rem = full_len - full_ext

        max_states = 2 * full_len + 4

        p_len_arr = np.zeros(max_states, dtype=np.int32)
        p_link_arr = np.full(max_states, -1, dtype=np.int32)
        p_t0_arr = np.full(max_states, -1, dtype=np.int32)
        p_t1_arr = np.full(max_states, -1, dtype=np.int32)
        p_size = np.int32(1)
        p_last = np.int32(0)
        p_has_zero = False

        for i in range(k):
            c = prefix_bits[i]
            if c == 0:
                p_has_zero = True
            cur = p_size
            p_len_arr[cur] = p_len_arr[p_last] + np.int32(1)
            p_link_arr[cur] = np.int32(-1)
            p_t0_arr[cur] = np.int32(-1)
            p_t1_arr[cur] = np.int32(-1)
            p_size += np.int32(1)
            p = p_last
            while p != -1:
                if c == 0:
                    if p_t0_arr[p] >= 0:
                        break
                    p_t0_arr[p] = cur
                else:
                    if p_t1_arr[p] >= 0:
                        break
                    p_t1_arr[p] = cur
                p = p_link_arr[p]
            if p == -1:
                p_link_arr[cur] = np.int32(0)
            else:
                q = p_t0_arr[p] if c == 0 else p_t1_arr[p]
                if p_len_arr[p] + 1 == p_len_arr[q]:
                    p_link_arr[cur] = q
                else:
                    clone = p_size
                    p_len_arr[clone] = p_len_arr[p] + np.int32(1)
                    p_link_arr[clone] = p_link_arr[q]
                    p_t0_arr[clone] = p_t0_arr[q]
                    p_t1_arr[clone] = p_t1_arr[q]
                    p_size += np.int32(1)
                    while p != -1:
                        if c == 0:
                            if p_t0_arr[p] == q:
                                p_t0_arr[p] = clone
                            else:
                                break
                        else:
                            if p_t1_arr[p] == q:
                                p_t1_arr[p] = clone
                            else:
                                break
                        p = p_link_arr[p]
                    p_link_arr[q] = clone
                    p_link_arr[cur] = clone
            p_last = cur

        s_len = np.empty(max_states, dtype=np.int32)
        s_link = np.empty(max_states, dtype=np.int32)
        s_t0 = np.empty(max_states, dtype=np.int32)
        s_t1 = np.empty(max_states, dtype=np.int32)
        state_cov = np.empty(max_states, dtype=np.int32)
        order_buf = np.empty(max_states, dtype=np.int32)
        cnt_buf = np.empty(max_states, dtype=np.int32)
        bucket_buf = np.empty(full_len + 2, dtype=np.int32)

        prefix_ones = np.int32(0)
        for i in range(k):
            if prefix_bits[i] == 1:
                prefix_ones += np.int32(1)

        best_zero = np.float64(-1.0)
        best_one = np.float64(-1.0)

        for eidx in range(L):
            sz = p_size
            for s in range(sz):
                s_len[s] = p_len_arr[s]
                s_link[s] = p_link_arr[s]
                s_t0[s] = p_t0_arr[s]
                s_t1[s] = p_t1_arr[s]
            s_size = sz
            s_last = p_last
            s_has_zero = p_has_zero

            ext_ones = prefix_ones
            for j in range(ext_len):
                c = extensions[eidx, j]
                if c == 0:
                    s_has_zero = True
                else:
                    ext_ones += np.int32(1)
                cur = s_size
                s_len[cur] = s_len[s_last] + np.int32(1)
                s_link[cur] = np.int32(-1)
                s_t0[cur] = np.int32(-1)
                s_t1[cur] = np.int32(-1)
                s_size += np.int32(1)
                p = s_last
                while p != -1:
                    if c == 0:
                        if s_t0[p] >= 0:
                            break
                        s_t0[p] = cur
                    else:
                        if s_t1[p] >= 0:
                            break
                        s_t1[p] = cur
                    p = s_link[p]
                if p == -1:
                    s_link[cur] = np.int32(0)
                else:
                    q = s_t0[p] if c == 0 else s_t1[p]
                    if s_len[p] + 1 == s_len[q]:
                        s_link[cur] = q
                    else:
                        clone = s_size
                        s_len[clone] = s_len[p] + np.int32(1)
                        s_link[clone] = s_link[q]
                        s_t0[clone] = s_t0[q]
                        s_t1[clone] = s_t1[q]
                        s_size += np.int32(1)
                        while p != -1:
                            if c == 0:
                                if s_t0[p] == q:
                                    s_t0[p] = clone
                                else:
                                    break
                            else:
                                if s_t1[p] == q:
                                    s_t1[p] = clone
                                else:
                                    break
                            p = s_link[p]
                        s_link[q] = clone
                        s_link[cur] = clone
                s_last = cur

            max_l = s_len[s_last]
            for i in range(max_l + 1):
                bucket_buf[i] = 0
            for s in range(s_size):
                bucket_buf[s_len[s]] += np.int32(1)
            total_tmp = np.int32(0)
            for i in range(max_l, -1, -1):
                c2 = bucket_buf[i]
                bucket_buf[i] = total_tmp
                total_tmp += c2
            for s in range(s_size):
                order_buf[bucket_buf[s_len[s]]] = np.int32(s)
                bucket_buf[s_len[s]] += np.int32(1)
            for s in range(s_size):
                cnt_buf[s] = 0
            for oi in range(s_size):
                s = order_buf[oi]
                val = np.int32(0)
                if s_t0[s] >= 0:
                    val += np.int32(1) + cnt_buf[s_t0[s]]
                if s_t1[s] >= 0:
                    val += np.int32(1) + cnt_buf[s_t1[s]]
                cnt_buf[s] = val
            if s_t1[0] >= 0:
                score = np.int32(1) + cnt_buf[s_t1[0]]
            else:
                score = np.int32(0)
            if s_has_zero:
                score += np.int32(1)

            if m_rem <= 0:
                rem = np.int64(0)
            elif m_rem <= max_known and known_a[m_rem] > 0:
                # ── Step 1: Build suffix-link chain ──────────────────────
                MAX_CHAIN = full_ext + np.int32(2)
                chain_nodes = np.empty(MAX_CHAIN, dtype=np.int32)
                chain_link_lens = np.empty(MAX_CHAIN, dtype=np.int32)
                chain_len = np.int32(0)
                v_ch = s_last
                while v_ch != np.int32(0):
                    chain_nodes[chain_len] = v_ch
                    lnk = s_link[v_ch]
                    chain_link_lens[chain_len] = s_len[lnk] if lnk >= 0 else np.int32(0)
                    chain_len += np.int32(1)
                    v_ch = lnk

                JMAX = np.int32(sam_jmax)
                MAX_FRONT = np.int32(1 << sam_jmax)
                savings = np.int64(0)
                frontier_cur = np.full(MAX_FRONT, np.int32(-1), dtype=np.int32)
                frontier_nxt = np.full(MAX_FRONT, np.int32(-1), dtype=np.int32)
                for s in range(s_size):
                    state_cov[s] = np.int32(-1)

                chain_pos = np.int32(0)
                for i in range(full_ext):
                    if i < k:
                        bit_i = prefix_bits[i]
                    else:
                        bit_i = extensions[eidx, i - k]
                    if bit_i != np.int8(1):
                        continue

                    suf_len = np.int32(full_ext - i)
                    while chain_pos < chain_len and chain_link_lens[chain_pos] >= suf_len:
                        chain_pos += np.int32(1)
                    if chain_pos >= chain_len:
                        continue
                    v_i = chain_nodes[chain_pos]

                    jmax_use = np.int32(JMAX if JMAX < m_rem else m_rem)
                    cached_cov = state_cov[v_i]
                    if cached_cov >= np.int32(0):
                        covered_depth = cached_cov
                    else:
                        frontier_cur[0] = v_i
                        front_size = np.int32(1)
                        covered_depth = np.int32(0)
                        for j in range(np.int32(1), jmax_use + np.int32(1)):
                            all_covered = True
                            nxt_size = np.int32(0)
                            for fi in range(front_size):
                                v = frontier_cur[fi]
                                for c2 in range(np.int32(2)):
                                    nxt2 = s_t0[v] if c2 == 0 else s_t1[v]
                                    if nxt2 < np.int32(0):
                                        all_covered = False
                                        break
                                    found_nxt = False
                                    for ni in range(nxt_size):
                                        if frontier_nxt[ni] == nxt2:
                                            found_nxt = True
                                            break
                                    if not found_nxt:
                                        if nxt_size < MAX_FRONT:
                                            frontier_nxt[nxt_size] = nxt2
                                            nxt_size += np.int32(1)
                                        else:
                                            all_covered = False
                                            break
                                if not all_covered:
                                    break
                            if not all_covered:
                                break
                            covered_depth += np.int32(1)
                            for fi in range(nxt_size):
                                frontier_cur[fi] = frontier_nxt[fi]
                                frontier_nxt[fi] = np.int32(-1)
                            front_size = nxt_size
                        state_cov[v_i] = covered_depth
                    savings += np.int64(covered_depth)
                    if covered_depth <= np.int32(ext_len):
                        _sp = v_i
                        for _sj in range(covered_depth):
                            _sc = extensions[eidx, _sj]
                            _sp = s_t0[_sp] if _sc == np.int8(0) else s_t1[_sp]
                        _sp_jend = (np.int32(ext_len) if np.int32(ext_len) < jmax_use
                                    else jmax_use)
                        for _sj2 in range(covered_depth, _sp_jend):
                            _sc2 = extensions[eidx, _sj2]
                            _sn2 = s_t0[_sp] if _sc2 == np.int8(0) else s_t1[_sp]
                            if _sn2 < np.int32(0):
                                break
                            savings += np.int64(1)
                            _sp = _sn2

                # Cross-boundary UB: ext_ones * m_rem bounds cross values,
                # a(m_rem) bounds suffix-only values, SAM savings subtracted.
                cross_new = np.int64(ext_ones) * np.int64(m_rem) - savings
                if cross_new < np.int64(0):
                    cross_new = np.int64(0)
                rem = known_a[m_rem] + cross_new
            else:
                rem = (np.int64(full_len) * (full_len + 1) // 2
                       - np.int64(full_ext) * (full_ext + 1) // 2)

            ub_val = np.float64(score) + np.float64(rem)
            if extensions[eidx, 0] == 0:
                if ub_val > best_zero:
                    best_zero = ub_val
            else:
                if ub_val > best_one:
                    best_one = ub_val

        return best_zero, best_one


    @_numba.njit(cache=True)
    def _score_bits_numba(bits):
        """Exact distinct-substring-value count. O(n^2) brute, JIT'd."""
        n = bits.shape[0]
        max_states = 2 * n + 4
        s_len = np.empty(max_states, dtype=np.int32)
        s_link = np.empty(max_states, dtype=np.int32)
        s_t0 = np.empty(max_states, dtype=np.int32)
        s_t1 = np.empty(max_states, dtype=np.int32)
        s_len[0] = np.int32(0)
        s_link[0] = np.int32(-1)
        s_t0[0] = np.int32(-1)
        s_t1[0] = np.int32(-1)
        s_size = np.int32(1)
        s_last = np.int32(0)
        has_zero = False
        for j in range(n):
            c = bits[j]
            if c == 0:
                has_zero = True
            cur = s_size
            s_len[cur] = s_len[s_last] + np.int32(1)
            s_link[cur] = np.int32(-1)
            s_t0[cur] = np.int32(-1)
            s_t1[cur] = np.int32(-1)
            s_size += np.int32(1)
            p = s_last
            while p != -1:
                if c == 0:
                    if s_t0[p] >= 0:
                        break
                    s_t0[p] = cur
                else:
                    if s_t1[p] >= 0:
                        break
                    s_t1[p] = cur
                p = s_link[p]
            if p == -1:
                s_link[cur] = np.int32(0)
            else:
                q = s_t0[p] if c == 0 else s_t1[p]
                if s_len[p] + 1 == s_len[q]:
                    s_link[cur] = q
                else:
                    clone = s_size
                    s_len[clone] = s_len[p] + np.int32(1)
                    s_link[clone] = s_link[q]
                    s_t0[clone] = s_t0[q]
                    s_t1[clone] = s_t1[q]
                    s_size += np.int32(1)
                    while p != -1:
                        if c == 0:
                            if s_t0[p] == q:
                                s_t0[p] = clone
                            else:
                                break
                        else:
                            if s_t1[p] == q:
                                s_t1[p] = clone
                            else:
                                break
                        p = s_link[p]
                    s_link[q] = clone
                    s_link[cur] = clone
            s_last = cur

        # Count distinct positive (1-leading) values via subtree-sum.
        max_l = s_len[s_last]
        bucket = np.zeros(max_l + 2, dtype=np.int32)
        order = np.empty(s_size, dtype=np.int32)
        cnt = np.zeros(s_size, dtype=np.int32)
        for s in range(s_size):
            bucket[s_len[s]] += np.int32(1)
        total = np.int32(0)
        for i in range(max_l, -1, -1):
            c2 = bucket[i]
            bucket[i] = total
            total += c2
        for s in range(s_size):
            order[bucket[s_len[s]]] = np.int32(s)
            bucket[s_len[s]] += np.int32(1)
        for oi in range(s_size):
            s = order[oi]
            v = np.int32(0)
            if s_t0[s] >= 0:
                v += np.int32(1) + cnt[s_t0[s]]
            if s_t1[s] >= 0:
                v += np.int32(1) + cnt[s_t1[s]]
            cnt[s] = v
        score = np.int32(0)
        if s_t1[0] >= 0:
            score += np.int32(1) + cnt[s_t1[0]]
        if has_zero:
            score += np.int32(1)
        return score

    # ──────────────────────────────────────────────────────────────────────
    # Production Numba DFS with collect_all + exact-tail + refine.
    #
    # Drop-in replacement for `_dfs_combined_worker` (single-process variant).
    # Returns best score AND the full set of optimal completions (capped by
    # `optimals_capacity` to bound memory; overflow_flag set if hit).
    #
    # Features:
    #   * Iterative DFS with stack-array representation
    #   * Lookahead UB at L = `lookahead` via `_sam_score_extensions`
    #   * Refine-path: when child UB ∈ [best, best+refine_margin] AND
    #     refine_lookahead > lookahead, recompute UB at refine_lookahead
    #   * Exact-tail solver: when n - k <= exact_tail_bits, enumerate all
    #     2^(n-k) completions exactly via `_score_bits_numba`
    #   * Parent-UB clamp (provably safe)
    #   * collect_all: stores prefix bits of every optimal full string into
    #     `optimals_buf`; on best-improvement, count resets and the new
    #     winner is stored
    # ──────────────────────────────────────────────────────────────────────

    @_numba.njit(cache=True, nogil=True)
    def _dfs_numba_collect(
        start_bits, start_len,
        n, initial_best, lookahead, global_max,
        ext_arr_L, ext_arr_Lp1, ext_arr_refine,
        refine_lookahead, refine_margin,
        known_a, max_known, sam_jmax,
        exact_tail_bits,
        optimals_buf, optimals_count_arr,
        shared_best_arr,        # int64[1] -- cross-thread incumbent. Pass a
                                # local zeros(1) for single-threaded use.
        max_nodes,              # int64: 0=unlimited; >0 cooperative chunk limit
        snap_bits_out,          # int8  [max_depth, n]  output: remaining stack bits
        snap_len_out,           # int32 [max_depth]     output: remaining stack lengths
        snap_ub_out,            # int64 [max_depth]     output: remaining stack UBs
        snap_count_arr,         # int32 [1]             output: number of snapshot entries
    ):
        """Iterative DFS B&B; collects all optimal full strings.

        Thread-safe via `shared_best_arr` polling: each thread reads the
        global best every 1024 expansions and bumps its local best if higher.
        Writes to `shared_best_arr[0]` are word-aligned int64 stores; the
        only invariant required is that the value monotonically ratchets up,
        which holds because we only ever store `max(local, shared)`.

        When max_nodes > 0, the loop exits early once nodes_expanded reaches
        the budget.  The remaining stack is written to snap_*_out and
        snap_count_arr[0] is set to the remaining stack pointer (sp).  The
        caller can then requeue those prefix/UB pairs as fresh subtasks.

        Returns (best, nodes_expanded, nodes_pruned, overflow_flag).
        """
        capacity = optimals_buf.shape[0]
        max_depth = 2 * n + 4
        stack_bits = np.empty((max_depth, n), dtype=np.int8)
        stack_len = np.empty(max_depth, dtype=np.int32)
        stack_ub = np.empty(max_depth, dtype=np.int64)

        child_pfx = np.empty(n + 1, dtype=np.int8)
        full = np.empty(n, dtype=np.int8)

        for i in range(start_len):
            stack_bits[0, i] = start_bits[i]
        stack_len[0] = start_len
        # Initial UB: we conservatively set to global_max; the popped node will
        # do the proper compute on its children (or score the leaf if k==n).
        stack_ub[0] = np.int64(global_max)

        sp = np.int32(1)
        best = np.int64(initial_best)
        # Adopt shared incumbent immediately if it's higher.
        if shared_best_arr[0] > best:
            best = shared_best_arr[0]
        nodes_expanded = np.int64(0)
        nodes_pruned = np.int64(0)
        overflow = np.int32(0)
        optimals_count_arr[0] = np.int32(0)

        ref_active = (refine_lookahead > lookahead) and (refine_margin >= 0)
        snap_count_arr[0] = np.int32(0)

        while sp > 0:
            # Cooperative chunk checkpoint: yield when budget exhausted.
            if max_nodes > np.int64(0) and nodes_expanded >= max_nodes:
                break

            # Poll shared incumbent every 1024 expansions.
            if (nodes_expanded & np.int64(1023)) == np.int64(0):
                sb = shared_best_arr[0]
                if sb > best:
                    best = sb
                    # Our locally collected optimals were valid for a lower
                    # threshold; clear them.
                    optimals_count_arr[0] = np.int32(0)

            sp -= np.int32(1)
            k = stack_len[sp]
            ub = stack_ub[sp]
            for i in range(k):
                child_pfx[i] = stack_bits[sp, i]

            if ub < best:
                nodes_pruned += np.int64(1)
                continue
            nodes_expanded += np.int64(1)

            # ── Exact-tail shortcut ─────────────────────────────────────
            m_rem_top = n - k
            if exact_tail_bits > 0 and m_rem_top > 0 and m_rem_top <= exact_tail_bits:
                # Enumerate all 2^m_rem_top suffixes; score each.
                num_suf = np.int64(1) << np.int64(m_rem_top)
                for i in range(k):
                    full[i] = child_pfx[i]
                for sfx in range(num_suf):
                    # Decode suffix bits MSB-first into full[k:].
                    bb = sfx
                    for j in range(m_rem_top - 1, -1, -1):
                        full[k + j] = np.int8(bb & 1)
                        bb >>= 1
                    sc = np.int64(_score_bits_numba(full))
                    if sc > best:
                        best = sc
                        if best > shared_best_arr[0]:
                            shared_best_arr[0] = best
                        # Reset optimals.
                        cnt = np.int32(1)
                        for i in range(n):
                            optimals_buf[0, i] = full[i]
                        optimals_count_arr[0] = cnt
                    elif sc == best:
                        cnt = optimals_count_arr[0]
                        if cnt < capacity:
                            for i in range(n):
                                optimals_buf[cnt, i] = full[i]
                            optimals_count_arr[0] = cnt + np.int32(1)
                        else:
                            overflow = np.int32(1)
                continue

            if k == n:
                # Internal leaf (only reached when exact_tail_bits == 0).
                sc = np.int64(_score_bits_numba(child_pfx[:n]))
                if sc > best:
                    best = sc
                    if best > shared_best_arr[0]:
                        shared_best_arr[0] = best
                    optimals_count_arr[0] = np.int32(1)
                    for i in range(n):
                        optimals_buf[0, i] = child_pfx[i]
                elif sc == best:
                    cnt = optimals_count_arr[0]
                    if cnt < capacity:
                        for i in range(n):
                            optimals_buf[cnt, i] = child_pfx[i]
                        optimals_count_arr[0] = cnt + np.int32(1)
                    else:
                        overflow = np.int32(1)
                continue

            # ── Expand: compute UBs for both children ───────────────────
            # When m_rem_child >= lookahead+1, one call to the split function
            # evaluates both children's UBs from the PARENT prefix.  This saves
            # one O(k) prefix-SAM rebuild compared to two per-child calls.
            # Correctness: _sam_score_extensions_split_by_first_bit with
            # (L+1)-bit extensions returns the same UB quality as two separate
            # _sam_score_extensions calls with L-bit extensions on each child.
            child_ub_0 = np.int64(-1)
            child_ub_1 = np.int64(-1)
            child_len = k + np.int32(1)
            m_rem_child = n - child_len
            needs_exact_tail = (m_rem_child == 0
                                or (exact_tail_bits > 0 and m_rem_child <= exact_tail_bits))
            use_main_L = (m_rem_child >= lookahead)
            # use_split_L: parent-prefix single call is valid only when:
            #   (a) enough remaining bits for L+1-step lookahead, AND
            #   (b) ext_arr_Lp1 actually has L+1 bits per extension
            #       (shape check guards the fallback case where Lp1 wasn't
            #        in the cache and ext_arr_L was substituted).
            use_split_L = (use_main_L
                           and (m_rem_child >= lookahead + np.int32(1))
                           and (ext_arr_Lp1.shape[1] > ext_arr_L.shape[1]))
            dyn_refine_margin = refine_margin
            if ref_active and m_rem_child <= refine_lookahead + np.int32(6):
                # Be more aggressive near the tail where tighter UB is most
                # likely to prune sibling work immediately.
                dyn_refine_margin = refine_margin + np.int32(1)
                if m_rem_child <= refine_lookahead + np.int32(2):
                    dyn_refine_margin = refine_margin + np.int32(2)

            if needs_exact_tail:
                # Near-leaf extraction: solve both child subtrees exactly now
                # instead of pushing children with loose global_max UBs.
                num_suf_child = np.int64(1) << np.int64(m_rem_child)
                for bit in range(2):
                    child_pfx[k] = np.int8(bit)
                    for i in range(child_len):
                        full[i] = child_pfx[i]
                    for sfx in range(num_suf_child):
                        bb = sfx
                        for j in range(m_rem_child - 1, -1, -1):
                            full[child_len + j] = np.int8(bb & 1)
                            bb >>= 1
                        sc = np.int64(_score_bits_numba(full))
                        if sc > best:
                            best = sc
                            if best > shared_best_arr[0]:
                                shared_best_arr[0] = best
                            cnt = np.int32(1)
                            for t in range(n):
                                optimals_buf[0, t] = full[t]
                            optimals_count_arr[0] = cnt
                        elif sc == best:
                            cnt = optimals_count_arr[0]
                            if cnt < capacity:
                                for t in range(n):
                                    optimals_buf[cnt, t] = full[t]
                                optimals_count_arr[0] = cnt + np.int32(1)
                            else:
                                overflow = np.int32(1)
                continue
            elif use_split_L:
                # Single call: build parent SAM once, split by first-extension
                # bit to get both children's UBs simultaneously.
                _pr = _sam_score_extensions_split_by_first_bit(
                    child_pfx[:k], ext_arr_Lp1, n, known_a, max_known, sam_jmax)
                _ub0 = np.int64(_pr[0])
                _ub1 = np.int64(_pr[1])
                child_ub_0 = _ub0 if _ub0 < np.int64(global_max) else np.int64(global_max)
                child_ub_1 = _ub1 if _ub1 < np.int64(global_max) else np.int64(global_max)
                # Adaptive refine for near-threshold children (uncommon path).
                if ref_active and m_rem_child >= refine_lookahead:
                    for _rb in range(2):
                        _rcu = child_ub_0 if _rb == 0 else child_ub_1
                        if best <= _rcu <= best + dyn_refine_margin:
                            child_pfx[k] = np.int8(_rb)
                            _rubf = _sam_score_extensions(
                                child_pfx[:child_len], ext_arr_refine, n,
                                known_a, max_known, _REFINE_SAM_JMAX)
                            if _rb == 0:
                                child_ub_0 = np.int64(_rubf)
                            else:
                                child_ub_1 = np.int64(_rubf)
            elif use_main_L:
                for bit in range(2):
                    child_pfx[k] = np.int8(bit)
                    ubf = _sam_score_extensions(
                        child_pfx[:child_len], ext_arr_L, n,
                        known_a, max_known, sam_jmax)
                    child_ub = np.int64(ubf)
                    if (ref_active
                            and best <= child_ub <= best + dyn_refine_margin
                            and m_rem_child >= refine_lookahead):
                        ubf2 = _sam_score_extensions(
                            child_pfx[:child_len], ext_arr_refine, n,
                            known_a, max_known, _REFINE_SAM_JMAX)
                        child_ub = np.int64(ubf2)
                    if child_ub > global_max:
                        child_ub = global_max
                    if bit == 0:
                        child_ub_0 = child_ub
                    else:
                        child_ub_1 = child_ub
            else:
                child_ub_0 = np.int64(global_max)
                child_ub_1 = np.int64(global_max)

            # Push lower-UB first so higher-UB is popped first (LIFO).
            push_first = np.int32(0) if child_ub_0 <= child_ub_1 else np.int32(1)
            for order_idx in range(2):
                bit = push_first if order_idx == 0 else (np.int32(1) - push_first)
                cu = child_ub_0 if bit == 0 else child_ub_1
                if cu < best:
                    nodes_pruned += np.int64(1)
                    continue
                for i in range(k):
                    stack_bits[sp, i] = child_pfx[i]
                stack_bits[sp, k] = np.int8(bit)
                stack_len[sp] = k + np.int32(1)
                stack_ub[sp] = cu
                sp += np.int32(1)

        # Copy remaining stack to snapshot arrays (sp > 0 means interrupted).
        snap_count_arr[0] = sp
        for _si in range(sp):
            for _sj in range(n):
                snap_bits_out[_si, _sj] = stack_bits[_si, _sj]
            snap_len_out[_si] = stack_len[_si]
            snap_ub_out[_si] = stack_ub[_si]

        return best, nodes_expanded, nodes_pruned, overflow

    _INCREMENTAL_SAM_OK = True
except Exception:
    _INCREMENTAL_SAM_OK = False
    _BASE_SAM_JMAX = 2
    _REFINE_SAM_JMAX = 3

# Fallback pure-Python evaluator (for very small n or when numba absent)
def _score_py(bits: list[int]) -> int:
    """Distinct substring-integer count (pure Python, O(n^2))."""
    n = len(bits)
    seen: set[int] = set()
    for i in range(n):
        v = 0
        for j in range(i, n):
            v = (v << 1) | bits[j]
            seen.add(v)
    return len(seen)


def _score(bits: list[int]) -> int:
    """Evaluate score for a complete or partial bit string."""
    if not bits:
        return 0
    if _eval_available:
        arr = np.array(bits, dtype=np.int8)
        result = int(_numba_eval_filtered(arr, np.int32(-1)))
        if result >= 0:
            return result
    return _score_py(bits)


# ---------------------------------------------------------------------------
# Known a(m) table -- used for the tighter "crossing" upper bound
# ---------------------------------------------------------------------------
# For any prefix p of length k and any (n-k)-bit extension e:
#   new distinct values <= a(n-k) + k*(n-k)
# where a(n-k)  bounds values from pure-extension substrings and
#       k*(n-k) bounds values from boundary-crossing substrings.
# This is always <= the simple n*(n+1)/2 - k*(k+1)/2 bound (proven since
# a(m) <= m*(m+1)/2), and is often substantially tighter.

try:
    from data.reference.known_values import KNOWN_VALUES as _kv_list
    # 1-indexed: _KNOWN_A[m] = a(m)
    _KNOWN_A: dict[int, int] = {i + 1: v for i, v in enumerate(_kv_list)}
except ImportError:
    _KNOWN_A = {}

# Load any branch_bound_exact results and merge them in (they are proven values).
# This allows subsequent searches to use these as anchors for tightening bounds.
def _load_branch_bound_exact_values() -> None:
    """Load a(n) from results/branch_bound_exact and merge into _KNOWN_A."""
    import json
    import os
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)
    )))
    exact_dir = os.path.join(project_root, "results", "branch_bound_exact")
    if not os.path.isdir(exact_dir):
        return
    for fname in os.listdir(exact_dir):
        if fname.endswith("_results.json"):
            try:
                n_str = fname.split("_")[1]
                n = int(n_str)
                fpath = os.path.join(exact_dir, fname)
                with open(fpath) as f:
                    data = json.load(f)
                val = data.get("a(n)")
                if val is not None:
                    val = int(val)
                    if n not in _KNOWN_A or val > _KNOWN_A[n]:
                        _KNOWN_A[n] = val
            except Exception:
                pass

_load_branch_bound_exact_values()

# Pre-built numpy array for the numba incremental SAM: _KNOWN_A_NP[m] = a(m) or 0.
# Index 0 unused; max valid index is len(_KNOWN_A).
_MAX_KNOWN_M = max(_KNOWN_A.keys()) if _KNOWN_A else 0
_KNOWN_A_NP = np.zeros(_MAX_KNOWN_M + 1, dtype=np.int64)
for _m, _v in _KNOWN_A.items():
    if _m <= _MAX_KNOWN_M:
        _KNOWN_A_NP[_m] = _v

# ---------------------------------------------------------------------------
# MH incumbent loader -- valid lower bound from stored heuristic results
# ---------------------------------------------------------------------------

def load_mh_incumbent(n: int) -> int:
    """Return the best MH score for bit-length *n*, or 0 if none found.

    Reads from ``results/mh_bounded/`` and ``results/mh_unbounded/``.
    The returned value is always <= a(n) (it is the score of an actual string),
    so it is a valid lower-bound / incumbent hint for branch-and-bound.
    """
    import json
    import os

    project_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)
    )))
    filename = f"n_{n:04d}_results.json"
    best = 0
    for subdir in ("mh_bounded", "mh_unbounded"):
        path = os.path.join(project_root, "results", subdir, filename)
        if os.path.isfile(path):
            try:
                with open(path) as f:
                    data = json.load(f)
                bv = int(data.get("best_value", 0))
                if bv > best:
                    best = bv
            except Exception:
                pass
    return best


# ---------------------------------------------------------------------------
# Upper bound helpers
# ---------------------------------------------------------------------------

def _remaining(prefix: tuple[int, ...], n: int) -> int:
    """Upper bound on additional distinct values for any (n-k)-bit suffix.

    Uses cross = ones * m, suffix = a(m).  No SAM savings at this level.
    Falls back to the simple formula when a(m) is unknown.
    """
    k = len(prefix)
    m = n - k
    ones = sum(prefix)          # popcount -- number of '1' bits in prefix
    if m in _KNOWN_A:
        cross = ones * m
        return _KNOWN_A[m] + cross
    # Fallback: simple substring count (loose but always valid)
    return n * (n + 1) // 2 - k * (k + 1) // 2


# Pre-build sorted extension tables at module load to avoid rebuilding
# Python tuples or NumPy arrays in the hot loop.
_EXT_CACHE: dict[int, list[tuple[int, ...]]] = {}
_EXT_CACHE_NP: dict[int, np.ndarray] = {}
for _L in range(1, 11):
    # (1,0) order first so all-ones is evaluated first (sets best high early)
    _exts = [tuple(int(bit) for bit in ext) for ext in _product((1, 0), repeat=_L)]
    _EXT_CACHE[_L] = _exts
    _EXT_CACHE_NP[_L] = np.array(_exts, dtype=np.int8)


# Solve the final few bits exactly instead of paying the generic UB cost on
# every node in a tiny residual subtree.
_EXACT_TAIL_BITS = 8
_EXACT_BOUND_BITS = 9
_BATCH_TAIL_SOLVER = True


def _exact_completion_best_score(prefix: tuple[int, ...], n: int) -> int:
    """Return the exact best score over all completions of this prefix."""
    m_rem = n - len(prefix)
    if m_rem < 0:
        raise ValueError("prefix longer than target length")
    if m_rem == 0:
        return _score(list(prefix))

    if _INCREMENTAL_SAM_OK:
        prefix_arr = np.array(prefix, dtype=np.int8)
        ext_arr = _EXT_CACHE_NP.get(m_rem)
        if ext_arr is None:
            ext_arr = np.array(list(_product((1, 0), repeat=m_rem)), dtype=np.int8)
        return int(_sam_score_extensions(
            prefix_arr, ext_arr, n, _KNOWN_A_NP, _MAX_KNOWN_M, _BASE_SAM_JMAX,
        ))

    suffixes = _EXT_CACHE.get(m_rem)
    if suffixes is None:
        suffixes = [tuple(int(bit) for bit in ext) for ext in _product((1, 0), repeat=m_rem)]
    prefix_list = list(prefix)
    best = 0
    for suffix in suffixes:
        score = _score(prefix_list + list(suffix))
        if score > best:
            best = score
    return best


def _solve_small_suffix_exact(
    prefix: tuple[int, ...],
    n: int,
    incumbent: int,
) -> tuple[int, list[str]]:
    """Return the exact optimum over all completions of a short suffix."""
    m_rem = n - len(prefix)
    if m_rem < 0:
        raise ValueError("prefix longer than target length")

    suffixes = _EXT_CACHE.get(m_rem)
    if suffixes is None:
        suffixes = [tuple(int(bit) for bit in ext) for ext in _product((1, 0), repeat=m_rem)]

    best = _exact_completion_best_score(prefix, n)
    if best < incumbent:
        return best, []

    prefix_str = "".join(map(str, prefix))

    if (_BATCH_TAIL_SOLVER and _eval_available and _numba_eval_batch is not None
            and m_rem > 0):
        prefix_arr = np.array(prefix, dtype=np.int8)
        suffix_arr = _EXT_CACHE_NP.get(m_rem)
        if suffix_arr is None:
            suffix_arr = np.array(list(_product((1, 0), repeat=m_rem)), dtype=np.int8)
        bits_batch = np.empty((suffix_arr.shape[0], n), dtype=np.int8)
        bits_batch[:, :len(prefix)] = prefix_arr[None, :]
        bits_batch[:, len(prefix):] = suffix_arr
        counts = _numba_eval_batch(bits_batch, np.int32(-1), len(prefix))
        optimals = [
            prefix_str + "".join(str(int(bit)) for bit in suffixes[idx])
            for idx, score in enumerate(counts)
            if int(score) == best
        ]
        return best, optimals

    optimals: list[str] = []
    prefix_list = list(prefix)
    for suffix in suffixes:
        full = prefix_list + list(suffix)
        score = _score(full)
        if score == best:
            optimals.append(prefix_str + "".join(map(str, suffix)))
    return best, optimals


def _lookahead_ub(prefix: tuple[int, ...], n: int, L: int, global_max: int) -> int:
    """L-step lookahead upper bound for a prefix at depth k = len(prefix).

    When _INCREMENTAL_SAM_OK, builds the prefix SAM once and evaluates all
    2^L extensions by cloning it -- O(k + 2^L * L) instead of O(2^L * (k+L)).

    For L=0 this reduces to  min(global_max, score(prefix) + remaining(prefix, n)).
    """
    k = len(prefix)
    if 0 < _EXACT_BOUND_BITS and n - k <= _EXACT_BOUND_BITS:
        return min(_exact_completion_best_score(prefix, n), global_max)

    L_use = min(L, n - k)

    if L_use == 0:
        # Either L=0 or prefix is already full-length -> exact score.
        s = _score(list(prefix))
        return min(s, global_max)           # remaining == 0 when k==n

    # ── Fast path: incremental SAM batch evaluation ───────────────────────
    if _INCREMENTAL_SAM_OK:
        prefix_arr = np.array(prefix, dtype=np.int8)
        ext_arr = _EXT_CACHE_NP.get(L_use)
        if ext_arr is None:
            ext_arr = np.array(list(_product((1, 0), repeat=L_use)), dtype=np.int8)
        best = int(_sam_score_extensions(
            prefix_arr, ext_arr, n, _KNOWN_A_NP, _MAX_KNOWN_M, _BASE_SAM_JMAX,
        ))
        return min(best, global_max)

    # ── Fallback: per-extension _score calls (original path) ──────────────
    ext_len = k + L_use
    m_rem = n - ext_len
    a_ext_len = _KNOWN_A.get(ext_len, ext_len * (ext_len + 1) // 2)
    a_m_rem = (_KNOWN_A[m_rem] if m_rem in _KNOWN_A
               else n * (n + 1) // 2 - ext_len * (ext_len + 1) // 2) if m_rem > 0 else 0
    prefix_ones = sum(prefix)

    best = 0
    for ext in (_EXT_CACHE.get(L_use) or list(_product((1, 0), repeat=L_use))):
        ext_ones = prefix_ones + sum(ext)
        if a_ext_len + a_m_rem + ext_ones * m_rem <= best:
            continue
        extended = prefix + ext
        s = _score(list(extended))
        rem = _remaining(extended, n)
        ub = s + rem
        if ub > best:
            best = ub
            if best >= global_max:
                return global_max
    return min(best, global_max)


def _refined_lookahead_ub(
    prefix: tuple[int, ...],
    n: int,
    refine_lookahead: int,
    global_max: int,
) -> int:
    """Return the tighter refinement UB used near the incumbent threshold."""
    k = len(prefix)
    L_use = min(refine_lookahead, n - k)
    if L_use == 0:
        return _lookahead_ub(prefix, n, refine_lookahead, global_max)

    if _INCREMENTAL_SAM_OK:
        prefix_arr = np.array(prefix, dtype=np.int8)
        ext_arr = _EXT_CACHE_NP.get(L_use)
        if ext_arr is None:
            ext_arr = np.array(list(_product((1, 0), repeat=L_use)), dtype=np.int8)
        best = int(_sam_score_extensions(
            prefix_arr, ext_arr, n, _KNOWN_A_NP, _MAX_KNOWN_M, _REFINE_SAM_JMAX,
        ))
        return min(best, global_max)

    return _lookahead_ub(prefix, n, refine_lookahead, global_max)


def _adaptive_lookahead_ub(
    prefix: tuple[int, ...],
    n: int,
    lookahead: int,
    global_max: int,
    target: int,
    refine_lookahead: int = -1,
    refine_margin: int = 2,
) -> int:
    """Target-aware UB cascade for DFS branch-and-bound.

    Compute the base lookahead UB first. Only if that UB is close to the
    current target do we pay for a tighter refinement pass. The returned value
    is still a valid UB because both passes independently compute rigorous UBs.
    """
    ub = _lookahead_ub(prefix, n, lookahead, global_max)
    if (refine_lookahead <= lookahead or refine_margin < 0
            or ub < target or ub > target + refine_margin):
        return ub
    return _refined_lookahead_ub(prefix, n, refine_lookahead, global_max)


# ---------------------------------------------------------------------------
# branch-and-bound search
# ---------------------------------------------------------------------------

def branch_bound_a112509(
    n: int,
    known_sequence: Optional[dict] = None,
    incumbent_hint: int = 0,
    lookahead: int = 1,
    verbose: bool = True,
    _start_prefix: Optional[tuple] = None,
    collect_all: bool = False,
    max_nodes: int = 0,
    _initial_heap: Optional[list] = None,
) -> tuple[int, list[str], dict]:
    """
    Compute a(n) exactly via branch-and-bound-search with L-step lookahead UB.

    Parameters
    ----------
    n : int
        Bit-string length (1-indexed per OEIS).
    known_sequence : dict, optional
        {n: a(n)} for seeding the incumbent (lower bound).
    incumbent_hint : int
        Additional floor for the incumbent (e.g. from a prior heuristic run).
    lookahead : int
        L -- number of extra bits to evaluate when computing the UB for each
        child.  L=0 uses the plain score+remaining bound; L=1 evaluates both
        grandchildren; L=2 evaluates all 4 great-grandchildren; etc.
        Higher L gives a tighter (provably smaller) UB and fewer nodes, at
        the cost of 2^L evaluations per child push.
    verbose : bool
        Print progress information.
    _start_prefix : tuple, optional
        Internal parameter for parallel search.  When set, the heap is
        initialised with this single prefix instead of the default (1,).
        The caller is responsible for ensuring the prefix has leading bit 1
        and UB >= incumbent.
    collect_all : bool
        When False (default), Rule 2 prunes any branch whose UB equals the
        current best once at least one optimal string has been found -- this
        is fast but may miss optimal strings.  When True, Rule 2 is disabled
        so every string achieving a(n) is collected, at the cost of exploring
        all UB==a(n) branches exhaustively.

    Returns
    -------
    (value, optimal_strings, stats)
        value           - a(n)
        optimal_strings - all n-bit strings achieving a(n)
        stats           - dict with search statistics
    """
    t0 = time.perf_counter()

    if n == 0:
        return 0, [], {"elapsed_seconds": 0.0, "nodes_expanded": 0,
                       "nodes_pruned": 0, "lookahead": lookahead}

    # ── Incumbent lower bound ────────────────────────────────────────────────
    incumbent = max(0, incumbent_hint)
    optimal_strings: list[str] = []
    if known_sequence and n in known_sequence:
        val = known_sequence[n]
        if val > incumbent:
            incumbent = val
        if verbose:
            print(f"  [branch-and-bound-L{lookahead}] n={n}: incumbent seeded = {incumbent}")

    # ── Global UB cap ────────────────────────────────────────────────────────
    global_max = _compute_global_max(n)

    # ── Special case n=1 ────────────────────────────────────────────────────
    if n == 1:
        s0 = _score([0]); s1 = _score([1])
        best = max(s0, s1)
        sols = ["0"] * (s0 == best) + ["1"] * (s1 == best)
        elapsed = time.perf_counter() - t0
        return best, sols, {"elapsed_seconds": round(elapsed, 4),
                            "nodes_expanded": 2, "nodes_pruned": 0,
                            "lookahead": lookahead}

    # ── Initialise heap ──────────────────────────────────────────────────────
    # All n-bit numbers start with 1 (leading bit convention from brute_force).
    # _start_prefix overrides the default root for parallel sub-tree search.
    # _initial_heap overrides both: resume from a saved heap (work-stealing bulk).
    heap: list[tuple[int, int, tuple[int, ...]]] = []
    if _initial_heap is not None:
        # Resume from a bulk continuation returned by an interrupted worker.
        # Entries are (ub, prefix); convert to the heap's (-ub, tiebreak, prefix) form.
        for _ub, _pfx in _initial_heap:
            heapq.heappush(heap, (-_ub, 0, _pfx))
    else:
        init_prefix: tuple[int, ...] = _start_prefix if _start_prefix is not None else (1,)
        init_ub = _lookahead_ub(init_prefix, n, lookahead, global_max)
        heapq.heappush(heap, (-init_ub, 0, init_prefix))

    nodes_expanded = 0
    nodes_pruned = 0
    evals_for_ub = 0          # extra evaluations spent on lookahead UBs

    while heap:
        neg_ub, _, prefix = heapq.heappop(heap)
        ub = -neg_ub
        k = len(prefix)

        # Prune: cannot improve on the incumbent.
        # Rule 1: UB < incumbent -> no string in this subtree can beat or match
        #          incumbent, so prune unconditionally.
        # Rule 2: UB == incumbent AND we already have at least one optimal string
        #          -> prune to avoid collecting duplicate optimals (optional).
        #          Disabled when collect_all=True to enumerate every optimum.
        if ub < incumbent or (not collect_all and ub == incumbent and optimal_strings):
            nodes_pruned += 1
            continue

        nodes_expanded += 1

        # Interrupt for dynamic work-stealing (parallel use with max_nodes > 0).
        # Bounded fan-out: top-K entries become individual parallel tasks;
        # the remainder is packaged as ONE bulk continuation task.
        # This guarantees correctness (no paths dropped) with bounded overhead.
        if 0 < max_nodes <= nodes_expanded:
            _SPLIT_K = 50
            # Sort heap to separate top-K (highest priority) from the rest.
            # heap is a min-heap with negated UBs, so nsmallest = highest UB.
            all_heap_sorted = sorted(heap)  # (-ub, tiebreak, prefix)
            top_k = all_heap_sorted[:_SPLIT_K]
            bulk  = all_heap_sorted[_SPLIT_K:]
            heap_snapshot = [(-neg_ub, pfx) for neg_ub, _, pfx in top_k]
            heap_bulk     = [(-neg_ub, pfx) for neg_ub, _, pfx in bulk]
            elapsed = time.perf_counter() - t0
            stats = {
                "elapsed_seconds": round(elapsed, 4),
                "nodes_expanded": nodes_expanded,
                "nodes_pruned": nodes_pruned,
                "evals_for_ub": evals_for_ub,
                "incumbent": incumbent,
                "lookahead": lookahead,
                "interrupted": True,
                "heap_snapshot": heap_snapshot,
                "heap_bulk": heap_bulk,
            }
            return incumbent, optimal_strings, stats

        if k == n:
            # Full string: UB already equals exact score (remaining=0, L_use=0).
            score = ub          # no need to re-evaluate
            s = "".join(map(str, prefix))
            if score > incumbent:
                incumbent = score
                optimal_strings = [s]
                if verbose:
                    label = s[:40] + ("..." if n > 40 else "")
                    print(f"  [branch-and-bound-L{lookahead}] n={n}: new best {incumbent} '{label}'")
            elif score == incumbent:
                optimal_strings.append(s)
            continue

        # Expand: append 0 and 1, compute L-step lookahead UB for each child.
        m_child = n - k - 1
        a_child_k = _KNOWN_A.get(k + 1, (k + 1) * (k + 2) // 2)
        a_child_m = (_KNOWN_A[m_child] if m_child in _KNOWN_A
                     else n * (n + 1) // 2 - (k + 1) * (k + 2) // 2) if m_child > 0 else 0
        ones = sum(prefix)
        for bit in (0, 1):
            child_ones = ones + bit
            # Fast pre-check: a(k+1) bounds score(child); skip _lookahead_ub if hopeless.
            if a_child_k + a_child_m + child_ones * m_child < incumbent:
                nodes_pruned += 1
                continue
            child = prefix + (bit,)
            child_ub = _lookahead_ub(child, n, lookahead, global_max)
            evals_for_ub += 1 << min(lookahead, n - k - 1)
            if child_ub > ub:
                child_ub = ub  # parent-UB clamp (provably safe; almost free)
            if child_ub < incumbent or (child_ub == incumbent and optimal_strings):
                nodes_pruned += 1
                continue
            heapq.heappush(heap, (-child_ub, nodes_expanded, child))

    elapsed = time.perf_counter() - t0
    stats = {
        "elapsed_seconds": round(elapsed, 4),
        "nodes_expanded": nodes_expanded,
        "nodes_pruned": nodes_pruned,
        "evals_for_ub": evals_for_ub,
        "incumbent": incumbent,
        "lookahead": lookahead,
        "interrupted": False,
    }
    if verbose:
        print(f"  [branch-and-bound-L{lookahead}] n={n}: a(n)={incumbent}, "
              f"{nodes_expanded:,} expanded, {nodes_pruned:,} pruned, "
              f"{evals_for_ub:,} UB-evals, {elapsed:.3f}s")
    return incumbent, optimal_strings, stats


# ---------------------------------------------------------------------------
# Parallel branch-and-bound search
# ---------------------------------------------------------------------------

# Module-level shared incumbent for parallel workers.
# Set by _init_shared_worker before pool starts.
_worker_shared_incumbent = None


def _init_shared_worker(shared_val):
    """Pool initializer: bind the shared incumbent Value into this worker."""
    global _worker_shared_incumbent
    _worker_shared_incumbent = shared_val




def _dfs_combined_worker_numba(
    start_prefix: tuple[int, ...],
    n: int,
    incumbent: int,
    lookahead: int,
    refine_lookahead: int,
    refine_margin: int,
    max_nodes: int = 0,         # >0 enables cooperative chunking
) -> tuple[int, list[str], dict]:
    """Numba-JIT'd inner loop for `_dfs_combined_worker`.

    Drop-in replacement that:
      * Reads the multiprocessing shared incumbent (`_worker_shared_incumbent`)
        once before launching, and merges it into the seed `best`.
      * Runs the JIT'd DFS with GIL released (`@njit(nogil=True)`).
      * Concurrently runs a Python shuttle thread that bridges the
        multiprocessing `Value('l')` ↔ a numpy int64[1] every ~0.2s, so
        the Numba loop sees other workers' improvements via its periodic
        poll of `shared_best_arr`.

    Returns (best, optimals, stats) in the same shape as `_dfs_combined_worker`.
    """
    import threading as _threading

    global_max = _compute_global_max(n)

    # Resolve refine cache (may equal main cache if refine disabled).
    if refine_lookahead > lookahead and refine_lookahead in _EXT_CACHE_NP:
        ext_arr_refine = _EXT_CACHE_NP[refine_lookahead]
        eff_refine_L = refine_lookahead
        eff_refine_margin = refine_margin
    else:
        ext_arr_refine = _EXT_CACHE_NP[lookahead]
        eff_refine_L = lookahead
        eff_refine_margin = -1  # disable refine
    ext_arr_L = _EXT_CACHE_NP[lookahead]
    # ext_arr_Lp1: extension array for L+1 bits, used to evaluate both children
    # from the PARENT prefix in one SAM call (saves one O(k) SAM rebuild/node).
    ext_arr_Lp1 = _EXT_CACHE_NP.get(lookahead + 1, ext_arr_L)

    # Initial best: max(incumbent, shared_incumbent).
    seed_best = int(incumbent)
    if _worker_shared_incumbent is not None:
        sv = int(_worker_shared_incumbent.value)
        if sv > seed_best:
            seed_best = sv

    # Shared int64[1] for cross-thread visibility from Numba.
    shared_best_arr = np.zeros(1, dtype=np.int64)
    shared_best_arr[0] = np.int64(seed_best)

    # Optimals buffer (capacity matches `branch_bound_a112509_numba` default).
    optimals_capacity = 1 << 15
    optimals_buf = np.zeros((optimals_capacity, n), dtype=np.int8)
    optimals_count_arr = np.zeros(1, dtype=np.int32)

    # Snapshot buffers for cooperative chunking (used when max_nodes > 0).
    max_depth = 2 * n + 4
    snap_bits_out  = np.zeros((max_depth, n), dtype=np.int8)
    snap_len_out   = np.zeros(max_depth, dtype=np.int32)
    snap_ub_out    = np.zeros(max_depth, dtype=np.int64)
    snap_count_arr = np.zeros(1, dtype=np.int32)

    # Convert start_prefix to int8 array.
    start_arr = np.asarray(start_prefix, dtype=np.int8)
    start_len = np.int32(len(start_prefix))

    # ── Shuttle thread ────────────────────────────────────────────────────
    # Bridges the mp.Value('l') and shared_best_arr[0] periodically.
    stop_evt = _threading.Event()

    def _shuttle():
        if _worker_shared_incumbent is None:
            return
        while not stop_evt.wait(0.2):
            try:
                mp_val = int(_worker_shared_incumbent.value)
                arr_val = int(shared_best_arr[0])
                if mp_val > arr_val:
                    shared_best_arr[0] = np.int64(mp_val)
                elif arr_val > mp_val:
                    with _worker_shared_incumbent.get_lock():
                        if arr_val > _worker_shared_incumbent.value:
                            _worker_shared_incumbent.value = arr_val
            except Exception:
                # Process is shutting down or value detached; bail out.
                return

    shuttle_t = None
    if _worker_shared_incumbent is not None:
        shuttle_t = _threading.Thread(target=_shuttle, daemon=True)
        shuttle_t.start()

    try:
        best, nodes_expanded, nodes_pruned, overflow = _dfs_numba_collect(
            start_arr, start_len,
            np.int32(n), np.int64(seed_best),
            np.int32(lookahead), np.int64(global_max),
            ext_arr_L, ext_arr_Lp1, ext_arr_refine,
            np.int32(eff_refine_L), np.int32(eff_refine_margin),
            _KNOWN_A_NP, _MAX_KNOWN_M, _BASE_SAM_JMAX,
            np.int32(_EXACT_TAIL_BITS),
            optimals_buf, optimals_count_arr,
            shared_best_arr,
            np.int64(max_nodes),
            snap_bits_out, snap_len_out, snap_ub_out, snap_count_arr,
        )
    finally:
        stop_evt.set()
        if shuttle_t is not None:
            shuttle_t.join(timeout=1.0)

    # Final flush of best to shared incumbent.
    if _worker_shared_incumbent is not None:
        try:
            with _worker_shared_incumbent.get_lock():
                if int(best) > _worker_shared_incumbent.value:
                    _worker_shared_incumbent.value = int(best)
        except Exception:
            pass

    cnt = int(optimals_count_arr[0])
    optimal_strings = [
        "".join(str(int(b)) for b in optimals_buf[i])
        for i in range(cnt)
    ]

    # Build stack snapshot if interrupted (cooperative chunk return).
    snap_count = int(snap_count_arr[0])
    stack_snapshot: list[tuple[tuple[int, ...], int]] = []
    if snap_count > 0:
        cur_best = int(best)
        for _si in range(snap_count):
            _slen = int(snap_len_out[_si])
            _pfx  = tuple(int(snap_bits_out[_si, _j]) for _j in range(_slen))
            _ub   = int(snap_ub_out[_si])
            if _ub >= cur_best:
                stack_snapshot.append((_pfx, _ub))

    return int(best), optimal_strings, {
        "nodes_expanded": int(nodes_expanded),
        "nodes_pruned":   int(nodes_pruned),
        "interrupted":    snap_count > 0,
        "overflow":       bool(overflow),
        "stack_snapshot": stack_snapshot,
    }

def _dfs_combined_worker(args: tuple) -> tuple[int, list[str], dict]:
    """Multiprocessing worker: iterative DFS branch-and-bound.

    Finds the maximum score in this subtask AND collects all optimal strings
    achieving that maximum.  Uses O(n) stack depth — no heap — so it is
    immune to the heap-OOM that afflicts branch-and-bound on large hard subtasks.

    Each stack entry stores (prefix, precomputed_ub) so _lookahead_ub is
    called exactly ONCE per node (at push time), not twice (push + pop).
    This halves the dominant SAM evaluation cost across all expanded nodes.

    Correctness: `best` can only increase (never decrease) during the search.
    A stale UB stored at push time is always >= the true UB (it's the same
    function, same prefix — UB is immutable).  The pop-time check `ub < best`
    still correctly prunes nodes whose UB is below the current (raised) best.

    args: (start_prefix, n, incumbent, lookahead, refine_lookahead,
           refine_margin, max_nodes, initial_stack)
    Returns: (best_value, optimal_strings, stats_dict)
    """
    start_prefix, n, incumbent, lookahead, refine_lookahead, refine_margin, max_nodes, initial_stack = args

    # ── Fast path: JIT'd Numba DFS ────────────────────────────────────────
    # Used when not resuming from a stolen stack snapshot.  Empirical 11–15×
    # per-process speedup vs the Python loop below, with full feature parity
    # (collect_all, exact-tail, refine, cooperative chunking via max_nodes).
    # The Python path is used only when initial_stack is supplied (resumed
    # snapshot tasks, which are already small).
    if (_INCREMENTAL_SAM_OK
            and initial_stack is None
            and lookahead in _EXT_CACHE_NP
            and (refine_lookahead in _EXT_CACHE_NP if refine_lookahead > lookahead else True)):
        return _dfs_combined_worker_numba(
            start_prefix, n, incumbent,
            lookahead, refine_lookahead, refine_margin,
            max_nodes=max_nodes,
        )

    global_max = _compute_global_max(n)

    best = incumbent
    if _worker_shared_incumbent is not None:
        shared_val = _worker_shared_incumbent.value
        if shared_val > best:
            best = shared_val
    optimal_strings: list[str] = []
    nodes_expanded = 0
    nodes_pruned = 0
    # Cache UB evaluations inside this worker invocation. This is most useful
    # for resumed/chunked tasks where local frontier prefixes can reappear.
    ub_cache: dict[tuple[tuple[int, ...], int], int] = {}
    _UB_CACHE_MAX = 500_000

    def _cached_adaptive_ub(prefix_bits: tuple[int, ...], best_floor: int) -> int:
        cache_key = (prefix_bits, best_floor)
        cached = ub_cache.get(cache_key)
        if cached is not None:
            return cached
        val = _adaptive_lookahead_ub(
            prefix_bits, n, lookahead, global_max, best_floor,
            refine_lookahead=refine_lookahead,
            refine_margin=refine_margin,
        )
        if len(ub_cache) >= _UB_CACHE_MAX:
            ub_cache.clear()
        ub_cache[cache_key] = val
        return val

    # Stack stores (prefix, ub) — UB precomputed at push time, not recomputed on pop.
    if initial_stack is not None:
        stack: list[tuple[tuple[int, ...], int]] = [
            (prefix, ub) for prefix, ub in initial_stack if ub >= best
        ]
    else:
        init_ub = _cached_adaptive_ub(start_prefix, best)
        stack = [(start_prefix, init_ub)]

    while stack:
        if _worker_shared_incumbent is not None and (nodes_expanded & 255) == 0:
            shared_val = _worker_shared_incumbent.value
            if shared_val > best:
                best = shared_val
                optimal_strings = []

        prefix, ub = stack.pop()

        if ub < best:
            nodes_pruned += 1
            continue

        k = len(prefix)
        nodes_expanded += 1

        if k == n:
            # Leaf: UB == exact score (no remaining bits, m_rem == 0).
            if ub > best:
                best = ub
                if _worker_shared_incumbent is not None:
                    with _worker_shared_incumbent.get_lock():
                        if best > _worker_shared_incumbent.value:
                            _worker_shared_incumbent.value = best
                optimal_strings = ["".join(map(str, prefix))]
            elif ub == best:
                optimal_strings.append("".join(map(str, prefix)))
            continue

        if 0 < _EXACT_TAIL_BITS and n - k <= _EXACT_TAIL_BITS:
            exact_best, exact_optimals = _solve_small_suffix_exact(prefix, n, best)
            if exact_best > best:
                best = exact_best
                if _worker_shared_incumbent is not None:
                    with _worker_shared_incumbent.get_lock():
                        if best > _worker_shared_incumbent.value:
                            _worker_shared_incumbent.value = best
                optimal_strings = exact_optimals
            elif exact_best == best:
                optimal_strings.extend(exact_optimals)
            continue

        # Compute both child UBs once, then push lower-UB children first so
        # the highest-UB child is popped and explored next (LIFO). This tends
        # to raise `best` earlier and improves pruning in the remaining stack.
        #
        # Parent-UB clamp: any descendant of `prefix` is also a descendant of
        # `prefix`'s parent, so the parent's UB (== `ub` here) is a valid
        # upper bound for every child. We clamp child_ub to ub. This is
        # provably tight: the lookahead UB recomputed at the child can be
        # *higher* than the parent's UB when the lookahead's longer-horizon
        # split sees something the parent's split missed; in that case the
        # parent's UB is the better (tighter) bound. Always-safe, almost free.
        children: list[tuple[int, tuple[int, ...]]] = []
        # Dynamic refine trigger: tighter toward the tail where near-threshold
        # nodes are more likely to be cut by a one-step deeper bound.
        dyn_refine_margin = refine_margin
        m_rem_child = n - (k + 1)
        if refine_lookahead > lookahead and refine_margin >= 0:
            if m_rem_child <= refine_lookahead + 2:
                dyn_refine_margin = refine_margin + 2
            elif m_rem_child <= refine_lookahead + 6:
                dyn_refine_margin = refine_margin + 1
        for bit in (0, 1):
            child = prefix + (bit,)
            child_ub = _cached_adaptive_ub(child, best)
            if child_ub >= best:
                children.append((child_ub, child))
            else:
                nodes_pruned += 1
        children.sort(key=lambda item: item[0])
        for child_ub, child in children:
            stack.append((child, child_ub))

        # Cooperative deferred split: return pending sibling subtrees so idle
        # workers can steal them from the external queue.
        if 0 < max_nodes <= nodes_expanded and stack:
            if _worker_shared_incumbent is not None:
                shared_val = _worker_shared_incumbent.value
                if shared_val > best:
                    best = shared_val
                    optimal_strings = []
            stack_snapshot = [(pfx, child_ub) for pfx, child_ub in stack if child_ub >= best]
            stack_snapshot.sort(key=lambda item: item[1], reverse=True)
            return best, optimal_strings, {
                "nodes_expanded": nodes_expanded,
                "nodes_pruned": nodes_pruned,
                "interrupted": True,
                "stack_snapshot": stack_snapshot,
            }

    return best, optimal_strings, {
        "nodes_expanded": nodes_expanded,
        "nodes_pruned": nodes_pruned,
        "interrupted": False,
    }


def _seed_incumbent_from_known(n: int) -> int:
    """Return a valid lower bound for a(n) minus 1.

    Uses a(n-1)-1 when available. The minus-1 ensures that if a(n)=a(n-1)
    (would be unusual; empirically doesn't happen for n>=1) we still detect
    the optimum -- the DFS uses `score > best` so we must seed strictly below
    the true optimum to collect any optima at all.
    """
    try:
        from data.reference.known_values import KNOWN_VALUES as _kv
        if n - 2 >= 0 and (n - 2) < len(_kv):
            return max(0, int(_kv[n - 2]) - 1)
    except Exception:
        pass
    return 0


def _resolve_incumbent(
    n: int,
    known_sequence: dict | None,
    incumbent_hint: int,
    use_mh_hint: bool = True,
    verbose: bool = False,
    label: str = "",
) -> int:
    """Compute the best safe incumbent lower bound from all available sources.

    Sources (applied in order, each can only raise the floor):
      1. ``incumbent_hint`` -- caller-supplied floor
      2. ``known_sequence[n]`` -- exact a(n) if available
      3. MH heuristic results -- best score from MH runs (always <= a(n))
      4. ``a(n-1) - 1`` -- monotonicity of a(n) guarantees this is safe
    """
    incumbent = max(0, incumbent_hint)
    if known_sequence and n in known_sequence:
        val = known_sequence[n]
        if val > incumbent:
            incumbent = val
    if use_mh_hint and (not known_sequence or n not in (known_sequence or {})):
        mh_val = load_mh_incumbent(n)
        if mh_val > incumbent:
            incumbent = mh_val
            if verbose:
                print(f"[{label}] n={n}: MH incumbent hint = {incumbent}")
    if (not known_sequence) or (n not in (known_sequence or {})):
        seed = _seed_incumbent_from_known(n)
        if seed > incumbent:
            incumbent = seed
            if verbose:
                print(f"[{label}] n={n}: a(n-1)-1 seed incumbent = {incumbent}")
    return incumbent


def _branch_bound_from_prefix(args: tuple) -> tuple[int, list[str], dict]:
    """Multiprocessing worker: run branch-and-bound from a fixed starting prefix.

    Used when collect_all=False (standard value-finding path).
    Not used when collect_all=True — _dfs_combined_worker handles that case.

    args is an 8-tuple:
        (start_prefix, n, incumbent, known_seq, lookahead, collect_all,
         max_nodes, initial_heap)
    start_prefix is None when initial_heap is provided (bulk continuation).
    initial_heap is None for normal prefix-rooted tasks.
    """
    start_prefix, n, incumbent, known_seq, lookahead, collect_all, max_nodes, initial_heap = args
    # If a shared incumbent is available, use its current value as the initial
    # floor -- another worker may have already raised the bar.
    if _worker_shared_incumbent is not None:
        shared_val = _worker_shared_incumbent.value
        if shared_val > incumbent:
            incumbent = shared_val
    result = branch_bound_a112509(
        n,
        known_sequence=known_seq,
        incumbent_hint=incumbent,
        lookahead=lookahead,
        verbose=False,
        _start_prefix=start_prefix,
        collect_all=collect_all,
        max_nodes=max_nodes,
        _initial_heap=initial_heap,
    )
    # Publish our best found value back into the shared counter.
    if _worker_shared_incumbent is not None:
        found_val = result[0]
        with _worker_shared_incumbent.get_lock():
            if found_val > _worker_shared_incumbent.value:
                _worker_shared_incumbent.value = found_val
    return result


def _enum_split_nodes(
    n: int,
    split_depth: int,
    incumbent: int,
    known_sequence: dict,
    lookahead: int,
    global_max: int,
) -> list[tuple[tuple[int, ...], int]]:
    """BFS to enumerate all prefixes at `split_depth` with UB >= incumbent.

    Returns (prefix, ub) pairs. Branches whose UB falls below the incumbent
    are pruned immediately, so the frontier stays small.
    """
    split_depth = min(split_depth, n - 1)
    frontier: list[tuple[int, ...]] = [(1,)]
    for _ in range(split_depth - 1):
        new_frontier: list[tuple[int, ...]] = []
        for prefix in frontier:
            for bit in (0, 1):
                child = prefix + (bit,)
                child_ub = _lookahead_ub(child, n, lookahead, global_max)
                if child_ub >= incumbent:
                    new_frontier.append(child)
        frontier = new_frontier

    # Pair each frontier node with its UB for adaptive re-splitting downstream.
    return [
        (prefix, _lookahead_ub(prefix, n, lookahead, global_max))
        for prefix in frontier
    ]


def _bfs_expand_prefix(
    prefix: tuple[int, ...],
    n: int,
    incumbent: int,
    lookahead: int,
    extra_levels: int,
    global_max: int,
) -> list[tuple[int, ...]]:
    """BFS-expand a prefix by extra_levels more levels, pruning by UB < incumbent.

    Used for tail re-splitting: converts one coarse DFS task into several
    finer ones so all workers stay busy through the expensive run tail.
    """
    frontier: list[tuple[int, ...]] = [prefix]
    for _ in range(extra_levels):
        new_frontier: list[tuple[int, ...]] = []
        for pfx in frontier:
            if len(pfx) >= n - 1:
                new_frontier.append(pfx)
                continue
            for bit in (0, 1):
                child = pfx + (bit,)
                if _lookahead_ub(child, n, lookahead, global_max) >= incumbent:
                    new_frontier.append(child)
        if new_frontier:
            frontier = new_frontier
    return frontier


def _probe_subtree_hardness(
    prefix: tuple[int, ...],
    n: int,
    incumbent: int,
    lookahead: int,
    global_max: int,
    refine_lookahead: int = -1,
    refine_margin: int = 2,
    max_pops: int = 200,
    near_margin: int = 2,
) -> dict[str, float | int | bool]:
    """Run a tiny bounded DFS probe to estimate subtree hardness.

    The probe intentionally mirrors the real DFS worker's expansion order so
    its residual stack size and child-pruning statistics correlate with true
    remaining work better than raw UB proximity alone.
    """
    if max_pops <= 0:
        return {
            "pops": 0,
            "expanded": 0,
            "residual_stack": 0,
            "avg_margin": 0.0,
            "child_prune_rate": 1.0,
            "near_child_rate": 0.0,
            "best_gain": 0,
            "exhausted_budget": False,
        }

    best = incumbent
    init_ub = _adaptive_lookahead_ub(
        prefix, n, lookahead, global_max, best,
        refine_lookahead=refine_lookahead,
        refine_margin=refine_margin,
    )
    stack: list[tuple[tuple[int, ...], int]] = [(prefix, init_ub)]

    pops = 0
    expanded = 0
    popped_pruned = 0
    child_pruned = 0
    child_count = 0
    near_children = 0
    margin_sum = 0.0
    margin_count = 0

    while stack and pops < max_pops:
        node, ub = stack.pop()
        pops += 1
        if ub < best:
            popped_pruned += 1
            continue

        margin_sum += float(ub - best)
        margin_count += 1
        expanded += 1

        k = len(node)
        if k == n:
            if ub > best:
                best = ub
            continue

        children: list[tuple[int, tuple[int, ...]]] = []
        for bit in (0, 1):
            child = node + (bit,)
            child_ub = _adaptive_lookahead_ub(
                child, n, lookahead, global_max, best,
                refine_lookahead=refine_lookahead,
                refine_margin=refine_margin,
            )
            child_count += 1
            if child_ub <= best + near_margin:
                near_children += 1
            if child_ub >= best:
                children.append((child_ub, child))
            else:
                child_pruned += 1

        children.sort(key=lambda item: item[0])
        for child_ub, child in children:
            stack.append((child, child_ub))

    return {
        "pops": pops,
        "expanded": expanded,
        "residual_stack": len(stack),
        "avg_margin": margin_sum / margin_count if margin_count else 0.0,
        "child_prune_rate": child_pruned / child_count if child_count else 1.0,
        "near_child_rate": near_children / child_count if child_count else 0.0,
        "best_gain": best - incumbent,
        "exhausted_budget": bool(stack) and pops >= max_pops,
        "popped_pruned": popped_pruned,
    }


def _probe_is_hard(probe: dict[str, float | int | bool]) -> bool:
    """Classify a task as hard using probe behavior, not just UB gap."""
    return bool(probe["exhausted_budget"]) and int(probe["residual_stack"]) >= 8


def _probe_priority_score(
    ub: int,
    incumbent: int,
    probe: dict[str, float | int | bool],
) -> float:
    """Priority score for task ordering after probing.

    Higher scores are scheduled earlier so likely-hard or likely-fruitful tasks
    start sooner, reducing the chance that they become the only remaining tail.
    """
    score = float(ub - incumbent)
    if bool(probe["exhausted_budget"]):
        score += 128.0
    score += 6.0 * float(probe["residual_stack"])
    score += 64.0 * (1.0 - float(probe["child_prune_rate"]))
    score += 32.0 * float(probe["near_child_rate"])
    score += 2.0 * float(probe["best_gain"])
    return score


def branch_bound_a112509_parallel(
    n: int,
    known_sequence: Optional[dict] = None,
    incumbent_hint: int = 0,
    lookahead: int = 1,
    split_depth: int = 0,
    num_workers: int = 0,
    verbose: bool = True,
    use_mh_hint: bool = True,
    collect_all: bool = False,
    hard_split_gap: int = 5,
    hard_split_extra: int = 7,
    max_nodes_per_task: int = 0,
    dfs_lookahead: int = 0,
    dfs_refine_lookahead: int = 0,
    dfs_refine_margin: int = 2,
    probe_nodes: int = 0,
) -> tuple[int, list[str], dict]:
    """Parallel branch-and-bound-search for a(n) via BFS-split work distribution.

    Algorithm
    ---------
    Phase 1 (enumeration): BFS from (1,) to depth ``split_depth``, collecting
    every prefix whose UB >= incumbent.  This takes O(2^split_depth) evaluations
    but prunes aggressively (typically yields 100-10 000 prefixes).

    Phase 2 (search): Each surviving prefix is a sub-tree root sent to an
    independent worker process.  Workers run the full branch-and-bound algorithm on their
    sub-tree using the exact (known) incumbent as the lower bound.

    Because the incumbent is exact, sub-trees are truly independent -- no
    communication is needed during search.

    Parameters
    ----------
    n : int
        Bit-string length.
    known_sequence : dict, optional
        {n: a(n)}.  Must contain key ``n`` for the search to be useful.
    incumbent_hint : int
        Extra floor on the incumbent (see branch_bound_a112509).
    lookahead : int
        Lookahead depth (passed to each worker).
    split_depth : int
        BFS expansion depth before distributing.  Should be chosen so that
        the number of surviving prefixes is >> num_workers.  Default 15 works
        well for n = 50-100.
    num_workers : int
        Number of parallel worker processes.  0 -> ``os.cpu_count()`` - 1.
    verbose : bool
        Progress reporting.
    use_mh_hint : bool
        When True (default), automatically load the best MH score for *n*
        from ``results/mh_bounded/`` and ``results/mh_unbounded/`` and use
        it as an additional incumbent lower bound.  This is always safe
        because an MH score is the score of a real string (<= a(n)).
        Set to False to disable auto-loading.
    collect_all : bool
        Passed to each worker.  When True, disables Rule 2 so every optimal
        string is collected (slower but exhaustive).  Default False.
    hard_split_gap : int
        Subtrees with UB <= incumbent + hard_split_gap are considered "hard"
        (near-optimal, expensive to prove).  They are re-split ``hard_split_extra``
        extra BFS levels to improve load balance.  Set to -1 to disable.
    hard_split_extra : int
        Number of extra BFS expansion levels applied to hard subtrees.
        Each level roughly doubles the number of sub-tasks for that region
        while halving their individual depth.  Default 7.
    max_nodes_per_task : int
        Work budget per worker task. For collect_all=False this bounds the
        best-first worker before returning a resumable heap snapshot. For
        collect_all=True, a positive value enables cooperative DFS chunking;
        when the budget is hit, the worker returns its pending sibling stack
        so those subtrees can be redistributed. 0 keeps the baseline behavior
        unless probe-driven chunking is enabled; a negative value disables
        DFS chunking even when probing is on.
    dfs_lookahead : int
        Lookahead depth used by DFS combined workers when collect_all=True.
        0 (default) = same as ``lookahead``.  Setting this to ``lookahead+1``
        (e.g. 3 when lookahead=2) tightens the per-node UB in DFS workers,
        pruning more aggressively at 2x cost per node — net win for hard
        subtasks where the L=2 bound is nearly tight.
    dfs_refine_lookahead : int
        Optional adaptive refinement depth for DFS combined workers.
        The base DFS UB is computed at ``dfs_lookahead`` (or ``lookahead`` when
        that is 0). If the base UB falls within ``dfs_refine_margin`` of the
        current best, it is recomputed at this tighter lookahead. 0 (default)
        auto-selects ``base_lookahead + 1`` when the base lookahead is at least 2;
        a negative value disables adaptive refinement.
    dfs_refine_margin : int
        Margin around the current best that triggers adaptive DFS refinement.
        Default 2; negative disables the near-threshold refinement pass.
    probe_nodes : int
        Number of DFS pops used for a shallow per-task hardness probe before
        scheduling collect_all subtrees. 0 disables probing and falls back to
        UB-gap-based splitting only. This is intended as an opt-in tuning knob
        for larger collect_all runs. Default 0.
    Returns
    -------
    (value, optimal_strings, stats)
    """
    import os

    t0 = time.perf_counter()
    if num_workers <= 0:
        num_workers = max(1, (os.cpu_count() or 2) - 1)
    # Auto-tune split_depth: target ~50x more subtrees than workers.
    # Rule of thumb: split_depth = n//3 clamped to [10, 25].
    if split_depth <= 0:
        split_depth = max(10, min(n // 3, 25))

    # ── Incumbent ─────────────────────────────────────────────────────────
    incumbent = _resolve_incumbent(
        n, known_sequence, incumbent_hint,
        use_mh_hint=use_mh_hint, verbose=verbose,
        label="parallel-branch-and-bound",
    )

    # ── Global UB cap ─────────────────────────────────────────────────────
    global_max = _compute_global_max(n)

    # ── Phase 1: enumerate split nodes ────────────────────────────────────
    t_enum = time.perf_counter()
    split_pairs: list[tuple[tuple[int, ...], int]] = _enum_split_nodes(
        n, split_depth, incumbent, known_sequence or {}, lookahead, global_max
    )

    _dfs_L = dfs_lookahead if dfs_lookahead > 0 else lookahead
    if dfs_refine_lookahead == 0:
        _dfs_refine_L = _dfs_L + 1 if _dfs_L >= 2 else -1
    elif dfs_refine_lookahead < 0:
        _dfs_refine_L = -1
    else:
        _dfs_refine_L = dfs_refine_lookahead

    probe_pairs: list[tuple[tuple[int, ...], int, float, bool]] | None = None
    probe_hard_prefixes: set[tuple[int, ...]] | None = None
    if collect_all and probe_nodes > 0 and split_pairs and n >= 40:
        t_probe = time.perf_counter()
        probe_pairs = []
        for prefix, ub in split_pairs:
            probe = _probe_subtree_hardness(
                prefix, n, incumbent, _dfs_L, global_max,
                refine_lookahead=_dfs_refine_L,
                refine_margin=dfs_refine_margin,
                max_pops=probe_nodes,
                near_margin=max(2, dfs_refine_margin),
            )
            probe_pairs.append(
                (prefix, ub, _probe_priority_score(ub, incumbent, probe), _probe_is_hard(probe))
            )
        if verbose:
            hard_count = sum(1 for _, _, _, is_hard in probe_pairs if is_hard)
            print(f"[parallel-branch-and-bound] n={n}: probe-ranked {len(probe_pairs)} subtrees in "
                  f"{time.perf_counter() - t_probe:.2f}s ({hard_count} predicted hard)")

    # Adaptive re-split: expand hard subtrees extra BFS levels so no single
    # task dominates wall-clock time. When probing is enabled, hardness comes
    # from probe behavior rather than UB proximity alone.
    if hard_split_extra > 0 and split_pairs:
        if probe_pairs is not None:
            easy_pairs = [(p, u, score) for p, u, score, is_hard in probe_pairs if not is_hard]
            hard_rows = [(p, score) for p, _, score, is_hard in probe_pairs if is_hard]
            # Estimate per-prefix remaining work from probe score and split until
            # each predicted task is near target size. This avoids over-splitting
            # easy hard-ish nodes while aggressively slicing true stragglers.
            target_cost = 180.0
            expanded_hard: list[tuple[tuple[int, ...], int, float]] = []
            for prefix, score in hard_rows:
                est_cost = max(1.0, score)
                extra_depth = 0
                while est_cost > target_cost and extra_depth < hard_split_extra:
                    est_cost *= 0.5
                    extra_depth += 1
                if extra_depth == 0:
                    extra_depth = 1

                hard_frontier = [prefix]
                for _ in range(extra_depth):
                    new_hard: list[tuple[int, ...]] = []
                    for hard_prefix in hard_frontier:
                        if len(hard_prefix) >= n - 1:
                            new_hard.append(hard_prefix)
                            continue
                        for bit in (0, 1):
                            child = hard_prefix + (bit,)
                            child_ub = _lookahead_ub(child, n, lookahead, global_max)
                            if child_ub >= incumbent:
                                new_hard.append(child)
                    hard_frontier = new_hard
                expanded_hard.extend((p, incumbent, score) for p in hard_frontier)
            probe_hard_prefixes = {p for p, _, _ in expanded_hard}
            split_ranked = expanded_hard + easy_pairs
            split_ranked.sort(key=lambda item: item[2], reverse=True)
            split_pairs = [(p, u) for p, u, _ in split_ranked]
        elif hard_split_gap >= 0:
            easy_pairs = [(p, u) for p, u in split_pairs
                          if u > incumbent + hard_split_gap]
            hard_prefixes = [p for p, u in split_pairs
                             if u <= incumbent + hard_split_gap]
            # BFS-expand hard prefixes hard_split_extra more levels.
            for _ in range(hard_split_extra):
                new_hard: list[tuple[int, ...]] = []
                for prefix in hard_prefixes:
                    if len(prefix) >= n - 1:
                        new_hard.append(prefix)   # already near leaf -- keep as-is
                        continue
                    for bit in (0, 1):
                        child = prefix + (bit,)
                        child_ub = _lookahead_ub(child, n, lookahead, global_max)
                        if child_ub >= incumbent:
                            new_hard.append(child)
                hard_prefixes = new_hard
            split_pairs = easy_pairs + [(p, incumbent) for p in hard_prefixes]

    t_enum = time.perf_counter() - t_enum

    # Run higher-UB subtrees first so a stronger incumbent is found earlier.
    # This improves global pruning for both local multiprocessing and sharded runs.
    split_pairs.sort(key=lambda item: item[1], reverse=True)

    if verbose:
        print(f"[parallel-branch-and-bound] n={n}: {len(split_pairs)} subtrees found at "
              f"depth {min(split_depth, n-1)} in {t_enum:.2f}s "
              f"({num_workers} workers)")

    if not split_pairs:
        # No candidates at all -- incumbent is proven tight with zero expansions.
        elapsed = time.perf_counter() - t0
        return incumbent, [], {"elapsed_seconds": round(elapsed, 4),
                               "nodes_expanded": 0, "nodes_pruned": 0,
                               "num_workers": num_workers,
                               "split_depth": split_depth}

    # ── Phase 2: parallel search ─────────────────────────────────────────
    # collect_all=True  → DFS combined workers: O(n) stack, no heap, no OOM.
    #                     Finds a(n) AND enumerates ALL optimal strings in one pass.
    # collect_all=False → branch-and-bound workers: best-first, standard value-finding path.
    total_nodes = 0
    total_pruned = 0
    worker_results: list[tuple[int, list[str]]] = []
    best_found = incumbent
    done = 0
    submitted = len(split_pairs)
    t_last = t0

    def _print_progress():
        nonlocal t_last
        now = time.perf_counter()
        if (done == 1 or done == submitted or now - t_last >= 30.0
                or (done > 0 and done % max(1, submitted // 20) == 0)):
            pct = 100 * done / submitted
            elapsed_so_far = now - t0
            print(
                f"  [{done}/{submitted}] "
                f"{pct:5.1f}%  nodes={total_nodes:,}  "
                f"elapsed={elapsed_so_far:.0f}s",
                flush=True,
            )
            t_last = now

    if collect_all:
        # ── DFS combined pool with iterative tail re-split ───────────────
        # Uses apply_async + external task queue instead of imap_unordered.
        # When (queue + active) drops below a threshold, every QUEUED task
        # is BFS-expanded `_tail_extra` more levels. This repeats every
        # time the queue depletes, until no further sub-division is
        # possible (each prefix already at depth >= n - 1) or we hit
        # `_max_resplits` rounds. Each round multiplies queued tasks by
        # ~2^_tail_extra, exponentially refining the tail to keep all
        # workers busy through the slow finish.
        #
        # NB: re-split only helps QUEUED tasks. In-flight tasks already
        # running on a worker cannot be preempted; if a few tasks are
        # individually slow (no further sub-division possible), they
        # will still serialise the final stretch. Cooperative checkpoint
        # in the JIT'd worker would be the next-level fix.
        # Trigger re-split at ~10% remaining, clamped to [3W, 30W] tasks.
        _tail_thresh = max(num_workers * 3,
                           min(len(split_pairs) // 10, num_workers * 30))
        _tail_extra = 3        # BFS levels added per re-split round (~8x finer)
        _max_resplits = 5      # cap: 8^5 = 32k-fold refinement of the tail
        _resplits_done = 0     # round counter
        if max_nodes_per_task < 0:
            _dfs_chunk_nodes = 0
        elif max_nodes_per_task > 0:
            _dfs_chunk_nodes = max_nodes_per_task
        else:
            # Always-on cooperative chunking for collect_all runs.
            # Default budget ~5s/chunk at Numba speed (~10-20 M nodes/s).
            # Small tasks complete in one chunk; only large tasks re-enqueue.
            _dfs_chunk_nodes = max(5_000_000, 50_000 * n)

        def _adaptive_tail_chunk(remaining_tasks: int, resplit_round: int) -> int:
            """Return chunk budget tuned for fairness as the queue thins."""
            if _dfs_chunk_nodes <= 0:
                return 0
            if remaining_tasks <= num_workers:
                return max(200_000, _dfs_chunk_nodes // 16)
            if remaining_tasks <= num_workers * 4:
                return max(500_000, _dfs_chunk_nodes // 8)
            if resplit_round > 0:
                return max(1_000_000, _dfs_chunk_nodes // 4)
            return _dfs_chunk_nodes

        if verbose:
            print(f"[parallel-branch-and-bound] n={n}: DFS mode (collect_all) - "
                  f"{len(split_pairs)} subtrees, {num_workers} workers "
                  f"(tail re-split at <={_tail_thresh} remaining, "
                  f"refine L={_dfs_refine_L if _dfs_refine_L > _dfs_L else 'off'}, "
                  f"chunk={_dfs_chunk_nodes if _dfs_chunk_nodes > 0 else 'off'})")

        task_queue: deque = deque(
            (
                sp,
                n,
                incumbent,
                _dfs_L,
                _dfs_refine_L,
                dfs_refine_margin,
                _dfs_chunk_nodes,
                None,
            )
            for sp, _ in split_pairs
        )
        ctx = multiprocessing.get_context("spawn")
        shared_incumbent = ctx.Value("l", incumbent)
        with ctx.Pool(processes=num_workers,
                      initializer=_init_shared_worker,
                      initargs=(shared_incumbent,)) as pool:
            active: list = []

            def _refill(n_submit: int) -> None:
                while task_queue and n_submit > 0:
                    t = list(task_queue.popleft())
                    if t[6] > 0:
                        rem_est = len(task_queue) + len(active) + 1
                        t[6] = _adaptive_tail_chunk(rem_est, _resplits_done)
                    active.append(pool.apply_async(_dfs_combined_worker, (t,)))
                    n_submit -= 1

            _refill(num_workers * 2)   # initial fill: one task queued per worker

            while active or task_queue:
                still_active = [f for f in active if not f.ready()]
                newly_done   = [f for f in active if f.ready()]

                for fut in newly_done:
                    best_val, opt_strs, dfs_stats = fut.get()
                    total_nodes  += dfs_stats["nodes_expanded"]
                    total_pruned += dfs_stats["nodes_pruned"]
                    if best_val > best_found:
                        best_found = best_val
                    worker_results.append((best_val, opt_strs))
                    done += 1
                    if dfs_stats.get("interrupted"):
                        cur_inc = max(best_found, shared_incumbent.value)
                        rem_est = len(task_queue) + len(active)
                        requeue_chunk = _adaptive_tail_chunk(rem_est, _resplits_done)
                        requeued: list[tuple] = []
                        for prefix, prefix_ub in dfs_stats.get("stack_snapshot", []):
                            if prefix_ub >= cur_inc:
                                requeued.append((
                                    prefix,
                                    n,
                                    cur_inc,
                                    _dfs_L,
                                    _dfs_refine_L,
                                    dfs_refine_margin,
                                    requeue_chunk,
                                    None,
                                ))
                        submitted += len(requeued)
                        for task in reversed(requeued):
                            task_queue.appendleft(task)
                    if verbose:
                        _print_progress()

                active[:] = still_active

                # Iterative tail re-split: every time (queue + active) drops
                # below the threshold AND there is still queued work, expand
                # each queued task `_tail_extra` BFS levels deeper. Skip
                # prefixes already at depth >= n - 1 (cannot expand further).
                # Stop after `_max_resplits` rounds.
                remaining = len(task_queue) + len(active)
                if (_resplits_done < _max_resplits
                        and remaining <= _tail_thresh
                        and task_queue):
                    _resplits_done += 1
                    to_expand = list(task_queue)
                    task_queue.clear()
                    expanded: list = []
                    n_unsplittable = 0
                    for sp, _n, _inc, _L, _refine_L, _margin, _chunk_nodes, _stack in to_expand:
                        if len(sp) >= _n - 1:
                            # Cannot expand further -- keep as-is.
                            _next_chunk = _adaptive_tail_chunk(
                                len(to_expand) + len(active), _resplits_done
                            )
                            expanded.append((sp, _n, _inc, _L, _refine_L,
                                             _margin, _next_chunk, None))
                            n_unsplittable += 1
                            continue
                        new_pfxs = _bfs_expand_prefix(
                            sp, _n, _inc, _L, _tail_extra, global_max)
                        _next_chunk = _adaptive_tail_chunk(
                            len(new_pfxs) + len(active), _resplits_done
                        )
                        expanded.extend(
                            (_p, _n, _inc, _L, _refine_L, _margin, _next_chunk, None)
                            for _p in new_pfxs
                        )
                    submitted += len(expanded) - len(to_expand)
                    task_queue.extend(expanded)
                    # Tighten the threshold for the next round so we don't
                    # re-split too eagerly when only a handful remain.
                    _tail_thresh = max(num_workers, _tail_thresh // 2)
                    if verbose:
                        print(
                            f"[parallel-branch-and-bound] n={n}: tail re-split #"
                            f"{_resplits_done} - {len(to_expand)} -> "
                            f"{len(expanded)} tasks "
                            f"(unsplittable={n_unsplittable})",
                            flush=True,
                        )
                    # If no progress was made (every prefix unsplittable),
                    # stop trying.
                    if len(expanded) == len(to_expand):
                        _resplits_done = _max_resplits

                _refill(max(0, num_workers * 2 - len(active)))

                if not newly_done:
                    if verbose and (time.perf_counter() - t_last) >= 30.0:
                        _print_progress()
                    time.sleep(0.05)

    else:
        # ── branch-and-bound pool (collect_all=False) ──────────────────────────────────
        tasks = [
            (sp, n, incumbent, known_sequence or {}, lookahead, False, max_nodes_per_task, None)
            for sp, _ in split_pairs
        ]
        ctx = multiprocessing.get_context("spawn")
        shared_incumbent = ctx.Value("l", incumbent)
        with ctx.Pool(processes=num_workers,
                      initializer=_init_shared_worker,
                      initargs=(shared_incumbent,)) as pool:
            pending = [pool.apply_async(_branch_bound_from_prefix, (task,)) for task in tasks]
            while pending:
                still_pending: list = []
                newly_ready = []
                for fut in pending:
                    if fut.ready():
                        newly_ready.append(fut)
                    else:
                        still_pending.append(fut)
                for fut in newly_ready:
                    val, optimals, stats = fut.get()
                    total_nodes += stats["nodes_expanded"]
                    total_pruned += stats.get("nodes_pruned", 0)
                    if val > best_found:
                        best_found = val
                    worker_results.append((val, optimals))
                    done += 1
                    if stats.get("interrupted"):
                        sv = shared_incumbent.value
                        cur_inc = max(best_found, sv)
                        for ub, prefix in stats.get("heap_snapshot", []):
                            if ub >= cur_inc:
                                new_task = (
                                    prefix, n, cur_inc,
                                    known_sequence or {}, lookahead,
                                    False, max_nodes_per_task, None,
                                )
                                submitted += 1
                                still_pending.append(
                                    pool.apply_async(_branch_bound_from_prefix, (new_task,))
                                )
                        heap_bulk = stats.get("heap_bulk", [])
                        if heap_bulk:
                            filtered_bulk = [(ub, pfx) for ub, pfx in heap_bulk
                                             if ub >= cur_inc]
                            if filtered_bulk:
                                bulk_task = (
                                    None, n, cur_inc,
                                    known_sequence or {}, lookahead,
                                    False, max_nodes_per_task, filtered_bulk,
                                )
                                submitted += 1
                                still_pending.append(
                                    pool.apply_async(_branch_bound_from_prefix, (bulk_task,))
                                )
                    if verbose:
                        _print_progress()
                pending = still_pending
                if not newly_ready:
                    if verbose and (time.perf_counter() - t_last) >= 30.0:
                        _print_progress()
                    time.sleep(0.05)

    all_optimals = [s for val, opts in worker_results if val == best_found for s in opts]

    elapsed = time.perf_counter() - t0
    result_stats = {
        "elapsed_seconds": round(elapsed, 4),
        "nodes_expanded": total_nodes,
        "nodes_pruned": total_pruned,
        "num_workers": num_workers,
        "split_depth": split_depth,
        "num_subtrees": submitted,
        "lookahead": lookahead,
    }
    if verbose:
        print(f"[parallel-branch-and-bound] n={n}: a(n)={best_found}, "
              f"{total_nodes:,} expanded, {total_pruned:,} pruned, "
              f"{elapsed:.3f}s ({num_workers} workers)")

    return best_found, all_optimals, result_stats


# ---------------------------------------------------------------------------
# Multi-machine sharding helpers
# ---------------------------------------------------------------------------

def enumerate_shard_tasks(
    n: int,
    known_sequence: Optional[dict] = None,
    incumbent_hint: int = 0,
    lookahead: int = 2,
    split_depth: int = 0,
    hard_split_gap: int = 5,
    hard_split_extra: int = 7,
    dfs_lookahead: int = 0,
    dfs_refine_lookahead: int = 0,
    dfs_refine_margin: int = 2,
    use_mh_hint: bool = True,
    verbose: bool = True,
) -> list[dict]:
    """Run Phase 1 (BFS enumeration) only and return task dicts for sharding.

    Returns a list of task dicts, each with keys:
      id, n, prefix, ub, incumbent, lookahead, refine_lookahead, refine_margin

    Each dict can be saved to a JSON file and claimed by a shard worker.
    """
    import os as _os

    if split_depth <= 0:
        split_depth = max(10, min(n // 3, 25))

    incumbent = _resolve_incumbent(
        n, known_sequence, incumbent_hint,
        use_mh_hint=use_mh_hint, verbose=verbose,
        label="shard-enum",
    )

    global_max = _compute_global_max(n)

    split_pairs: list[tuple[tuple[int, ...], int]] = _enum_split_nodes(
        n, split_depth, incumbent, known_sequence or {}, lookahead, global_max
    )

    _dfs_L = dfs_lookahead if dfs_lookahead > 0 else lookahead
    if dfs_refine_lookahead == 0:
        _dfs_refine_L = _dfs_L + 1 if _dfs_L >= 2 else -1
    elif dfs_refine_lookahead < 0:
        _dfs_refine_L = -1
    else:
        _dfs_refine_L = dfs_refine_lookahead

    if hard_split_extra > 0 and hard_split_gap >= 0 and split_pairs:
        easy_pairs = [(p, u) for p, u in split_pairs if u > incumbent + hard_split_gap]
        hard_prefixes = [p for p, u in split_pairs if u <= incumbent + hard_split_gap]
        for _ in range(hard_split_extra):
            new_hard: list[tuple[int, ...]] = []
            for prefix in hard_prefixes:
                if len(prefix) >= n - 1:
                    new_hard.append(prefix)
                    continue
                for bit in (0, 1):
                    child = prefix + (bit,)
                    if _lookahead_ub(child, n, lookahead, global_max) >= incumbent:
                        new_hard.append(child)
            hard_prefixes = new_hard
        split_pairs = easy_pairs + [(p, incumbent) for p in hard_prefixes]

    if verbose:
        print(f"[shard-enum] n={n}: {len(split_pairs)} subtrees, incumbent={incumbent}")

    tasks = []
    for i, (prefix, ub) in enumerate(split_pairs):
        tasks.append({
            "id": i,
            "n": n,
            "prefix": list(prefix),
            "ub": ub,
            "incumbent": incumbent,
            "lookahead": _dfs_L,
            "refine_lookahead": _dfs_refine_L,
            "refine_margin": dfs_refine_margin,
        })
    return tasks


def _shard_claim_task(task_path: str, claim_dir: str) -> bool:
    """Attempt to exclusively claim a task file by creating a .pid sentinel.

    Returns True if claim succeeded (this worker owns the task).
    Uses O_CREAT|O_EXCL for atomic creation — safe on NTFS and POSIX.
    """
    import socket as _socket

    task_id = os.path.splitext(os.path.basename(task_path))[0]
    claim_path = os.path.join(claim_dir, task_id + ".pid")
    try:
        fd = os.open(claim_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        with os.fdopen(fd, "w") as f:
            f.write(f"{_socket.gethostname()}:{os.getpid()}\n")
        return True
    except FileExistsError:
        return False
    except OSError:
        return False


def _read_shard_incumbent(incumbent_path: str, fallback: int) -> int:
    """Read shard-wide incumbent (if present), else return fallback."""
    try:
        with open(incumbent_path, encoding="utf-8") as f:
            data = json.load(f)
        val = int(data.get("incumbent", fallback))
        return val if val >= fallback else fallback
    except Exception:
        return fallback


def _ratchet_shard_incumbent(incumbent_path: str, candidate: int) -> int:
    """Monotonically raise shard-wide incumbent via lockfile + atomic replace."""
    lock_path = incumbent_path + ".lock"
    tmp_path = incumbent_path + f".{os.getpid()}.tmp"
    locked = False
    try:
        while not locked:
            try:
                fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.close(fd)
                locked = True
            except FileExistsError:
                time.sleep(0.01)

        cur = _read_shard_incumbent(incumbent_path, 0)
        new_val = candidate if candidate > cur else cur
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump({"incumbent": int(new_val), "pid": os.getpid()}, f)
        os.replace(tmp_path, incumbent_path)
        return int(new_val)
    finally:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        if locked:
            try:
                os.remove(lock_path)
            except OSError:
                pass


def run_shard_worker_loop(
    shard_dir: str,
    n: int,
    chunk_nodes: int = 0,
    poll_interval: float = 1.0,
    verbose: bool = True,
) -> int:
    """Claim and process tasks from a shard work directory.

    Directory layout expected under ``shard_dir/n_{n:04d}/``:
      tasks/task_{id:06d}.json   — task descriptors written by coordinator
      claims/task_{id:06d}.pid   — ownership sentinels (atomic O_CREAT|O_EXCL)
      done/task_{id:06d}.json    — result files written by workers

    Returns the number of tasks processed by this worker.
    """
    import glob as _glob

    tasks_dir  = os.path.join(shard_dir, f"n_{n:04d}", "tasks")
    claims_dir = os.path.join(shard_dir, f"n_{n:04d}", "claims")
    done_dir   = os.path.join(shard_dir, f"n_{n:04d}", "done")
    incumbent_path = os.path.join(shard_dir, f"n_{n:04d}", "incumbent.json")
    for d in (claims_dir, done_dir):
        os.makedirs(d, exist_ok=True)

    processed = 0
    while True:
        pending = sorted(_glob.glob(os.path.join(tasks_dir, "task_*.json")))
        if not pending:
            if verbose:
                print(f"[shard-worker] n={n}: no task files found in {tasks_dir}")
            break

        claimed_any = False
        for task_path in pending:
            task_id = os.path.splitext(os.path.basename(task_path))[0]
            done_path = os.path.join(done_dir, task_id + ".json")
            if os.path.exists(done_path):
                continue  # already completed by another worker
            if not _shard_claim_task(task_path, claims_dir):
                continue  # claimed by another worker

            # Load task.
            with open(task_path, encoding="utf-8") as f:
                task = json.load(f)

            prefix   = tuple(task["prefix"])
            base_inc = int(task["incumbent"])
            live_inc = _read_shard_incumbent(incumbent_path, base_inc)
            inc      = live_inc if live_inc > base_inc else base_inc
            L        = task["lookahead"]
            refine_L = task["refine_lookahead"]
            margin   = task["refine_margin"]
            mn       = task.get("n", n)
            c_nodes  = chunk_nodes if chunk_nodes > 0 else max(5_000_000, 50_000 * mn)

            if verbose:
                print(f"[shard-worker] claiming {task_id} (prefix len={len(prefix)}, ub={task['ub']})",
                      flush=True)

            t0 = time.perf_counter()
            best_val, opt_strs, stats = _dfs_combined_worker(
                (prefix, mn, inc, L, refine_L, margin, c_nodes, None)
            )
            elapsed = time.perf_counter() - t0

            # Publish improvements promptly for other workers.
            live_after = _ratchet_shard_incumbent(incumbent_path, int(best_val))

            result = {
                "task_id": task_id,
                "n": mn,
                "prefix": list(prefix),
                "best": best_val,
                "num_optimal": len(opt_strs),
                "optimal_strings": opt_strs,
                "stats": {k: v for k, v in stats.items() if k != "stack_snapshot"},
                "elapsed_s": round(elapsed, 3),
            }
            with open(done_path, "w", encoding="utf-8") as f:
                json.dump(result, f)

            if verbose:
                print(f"[shard-worker] done  {task_id}: best={best_val}, live_inc={live_after}, "
                      f"nodes={stats['nodes_expanded']:,}, {elapsed:.1f}s", flush=True)

            processed += 1
            claimed_any = True

        if not claimed_any:
            # All tasks claimed; check if done dir has everything.
            n_done = len(_glob.glob(os.path.join(done_dir, "task_*.json")))
            if n_done >= len(pending):
                break
            if verbose:
                print(f"[shard-worker] waiting for other workers ({n_done}/{len(pending)} done) ...",
                      flush=True)
            time.sleep(poll_interval)

    return processed


def merge_shard_results(
    shard_dir: str,
    n: int,
    verbose: bool = True,
) -> tuple[int, list[str], dict]:
    """Merge all done-task result files into a final (best, optimals, stats).

    Returns (best_value, all_optimal_strings, aggregate_stats).
    """
    import glob as _glob

    done_dir = os.path.join(shard_dir, f"n_{n:04d}", "done")
    done_files = sorted(_glob.glob(os.path.join(done_dir, "task_*.json")))
    if not done_files:
        raise FileNotFoundError(f"No done files found in {done_dir}")

    best_found = 0
    all_optimals: list[str] = []
    total_nodes = 0
    total_pruned = 0

    for path in done_files:
        with open(path, encoding="utf-8") as f:
            r = json.load(f)
        bv = r["best"]
        if bv > best_found:
            best_found = bv
            all_optimals = list(r["optimal_strings"])
        elif bv == best_found:
            all_optimals.extend(r["optimal_strings"])
        total_nodes  += r["stats"].get("nodes_expanded", 0)
        total_pruned += r["stats"].get("nodes_pruned", 0)

    # Deduplicate.
    all_optimals = list(dict.fromkeys(all_optimals))

    if verbose:
        print(f"[shard-merge] n={n}: a(n)={best_found}, "
              f"{len(all_optimals)} optimal strings, "
              f"{total_nodes:,} nodes from {len(done_files)} tasks")

    return best_found, all_optimals, {
        "nodes_expanded": total_nodes,
        "nodes_pruned": total_pruned,
        "num_tasks": len(done_files),
    }
