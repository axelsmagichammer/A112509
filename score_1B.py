"""Compute a(1,000,000,000) lower bound: build a counter-pattern seed and
score it exactly via suffix array + LCP array.

Memory: ~19 GB peak (fits in 27 GB available)
Time:   ~20-40 minutes total
"""
import numpy as np
import time
import sys
import gc
import json
from pathlib import Path
from datetime import datetime

sys.path.insert(0, "src/algorithms")
from surgical_nudge import (
    _build_from_pattern,
    build_sa_lcp,
    score_from_sa_lcp,
    _DIVSUFSORT_AVAILABLE,
)

assert _DIVSUFSORT_AVAILABLE, "pydivsufsort is required for n=1B"

N = 2_000_000_000

print(f"{'='*60}")
print(f"  Computing a({N:,}) lower bound")
print(f"  (seed-only, no optimization)")
print(f"{'='*60}")
print()

# Step 1: Generate seed
print("[1/3] Generating counter-pattern seed...")
t0 = time.perf_counter()
bits = _build_from_pattern(N)
t1 = time.perf_counter()
print(f"  Done in {t1-t0:.1f}s  ({bits.nbytes/1e9:.1f} GB, dtype={bits.dtype})")
print()

# Step 2: Build suffix array + LCP
print("[2/3] Building suffix array + LCP (this is the slow part)...")
print("  SA construction (divsufsort, O(n))...")
t_sa_start = time.perf_counter()
sa, lcp = build_sa_lcp(bits)
t_sa_end = time.perf_counter()
gc.collect()
print(f"  Done in {t_sa_end-t_sa_start:.1f}s")
print(f"  sa: {sa.nbytes/1e9:.1f} GB, lcp: {lcp.nbytes/1e9:.1f} GB")
print()

# Step 3: Score
print("[3/3] Scoring...")
t0 = time.perf_counter()
score = score_from_sa_lcp(bits, sa, lcp)
t1 = time.perf_counter()
print(f"  Done in {t1-t0:.1f}s")
print()

total_time = t1 - t_sa_start + (t_sa_start - t0)  # approximate
n_sq = N * N

print(f"{'='*60}")
print(f"  n          = {N:,}")
print(f"  a(n)      >= {score:,}")
print(f"  a(n) / n²  = {score / n_sq:.12f}")
print(f"  max LCP    = {int(lcp.max()):,}")
print(f"{'='*60}")

# Save results
results_dir = Path("results/large_n")
results_dir.mkdir(parents=True, exist_ok=True)
seeds_dir = Path("seeds")
seeds_dir.mkdir(parents=True, exist_ok=True)

# Save seed as .npy binary (4 GB, much faster than text)
seed_path = seeds_dir / f"seed_n{N}.npy"
np.save(seed_path, bits)
print(f"\nSeed saved: {seed_path} ({seed_path.stat().st_size / 1e9:.1f} GB)")

# Save result JSON (no bitstring — way too large for JSON)
result = {
    "n": N,
    "score": score,
    "a_over_n_sq": score / n_sq,
    "max_lcp": int(lcp.max()),
    "ones": int(np.sum(bits == 1)),
    "zeros": int(np.sum(bits == 0)),
    "density": float(np.sum(bits == 1)) / N,
    "method": "counter-pattern seed (no optimization)",
    "timestamp": datetime.now().isoformat(),
    "sa_build_seconds": t_sa_end - t_sa_start,
}
result_path = results_dir / f"n_{N:010d}_results.json"
result_path.write_text(json.dumps(result, indent=2))
print(f"Result saved: {result_path}")
