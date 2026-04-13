"""
populate_large_n_solutions.py
==============================
For each large-n result file that has a seed file but no stored solution string,
this script:
  1. Loads the seed (npy or txt)
  2. Re-scores it exactly using the suffix-array method
  3. Computes the run-length encoding (RLE) of the seed
  4. Updates the result JSON with:
       - solution_rle:  comma-separated run lengths (reconstruct with np.repeat)
       - best_value / score:  the exact verified score of the seed
       - method:  updated to reflect this script
  5. Writes the RLE to  seeds/<tag>.rle.txt  for GitHub storage

Usage:
    python -m src.tools.populate_large_n_solutions          # all entries
    python -m src.tools.populate_large_n_solutions --n 10000000
    python -m src.tools.populate_large_n_solutions --dry-run

Timing (rough):
    n=10M   ~30s
    n=100M  ~5 min
    n=1B    ~2 h
    n=2B    ~4 h
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

# ── repo root ──────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(ROOT))
from src.algorithms.surgical_nudge import build_sa_lcp, score_from_sa_lcp

# ── entries to process ─────────────────────────────────────────────────────
# (result_json, seed_file)
ENTRIES = [
    (
        ROOT / "results/large_n/n_10000000_results.json",
        ROOT / "seeds/n_10000000_seed.npy",
    ),
    (
        ROOT / "results/large_n/n_100000000_results.json",
        ROOT / "seeds/n_100000000_seed.npy",
    ),
    (
        ROOT / "results/large_n/n_1000000000_results.json",
        ROOT / "seeds/seed_n1000000000.npy",
    ),
    (
        ROOT / "results/large_n/n_2000000000_results.json",
        ROOT / "seeds/seed_n2000000000.npy",
    ),
]


def bits_to_rle(bits: np.ndarray) -> list[int]:
    """Return run-length encoding as a list of integers (alternating 1-runs, 0-runs)."""
    runs = []
    cur = int(bits[0])
    cnt = 1
    for b in bits[1:]:
        if int(b) == cur:
            cnt += 1
        else:
            runs.append(cnt)
            cur = int(b)
            cnt = 1
    runs.append(cnt)
    return runs


def rle_to_bits(runs: list[int], n: int) -> np.ndarray:
    """Reconstruct bit array from RLE produced by bits_to_rle().
    Starts with a 1-run. Asserts total length == n."""
    vals = np.zeros(len(runs), dtype=np.int8)
    vals[::2] = 1  # even indices are 1-runs
    bits = np.repeat(vals, runs)
    assert len(bits) == n, f"RLE length {len(bits)} != n={n}"
    return bits.astype(np.int32)


def load_seed(seed_path: Path) -> np.ndarray:
    """Load a seed as an int32 numpy array regardless of file format."""
    if seed_path.suffix == ".npy":
        return np.load(seed_path).astype(np.int32)
    # Plain text: one character per bit
    text = seed_path.read_text().strip()
    return np.frombuffer(text.encode(), dtype=np.uint8).astype(np.int32) - ord("0")


def process_entry(result_path: Path, seed_path: Path, dry_run: bool = False) -> None:
    if not seed_path.exists():
        print(f"  SKIP: seed not found: {seed_path}")
        return

    with open(result_path) as f:
        result = json.load(f)

    n = result.get("n", 0)
    if not n:
        print(f"  SKIP: no n in {result_path.name}")
        return

    stored_score = result.get("best_value") or result.get("score") or 0
    already_has_rle = "solution_rle" in result

    print(f"\n{'='*60}")
    print(f"  n = {n:,}")
    print(f"  result: {result_path.name}")
    print(f"  seed:   {seed_path.name}")
    print(f"  stored score: {stored_score:,}")
    print(f"  has RLE already: {already_has_rle}")

    print("  Loading seed...", end=" ", flush=True)
    t0 = time.time()
    bits = load_seed(seed_path)
    print(f"{time.time()-t0:.1f}s  ({len(bits):,} bits)")

    if len(bits) != n:
        print(f"  WARNING: seed length {len(bits):,} != n={n:,}. Skipping.")
        return

    # Compute RLE
    print("  Computing RLE...", end=" ", flush=True)
    t0 = time.time()
    runs = bits_to_rle(bits)
    rle_str = ",".join(map(str, runs))
    print(f"{time.time()-t0:.2f}s  ({len(runs):,} runs, {len(rle_str):,} bytes)")

    # Re-score the seed
    print("  Scoring seed (SA+LCP)...", end=" ", flush=True)
    t0 = time.time()
    sa, lcp = build_sa_lcp(bits)
    seed_score = score_from_sa_lcp(bits, sa, lcp)
    elapsed = time.time() - t0
    print(f"{elapsed:.1f}s  =>  {seed_score:,}")

    ratio = seed_score / n**2
    print(f"  ratio a(n)/n^2 = {ratio:.10f}")

    # Compare
    if seed_score == stored_score:
        print("  ✓ Seed score matches stored value exactly.")
    elif seed_score > stored_score:
        print(f"  ↑ Seed score EXCEEDS stored value by {seed_score - stored_score:,}")
        print(f"    (stored was from: {result.get('method', 'unknown')})")
    else:
        print(f"  ↓ Seed score is LOWER than stored value by {stored_score - seed_score:,}")
        print("    Keeping stored value; not updating score.")

    if dry_run:
        print("  [dry-run] No files written.")
        return

    # Write RLE file to seeds/
    rle_path = seed_path.parent / (seed_path.stem + ".rle.txt")
    print(f"  Writing RLE to {rle_path.name}...", end=" ", flush=True)
    rle_path.write_text(rle_str)
    print("done")

    # Update result JSON
    # Only update score if seed scored higher (never reduce a stored lower bound)
    new_score = max(seed_score, stored_score)
    result["solution_rle"] = rle_str
    result["solution_rle_path"] = str(rle_path.relative_to(ROOT))
    result["solution_rle_runs"] = len(runs)
    if seed_score >= stored_score:
        # Update whichever field held the score
        if "best_value" in result:
            result["best_value"] = new_score
        else:
            result["score"] = new_score
        result["a_over_n_sq"] = ratio
        result["method"] = "seed (re-scored and verified)"
    result["seed_score_verified"] = seed_score
    result["seed_score_matches_stored"] = (seed_score == stored_score)
    result["verification_timestamp"] = datetime.now().isoformat()

    # Atomic write
    tmp = result_path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(result, f, indent=2)
    tmp.replace(result_path)
    print(f"  Updated {result_path.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=None, help="Process only this n value")
    parser.add_argument("--dry-run", action="store_true", help="Score and report but don't write files")
    args = parser.parse_args()

    entries = ENTRIES
    if args.n is not None:
        entries = [(r, s) for r, s in ENTRIES if str(args.n) in r.name]
        if not entries:
            print(f"No entry found for n={args.n}")
            sys.exit(1)

    for result_path, seed_path in entries:
        try:
            process_entry(result_path, seed_path, dry_run=args.dry_run)
        except Exception as e:
            print(f"\nERROR processing {result_path.name}: {e}")
            import traceback
            traceback.print_exc()

    print("\nDone.")


if __name__ == "__main__":
    main()
