"""
export_large_n_rle.py
=====================
For each large-n seed (.npy file), this script:

  1. Loads the seed and computes its run-length encoding (RLE).
  2. Writes an .rle.txt file (single comma-separated line of run lengths)
     alongside the .npy file in seeds/. These small text files are
     committed to GitHub while the bulky .npy files are gitignored.
  3. Scores the seed using the suffix-array exact evaluator and compares
     it against the stored value in the corresponding result JSON.
  4. Updates the result JSON to add a `solution_rle` field (the RLE
     string) and, if the seed score is higher, updates the score field.

Usage:
    python -m src.tools.export_large_n_rle          # process all seeds
    python -m src.tools.export_large_n_rle 10000000 # process one n value

Reconstruct a seed from its .rle.txt:
    import numpy as np
    runs = list(map(int, open('seeds/n_10000000_seed.rle.txt').read().split(',')))
    bits = np.repeat(np.tile([1, 0], (len(runs) + 1) // 2)[:len(runs)], runs).astype(np.int32)
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Locate repo root regardless of cwd
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent.parent
SEEDS_DIR = ROOT / "seeds"
RESULTS_DIR = ROOT / "results" / "large_n"

# ---------------------------------------------------------------------------
# Mapping: n → (npy filename, result json filename)
# ---------------------------------------------------------------------------
TARGETS = {
    10_000_000:        ("n_10000000_seed.npy",      "n_10000000_results.json"),
    100_000_000:       ("n_100000000_seed.npy",     "n_100000000_results.json"),
    1_000_000_000:     ("seed_n1000000000.npy",     "n_1000000000_results.json"),
    2_000_000_000:     ("seed_n2000000000.npy",     "n_2000000000_results.json"),
}


def compute_rle(bits: np.ndarray) -> list:
    """Return run-length encoding of a binary array as a list of integers.
    The first element is always the length of the leading 1-block.
    Alternates 1-runs / 0-runs / 1-runs / ...
    Uses numpy vectorisation — fast even for n = 2×10⁹.
    """
    bits8 = bits.astype(np.uint8)
    # Find positions where the value changes
    changes = np.flatnonzero(np.diff(bits8)) + 1   # indices of run starts (excl. 0)
    boundaries = np.concatenate(([0], changes, [len(bits8)]))
    runs = np.diff(boundaries).tolist()
    return runs


def rle_to_bits(runs: list) -> np.ndarray:
    """Reconstruct a binary int32 array from a run-length list."""
    vals = np.zeros(len(runs), dtype=np.int32)
    vals[::2] = 1          # even indices are 1-runs, odd are 0-runs
    return np.repeat(vals, runs)


def score_bits(bits: np.ndarray) -> int:
    """Compute a(n) exactly via suffix array + Kasai LCP."""
    try:
        from src.algorithms.surgical_nudge import build_sa_lcp, score_from_sa_lcp
    except ImportError:
        sys.path.insert(0, str(ROOT))
        from src.algorithms.surgical_nudge import build_sa_lcp, score_from_sa_lcp
    sa, lcp = build_sa_lcp(bits)
    return int(score_from_sa_lcp(bits, sa, lcp))


def get_stored_score(result: dict) -> int:
    """Return the stored score from a result dict regardless of key name."""
    return result.get("score") or result.get("best_value") or 0


def process(n: int, npy_name: str, json_name: str, skip_score: bool = False) -> None:
    npy_path  = SEEDS_DIR / npy_name
    rle_path  = SEEDS_DIR / (npy_path.stem + ".rle.txt")
    json_path = RESULTS_DIR / json_name

    if not npy_path.exists():
        print(f"  SKIP  seed not found: {npy_path}")
        return

    # ------------------------------------------------------------------
    # 1. Load seed
    # ------------------------------------------------------------------
    print(f"\nn = {n:,}")
    print(f"  Loading {npy_path.name} ...", flush=True)
    t0 = time.time()
    bits = np.load(str(npy_path)).astype(np.int32)
    assert len(bits) == n, f"Expected {n} bits, got {len(bits)}"
    print(f"  Loaded in {time.time()-t0:.1f}s")

    # ------------------------------------------------------------------
    # 2. Compute + write RLE
    # ------------------------------------------------------------------
    print("  Computing RLE ...", flush=True)
    t0 = time.time()
    runs = compute_rle(bits)
    rle_str = ",".join(map(str, runs))
    print(f"  {len(runs):,} runs  ({len(rle_str):,} bytes)  [{time.time()-t0:.1f}s]")

    # Verify round-trip
    recovered = rle_to_bits(runs)
    assert np.array_equal(bits, recovered), "RLE round-trip mismatch!"
    print("  Round-trip check: OK")

    rle_path.write_text(rle_str, encoding="ascii")
    print(f"  Written: {rle_path.name}")

    # ------------------------------------------------------------------
    # 3. Score the seed
    # ------------------------------------------------------------------
    if skip_score:
        seed_score = None
        print("  Scoring skipped (--skip-score)")
    else:
        print("  Scoring (suffix array) ...", flush=True)
        t0 = time.time()
        seed_score = score_bits(bits)
        elapsed = time.time() - t0
        print(f"  Seed score:   {seed_score:,}  [{elapsed:.1f}s]")

    # ------------------------------------------------------------------
    # 4. Update result JSON
    # ------------------------------------------------------------------
    if not json_path.exists():
        print(f"  WARNING: result JSON not found: {json_path}")
        return

    with open(json_path) as f:
        result = json.load(f)

    stored_score = get_stored_score(result)
    print(f"  Stored score: {stored_score:,}")

    changed = False

    # Add / update RLE field
    if result.get("solution_rle") != rle_str:
        result["solution_rle"] = rle_str
        changed = True
        print("  Added solution_rle field.")

    # Add / update score if seed scored higher (or score was missing)
    if seed_score is not None:
        if seed_score > stored_score:
            print(f"  Seed IMPROVES stored score by {seed_score - stored_score:,}!")
            if "score" in result:
                result["score"] = seed_score
            else:
                result["best_value"] = seed_score
            result["a_over_n_sq"] = seed_score / n**2
            changed = True
        elif seed_score == stored_score:
            print("  Seed score matches stored score. ✓")
        else:
            print(f"  WARNING: seed score ({seed_score:,}) is LOWER than stored "
                  f"({stored_score:,}) by {stored_score - seed_score:,}. "
                  "JSON left unchanged.")

    if changed:
        tmp = json_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(result, indent=2), encoding="utf-8")
        tmp.replace(json_path)
        print(f"  Updated: {json_path.name}")
    else:
        print("  No changes to JSON.")


def main():
    targets = TARGETS

    # Optional: filter to a specific n passed as command-line argument
    if len(sys.argv) > 1:
        try:
            requested_n = int(sys.argv[1].replace("_", "").replace(",", ""))
        except ValueError:
            print(f"Usage: python -m src.tools.export_large_n_rle [n]")
            sys.exit(1)
        if requested_n not in TARGETS:
            print(f"Unknown n={requested_n}. Valid values: {list(TARGETS)}")
            sys.exit(1)
        targets = {requested_n: TARGETS[requested_n]}

    skip_score = "--skip-score" in sys.argv

    for n, (npy_name, json_name) in sorted(targets.items()):
        process(n, npy_name, json_name, skip_score=skip_score)

    print("\nDone.")


if __name__ == "__main__":
    main()
