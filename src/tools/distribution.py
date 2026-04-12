"""
distribution.py — Full distribution of b_n(S) across all 2^n binary strings.

b_n(S) = number of distinct integers representable by contiguous substrings
         of the n-bit binary string S.

For each n we iterate all 2^n bit strings, compute b_n(S), and record
the frequency of each value — giving a complete tally chart.

This is brute-force O(2^n × n²), so it is only practical for small n.
Rule-of-thumb feasibility on a modern laptop (with multiprocessing):
    n ≤ 18  →  seconds
    n ≤ 22  →  minutes
    n ≤ 25  →  tens of minutes
    n ≥ 26  →  hours / not recommended

Usage
-----
    python -m src.distribution 18
    python -m src.distribution 20 --plot --workers 8
    python -m src.distribution 22 --no-leading-zero-skip

API
---
    from src.distribution import compute_distribution, load_distribution, plot_distribution
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Force non-interactive matplotlib backend before any matplotlib import.
# This prevents a macOS hang caused by `system_profiler` font scanning.
os.environ.setdefault("MPLBACKEND", "Agg")


# ---------------------------------------------------------------------------
# Fast per-string evaluator (pure Python, no numpy overhead per call).
# For n ≤ 62 all substring values fit in a 64-bit int, so no bignum is made.
# ---------------------------------------------------------------------------

def _coverage_count_int_small(num: int, n: int) -> int:
    """Count distinct integer values from all substrings of the n-bit number `num`.

    Optimised for small n (≤ 62): uses direct bit arithmetic, no numpy overhead.
    Complexity: O(n²).
    """
    masks = [0] * (n + 1)
    for length in range(1, n + 1):
        masks[length] = (1 << length) - 1

    seen: set = set()
    for start in range(n):
        first_bit = (num >> (n - 1 - start)) & 1
        if first_bit == 0:
            seen.add(0)
            continue
        for length in range(1, n - start + 1):
            shift = n - start - length
            seen.add((num >> shift) & masks[length])
    return len(seen)


# ---------------------------------------------------------------------------
# Worker function (module-level so it is picklable for multiprocessing)
# ---------------------------------------------------------------------------

def _worker_chunk(args: Tuple[int, int, int]) -> List[int]:
    """Evaluate b_n(S) for all integers in [start, end).

    Returns a list of coverage counts, one per integer.
    """
    start_int, end_int, n = args
    return [_coverage_count_int_small(num, n) for num in range(start_int, end_int)]


# ---------------------------------------------------------------------------
# Main distribution computation
# ---------------------------------------------------------------------------

def compute_distribution(
    n: int,
    num_workers: Optional[int] = None,
    chunk_size: int = 4096,
    verbose: bool = True,
) -> Dict[int, int]:
    """Compute the full distribution of b_n(S) over all 2^n bit strings.

    Parameters
    ----------
    n:
        Bit length.  Recommended n ≤ 22; above 25 takes a very long time.
    num_workers:
        Parallel worker processes.  Defaults to cpu_count.
    chunk_size:
        Strings per work unit.  Larger = less IPC overhead; smaller = better
        progress granularity.
    verbose:
        Print progress to stdout.

    Returns
    -------
    dict mapping coverage_count → frequency (number of n-bit strings
    with that coverage).
    """
    if n > 26:
        print(
            f"WARNING: n={n} requires evaluating {2**n:,} strings.  "
            "This may take many hours.",
            flush=True,
        )

    total = 2 ** n
    workers = num_workers or (os.cpu_count() or 1)
    workers = min(workers, total // chunk_size + 1)

    # Build work chunks: [0, total) split into blocks of chunk_size.
    chunks = []
    for start_int in range(0, total, chunk_size):
        end_int = min(start_int + chunk_size, total)
        chunks.append((start_int, end_int, n))

    if verbose:
        print(
            f"n={n}:  {total:,} strings  |  "
            f"{len(chunks)} chunks  |  {workers} workers",
            flush=True,
        )

    freq: Counter = Counter()
    t0 = time.time()
    done = 0

    with multiprocessing.Pool(processes=workers) as pool:
        try:
            for counts in pool.imap_unordered(_worker_chunk, chunks,
                                              chunksize=max(1, len(chunks) // (workers * 8))):
                for c in counts:
                    freq[c] += 1
                done += len(counts)
                if verbose:
                    elapsed = time.time() - t0
                    pct = 100 * done / total
                    rate = done / max(elapsed, 1e-9)
                    eta = (total - done) / rate if rate > 0 else 0
                    print(
                        f"  {done:>{len(str(total))},}/{total:,}  "
                        f"({pct:5.1f}%)  "
                        f"elapsed={elapsed:.1f}s  eta={eta:.0f}s    ",
                        end="\r",
                        flush=True,
                    )
        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            if verbose:
                print(f"\n[Interrupted at {done:,}/{total:,} strings]", flush=True)
            if not freq:
                return {}

    if verbose:
        print(f"\n  Done in {time.time() - t0:.1f}s", flush=True)

    return dict(freq)


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------

def save_distribution(
    n: int,
    dist: Dict[int, int],
    results_dir: str = "results",
) -> Path:
    """Save distribution to results/distribution_n_XXXX.json.

    Returns the path written.
    """
    results_dir_path = Path(results_dir)
    results_dir_path.mkdir(exist_ok=True)
    path = results_dir_path / f"distribution_n_{n:04d}.json"

    total = sum(dist.values())
    if total == 0:
        return path

    values = sorted(dist.keys())
    mean = sum(k * v for k, v in dist.items()) / total
    mode = max(dist, key=lambda k: dist[k])
    variance = sum(dist[k] * (k - mean) ** 2 for k in dist) / total

    data = {
        "n": n,
        "total_strings": total,
        "min_coverage": min(values),
        "max_coverage": max(values),
        "mean_coverage": round(mean, 6),
        "mode_coverage": mode,
        "std_coverage": round(math.sqrt(variance), 6),
        "distribution": {str(k): dist[k] for k in values},
        "computed_at": datetime.now().isoformat(),
    }

    with open(path, "w") as f:
        json.dump(data, f, indent=2)

    return path


def load_distribution(
    n: int,
    results_dir: str = "results",
) -> Optional[Dict]:
    """Load a previously saved distribution file.

    Returns the full JSON dict, or None if not found.
    """
    path = Path(results_dir) / f"distribution_n_{n:04d}.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def print_summary(n: int, dist: Dict[int, int]) -> None:
    """Print a human-readable summary of the distribution."""
    if not dist:
        print("(empty distribution)")
        return

    total = sum(dist.values())
    values = sorted(dist.keys())
    mean = sum(k * v for k, v in dist.items()) / total
    variance = sum(dist[k] * (k - mean) ** 2 for k in dist) / total
    std = math.sqrt(variance)
    mode = max(dist, key=lambda k: dist[k])

    # Percentiles
    cum = 0
    p25 = p50 = p75 = p95 = None
    for v in values:
        cum += dist[v]
        pct = 100 * cum / total
        if p25 is None and pct >= 25:
            p25 = v
        if p50 is None and pct >= 50:
            p50 = v
        if p75 is None and pct >= 75:
            p75 = v
        if p95 is None and pct >= 95:
            p95 = v

    print(f"\n{'='*60}")
    print(f"  Distribution of b_{n}(S)  —  all {total:,} n-bit strings")
    print(f"{'='*60}")
    print(f"  min         : {min(values)}")
    print(f"  max (a({n}))  : {max(values)}")
    print(f"  mean        : {mean:.3f}")
    print(f"  std dev     : {std:.3f}")
    print(f"  mode        : {mode}  (freq={dist[mode]:,},  {100*dist[mode]/total:.1f}%)")
    print(f"  p25 / p50 / p75 / p95 : {p25} / {p50} / {p75} / {p95}")
    print(f"  distinct coverage values : {len(values)}")
    print(f"{'='*60}")

    # Top 10 most frequent values
    top = sorted(dist.items(), key=lambda x: -x[1])[:10]
    print(f"\n  Top 10 most frequent b_{n} values:")
    print(f"  {'coverage':>10}  {'count':>10}  {'%':>6}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*6}")
    for v, cnt in top:
        print(f"  {v:>10}  {cnt:>10,}  {100*cnt/total:>5.1f}%")
    print()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_distribution(
    n: int,
    dist: Dict[int, int],
    output_dir: str = "output",
    show: bool = False,
) -> Optional[Path]:
    """Plot the distribution as a bar chart + save to output/distribution_n_XXXX.png.

    Returns the saved path, or None if matplotlib is unavailable.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        import numpy as np
    except ImportError:
        print("matplotlib / numpy not available — skipping plot.", flush=True)
        return None

    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    save_path = output_path / f"distribution_n_{n:04d}.png"

    total = sum(dist.values())
    values = sorted(dist.keys())
    freqs = [dist[v] for v in values]
    mean = sum(k * v for k, v in dist.items()) / total
    variance = sum(dist[k] * (k - mean) ** 2 for k in dist) / total
    std = math.sqrt(variance)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f"Distribution of $b_{{{n}}}(S)$ over all $2^{{{n}}}$ binary strings",
        fontsize=14,
    )

    # --- Left: bar chart (frequency) ---
    ax = axes[0]
    bar_colors = ["#e74c3c" if v == max(values) else "#3498db" for v in values]
    ax.bar(values, freqs, color=bar_colors, alpha=0.85, width=max(1, (max(values)-min(values))//200))
    ax.axvline(mean, color="orange", linestyle="--", linewidth=1.5, label=f"mean={mean:.1f}")
    ax.set_xlabel(f"$b_{{{n}}}(S)$  (distinct integers covered)", fontsize=11)
    ax.set_ylabel("Number of strings", fontsize=11)
    ax.set_title("Frequency distribution")
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.annotate(
        f"max = {max(values)}\n(a({n}))",
        xy=(max(values), dist[max(values)]),
        xytext=(-60, 20),
        textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color="red"),
        color="red",
        fontsize=9,
    )

    # --- Right: normalised density + cumulative ---
    ax2 = axes[1]
    density = np.array(freqs) / total
    cumulative = np.cumsum(density)

    ax2.bar(values, density, color="#3498db", alpha=0.7, width=max(1, (max(values)-min(values))//200), label="density")
    ax2_r = ax2.twinx()
    ax2_r.plot(values, cumulative, color="darkorange", linewidth=2, label="CDF")
    ax2_r.set_ylabel("Cumulative fraction", fontsize=10)
    ax2_r.set_ylim(0, 1.05)
    ax2_r.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))

    ax2.set_xlabel(f"$b_{{{n}}}(S)$", fontsize=11)
    ax2.set_ylabel("Fraction of strings", fontsize=11)
    ax2.set_title("Density + CDF")

    # Annotation box
    textstr = (
        f"n = {n}\n"
        f"strings = $2^{{{n}}}$ = {total:,}\n"
        f"min = {min(values)}\n"
        f"max = {max(values)}\n"
        f"mean = {mean:.1f}\n"
        f"std = {std:.1f}"
    )
    props = dict(boxstyle="round", facecolor="lightyellow", alpha=0.7)
    ax2.text(
        0.03, 0.97, textstr, transform=ax2.transAxes,
        fontsize=8.5, verticalalignment="top", bbox=props,
    )

    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_r.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"  Plot saved: {save_path}", flush=True)

    if show:
        try:
            plt.switch_backend("TkAgg")  # try interactive; falls back gracefully
            plt.show()
        except Exception:
            pass
    plt.close(fig)

    return save_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute full distribution of b_n(S) for all 2^n binary strings"
    )
    p.add_argument("n", type=int, help="Bit length")
    p.add_argument("--workers", type=int, default=None,
                   help="Parallel workers (default: cpu_count)")
    p.add_argument("--results_dir", type=str, default="results")
    p.add_argument("--output_dir", type=str, default="output")
    p.add_argument("--plot", action="store_true", default=True,
                   help="Save a distribution plot (default: on)")
    p.add_argument("--no-plot", dest="plot", action="store_false")
    p.add_argument("--show", action="store_true", help="Display plot interactively")
    p.add_argument("--quiet", action="store_true")
    return p


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    n = args.n
    verbose = not args.quiet

    if n < 1:
        print("n must be ≥ 1")
        sys.exit(1)
    if n > 28:
        print(f"n={n} would require evaluating {2**n:,} strings — not recommended.")
        print("Aborting. Use --force (not implemented) to override.")
        sys.exit(1)

    # Check for existing file
    existing = load_distribution(n, args.results_dir)
    if existing is not None:
        print(f"Loaded existing distribution from results/distribution_n_{n:04d}.json")
        dist = {int(k): v for k, v in existing["distribution"].items()}
    else:
        dist = compute_distribution(
            n,
            num_workers=args.workers,
            verbose=verbose,
        )
        if dist:
            path = save_distribution(n, dist, args.results_dir)
            print(f"  Saved: {path}")

    if dist:
        print_summary(n, dist)
        if args.plot:
            plot_distribution(n, dist, args.output_dir, show=args.show)


if __name__ == "__main__":
    # Required on macOS/Windows for multiprocessing spawn safety
    multiprocessing.freeze_support()
    main()
