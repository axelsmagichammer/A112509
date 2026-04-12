"""Generate OEIS-formatted b-files in data/oeis/ from cached_results.json and mh_bounded results."""

import json
import os
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent.parent
DATA = BASE / "data"
OEIS_DIR = DATA / "oeis"
CACHE_FILE = DATA / "cached_results.json"
MH_DIR = BASE / "results" / "mh_bounded"

OEIS_DIR.mkdir(exist_ok=True)

# Load cached results (n=1-150, has all sequences)
with open(CACHE_FILE) as f:
    cache = json.load(f)

# Load mh_bounded results for n=151-200
mh_data = {}
for n in range(151, 201):
    path = MH_DIR / f"n_{n:04d}_results.json"
    if path.exists():
        with open(path) as f:
            d = json.load(f)
        mh_data[n] = d

print(f"Loaded {len(cache)} cached entries, {len(mh_data)} mh_bounded entries")


def write_bfile(filename, seq_id, description, values, offset=1):
    """Write a single OEIS b-file."""
    path = OEIS_DIR / filename
    with open(path, "w", newline="\n") as f:
        f.write(f"# OEIS {seq_id}\n")
        f.write(f"# {description}\n")
        f.write(f"# Source: Computed values from A112509 research project\n")
        f.write(f"#\n")
        f.write(f"# Format: n, a(n)\n")
        for n, val in sorted(values.items()):
            f.write(f"{n} {val}\n")
    print(f"  Wrote {path.name}: n={min(values)}-{max(values)} ({len(values)} values)")


# ── A112509: max distinct integers from n-bit binary substrings ──
a112509 = {}
for n_str, entry in cache.items():
    a112509[int(n_str)] = entry["a(n)"]
for n, d in mh_data.items():
    a112509[n] = d["best_value"]

write_bfile("b112509.txt", "A112509",
            "Maximum number of distinct integers represented by substrings of an n-bit binary number",
            a112509)

# ── A112510: smallest n-bit number achieving A112509(n) ──
# Only from cache (needs all optimal strings to determine min)
a112510 = {}
for n_str, entry in cache.items():
    a112510[int(n_str)] = entry["A112510_min_decimal"]

write_bfile("b112510.txt", "A112510",
            "Smallest n-bit number achieving A112509(n) distinct substring values",
            a112510)

# ── A112511: largest n-bit number achieving A112509(n) ──
# Only from cache (needs all optimal strings to determine max)
a112511 = {}
for n_str, entry in cache.items():
    a112511[int(n_str)] = entry["A112511_max_decimal"]

write_bfile("b112511.txt", "A112511",
            "Largest n-bit number achieving A112509(n) distinct substring values",
            a112511)

# ── A156022: max positive integers = A112509(n)-1 for n>=2, 1 for n=1 ──
a156022 = {}
for n_str, entry in cache.items():
    a156022[int(n_str)] = entry["A156022_positive_only"]
for n, d in mh_data.items():
    a156022[n] = 1 if n == 1 else d["best_value"] - 1

write_bfile("b156022.txt", "A156022",
            "Maximum number of distinct positive integers from substrings of an n-bit binary number",
            a156022)

# ── A156023: n*(n+1)/2 - A112509(n) ──
a156023 = {}
for n_str, entry in cache.items():
    a156023[int(n_str)] = entry["A156023_complement"]
for n, d in mh_data.items():
    a156023[n] = n * (n + 1) // 2 - d["best_value"]

write_bfile("b156023.txt", "A156023",
            "n*(n+1)/2 - A112509(n): complement count of binary substring values",
            a156023)

# ── A156024: n*(n+1)/2 - A156022(n) ──
a156024 = {}
for n_str, entry in cache.items():
    a156024[int(n_str)] = entry["A156024_complement_positive"]
for n, d in mh_data.items():
    a156022_val = 1 if n == 1 else d["best_value"] - 1
    a156024[n] = n * (n + 1) // 2 - a156022_val

write_bfile("b156024.txt", "A156024",
            "n*(n+1)/2 - A156022(n): complement count of positive binary substring values",
            a156024)

# ── A156025: number of n-bit strings achieving A112509(n) ──
# Only from cache (needs complete enumeration of all optimal strings)
a156025 = {}
for n_str, entry in cache.items():
    a156025[int(n_str)] = entry["num_optimal"]

write_bfile("b156025.txt", "A156025",
            "Number of n-bit binary strings achieving A112509(n) distinct substring values",
            a156025)

print(f"\nAll files written to {OEIS_DIR}")
