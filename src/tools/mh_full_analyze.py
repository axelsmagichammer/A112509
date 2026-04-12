import json
import os
from collections import defaultdict

def find_runs(s):
    if not s: return []
    runs = []; cur = s[0]; cnt = 1
    for c in s[1:]:
        if c == cur:
            cnt += 1
        else:
            runs.append((cur, cnt)); cur = c; cnt = 1
    runs.append((cur, cnt))
    return runs

def extract(s):
    runs = find_runs(s)
    blocks = []; seps = []
    i = 0
    while i < len(runs):
        ch, ln = runs[i]
        if ch == '1':
            blocks.append(ln); i += 1
            if i < len(runs) and runs[i][0] == '0':
                seps.append(runs[i][1]); i += 1
        else:
            break
    return blocks, seps

_BASE = os.path.join(os.path.dirname(__file__), "..", "..")
with open(os.path.join(_BASE, 'results', 'mh_unbounded', 'n_0198_results.json')) as f:
    data = json.load(f)

structures = [extract(s) for s in data['solutions']]
n_sols = len(structures)
n = data['n']

max_blocks = max(len(b) for b, _ in structures)
min_blocks = min(len(b) for b, _ in structures)
print(f"n={n}: {n_sols} solutions  (best_value={data['best_value']})")
print(f"1-block counts: {min_blocks} to {max_blocks}")
print()

print("ALL BLOCK POSITIONS (full observed range):")
for i in range(max_blocks):
    sizes = [blocks[i] for blocks, _ in structures if len(blocks) > i]
    cov = len(sizes)
    freq = defaultdict(int)
    for s in sizes: freq[s] += 1
    lo, hi = min(sizes), max(sizes)
    pct = 100 * cov // n_sols
    pct_str = f"{pct}%" if pct > 0 or cov == 0 else "<1%"
    if lo == hi:
        print(f"  Block {i+1:2d}: {lo}           (fixed, {pct_str} coverage)")
    else:
        print(f"  Block {i+1:2d}: {lo}-{hi}  freq={dict(sorted(freq.items()))}  coverage={pct_str}")

print()
print("ALL SEPARATOR POSITIONS (full observed range):")
max_seps = max(len(s) for _, s in structures)
for i in range(max_seps):
    sizes = [seps[i] for _, seps in structures if len(seps) > i]
    cov = len(sizes)
    freq = defaultdict(int)
    for s in sizes: freq[s] += 1
    lo, hi = min(sizes), max(sizes)
    pct = 100 * cov // n_sols
    pct_str = f"{pct}%" if pct > 0 or cov == 0 else "<1%"
    if lo == hi:
        print(f"  Sep {i+1:2d}: {lo}           (fixed, {pct_str} coverage)")
    else:
        print(f"  Sep {i+1:2d}: {lo}-{hi}  freq={dict(sorted(freq.items()))}  coverage={pct_str}")

print()
total_ones = [sum(b) for b, _ in structures]
print(f"Total 1s range: {min(total_ones)} to {max(total_ones)}")
print(f"Density range:  {min(total_ones)/n:.4f} to {max(total_ones)/n:.4f}")
print(f"1-block counts: {min_blocks} to {max_blocks}")
print(f"0-block counts: {min(len(s) for _,s in structures)} to {max(len(s) for _,s in structures)}")

# Find longest separator prefix common to 100% of solutions
all_seps = [s for _, s in structures]
K_common = 0
for k in range(1, max_seps + 1):
    prefix = all_seps[0][:k]
    if all(len(s) >= k and s[:k] == prefix for s in all_seps):
        K_common = k
    else:
        break
common_seps = all_seps[0][:K_common]
print(f"\nK_common (separator prefix length): {K_common}")
print(f"common_seps: {common_seps}")

print()
print("BLOCK RANGES for common prefix (first K_common blocks of ones):")
block_ranges = []
for i in range(K_common):
    sizes = [blocks[i] for blocks, _ in structures if len(blocks) > i]
    lo, hi = min(sizes), max(sizes)
    block_ranges.append([lo, hi])
    print(f"  Block {i+1:2d}: [{lo}, {hi}]")

print()
print("LEARNED BOUNDS JSON entry:")
max_zeros_block = max(max(s) for _, s in structures if s)
entry = {
    "target_n": n,
    "K_common": K_common,
    "common_seps": common_seps,
    "block_ranges": block_ranges,
    "min_total_1s": min(total_ones),
    "max_total_1s": max(total_ones),
    "max_zeros_block": max_zeros_block,
    "min_one_blocks": min_blocks,
    "max_one_blocks": max_blocks,
    "note": ""
}
print(json.dumps(entry, indent=2))
