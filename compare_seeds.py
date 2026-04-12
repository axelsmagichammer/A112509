"""Compare the RLE structure of the good n=100M seed vs the bad n=1B seed."""
import numpy as np
from collections import Counter

def rle(bits):
    runs = []
    current = int(bits[0])
    count = 1
    for i in range(1, len(bits)):
        if bits[i] == current:
            count += 1
        else:
            runs.append((current, count))
            current = int(bits[i])
            count = 1
    runs.append((current, count))
    return runs

def analyze(name, bits, known_score=None):
    n = len(bits)
    ones = int(np.sum(bits == 1))
    zeros = int(np.sum(bits == 0))
    print(f"\n{'='*60}")
    print(f"  {name}: n={n:,}, ones={ones:,}, zeros={zeros:,}")
    print(f"  density={ones/n:.8f}")
    if known_score:
        print(f"  a(n)/n^2 = {known_score/(n*n):.8f}")
    print(f"{'='*60}")

    runs = rle(bits)
    ones_blocks = [c for v, c in runs if v == 1]
    zero_blocks = [c for v, c in runs if v == 0]
    print(f"Ones blocks: {len(ones_blocks)}")
    print(f"  min={min(ones_blocks):,}, max={max(ones_blocks):,}, "
          f"median={sorted(ones_blocks)[len(ones_blocks)//2]:,}")
    print(f"Zero blocks: {len(zero_blocks)}")
    print(f"  sizes: {dict(sorted(Counter(zero_blocks).items()))}")
    
    print(f"\nFirst 30 runs:")
    for i, (v, c) in enumerate(runs[:30]):
        label = "1s" if v == 1 else "0s"
        print(f"  [{i:2d}] {label}: {c:,}")
    
    # Check for counter pattern
    print(f"\nOnes blocks (first 20): {ones_blocks[:20]}")
    
    return runs

# Load and analyze n=100M (the good seed)
bits100m = np.load("seeds/n_100000000_seed.npy")
analyze("n=100M seed (good)", bits100m, known_score=4997781203807648)

# Load and analyze n=1B (the bad seed)
bits1b = np.load("seeds/seed_n1000000000.npy")
analyze("n=1B seed (bad)", bits1b, known_score=457255790442423061)
