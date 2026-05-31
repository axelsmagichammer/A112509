# OEIS A112509: Maximum Distinct Integers from Binary Substrings

> **A112509(n)** = the maximum number of distinct nonnegative integers whose binary representations occur as contiguous substrings of some binary string of length n.

This repository contains algorithms, analysis tools, and precomputed data for exploring [OEIS sequence A112509](https://oeis.org/A112509) and six related sequences. Exact proven values are certified through **n = 111**: Martin Fuller's A* computation certifies n ≤ 80, and this project's branch and bound solver certifies n = 81 to 111. Beyond n = 111 the stored values are best known heuristic records or certified lower bounds for explicit strings, not proofs of optimality.

## The Problem

Given a binary string of length n, every contiguous substring can be read as a binary number (with leading zeros collapsing to 0). For example, the length 4 string `1101` contains substrings:

| Length | Substrings | Integer values |
| -------- | ----------- | ---------------- |
| 1 | `1`, `1`, `0`, `1` | 1, 0 |
| 2 | `11`, `10`, `01` | 3, 2, 1 |
| 3 | `110`, `101` | 6, 5 |
| 4 | `1101` | 13 |

Distinct values: {0, 1, 2, 3, 5, 6, 13} → **7 distinct integers**, which happens to be the maximum achievable for n = 4.

A112509 asks: for each n, what is the **maximum** count of distinct integers achievable by any binary string of length n?

## First 20 Values

```text
n:    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
a(n): 1   3   5   7  10  13  17  22  27  33  40  47  55  64  73  83  94 106 118 131
```

## Related OEIS Sequences

| Sequence | Description | Proven (a-file) | Best known (b-file) |
| ---------- | ------------- | :-: | :-: |
| **A112509** | Max distinct integers, including 0, from length n binary substrings | n = 1 to 111 | n = 1 to 200 |
| **A112510** | Smallest length n binary number achieving A112509(n) | n = 1 to 111 | n = 1 to 150 |
| **A112511** | Largest length n binary number achieving A112509(n) | n = 1 to 111 | n = 1 to 150 |
| **A156022** | Max distinct *positive* integers from length n binary substrings | n = 1 to 111 | n = 1 to 200 |
| **A156023** | n(n+1)/2 − A112509(n), the missing count | n = 1 to 111 | n = 1 to 200 |
| **A156024** | n(n+1)/2 − A156022(n), the missing positive count | n = 1 to 111 | n = 1 to 200 |
| **A156025** | Number of length n binary strings achieving A112509(n) | n = 1 to 111 | n = 1 to 150 |

The **a files** (`data/oeis/a*.txt`) contain only rigorously proven values (Fuller n ≤ 80, branch and bound n = 81 to 111). The **b files** (`data/oeis/b*.txt`) additionally include best known heuristic values beyond the proven range. The b file format permits best known values.

## Key Findings

- **Asymptotic growth**: lim a(n)/n² = 1/2 (proved; see [Clare (2026)](paper/A112509_asymptotic.pdf))
- **Quantitative lower bound**: a(n) ≥ n²/2 − 15n^{3/2} for all n ≥ 64 (proved, explicit constant)
- **Ones density for near optimal strings**: explicit near optimal strings have ones density 1 − O(n^{−1/2}) (proved)
- **Ones density for optimal strings**: every optimal string has ones density 1 − O(n^{−1/4}) (proved)
- **Zero count penalty**: any z zeros impose a penalty ≥ z(z+1)/2 on f(s) relative to the upper bound (proved)
- **Fragmentation penalty**: if the z zeros form m separate runs, the additional penalty is ≥ C(m, 2) (proved)
- **Structure in certified optima**: the optimal strings certified through n = 111 start with a long block of 1s, followed by decreasing 1 blocks separated by short 0 blocks. Their common zero separator prefix begins `[1, 2, 1, ...]` and grows with n. Best known heuristic strings beyond n = 111 continue the same pattern.
- **De Bruijn connection**: for certified optima through n = 111, the minimum de Bruijn order k matches the leading 1 block length. Best known heuristic records beyond n = 111 continue to support this pattern, but do not prove it.
- **Lower bounds at scale**: a(1,000,000,000) ≥ 499,943,911,277,037,650
- **Second differences**: Δ²a(n) = Δa(n) − Δa(n−1) ∈ {0, 1, 2} for the certified exact values through n = 111 and for the best known continuation through n = 200. In that continuation, the value 2 occurs at n = 52, 71, 117, 157, 178, 191, and 194.
- **Heuristic structure detail**: leading 1 block lengths and later 1 block sizes in the best known strings scale approximately linearly with n, and the overall 1 density tends to 1.

## Elementary Bounds

Two simple bounds bracket the asymptotic growth rate of a(n).

### Upper bound

**Theorem.** For all n ≥ 1, a(n) ≤ n(n+1)/2, and therefore lim sup a(n)/n² ≤ 1/2.

**Proof.** Two substrings represent the same integer if and only if they are identical after stripping leading zeros. Consequently, every distinct positive integer value corresponds to a unique substring beginning with `1`. There are two cases. If the string contains at least one `0` bit, then 0 is also represented, but by substrings beginning with `0`; at least one of the n(n+1)/2 substrings begins with `0`, so the number of leading-`1` substrings is at most n(n+1)/2 − 1, giving at most (n(n+1)/2 − 1) + 1 = n(n+1)/2 distinct values. If the string is all ones, there is no value 0, so distinct values equals distinct positive values ≤ n(n+1)/2. In either case a(n) ≤ n(n+1)/2. Dividing by n² and taking the limit gives the stated bound.

### Lower bound

**Theorem.** For all n ≥ 2, a(n) ≥ ⌊n/2⌋·⌈n/2⌉ + ⌊n/2⌋ + 1 ≥ n²/4 + n/2, and therefore lim inf a(n)/n² ≥ 1/4.

**Proof.** Consider the explicit string s = 1^x 0^y for positive integers x, y with x + y = n (x ones followed by y zeros). The substrings of s fall into three disjoint types:

1. **All-ones:** the x strings `1`, `11`, ..., `1^x`, taking values 1, 3, 7, ..., 2^x − 1. These give x distinct odd positive integers.
2. **All-zeros:** all such substrings have value 0, contributing exactly one distinct value.
3. **Mixed:** substrings of the form `1^a 0^b` with a ∈ {1,...,x} and b ∈ {1,...,y}, taking value (2^a − 1)·2^b. These xy values are pairwise distinct (different b gives different 2-adic valuation; different a with the same b gives different odd part). Mixed values are even and positive, so they are disjoint from both the all-ones values and 0.

The total is f(s) = x + 1 + xy = x(n − x + 1) + 1. Setting x = ⌊n/2⌋ gives a lower bound of n²/4 + n/2.

### Combining the bounds

Together these give 1/4 + O(1/n) ≤ a(n)/n² ≤ 1/2 + O(1/n). The elementary lower bound is far from tight in general. The tight asymptotic is established in the companion paper:

> **Theorem** (Clare, 2026). lim a(n)/n² = 1/2.

The proof constructs an explicit family of length n binary strings with strictly decreasing one block lengths, achieving a(n) ≥ n²/2 − 15n^{3/2} for all n ≥ 64. The paper also proves that every optimal string has ones density 1 − O(n^{−1/4}), and derives zero count and zero fragmentation penalties on f(s). See [A112509_asymptotic.pdf](paper/A112509_asymptotic.pdf) for the full proofs.

## Asymptotic Numerical Evidence

The following certified lower bounds, computed using the surgical nudge algorithm on optimised seeds, confirm the proved limit and illustrate the rate of convergence:

| n | Lower bound on a(n) | a(n)/n² ≥ |
| --- | --- | --- |
| 3,000 | 4,261,266 | 0.4735 |
| 10,000 | 48,490,001 | 0.4849 |
| 1,000,000 | 498,396,662,200 | 0.4984 |
| 1,000,000,000 | 499,943,911,277,037,650 | 0.49994 |
| 2,000,000,000 | 1,999,841,243,789,594,574 | 0.49996 |

The computation at n = 2 × 10⁹ required scoring a two billion bit string (approximately 238 MB as a NumPy boolean array) using the `pydivsufsort` C library. The ratio a(n)/n² converges rapidly to the proved limit of 1/2 and is within 0.004% at n = 2 × 10⁹.

For comparison, certified and best known medium n values give these ratios:

```text
n=100: a(100)=3875,  a(n)/n²=0.3875
n=120: a(120)=5684,  a(n)/n²=0.3947
n=150: a(150)=9070,  a(n)/n²=0.4031
n=200: a(200)=16536, a(n)/n²=0.4134
```

Only the entries through n = 111 are exact. The moderate n values approach 0.5 more slowly than the large n lower bounds because the subleading correction terms are significant. Fitting a constrained log₂ tail model `a(n)/n² ≈ L − c/(log₂ n)^α` to the certified and best known data gives L ≈ 0.4997 with very tight confidence across multiple tail cutoffs, consistent with the limit being exactly 1/2.

## Installation

```bash
git clone <repo-url>
cd A112509
python -m venv .venv

# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt
```

### Requirements

```text
numpy>=1.24.0
matplotlib>=3.7.0
scipy>=1.10.0
sympy>=1.12
jupyter>=1.0.0
notebook>=6.5.0
ipykernel>=6.20.0
numba>=0.58.0
```

**Optional** (for n ≥ 1 billion scoring):

```bash
pip install pydivsufsort
```

## Repository Structure

```text
A112509/
├── README.md                     ← you are here
├── requirements.txt              ← Python dependencies
├── run_branch_bound.py           ← entry point: exact branch and bound solver
│
├── src/
│   ├── algorithms/               ← core search algorithms
│   │   ├── branch_bound_search.py    ← exact parallel branch and bound with SAM based UB (~2900 lines)
│   │   ├── brute_force.py            ← exact search for small n (≤ 25)
│   │   ├── evaluate.py               ← O(n log²n) distinct value counter
│   │   ├── structured_search.py      ← heuristic structural search, exact agreement for n < 112
│   │   ├── MH_algorithm.py           ← probabilistic Metropolis Hastings MCMC search (n ≤ 1000)
│   │   ├── greedy_search.py          ← template based greedy with learned prefixes
│   │   ├── large_n_lower_bound.py    ← certified lower bounds for n ≥ 3,000
│   │   └── surgical_nudge.py         ← LCP targeted local refinement
│   │
│   └── tools/                    ← analysis and utility scripts
│       ├── generate_oeis_bfiles.py   ← export to OEIS b-file format
│       ├── enhance_cached_results.py ← add derived sequences to cache
│       ├── debruijn_analysis.py      ← de Bruijn embedding analysis
│       ├── compute_hamming.py        ← pairwise Hamming distance stats
│       ├── distribution.py           ← full distribution of b_n(S) for small n
│       └── ...
│
├── config/
│   ├── learned_bounds.json       ← structural constraints for block search
│   └── search_constraints.json   ← per-n constraints for extended search
│
├── data/
│   ├── cached_results.json       ← exact through n = 111, best known after that
│   ├── oeis/                     ← OEIS a files (proven) and b files (best known)
│   └── reference/
│       ├── known_values.py       ← proven exact values used directly in the solver (n = 1 to 110)
│       └── A156025_values.py     ← optimal string counts (n = 1 to 80)
│
├── results/
│   ├── branch_bound_exact/       ← proven branch and bound results (n = 5 to 111+)
│   ├── mh_bounded/               ← constrained MH results (n = 81 to 200)
│   ├── mh_unbounded/             ← unconstrained MH results (n ≥ 80)
│   ├── large_n/                  ← lower bounds (n = 3K to 2B)
│   ├── distributions/            ← full distributions for n = 12, 16, 20
│   ├── debruijn_analysis.json    ← de Bruijn embedding metadata
│   └── hamming_distances.json    ← Hamming distance statistics
│
├── seeds/                        ← best known bit strings for large n
│
├── paper/                        ← LaTeX: asymptotic proof & analysis
│
├── notebooks/
│   ├── A112509_Extended_Results.ipynb   ← tables, charts, structural analysis
│   └── A112509_Results_Comparison.ipynb ← agreement between methods
│
└── tests/
```

---

## Branch and Bound Exact Solver

The branch and bound (B&B) solver is the primary tool for proving exact values of a(n). It exhaustively searches the space of all length n binary strings beginning with `1`, pruning branches whose upper bound cannot exceed the current best known value. The result is a mathematically rigorous proof that the returned value is optimal, not merely a heuristic estimate.

### Quick Start

```bash
# Prove a(n) and find all optimal strings:
python run_branch_bound.py --n 112

# Faster run without collecting all optimal strings:
python run_branch_bound.py --n 112 --no-collect-all

# Override worker count:
python run_branch_bound.py --n 112 --workers 8
```

Results are saved to `results/branch_bound_exact/n_NNNN_results.json`.

### Multi Machine Sharding

For very large n, the search can be distributed across multiple machines via filesystem based sharding:

```bash
# Machine 1: enumerate tasks
python run_branch_bound.py --n 120 --shard-mode coordinator

# Machines 1 to N: run workers (each claims tasks from shared directory)
python run_branch_bound.py --n 120 --shard-mode worker --shard-dir /shared/path/

# After all workers finish: merge results
python run_branch_bound.py --n 120 --shard-mode merge --shard-dir /shared/path/
```

### Algorithm Overview

The solver operates in three phases:

```text
Phase 1: BFS Split        Phase 2: Parallel DFS         Phase 2b: Tail Split
┌──────────────────┐      ┌───────────────────────┐      ┌─────────────────────┐
│ Expand root to   │      │ Workers claim subtrees │      │ When few tasks left,│
│ depth d, prune   │─────▶│ from queue and run     │─────▶│ split remaining     │
│ unpromising       │      │ iterative DFS with     │      │ tasks for load      │
│ prefixes, split  │      │ SAM based UB pruning   │      │ balancing            │
│ hard subtrees    │      │ + shared incumbent     │      │                      │
└──────────────────┘      └───────────────────────┘      └─────────────────────┘
```

#### Phase 1: BFS Enumeration and Task Creation

Starting from the root prefix `(1,)`, BFS expands the search tree to a configurable split depth (auto-selected as approximately n/3, clamped to [10, 25]). At each level, children whose upper bound is below the current incumbent are pruned. The surviving leaf prefixes become independent DFS tasks.

**Hard subtree splitting:** Subtrees with upper bounds only slightly above the incumbent (within `hard_split_gap`, default 5) are classified as hard. They are unlikely to contain the optimum but cannot be pruned. These are expanded an additional `hard_split_extra` levels (default 7), breaking them into many fine grained tasks that finish quickly and keep workers busy.

**Probe based adaptive splitting** (optional): A shallow DFS probe can be run on each subtree to classify it as easy or hard based on behavioural metrics (prune rate, residual stack size). Hard tasks get deeper splitting; all tasks are ranked by priority so the most promising subtrees execute first.

#### Phase 2: Parallel DFS Workers

Each task is an independent iterative DFS rooted at a BFS frontier prefix. Workers run in parallel across all available CPU cores (default: `os.cpu_count() - 1`).

**Worker architecture:**

- Each worker uses an **explicit LIFO stack** (not a heap), giving O(n²) memory and no risk of OOM.
- Stack entries store `(prefix_bits, precomputed_ub)`. The upper bound is computed once at push time and rechecked against the monotonically rising incumbent at pop time.
- The entire DFS inner loop is compiled to native machine code via **Numba JIT** (`@njit(nogil=True, cache=True)`).
- A **shuttle thread** bridges `multiprocessing.Value` ↔ `numpy.int64[1]` every 0.2 seconds, allowing Numba code that is free of the GIL to poll the cross process shared incumbent without Python overhead.

**Shared incumbent:** A single `multiprocessing.Value('l')` holds the best score found by any worker. When one worker finds a better solution, all other workers see the update within 0.2 seconds and can prune more aggressively. This is critical for `collect_all=True` mode where workers must know the exact optimum to decide which strings to keep.

**Cooperative chunking** (optional): If `max_nodes_per_task > 0`, a DFS worker that exceeds its node budget snapshots its remaining stack and returns. The orchestrator re-enqueues the snapshot entries as fresh tasks, ensuring no single subtree monopolises a worker for hours.

#### Phase 2b: Tail Split

When the number of remaining tasks drops below a threshold (~10% of the original count), every queued task is expanded by BFS for 3 additional levels. This breaks stragglers into fine grained work to keep all cores busy through the long tail. Up to 5 rounds of additional splitting can occur.

### Upper Bound Formula

The upper bound is the heart of B&B pruning. For a prefix p of length k, extended by L lookahead bits to form q = p || e of length k + L, with m = n − k − L bits remaining:

```text
UB_L(p) = max over all e in {0,1}^L of [ score(q) + a(m) + ones(q) × m − SAM_savings ]
```

**Components:**

| Term | Meaning |
| --- | --- |
| score(q) | Exact distinct value count of the (k+L) bit prefix q, computed via the suffix automaton |
| a(m) | Proven exact value for a length m binary string (from `known_values.py` + B&B results). When unknown, falls back to the loose bound n(n+1)/2 − (k+L)(k+L+1)/2 |
| ones(q) × m | Upper bound on cross boundary values: each `1` bit in q can contribute at most m new distinct values via substrings that start in q and extend into the unknown suffix |
| SAM_savings | Cross boundary values provably already present in the suffix automaton of q (see below) |

**Additional tightening mechanisms:**

- **Global UB cap**: a sequence level upper bound derived from the maximum possible distinct values across all substring lengths.
- **Parent UB clamp**: `child_ub = min(child_ub, parent_ub)`. Any descendant of a child is also a descendant of the parent, so the child's UB cannot exceed the parent's.
- **Exact tail**: when m ≤ 8 remaining bits, all 2^m completions are scored exactly instead of using the UB estimate.
- **Structural pruning**: prefixes containing ≥2 consecutive zeros after the leading 1 block are pruned (provably suboptimal for n ≥ 8).

### Suffix Automaton (SAM): Dual Role

The suffix automaton is used in two distinct ways:

#### 1. Exact Scoring

For a complete length n binary string, the SAM is built online (one character at a time) and the distinct value count is computed via a topological subtree sum on the `1` child of the root:

```text
f(s) = subtree_count(t1[root]) + 𝟙[0 ∈ s]
```

This counts all distinct substrings starting with `1` (each corresponding to a unique positive integer), plus 1 for the value 0 if any `0` bit exists.

#### 2. Upper Bound Computation (SAM Savings)

When computing the UB for a prefix q, the SAM is used to tighten the cross boundary estimate. For each `1` bit at position i in q, the suffix starting there can produce substrings that extend into the unknown suffix. Some of these extensions are guaranteed to already exist in the SAM of q regardless of what the suffix contains:

1. Walk the suffix link chain to find the SAM state representing the suffix of q starting at position i.
2. BFS from that state up to JMAX steps (2 for base UB, 3 for refined UB): check whether all 2^j possible j-step continuations exist as transitions in the SAM.
3. Each fully covered step means one fewer new distinct value can be produced by extending into the suffix. Subtract 1 from the cross boundary budget.
4. **Caching**: `state_cov[v]` memoises the coverage depth per SAM state within each extension evaluation, avoiding redundant BFS walks.

The net effect is that the SAM savings term subtracts a provably correct amount from the naive `ones(q) × m` cross boundary estimate, making the UB significantly tighter. This is the key innovation that makes the B&B solver tractable for n > 80.

### Lookahead

The solver uses a multi tier lookahead system:

| Tier | When Used | Lookahead L | SAM JMAX |
| --- | --- | --- | --- |
| Base | Default UB | 2 | 2 |
| Refined | Near threshold: incumbent ≤ UB ≤ incumbent + margin | L + 1 | 3 |
| Adaptive | Cascade: base first, refine only if UB is within margin of incumbent | 2 → 3 | 2 → 3 |

**Base lookahead** (L = 2): enumerate all 2^L = 4 two bit extensions of the prefix, build the SAM for each, compute the UB, and return the maximum. This is the default for every node.

**Refined lookahead** (L = 3): used only when the base UB is tantalisingly close to the incumbent (within `refine_margin`, default 2). The extra lookahead bit and deeper SAM BFS (JMAX = 3 vs 2) often prove that the branch is actually suboptimal, saving the cost of exploring it. The dynamic margin widens near leaf depth to catch more near misses.

### Provability

The B&B solver produces a mathematically rigorous proof that the returned a(n) is optimal:

1. **Completeness**: every length n binary string beginning with `1` is either explicitly scored or its entire subtree is pruned by an upper bound that is provably ≥ the true maximum achievable by any string in that subtree.
2. **UB correctness**: the upper bound formula uses only (a) exact scoring of known prefixes via the SAM, (b) proven values a(m) for smaller m, and (c) combinatorial arguments that overcount (never undercount) the number of new distinct values that the unknown suffix can add.
3. **a(m) integrity**: the `known_values.py` file contains only values proven by either Fuller (n ≤ 80) or prior B&B runs. Heuristic values from MH or structured search are **never** used in the UB formula. Using an underestimate of a(m) could cause incorrect pruning.
4. **Incumbent seeding**: the initial lower bound is seeded from MH heuristic results when available. This is safe because it only makes pruning more aggressive and never causes incorrect pruning. If no heuristic is available, the search still finds the correct answer; it just explores more nodes.
5. **`collect_all=True`**: when enabled, the solver continues searching even after finding the optimum, collecting every string that achieves it. The returned `num_optimal` count is exact (not a lower bound).

### Result Format

Each run saves a JSON file to `results/branch_bound_exact/`:

```json
{
  "a(n)": 4826,
  "num_optimal": 80,
  "optimal_strings": ["1111111111111111111101111111111...", "..."],
  "K_common": 10,
  "common_seps": [1, 2, 1, 1, 1, 1, 1, 1, 1, 1],
  "meta": {
    "n": 111,
    "preset": "stable",
    "lookahead": 2,
    "collect_all": true,
    "resolved_workers": 15,
    "nodes_expanded": 36956030926,
    "nodes_pruned": 20451182412,
    "num_subtrees": 2896,
    "elapsed_s": 212429.27
  }
}
```

### Performance

Representative timings on a 16-core machine (15 workers, `collect_all=True`, L=2):

| n | a(n) | # optimal | Nodes expanded | Time |
| --- | --- | --- | --- | --- |
| 20 | 131 | 19 | 58 | 8 s |
| 50 | 899 | not recorded | ~10⁵ | ~30 s |
| 80 | 2,423 | 10 | ~10⁷ | ~20 min |
| 100 | 3,875 | 535 | ~10⁹ | ~12 h |
| 111 | 4,826 | 80 | ~3.7×10¹⁰ | ~59 h |

The node count grows roughly 3 to 5× per increment of n, so each additional value beyond n = 111 requires substantially more compute.

---

## Scoring

All algorithms share the same inner loop. For small n (≲ 30) a brute force hash set over all O(n²) substrings is used. For larger n, the count is computed exactly in O(n log n) via the suffix array and LCP identity.

### Suffix Array (SA)

The **suffix array** SA of a string s of length n is a permutation of {0, 1, ..., n−1} such that the suffixes s[SA[0]..], s[SA[1]..], ..., s[SA[n−1]..] are in lexicographic order. This implementation uses `pydivsufsort`, a Python wrapper around the C library libdivsufsort, which achieves O(n) time with small constants.

### LCP Array and Kasai's Algorithm

The **LCP array** LCP has LCP[i] = length of the longest common prefix between the suffixes at consecutive SA positions. **Kasai's algorithm** computes the full LCP array in O(n) time and O(n) space from the suffix array. In this codebase Kasai's algorithm is compiled to native machine code with Numba's JIT compiler, giving a further 10 to 50× speedup over a pure Python implementation.

### The Scoring Identity

The key observation is that two substrings represent the same nonnegative integer if and only if they are identical after stripping leading zeros. It follows that:

```text
distinct values = (1 if any zero bit is present) + (# distinct substrings starting with '1')
```

Summing only over positions where the suffix starts with `1`:

```text
f(s) = Σ_{i: s[SA[i]] = '1'}  (n − SA[i] − LCP[i])  +  𝟙[0 ∈ s]
```

Each term `(n − SA[i] − LCP[i])` counts the substrings that begin at SA[i], start with `1`, and are not a prefix of the lexicographically preceding leading 1 suffix. The entire computation runs in O(n) time once the SA and LCP arrays are available, making the full scoring pipeline O(n log n).

## Other Algorithms

### 1. Brute Force (`src/algorithms/brute_force.py`)

Exhaustive enumeration of all 2^(n−1) candidate length n binary strings (leading bit is always 1). For each string, computes the set of distinct integer values from all O(n²) substrings.

**Use for**: n ≤ 20 (exact, complete enumeration).

```bash
python -m src.algorithms.brute_force
```

### 2. Evaluate (`src/algorithms/evaluate.py`)

Counts distinct substring values for a single bit string using a suffix array approach.

```bash
# Score a specific bit string
python -m src.algorithms.evaluate 11111101111111001111101011110110

# Pipe from stdin
echo 1101 | python -m src.algorithms.evaluate
```

### 3. Structured Search (`src/algorithms/structured_search.py`)

**Use for**: n ≤ 150 as a heuristic structural search. It agrees with every certified exact value for n < 112, but it is not a proof method for larger n.

Exploits an observed structural template: certified optima have a long leading 1 block followed by alternating zero separator blocks and decreasing 1 blocks. The block size bounds are inferred from observed optimal structures, reducing the search to enumerating integer tuples within tight windows. This makes the method extremely useful, and it agrees with all exact values through n = 111. For n ≥ 112 its results remain heuristic evidence unless branch and bound later proves them.

```bash
python -m src.algorithms.structured_search
```

### 4. Metropolis Hastings MCMC (`src/algorithms/MH_algorithm.py`)

**Use for**: n = 1 to ~1000 as a probabilistic heuristic. It agrees with every certified exact value through n = 111, but it does not certify optimality.

Operates on run length encodings with five move types (transfer, split, merge, swap, multi transfer). Geometric cooling with periodic reheating. Despite having no structural priors, unconstrained MH repeatedly finds the certified optima through n = 111 and agrees strongly with the structured search records beyond that. This is independent evidence for the block template, not a proof.

```bash
python -m src.algorithms.MH_algorithm
```

### 5. Greedy Search (`src/algorithms/greedy_search.py`)

Template based construction exploiting the linear scaling pattern. Fast initial lower bound for any n.

```bash
python -m src.algorithms.greedy_search
```

### 6. Large n Lower Bound (`src/algorithms/large_n_lower_bound.py`)

Computes mathematically valid lower bounds for n ≥ 3,000 by counting only a chosen subset of substrings exactly. Used as a library by other scripts (e.g. `score_1B.py`).

### 7. Surgical Nudge (`src/algorithms/surgical_nudge.py`)

Local refinement starting from a known good seed. Identifies the highest LCP collision in the suffix array and makes the smallest density preserving bit swap to break it. This is the main practical method for improving records at n ≥ 10⁸.

```bash
python -m src.algorithms.surgical_nudge
```

## Tools

| Script | Purpose |
| -------- | --------- |
| `generate_oeis_bfiles.py` | Export all sequences to OEIS b file format in `data/oeis/` |
| `enhance_cached_results.py` | Add A112510, A112511, and A156022 to A156024 to `cached_results.json` |
| `debruijn_analysis.py` | Compute minimum de Bruijn embedding order for each optimal string |
| `compute_hamming.py` | Compute min/max/mean pairwise Hamming distances between optimal solutions |
| `distribution.py` | Full distribution of distinct value counts across all 2^n strings |

Run any tool with:

```bash
python -m src.tools.<tool_name>
```

## Notebooks

### Extended Results (`notebooks/A112509_Extended_Results.ipynb`)

Interactive presentation of all computed results:

- Sequence values and all optimal strings or best known strings (tabbed by groups of 10)
- Structural analysis: block counts, 1 density, run lengths, K_common
- K_common vs n charts
- Hamming distance analysis between optimal solutions
- Growth analysis: Δa, ΔΔa, a(n)/n, a(n)/n²

### Results Comparison (`notebooks/A112509_Results_Comparison.ipynb`)

Compares three independent methods for n = 80 to 200:

- Unbounded Metropolis Hastings
- Bounded constrained Metropolis Hastings
- Heuristic structured search

## Data Files

### `data/oeis/`

Two classes of OEIS formatted data files:

- **a files** (`a*.txt`): contain only rigorously proven values (Fuller n ≤ 80, branch and bound n = 81 to 111). Suitable for OEIS submission.
- **b files** (`b*.txt`): contain best known values including heuristic estimates beyond the proven range. The b file format permits best known values.

### `results/branch_bound_exact/`

Proven branch and bound results for individual values of n. Each file contains `a(n)`, all optimal strings, structural metadata (`K_common`, `common_seps`), and solver statistics (nodes expanded, elapsed time, etc.).

### `data/cached_results.json`

Authoritative exact entries for n = 1 to 111, including structural metadata and all optimal strings. Entries after n = 111 are best known heuristic records or lower bound data and should not be cited as exact unless later certified by branch and bound.

### `data/reference/known_values.py`

Proven exact values used by the B&B solver's upper bound formula. Contains only values proven by Fuller (n ≤ 80) or prior branch and bound runs (n = 81 to 110). **Never** includes heuristic values. Using an underestimate of a(m) in the UB formula could cause incorrect pruning.

## Configuration

### `config/learned_bounds.json`

Structural constraints inferred from analysed optimal solutions and best known strings, used by the heuristic `structured_search.py`.

### `config/search_constraints.json`

Per n structural constraints for n ≥ 81, used by the bounded MH algorithm.

## Typical Workflows

### Prove the next exact value

```bash
# Run branch and bound for n = 112
python run_branch_bound.py --n 112

# Results saved to results/branch_bound_exact/n_0112_results.json
```

The solver automatically loads all previously proven values from `known_values.py` and `results/branch_bound_exact/` to use in its upper bound formula. Each new proven value strengthens the UB for subsequent runs.

### Chain multiple values

Run values sequentially:

```bash
python run_branch_bound.py --n 112
python run_branch_bound.py --n 113
python run_branch_bound.py --n 114
```

### Explore a new large n

1. Construct an initial seed using `greedy_search.py` or `large_n_lower_bound.py`
2. Refine with `surgical_nudge.py`
3. Score the final string with `evaluate.py`

### Check heuristic candidates

1. Run MH (bounded and/or unbounded) for the target n
2. Run branch and bound when a proof is required and computationally feasible
3. Compare in `A112509_Results_Comparison.ipynb`
4. Treat structured search and MH agreement as evidence, not certification

## Second Differences

Define the first difference Δa(n) = a(n) − a(n−1) and the second difference Δ²a(n) = Δa(n) − Δa(n−1). Since a(n) ≈ n²/2, one expects Δa(n) ≈ n and Δ²a(n) ≈ 1. The precise behaviour is striking.

> **Observation.** For all certified values n ≤ 111, and for the best known continuation through n ≤ 200, Δ²a(n) ∈ {0, 1, 2}. The values 0 and 1 alternate in a near regular pattern. The value 2 occurs at n = 52, 71, 117, 157, 178, 191, and 194 (bold in the tables below).

The predominance of Δ²a(n) ∈ {0, 1} reflects the fact that Δa(n) increases by exactly 0 or 1 at each step, a remarkable regularity for a sequence defined by a combinatorial optimisation. In the certified and best known data, the rare value 2 is associated with a simultaneous resolution of two de Bruijn order tiers at step n−1.

> **Conjecture.** Δ²a(n) ∈ {0, 1, 2} for all n ≥ 1.

### Table 1: a(n), Δa(n), Δ²a(n) for n = 1 to 50

| n | a(n) | Δa | Δ²a |
| --: | -----: | ---: | ----: |
| 1 | 1 | n/a | n/a |
| 2 | 3 | 2 | n/a |
| 3 | 5 | 2 | 0 |
| 4 | 7 | 2 | 0 |
| 5 | 10 | 3 | 1 |
| 6 | 13 | 3 | 0 |
| 7 | 17 | 4 | 1 |
| 8 | 22 | 5 | 1 |
| 9 | 27 | 5 | 0 |
| 10 | 33 | 6 | 1 |
| 11 | 40 | 7 | 1 |
| 12 | 47 | 7 | 0 |
| 13 | 55 | 8 | 1 |
| 14 | 64 | 9 | 1 |
| 15 | 73 | 9 | 0 |
| 16 | 83 | 10 | 1 |
| 17 | 94 | 11 | 1 |
| 18 | 106 | 12 | 1 |
| 19 | 118 | 12 | 0 |
| 20 | 131 | 13 | 1 |
| 21 | 145 | 14 | 1 |
| 22 | 160 | 15 | 1 |
| 23 | 176 | 16 | 1 |
| 24 | 192 | 16 | 0 |
| 25 | 209 | 17 | 1 |
| 26 | 227 | 18 | 1 |
| 27 | 246 | 19 | 1 |
| 28 | 265 | 19 | 0 |
| 29 | 285 | 20 | 1 |
| 30 | 306 | 21 | 1 |
| 31 | 328 | 22 | 1 |
| 32 | 351 | 23 | 1 |
| 33 | 375 | 24 | 1 |
| 34 | 399 | 24 | 0 |
| 35 | 424 | 25 | 1 |
| 36 | 450 | 26 | 1 |
| 37 | 477 | 27 | 1 |
| 38 | 504 | 27 | 0 |
| 39 | 532 | 28 | 1 |
| 40 | 561 | 29 | 1 |
| 41 | 591 | 30 | 1 |
| 42 | 622 | 31 | 1 |
| 43 | 654 | 32 | 1 |
| 44 | 687 | 33 | 1 |
| 45 | 720 | 33 | 0 |
| 46 | 754 | 34 | 1 |
| 47 | 789 | 35 | 1 |
| 48 | 825 | 36 | 1 |
| 49 | 862 | 37 | 1 |
| 50 | 899 | 37 | 0 |

Note: Δ²a(n) ∈ {0, 1} throughout this range.

### Table 2: a(n) and Δ²a(n) for n = 51 to 200

Bold entries mark occurrences of Δ²a(n) = 2. Values for n ≤ 111 are exact, with Fuller certifying n ≤ 80 and branch and bound certifying n = 81 to 111. Values for n = 112 to 200 are best known heuristic records from structured search and independent MH runs. They are not certified exact values.

| n | a(n) | Δ²a | n | a(n) | Δ²a | n | a(n) | Δ²a | n | a(n) | Δ²a |
| --: | -----: | ----: | --: | -----: | ----: | --: | -----: | ----: | --: | -----: | ----: |
| 51 | 937 | 1 | 89 | 3,034 | 0 | 127 | 6,402 | 1 | 165 | 11,070 | 1 |
| 52 | 977 | **2** | 90 | 3,106 | 1 | 128 | 6,508 | 1 | 166 | 11,211 | 1 |
| 53 | 1,017 | 0 | 91 | 3,179 | 1 | 129 | 6,615 | 1 | 167 | 11,353 | 1 |
| 54 | 1,058 | 1 | 92 | 3,253 | 1 | 130 | 6,723 | 1 | 168 | 11,496 | 1 |
| 55 | 1,100 | 1 | 93 | 3,328 | 1 | 131 | 6,832 | 1 | 169 | 11,639 | 0 |
| 56 | 1,143 | 1 | 94 | 3,404 | 1 | 132 | 6,942 | 1 | 170 | 11,783 | 1 |
| 57 | 1,186 | 0 | 95 | 3,480 | 0 | 133 | 7,052 | 0 | 171 | 11,928 | 1 |
| 58 | 1,230 | 1 | 96 | 3,557 | 1 | 134 | 7,163 | 1 | 172 | 12,074 | 1 |
| 59 | 1,275 | 1 | 97 | 3,635 | 1 | 135 | 7,275 | 1 | 173 | 12,221 | 1 |
| 60 | 1,321 | 1 | 98 | 3,714 | 1 | 136 | 7,388 | 1 | 174 | 12,369 | 1 |
| 61 | 1,368 | 1 | 99 | 3,794 | 1 | 137 | 7,502 | 1 | 175 | 12,518 | 1 |
| 62 | 1,416 | 1 | 100 | 3,875 | 1 | 138 | 7,617 | 1 | 176 | 12,668 | 1 |
| 63 | 1,465 | 1 | 101 | 3,957 | 1 | 139 | 7,733 | 1 | 177 | 12,818 | 0 |
| 64 | 1,514 | 0 | 102 | 4,040 | 1 | 140 | 7,850 | 1 | 178 | 12,970 | **2** |
| 65 | 1,564 | 1 | 103 | 4,124 | 1 | 141 | 7,968 | 1 | 179 | 13,122 | 0 |
| 66 | 1,615 | 1 | 104 | 4,209 | 1 | 142 | 8,087 | 1 | 180 | 13,275 | 1 |
| 67 | 1,667 | 1 | 105 | 4,295 | 1 | 143 | 8,207 | 1 | 181 | 13,429 | 1 |
| 68 | 1,720 | 1 | 106 | 4,381 | 0 | 144 | 8,328 | 1 | 182 | 13,584 | 1 |
| 69 | 1,774 | 1 | 107 | 4,468 | 1 | 145 | 8,450 | 1 | 183 | 13,740 | 1 |
| 70 | 1,828 | 0 | 108 | 4,556 | 1 | 146 | 8,572 | 0 | 184 | 13,897 | 1 |
| 71 | 1,884 | **2** | 109 | 4,645 | 1 | 147 | 8,695 | 1 | 185 | 14,055 | 1 |
| 72 | 1,941 | 1 | 110 | 4,735 | 1 | 148 | 8,819 | 1 | 186 | 14,214 | 1 |
| 73 | 1,998 | 0 | 111 | 4,826 | 1 | 149 | 8,944 | 1 | 187 | 14,374 | 1 |
| 74 | 2,056 | 1 | 112 | 4,918 | 1 | 150 | 9,070 | 1 | 188 | 14,535 | 1 |
| 75 | 2,115 | 1 | 113 | 5,011 | 1 | 151 | 9,197 | 1 | 189 | 14,697 | 1 |
| 76 | 2,175 | 1 | 114 | 5,104 | 0 | 152 | 9,325 | 1 | 190 | 14,859 | 0 |
| 77 | 2,236 | 1 | 115 | 5,198 | 1 | 153 | 9,454 | 1 | 191 | 15,023 | **2** |
| 78 | 2,298 | 1 | 116 | 5,293 | 1 | 154 | 9,584 | 1 | 192 | 15,187 | 0 |
| 79 | 2,360 | 0 | 117 | 5,390 | **2** | 155 | 9,714 | 0 | 193 | 15,352 | 1 |
| 80 | 2,423 | 1 | 118 | 5,487 | 0 | 156 | 9,845 | 1 | 194 | 15,519 | **2** |
| 81 | 2,487 | 1 | 119 | 5,585 | 1 | 157 | 9,978 | **2** | 195 | 15,686 | 0 |
| 82 | 2,552 | 1 | 120 | 5,684 | 1 | 158 | 10,111 | 0 | 196 | 15,854 | 1 |
| 83 | 2,618 | 1 | 121 | 5,784 | 1 | 159 | 10,245 | 1 | 197 | 16,023 | 1 |
| 84 | 2,685 | 1 | 122 | 5,885 | 1 | 160 | 10,380 | 1 | 198 | 16,193 | 1 |
| 85 | 2,753 | 1 | 123 | 5,987 | 1 | 161 | 10,516 | 1 | 199 | 16,364 | 1 |
| 86 | 2,822 | 1 | 124 | 6,090 | 1 | 162 | 10,653 | 1 | 200 | 16,536 | 1 |
| 87 | 2,892 | 1 | 125 | 6,193 | 0 | 163 | 10,791 | 1 | | | |
| 88 | 2,963 | 1 | 126 | 6,297 | 1 | 164 | 10,930 | 1 | | | |

## Open Problems

The following open problems are suggested by the computational findings.

**(i) Asymptotic rate.** The limit lim a(n)/n² = 1/2 is proved (Clare, 2026). The quantitative bound is a(n) ≥ n²/2 − 15n^{3/2} for n ≥ 64, but the constant 15 is not tight and the true subleading coefficient is unknown. Closing the gap between the lower bound O(n^{3/2}) and the upper bound O(n) would be progress.

**(ii) Separator structure.** Prove or disprove that the zero separator sequence (1, 2, 1, 1, ...) is universal for all sufficiently large n. Can the full separator sequence be determined analytically? The second separator being 2 while all others are 1 is a highly nontrivial constraint whose origin is not understood.

**(iii) Second differences.** Prove that Δ²a(n) ∈ {0, 1, 2} for all n. This would require understanding the growth rate of Δa(n) with precision.

**(iv) Exact values.** Extend the branch and bound computation beyond n = 111. Each additional value requires a substantial increase in the search budget; reaching n = 150 by B&B appears to require many CPU months with the current solver. Multi machine sharding may make this tractable.

**(v) One bit density.** It is now proved that every optimal string has ones density 1 − O(n^{−1/4}) (Clare, 2026). The empirical rate at large n is far tighter: the best known seeds have density 0.999950 at n = 10⁹ and 0.999964 at n = 2 × 10⁹, suggesting the true rate is 1 − O(n^{−1}) or faster. Proving this sharper rate remains open.

**(vi) Alphabet generalisation.** What is the analogous sequence for ternary (base 3) or higher base alphabets? Does the ratio a_q(n)/n² for base q approach (q−1)/2 or some other constant?

## Project History and AI Assistance

This project spanned several months of active research and consumed thousands of hours of CPU time across multiple machines. What began as a straightforward brute force exploration of the sequence grew into a layered computational investigation requiring increasingly sophisticated algorithms at each scale barrier.

The work was carried out with assistance from **Claude Sonnet 4.6** and **Claude Opus 4.6** (Anthropic). These AI systems contributed throughout the project with performance optimisation: identifying bottlenecks and profiling.

### Algorithm Evolution

The pipeline was not designed upfront, it evolved iteratively as each approach hit its limits:

1. **Brute force**: exhaustive enumeration of all 2^(n−1) strings worked to n ≈ 40 but was computationally intractable beyond that.

2. **Block structure discovery**: manual inspection of brute force solutions revealed the run length pattern. This motivated the structured search heuristic. It agrees with all certified exact values less than n = 112, but it does not by itself prove global optimality as the whole search space is not searched.

3. **Suffix array scoring**: as n grew past 30, the O(n²) substring hash set became the bottleneck. The key insight that distinct positive integer values equal distinct leading 1 suffixes in the suffix array, countable via the LCP array, reduced scoring from O(n²) to O(n log n).

4. **Metropolis Hastings MCMC**: to escape the constraint design bottleneck and provide independent evidence, an unconstrained MH search was built operating directly on run length encodings. Five move types (transfer, split, merge, swap, multi transfer) ensure ergodicity.

5. **Greedy and large n seeding**: pushing to n in the thousands required template based greedy construction using the observed linear scaling of block sizes.

6. **Surgical nudge**: at n ≥ 10⁶, targeted local improvement guided by the suffix array structure became the main practical method for making progress at n = 10⁸ and above.

7. **Branch and bound**: the definitive exact solver. Parallel DFS with SAM based upper bound pruning, Numba JIT inner loops, shared incumbent propagation, and adaptive work splitting. Extended proven values from n = 80 to n = 111 and counting.

## Licence

This project is for research and educational purposes. The OEIS reference values are from [A112509](https://oeis.org/A112509) (Martin Fuller).
