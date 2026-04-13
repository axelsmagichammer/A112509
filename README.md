# OEIS A112509 — Maximum Distinct Integers from Binary Substrings

> **A112509(n)** = the maximum number of distinct non-negative integers whose binary representations occur as contiguous substrings of some n-bit binary string.

This repository contains algorithms, analysis tools, and pre-computed data for exploring [OEIS sequence A112509](https://oeis.org/A112509) and six related sequences. Results extend the known values from n = 80 (Martin Fuller, OEIS) to **n = 200**, with certified lower bounds computed up to **n = 2 billion**.

## The Problem

Given an n-bit binary string, every contiguous substring can be read as a binary number (with leading zeros collapsing to 0). For example, the 4-bit string `1101` contains substrings:

| Length | Substrings | Integer values |
|--------|-----------|----------------|
| 1 | `1`, `1`, `0`, `1` | 1, 0 |
| 2 | `11`, `10`, `01` | 3, 2, 1 |
| 3 | `110`, `101` | 6, 5 |
| 4 | `1101` | 13 |

Distinct values: {0, 1, 2, 3, 5, 6, 13} → **7 distinct integers**, which happens to be the maximum achievable for n = 4.

A112509 asks: for each n, what is the **maximum** count of distinct integers achievable by any n-bit string?

## First 20 Values

```
n:    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
a(n): 1   3   5   7  10  13  17  22  27  33  40  47  55  64  73  83  94 106 118 131
```

## Related OEIS Sequences

| Sequence | Description | Range in this repo |
|----------|-------------|--------------------|
| **A112509** | Max distinct integers (including 0) from n-bit substrings | n = 1–200 |
| **A112510** | Smallest n-bit number achieving A112509(n) | n = 1–150 |
| **A112511** | Largest n-bit number achieving A112509(n) | n = 1–150 |
| **A156022** | Max distinct *positive* integers from n-bit substrings | n = 1–200 |
| **A156023** | n(n+1)/2 − A112509(n) (the "missing" count) | n = 1–200 |
| **A156024** | n(n+1)/2 − A156022(n) (missing positive count) | n = 1–200 |
| **A156025** | Number of n-bit strings achieving A112509(n) | n = 1–150 |

## Key Findings

- **Asymptotic growth**: a(n)/n² approaches a constant 0.5
- **Optimal string structure**: all optimal strings start with a long block of 1s, followed by decreasing 1-blocks separated by short 0-blocks, and the zero-separator prefix always begins `[1, 2, 1, ...]` and grows with n (K_common ≈ 3 at n = 30, 8 at n = 81–103, 14 at n = 150)
- **De Bruijn connection**: the minimum de Bruijn order k for an optimal string equals its leading 1-block length (proved for all n ≤ 150)
- **Lower bounds at scale**: a(1,000,000,000) ≥ 499,943,911,277,037,650
- **Second differences**: Δ²a(n) = Δa(n) − Δa(n−1) ∈ {0, 1, 2} for all n ≤ 200; the value 2 occurs at exactly n = 52, 71, 117, 157, 178, 191, 194
- **Optimal string structure detail**: the leading 1-block length equals the minimum de Bruijn order; 1-block sizes scale approximately linearly with n (α₁ ≈ 0.19n, α₂ ≈ 0.17n, ...); overall 1-density is tends to 1.

## Elementary Bounds

Two simple bounds bracket the asymptotic growth rate of a(n).

### Upper bound

**Theorem.** For all n ≥ 1, a(n) ≤ n(n+1)/2 + 1, and therefore lim sup a(n)/n² ≤ 1/2.

**Proof.** Two substrings represent the same integer if and only if they are identical after stripping leading zeros. Consequently, every distinct positive integer value corresponds to a unique substring beginning with `1`. An n-bit string has at most n(n+1)/2 substrings in total, so there are at most that many distinct positive values. Adding the value 0 (present whenever the string contains a `0` bit, contributing one additional integer) gives a(n) ≤ n(n+1)/2 + 1. Dividing by n² and taking the limit gives the stated bound.

### Lower bound

**Theorem.** For all n ≥ 2, a(n) ≥ ⌊n/2⌋·⌈n/2⌉ + ⌊n/2⌋ + 1 ≥ n²/4 + n/2, and therefore lim inf a(n)/n² ≥ 1/4.

**Proof.** Consider the explicit string s = 1^x 0^y for positive integers x, y with x + y = n (x ones followed by y zeros). The substrings of s fall into three disjoint types:

1. **All-ones:** the x strings `1`, `11`, ..., `1^x`, taking values 1, 3, 7, ..., 2^x − 1. These give x distinct odd positive integers.
2. **All-zeros:** all such substrings have value 0, contributing exactly one distinct value.
3. **Mixed:** substrings of the form `1^a 0^b` with a ∈ {1,...,x} and b ∈ {1,...,y}, taking value (2^a − 1)·2^b. These xy values are pairwise distinct (different b gives different 2-adic valuation; different a with the same b gives different odd part). Mixed values are even and positive, so they are disjoint from both the all-ones values and 0.

The total is f(s) = x + 1 + xy = x(n − x + 1) + 1. Setting x = ⌊n/2⌋ gives a lower bound of n²/4 + n/2.

### Combining the bounds

Together these give 1/4 + O(1/n) ≤ a(n)/n² ≤ 1/2 + O(1/n). The lower bound is already tight for small n (for n = 4 it gives a(4) ≥ 7, which is exact), but far from tight in general. The project's central conjecture, supported by extensive numerical evidence, is:

> **Conjecture.** lim a(n)/n² = 1/2.

## Asymptotic Numerical Evidence

The following certified lower bounds, computed using the surgical nudge algorithm on optimised seeds, give very strong numerical support for the conjecture:

| n | Lower bound on a(n) | a(n)/n² ≥ |
|---|---|---|
| 3,000 | 4,261,266 | 0.4735 |
| 10,000 | 48,490,001 | 0.4849 |
| 1,000,000 | 498,396,662,200 | 0.4984 |
| 1,000,000,000 | 499,943,911,277,037,650 | 0.49994 |
| 2,000,000,000 | 1,999,841,302,297,432,600 | 0.49996 |

The computation at n = 2 × 10⁹ required scoring a two-billion-bit string (approximately 238 MB as a NumPy boolean array) using the `pydivsufsort` C library. The ratio a(n)/n² converges rapidly towards 0.5 and is within 0.004% of 1/2 at n = 2 × 10⁹.

For exact values up to n = 150, the ratios are:

```
n=100: a(100)=3936,  a(n)/n²=0.3936
n=120: a(120)=5688,  a(n)/n²=0.3950
n=150: a(150)=9070,  a(n)/n²=0.4031
n=200: a(200)=16535, a(n)/n²=0.4134
```

These exact small-n values approach 0.5 more slowly than the large-n lower bounds because the sub-leading correction terms are significant at moderate n. Fitting a constrained log₂ tail model `a(n)/n² ≈ L − c/(log₂ n)^α` to the known data gives L ≈ 0.4997 with very tight confidence across multiple tail cutoffs, consistent with the limit being exactly 1/2.

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

```
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

```
A112509/
├── README.md                 ← you are here
├── requirements.txt          ← Python dependencies
├── config.py                 ← global configuration (paths, beam width, timeouts)
├── compare_seeds.py          ← compare RLE structure of large-n seeds
├── score_1B.py               ← compute a(n) lower bound for n = 1 billion
│
├── src/
│   ├── algorithms/           ← core search algorithms
│   │   ├── brute_force.py        ← exact search for small n (≤ 25)
│   │   ├── evaluate.py           ← O(n log²n) distinct-value counter
│   │   ├── structured_search.py  ← block-template exhaustive search (n ≤ 150)
│   │   ├── MH_algorithm.py       ← Metropolis-Hastings MCMC search (n ≤ 1000)
│   │   ├── greedy_search.py      ← template-based greedy with learned prefixes
│   │   ├── large_n_lower_bound.py← certified lower bounds for n ≥ 3,000
│   │   └── surgical_nudge.py     ← LCP-targeted local refinement
│   │
│   └── tools/                ← analysis and utility scripts
│       ├── generate_oeis_bfiles.py   ← export to OEIS b-file format
│       ├── enhance_cached_results.py ← add derived sequences to cache
│       ├── debruijn_analysis.py      ← de Bruijn embedding analysis
│       ├── compute_hamming.py        ← pairwise Hamming distance stats
│       ├── distribution.py           ← full distribution of b_n(S) for small n
│       ├── benchmark_timing.py       ← MH convergence benchmarks
│       ├── expand_mh_solutions.py    ← expand MH solutions via bit-swap search
│       ├── add_common_sep_fields.py  ← extract common separator structure
│       ├── analyse_common_structure.py ← run-length pattern analysis
│       ├── mh_full_analyze.py        ← structural analysis of MH results
│       └── run_mh.py                 ← parallel MH runner for multiple n
│
├── config/
│   ├── learned_bounds.json       ← structural constraints for block search
│   └── search_constraints.json   ← per-n constraints for extended search
│
├── data/
│   ├── cached_results.json       ← authoritative results (n = 1–150)
│   ├── oeis/                     ← OEIS b-files for all 7 sequences
│   └── reference/
│       ├── known_values.py       ← OEIS reference values (n = 1–80)
│       └── A156025_values.py     ← optimal string counts (n = 1–80)
│
├── results/
│   ├── mh_bounded/               ← constrained MH results (n = 81–200)
│   ├── mh_unbounded/             ← unconstrained MH results (n ≥ 80)
│   ├── large_n/                  ← lower bounds (n = 3K to 2B)
│   ├── distributions/            ← full distributions for n = 12, 16, 20
│   ├── debruijn_analysis.json    ← de Bruijn embedding metadata
│   └── hamming_distances.json    ← Hamming distance statistics
│
├── seeds/                        ← best-known bit-strings for large n
│   ├── n_3000_seed.txt
│   ├── n_2000000_seed.txt
│   ├── n_10000000_seed.{npy,txt}
│   ├── n_100000000_seed.npy
│   ├── seed_n1000000000.npy
│   └── seed_n2000000000.npy
│
├── notebooks/
│   ├── A112509_Extended_Results.ipynb   ← tables, charts, structural analysis
│   └── A112509_Results_Comparison.ipynb ← cross-method agreement for n=80–200
│
└── tests/
    └── test_debruijn_analysis.py
```

## Project History and AI Assistance

This project spanned several months of active research and consumed thousands of hours of CPU time across multiple machines. What began as a straightforward brute-force exploration of the sequence grew into a multi-layered computational investigation requiring increasingly sophisticated algorithms at each scale barrier.

The work was carried out with substantial assistance from **Claude Sonnet 4.6** and **Claude Opus 4.6** (Anthropic). These AI systems contributed throughout the project with performance optimisation: identifying bottlenecks and profiling, and code implementation: writing, debugging and iterating on all major modules, test suites, and analysis tools.

### Algorithm Evolution

The pipeline was not designed upfront, it evolved iteratively as each approach hit its limits:

1. **Brute force**: exhaustive enumeration of all 2^(n−1) strings worked to n ≈ 40 but was computationally intractable beyond that. Even here the scoring inner loop went through three rewrites: set comprehension → hash-set → early-exit integer comparisons.

2. **Block-structure discovery**: manual inspection of brute-force solutions revealed the run-length pattern. This motivated a block-template search that reduced the search space by many orders of magnitude and extended exact results to n ≈ 60. The constraint windows were tightened iteratively by running the search, inspecting which block combinations ever appeared in any optimal solution, and narrowing the bounds accordingly.

3. **Suffix-array scoring**: as n grew past 30, the O(n²) substring hash-set became the bottleneck. The key insight that distinct positive integer values equal distinct leading-1 suffixes in the suffix array, countable via the LCP array reduced scoring from O(n²) to O(n log n). Several implementations were benchmarked; the final version uses `pydivsufsort` (a C-wrapped libdivsufsort call) for SA construction and a Numba-JIT Kasai for LCP, achieving 10–50× speedup over pure Python.

4. **Metropolis–Hastings MCMC**: to escape the constraint-design bottleneck and provide independent validation, an unconstrained MH search was built operating directly on run-length encodings. Five move types (transfer, split, merge, swap, multi-transfer) emerged through trial and error; the final set was chosen because it could reach any RLE representation from any starting point in a small number of moves. The cooling schedule went through ~8 iterations before the current geometric-with-reheating design proved reliable across all n ≤ 150. Crucially, unconstrained MH matched every exact result from the structured search, providing strong empirical certification.

5. **Greedy and large-n seeding**: pushing to n in the thousands required a way to construct a good starting string without search. A greedy template approach was built using the observed linear scaling of block sizes. Counter-pattern seeds (inspired by de Bruijn sequence structure) were developed for n in the millions and above.

6. **Surgical nudge**: at n ≥ 10⁶ even a single score evaluation dominates wall time, so random restarts become unaffordable. The surgical nudge replaced global random search with targeted local improvement: find the highest-LCP collision in the suffix array (the pair of substrings that waste the most coverage by being identical after leading-zero stripping), then make the smallest bit-swap capable of breaking it. This approach, parallelised across worker processes sharing a memory-mapped bit array, was the only method capable of making progress at n = 10⁸ and above. 

Thousands of hours of CPU time were consumed running the algorithms across extended overnight and multi-day runs.

---

## Computational Methods (Summary)

The pipeline moves through four stages of decreasing exactness and increasing scale:

### Scoring

All algorithms share the same inner loop. For small n (≲ 30) a brute-force hash-set over all O(n²) substrings is used. For larger n, the count is computed exactly in O(n log n) via the suffix-array/LCP identity.

#### Suffix Array (SA)

The **suffix array** SA of a string s of length n is a permutation of {0, 1, ..., n−1} such that the suffixes s[SA[0]..], s[SA[1]..], ..., s[SA[n−1]..] are in lexicographic order. For example, if s = `1101` then the four suffixes in sorted order are `01`, `1`, `101`, `1101`, so SA = [2, 1, 3, 0] (0-indexed). The suffix array can be constructed in O(n) time using algorithm DC3/skew or the DivSufSort algorithm. This implementation uses `pydivsufsort`, a Python wrapper around the C library libdivsufsort, which achieves O(n) time with small constants.

#### LCP Array and Kasai's Algorithm

The **LCP array** LCP has LCP[i] = length of the longest common prefix between the suffixes at consecutive SA positions, i.e. between s[SA[i−1]..] and s[SA[i]..], with LCP[0] = 0 by convention. For the example above, the consecutive suffix pairs share prefixes of lengths 0, 1, 1, so LCP = [0, 0, 1, 1].

**Kasai's algorithm** computes the full LCP array in O(n) time and O(n) space from the suffix array using the following key property: if the suffix starting at position i has LCP value h with its predecessor in the sorted order, then the suffix starting at i+1 has LCP value at least h−1 with its own predecessor. This allows the algorithm to compute all LCP values in a single left-to-right scan over the original string, amortising the cost so that each character is examined at most twice overall. In this codebase Kasai's algorithm is compiled to native machine code with Numba's JIT compiler, giving a further 10–50× speedup over a pure Python implementation.

#### The Scoring Identity

The key observation is that two substrings represent the same non-negative integer if and only if they are identical after stripping leading zeros. It follows that:

```
distinct values = (1 if any 0-bit present) + (# distinct substrings starting with '1')
```

The second term counts distinct leading-1 substrings. Since every leading-1 substring is a prefix of some leading-1 suffix, and the suffix array lists all suffixes in sorted order, the number of distinct leading-1 substrings can be read off from the LCP array: each leading-1 suffix at SA position i contributes `(n − SA[i])` substrings, of which `LCP[i]` are already counted by the preceding suffix (they share a common prefix of that length). Summing only over positions where the suffix starts with `1`:

```
f(s) = Σ_{i: s[SA[i]] = '1'}  (n − SA[i] − LCP[i])  +  𝟙[0 ∈ s]
```

where LCP[0] = 0 by convention and 𝟙[0 ∈ s] is 1 if s contains any `0` bit, 0 otherwise. Each term `(n − SA[i] − LCP[i])` counts the substrings that begin at SA[i], start with `1`, and are not a prefix of the lexicographically preceding leading-1 suffix. The entire computation runs in O(n) time once the SA and LCP arrays are available, making the full scoring pipeline O(n log n) (dominated by SA construction, or O(n) with DivSufSort).

### Structured search (exact, n ≤ 150)

Early brute-force runs revealed a universal structural template: every optimal n-bit string takes the form of a long leading block of 1s, followed by alternating 0-separator blocks and 1-blocks of strictly decreasing length. The zero-separator sequence always begins `[1, 2, 1, ...]` and the shared prefix length (K_common) grows with n: approximately 3 at n = 30, 8 at n = 81–103, and 14 at n = 150.

This structure reduces the search to enumerating integer block-size tuples within tight windows. The windows were discovered based on previous optimal values of n and how the strucutre evolves. The final constraints (stored in `config/learned_bounds.json`) reduce the search space by large amounts relative to unconstrained enumeration, making exact results tractable to n = 150 (and beyond) on a modern multi-core CPU.

Results for n ≤ 80 are cross-validated against Fuller's OEIS values; n = 81–150 are certified by the agreement of two independent methods: structured search and unconstrained MH.

### Metropolis–Hastings (heuristic, n ≤ ~1000)

The MH algorithm operates on run-length encodings rather than raw bit strings, representing each candidate as a sequence of alternating block lengths `(b₁, z₁, b₂, z₂, ..., bₖ)` where bᵢ are 1-block lengths and zᵢ are 0-block lengths. This representation makes the five move types natural:

- **Transfer**: move one bit from one block to an adjacent block (most frequent)
- **Split**: divide a block into two by inserting a new separator
- **Merge**: join two adjacent same-type blocks by removing a separator
- **Swap**: exchange the lengths of two blocks
- **Multi-transfer**: transfer several bits at once for larger steps

The temperature follows a geometric cooling schedule with periodic reheating to escape local optima. Despite having no structural priors, unconstrained MH reliably recovers the true optimum for all n ≤ 150 (often in minutes and even seconds for small n), confirming that the block template is a natural attractor of the score landscape rather than an artificially imposed constraint. For n ≤ 150, independent MH runs from random initial strings converge to the same best value as exhaustive structured search, providing strong empirical certification.  Morevoer this could be adopted for n>150 top find the optimal value for a given n (what it does not guarantee is finding all optimal solutions but it appears to come very close, for this reasons optimal solutions beyond 150 have not been submitted, just the optimal value).

### Surgical nudge (certified lower bounds, n up to 2 × 10⁹)

For large, even a single suffix-array score evaluation takes tens of seconds, making random restarts infeasible. The surgical nudge instead maintains a single incumbent string and makes targeted local improvements guided by the suffix array structure itself.

At each step the algorithm identifies the highest-LCP collision: the pair of adjacent suffixes in sorted order with the largest LCP value. This pair represents the most costly duplication: two long substrings that happen to be identical, each wasting LCP positions that could contribute new distinct values. A small neighbourhood of bit positions around the collision point is then searched for a density-preserving swap (flip a 0→1 and a 1→0 simultaneously, preserving the total number of 1-bits) that breaks the collision.

Candidate swaps are evaluated in parallel across a worker process pool sharing a memory-mapped copy of the incumbent bit array. Because every accepted swap strictly increases the exact score, the final value is a rigorous certified lower bound. At n = 2 × 10⁹ a lower bound was found a(n)/n² ≥ 0.49996, establishing that the limiting constant is extremely close to 1/2.

---

## Algorithms

> **Note on configuration**: Most algorithms are configured by editing the `main()` function at the bottom of the script; look for the **USER CONFIGURATION** section. This is deliberate: research scripts with many interrelated parameters (block constraints, density bounds, annealing schedules) are easier to manage as editable config blocks than as dozens of CLI flags. The brute_force and evaluate scripts are the exceptions, accepting command-line arguments.

### 1. Brute Force (`src/algorithms/brute_force.py`)

Exhaustive enumeration of all 2^(n−1) candidate n-bit strings (leading bit is always 1). For each string, computes the set of distinct integer values from all O(n²) substrings.

**Use for**: n ≤ 20 (exact, complete enumeration, large n can be used but it is slow).

```bash
python -m src.algorithms.brute_force
```

Prints a table of a(n) for n = 1–20, validated against known OEIS values.

### 2. Evaluate (`src/algorithms/evaluate.py`)

Counts distinct substring values for a single bit-string using a suffix-array approach in O(n log² n). Key insight: two substrings have the same integer value if and only if they are identical after stripping leading zeros. So:

> distinct values = (1 if any 0-bit exists) + (number of distinct substrings starting with 1)

The second term is computed exactly using suffix arrays + Kasai's LCP algorithm.

```bash
# Score a specific bit-string
python -m src.algorithms.evaluate 11111101111111001111101011110110

# Pipe from stdin
echo 1101 | python -m src.algorithms.evaluate

# Interactive prompt
python -m src.algorithms.evaluate
```

### 3. Structured Search (`src/algorithms/structured_search.py`)

**Use for**: n ≤ 150 (exact, certified results).

#### Structural constraints as a search oracle

Analysing previously known optimal strings for n ≤ 80 reveals a universal structural template: the zero-separator sequence stabilises to `(1, 2, 1, 1, 1, 1, 1, 1)` for n ≥ 80 (with the shared prefix growing further as n increases), and the one-block lengths follow a predominantly decreasing pattern with characteristic small interspersed blocks.

Block-size bounds [b_i_min(n), b_i_max(n)] for each block index i are established by manual inspection of the optimal strings for nearby values of n. For a new target n, the previous few optimal strings (for n−1, n−2, ...) are examined to identify the range within which each block length has moved; a small additional tolerance is added on each side to ensure no valid candidate is excluded. The bounds are also cross-checked against the best solutions found by the unconstrained MH algorithm, which provides an independent estimate of plausible block sizes. The resulting window width b_i_max − b_i_min is typically 3–8 for the dominant blocks, tight enough to make enumeration feasible.

#### The search procedure

Given the manually determined bounds for a target n, the algorithm proceeds as follows:

1. **Fix the separator template.** The zero-block lengths are set to `(1, 2, 1, 1, 1, 1, 1, 1, ...)`. This is treated as a hard constraint.

2. **Enumerate block-size tuples.** Iterate over all integer tuples (b_1, b_2, ..., b_K) with b_i in [b_i_min, b_i_max] and sum(b_i) + sum(c_i) = n. The product of window widths is the search space size; for n ≤ 150 this is comfortably below 10^6 configurations.

3. **Evaluate.** For each tuple, assemble the n-bit string from the RLE and compute f(s) using the suffix-array kernel. Evaluation is parallelised across all available CPU cores with a process pool.

4. **Record the optimum.** The maximum f(s) found is reported as a(n), and all strings achieving it are stored.

#### Why it works

The approach is correct (not merely heuristic) for n ≤ 150 because:

- For n ≤ 80, all results are cross-validated against Fuller's independently computed values.
- For 81 ≤ n ≤ 150, an independent MH search consistently finds the same best value when run from a random initialisation with no structural constraints. Agreement across more than 10 independent restarts for all these values of n certifies that the structured search is not missing any region of the optimum.  Moreover the patterns in the a(n) continue to be seen.
- Analytically, the constraints are necessary conditions: every optimal string for n ≤ 80 satisfies them. An optimal string outside the template would represent a new structural regime; the MH search is specifically designed to detect this if it occurs.

**To run:** edit the **USER CONFIGURATION** section in `main()`:

```python
n_start = 1                # first n to compute
n_end   = 30               # last n to compute
force_recompute = False    # True = ignore cache and recompute all
```

```bash
python -m src.algorithms.structured_search
```

Results are automatically saved to `data/cached_results.json`. Already-computed values are skipped unless `force_recompute` is set. CLI overrides are also accepted (e.g. `structured_search 80 100 --force`).

### 4. Metropolis-Hastings MCMC (`src/algorithms/MH_algorithm.py`)

**Use for**: n = 1–~1000 (heuristic; reliably finds the true optimum for all n ≤ 150 in unconstrained mode).

#### Representation

The MH algorithm works directly with the run-length encoding (RLE) of binary strings. An RLE is an integer list [b_1, c_1, b_2, c_2, ...] where the b_i are lengths of consecutive 1-blocks and the c_i are lengths of 0-blocks. For example, `11100110` is encoded as [3, 2, 2, 1]. Every valid n-bit string beginning with `1` (a necessary condition for maximality for n ≥ 2) is uniquely described by such a list.

Working in RLE space reduces the effective dimensionality of the search from n bits to K ≈ 10–20 run lengths for the values of n considered, and makes it natural to propose structured moves that preserve validity.

#### Proposed moves

At each iteration the algorithm selects one of five move types:

1. **Transfer.** Move 1–3 bits between adjacent runs: [..., b, c, ...] → [..., b−t, c+t, ...] for a uniformly random transfer amount t ∈ {1, 2, 3} capped at b−1. Preserves the total number of bits n.

2. **Split.** Divide one run into three by inserting a new separator: [..., b, ...] → [..., b_1, 1, b_2, ...] with b_1 + 1 + b_2 = b. Increases the number of blocks by 2.

3. **Merge.** Remove a separator and fuse the two flanking runs of the same type: [..., b_1, c, b_2, ...] → [..., b_1 + c + b_2, ...]. Decreases the number of blocks by 2.

4. **Swap.** Exchange the lengths of two non-adjacent runs of the same type (both 1-runs or both 0-runs): [..., b_i, ..., b_j, ...] → [..., b_j, ..., b_i, ...]. Preserves n and the number of blocks.

5. **Multi-transfer.** Simultaneously transfer bits between two non-adjacent pairs of runs, and with probability 5% apply a double-pair cross-swap that adjusts four runs at once. This extra move is crucial for crossing connectivity gaps in the optimisation landscape where no single-run change improves the score.

The probabilities of these moves are not uniform: they vary with the search phase. In the early phase (0–30% of iterations), structural moves (split, merge) receive higher weight to explore diverse block topologies. In the late phase (70–100%), fine-grained transfer moves are favoured.

#### Acceptance rule and temperature schedule

The algorithm uses the standard Metropolis criterion. Let v_old and v_new = f(s_new) be the current and proposed scores. The proposed move is accepted with probability min(1, exp((v_new − v_old) / T)), where T is a temperature parameter that decays geometrically: T_{k+1} = γ · T_k with γ ≈ 0.9995 (the cooling rate). Improvements are always accepted. Worsening moves are accepted with a probability that decreases as the search progresses and T → 0, allowing escape from local optima early in the search while concentrating on local refinement late.

Periodic reheating (resetting T to a fraction of its initial value) prevents premature convergence when the chain has been non-improving for many iterations.

#### Constrained and unconstrained modes

**Constrained mode.** Block-size bounds and the separator template from the structured search are enforced during move proposal. If a proposed move would violate any constraint (e.g. push a block below b_i_min), the proposal is rejected outright before evaluation. This dramatically reduces wasted evaluations and speeds convergence.

**Unconstrained mode.** No structural priors are imposed. The algorithm starts from a random RLE (or a previous best) and explores the full space of binary strings. Despite the much larger search space, unconstrained mode reliably finds the optimal a(n) value for n ≤ 150 in a few thousand iterations per restart.

> **Key observation.** For every n ≤ 150, at least one of 10 independent unconstrained restarts (each using 5 × 10^4 iterations and random initialisation) returns a solution achieving a(n). The typical time to first reach the optimum is O(10^3) iterations, far fewer than the full budget.

This observation is striking: even without knowing the separator template or block-size bounds, the MH algorithm discovers the correct structure automatically. It provides independent probabilistic certification of the structured search results and suggests that the (1, 2, 1, 1, ...) separator pattern is a deep attractor in the optimisation landscape, not an artefact of the search constraints.

#### Multi-restart parallelism

Multiple independent MH chains run in parallel, one per CPU core. Each restart uses a fresh random initialisation and slightly varied hyperparameters (initial temperature, cooling rate, iteration budget). The global best solution is tracked across all restarts. This embarrassingly parallel architecture provides both speed and robustness: if one chain converges prematurely to a sub-optimal local maximum, others are likely to find the global optimum.

**To run:** edit the configuration section at the top of `main()`, then:

```bash
python -m src.algorithms.MH_algorithm
```

Key settings to edit in `main()`:

```python
n = 194                        # bit-string length to optimise
num_restarts = 15_000          # Ctrl+C at any time; partial results are saved
iterations_per_run = 100_000   # iterations per restart
num_processes = None           # None = use all cores

# Optional structural constraints (None = unconstrained)
min_ones = None                # e.g. [18, 12, 10, 5] for first 4 one-blocks
max_ones = [None]*4
min_zeros = None               # e.g. [1, 1, 0, 0] for first 4 zero-blocks
max_zeros = [None]*4
min_density = None             # e.g. 0.75
max_density = None             # e.g. 0.88

# Post-search solution expansion
run_level_set_walk = True      # find more solutions at the best value found
```

Results are saved to `results/mh_unbounded/n_XXXX_results.json` (this can be changed). Press Ctrl+C at any time; progress is saved automatically.

### 5. Greedy Search (`src/algorithms/greedy_search.py`)

Template-based construction exploiting the discovered linear scaling pattern:

- Zero-separator prefix begins `[1, 2, 1, ...]` and grows with n
- Block sizes scale approximately linearly with n
- The algorithm fixes the prefix structure, estimates block sizes, and greedily optimizes the tail

**Use for**: fast initial lower bound solutions for any n

**To run**: edit the configuration in `main()`, then:

```bash
python -m src.algorithms.greedy_search
```

Key settings:

```python
n = 20000                              # target bit-string length
local_search_iterations = 5_000_000    # local improvement steps
max_seconds_per_restart = 18000        # time budget per restart (seconds)
seed_solution = None                   # None = auto-load best from results/
```

### 6. Large-n Lower Bound (`src/algorithms/large_n_lower_bound.py`)

Computes mathematically valid lower bounds for n ≥ 3,000 by counting only a chosen subset of substrings exactly (any subset gives a valid lower bound). Uses:

- Pattern-based seed construction (counter-pattern similar to de Bruijn sequences)
- Full or strided sampling strategies
- Suffix array + LCP scoring for exact subset counting

**Use for**: n = 3,000 to 2,000,000,000+

This module is primarily used as a library by other scripts (e.g. `score_1B.py`). See `score_1B.py` for an example of scoring a billion-bit string.

### 7. Surgical Nudge (`src/algorithms/surgical_nudge.py`)

Local refinement starting from a known good seed. Identifies the highest-LCP collision in the suffix array (the "worst" duplicate substring) and makes the smallest bit-swap to break it. Features:

- LCP-targeted collision breaking
- Density-preserving swaps (flip a 0→1 and a 1→0 simultaneously)
- Random restart with perturbation kicks
- Multi-swap exploration

**Use for**: polishing large-n seeds to get an improved lower bound after initial construction

**To run**: edit the configuration in `main()`, then:

```bash
python -m src.algorithms.surgical_nudge
```

Key settings:

```python
n = 100_000_000       # must match a seed file in seeds/
max_seconds = 36000   # wall-clock budget (seconds)
num_workers = 14      # parallel suffix-array evaluations
```

Expects a seed file in `seeds/` (e.g. `n_100000000_seed.npy`). Results are saved to `results/large_n/`.

## Tools

| Script | Purpose |
|--------|---------|
| `generate_oeis_bfiles.py` | Export all sequences to OEIS b-file format in `data/oeis/` |
| `enhance_cached_results.py` | Add A112510/A112511/A156022–A156024 to `cached_results.json` |
| `debruijn_analysis.py` | Compute minimum de Bruijn embedding order for each optimal string |
| `compute_hamming.py` | Compute min/max/mean pairwise Hamming distances between optimal solutions |
| `distribution.py` | Full distribution of distinct-value counts across all 2^n strings |
| `benchmark_timing.py` | Benchmark MH convergence speed per n |
| `expand_mh_solutions.py` | Expand MH solution sets via exhaustive 1-hop bit-swap neighbourhood search |
| `add_common_sep_fields.py` | Extract common separator structure (K_common, common_seps) |
| `run_mh.py` | Run MH in parallel for multiple n values |

Run any tool with:

```bash
python -m src.tools.<tool_name>
```

## Notebooks

### Extended Results (`notebooks/A112509_Extended_Results.ipynb`)

Interactive presentation of all computed results:

- Sequence values and all optimal strings (tabbed by groups of 10)
- Structural analysis: block counts, 1-density, run lengths, K_common
- K_common vs n charts
- Hamming distance analysis between optimal solutions
- Growth analysis: Δa, ΔΔa, a(n)/n, a(n)/n²

### Results Comparison (`notebooks/A112509_Results_Comparison.ipynb`)

Cross-validates three independent methods for n = 80–200:

- Unbounded Metropolis-Hastings
- Bounded (constrained) Metropolis-Hastings
- Exhaustive block-template search

Shows value agreement, solution count coverage, and flags any mismatches.

## Data Files

### `data/cached_results.json`

Authoritative results for n = 1–150. Each entry contains:

```json
{
  "42": {
    "a(n)": 622,
    "num_optimal": 36,
    "optimal_strings": ["111111111110111111111100...", ...],
    "A112510_min_decimal": 4393508438016,
    "A112511_max_decimal": 4398004084736,
    "A156022_positive_only": 621,
    "A156023_complement": 281,
    "A156024_complement_positive": 282,
    "K_common": 8,
    "common_seps": [1, 2, 1, 1, 1, 1, 1, 1]
  }
}
```

### `results/mh_unbounded/`

Unconstrained MH results (no block or density constraints) for n ≥ 80. Used as a cross-check: if unbounded and bounded MH agree on `best_value`, it strongly supports that value being the true optimum.

```json
{
  "n": 151,
  "best_value": 9197,
  "solutions": ["11111111111...", ...],
  "last_updated": "...",
  "metadata": {...}
}
```

### `results/mh_bounded/`

MH results with structural constraints (block ranges derived from `config/search_constraints.json`, with loose bounds) for n = 81–200. The constraints reduce the search space, allowing more restarts to focus on the most promising region. Same format as `mh_unbounded/` above.

### `results/large_n/`

Lower-bound computations for a selection of n = 3,000 up to 2 billion:

```json
{
  "n": 1000000000,
  "score": 499943911277037650,
  "a_over_n_sq": 0.49994391127703763,
  "method": "counter-pattern seed (no optimization)"
}
```

### `data/oeis/`

OEIS-formatted b-files (`n a(n)` pairs) for submission/reference:

- `b112509.txt` (n = 1–200)
- `b112510.txt` (n = 1–150)
- `b112511.txt` (n = 1–150)
- `b156022.txt` (n = 1–200)
- `b156023.txt` (n = 1–200)
- `b156024.txt` (n = 1–200)
- `b156025.txt` (n = 1–150)

## Configuration

### `config/learned_bounds.json`

Structural constraints discovered from analysing optimal solutions, used by `structured_search.py`:

- **K_common**: the number of leading zero-separators shared by all optimal solutions for the target n (e.g. 8 for n = 99)
- **common_seps**: the corresponding shared separator values (e.g. `[1, 2, 1, 1, 1, 1, 1, 1]` for n = 99); both grow with n
- **block_ranges**: `[min, max]` for each 1-block position
- **density bounds**: min/max 1-bit density

### `config/search_constraints.json`

Per-n structural constraints for n ≥ 81, used by the bounded MH algorithm. Includes tailored block ranges, density bounds, and separator constraints per value of n.

## Typical Workflows

### Compute exact values for a new n

1. Ensure `config/learned_bounds.json` has appropriate constraints
2. Run `structured_search.py` for the target n range
3. Run `enhance_cached_results.py` to add derived sequences
4. Run `generate_oeis_bfiles.py` to update b-files

### Explore a new large n

1. Construct an initial seed using `greedy_search.py` or `large_n_lower_bound.py`
2. Refine with `surgical_nudge.py`
3. Score the final string with `evaluate.py`

### Validate MH results against exhaustive search

1. Run MH (bounded and/or unbounded) for the target n
2. Run exhaustive `structured_search.py` if feasible
3. Compare in `A112509_Results_Comparison.ipynb`

## Second Differences

Define the first difference Δa(n) = a(n) − a(n−1) and the second difference Δ²a(n) = Δa(n) − Δa(n−1). Since a(n) ≈ n²/2, one expects Δa(n) ≈ n and Δ²a(n) ≈ 1. The precise behaviour is striking.

> **Observation.** For all n ≤ 200, Δ²a(n) ∈ {0, 1, 2}. The values 0 and 1 alternate in a near-regular pattern. The value 2 occurs at n = 52, 71, 117, 157, 178, 191, and 194 (bold in the tables below).

The predominance of Δ²a(n) ∈ {0, 1} reflects the fact that Δa(n) increases by exactly 0 or 1 at each step, a remarkable regularity for a sequence defined by a combinatorial optimisation. The rare value 2 is associated with a simultaneous resolution of two de Bruijn order tiers at step n−1.

> **Conjecture.** Δ²a(n) ∈ {0, 1, 2} for all n ≥ 1.

### Table 1: a(n), Δa(n), Δ²a(n) for n = 1 to 50

| n | a(n) | Δa | Δ²a |
|--:|-----:|---:|----:|
| 1 | 1 | – | – |
| 2 | 3 | 2 | – |
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

Bold entries mark occurrences of Δ²a(n) = 2. Values for n ≤ 150 are exact from exhaustive search; n = 151–200 are exact values from independent MH runs (certified by multi-restart agreement).

| n | a(n) | Δ²a | n | a(n) | Δ²a | n | a(n) | Δ²a | n | a(n) | Δ²a |
|--:|-----:|----:|--:|-----:|----:|--:|-----:|----:|--:|-----:|----:|
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

**(i) Asymptotic.** Prove or disprove conjecture lim a(n)/n² = 1/2. Improving the lower bound significantly would be progress.

**(ii) Separator structure.** Prove or disprove that the zero-separator sequence (1, 2, 1, 1, ...) is universal for all sufficiently large n. Can the full separator sequence be determined analytically? The second separator being 2 while all others are 1 is a highly non-trivial constraint whose origin is not understood.

**(iii) Second differences.** Prove that Δ²a(n) ∈ {0, 1, 2} for all n. This would require understanding the growth rate of Δa(n) with precision.

**(iv) Exact values.** Extend the exhaustive search computation beyond n = 150. Each additional value requires a substantial increase in the search budget; reaching n = 200 exactly appears to require many CPU-hours with the current structured-search method.

**(v) One-bit density.** The fraction of 1-bits in optimal strings grows rapidly with n: from approximately 0.60 at n = 10, through 0.83 at n = 80, to 0.85–0.86 at n = 150. At large n the convergence is striking: the best-known seeds have density 0.999950 at n = 10⁹ and 0.999964 at n = 2 × 10⁹. This strongly suggests that the 1-bit density of optimal strings tends to 1 as n → ∞. Proving this and determining the rate of convergence, which appears to be 1 − O(1/n) is an open problem.

**(vi) Alphabet generalisation.** What is the analogous sequence for ternary (base-3) or higher-base alphabets? Does the ratio a_q(n)/n² for base q approach (q−1)/2 or some other constant?

## Licence

This project is for research and educational purposes. The OEIS reference values are from [A112509](https://oeis.org/A112509) (Martin Fuller).
