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
- **Optimal string structure**: all optimal strings start with a long block of 1s, followed by decreasing 1-blocks separated by short 0-blocks — the first 8 separators are always `[1, 2, 1, 1, 1, 1, 1, 1]` for n ≥ 30
- **De Bruijn connection**: the minimum de Bruijn order k for an optimal string equals its leading 1-block length (proved for all n ≤ 150)
- **Lower bounds at scale**: a(1,000,000,000) ≥ 499,943,911,277,037,650

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

## Algorithms

> **Note on configuration**: Most algorithms are configured by editing the `main()` function at the bottom of the script — look for the **USER CONFIGURATION** section. This is deliberate: research scripts with many interrelated parameters (block constraints, density bounds, annealing schedules) are easier to manage as editable config blocks than as dozens of CLI flags. The brute_force and evaluate scripts are the exceptions, accepting command-line arguments.

### 1. Brute Force (`src/algorithms/brute_force.py`)

Exhaustive enumeration of all 2^(n−1) candidate n-bit strings (leading bit is always 1). For each string, computes the set of distinct integer values from all O(n²) substrings.

**Use for**: n ≤ 20 (exact, complete enumeration)

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

Block-template exhaustive search exploiting discovered structural properties. Optimal strings follow a pattern: a leading 1-block, then alternating 0-separators and decreasing 1-blocks. The search uses pre-learned constraints from `config/learned_bounds.json`:

- **K_common**: number of initial separators shared across all optimal solutions
- **common_seps**: the shared separator values (seem to be at least `[1, 2, 1, 1, 1, 1, 1, 1]` for large enough n)
- **block_ranges**: min/max length for each 1-block
- **density bounds**: min/max fraction of 1-bits

The algorithm fixes the common prefix and enumerates all valid block-size permutations, then exhaustively searches tail combinations.

**Use for**: n ≤ 150 (exact results, near guaranteed optimal), it will work for larger n but slow due to search space.

Edit the **USER CONFIGURATION** section in `main()`:

```python
n_start = 1                # first n to compute
n_end   = 30               # last n to compute
force_recompute = False    # True = ignore cache and recompute all
```

Then run:

```bash
python -m src.algorithms.structured_search
```

Results are automatically saved to `data/cached_results.json`. Already-computed values are skipped unless `force_recompute` is set. CLI overrides are also accepted (e.g. `structured_search 80 100 --force`).

### 4. Metropolis-Hastings MCMC (`src/algorithms/MH_algorithm.py`)

Stochastic search using Markov Chain Monte Carlo with run-length encoding. The algorithm represents strings as alternating blocks of 1s and 0s, and proposes moves:

- **Transfer**: move bits between adjacent blocks
- **Split**: divide one block into two
- **Merge**: combine two adjacent blocks
- **Swap**: exchange positions of blocks

Uses simulated annealing with configurable temperature schedules and supports structural constraints (bounded mode uses `config/search_constraints.json`).

**Use for**: n = 1–~1000 (heuristic, may not find all the optimal solutions but appears to be very efficient and quick at finding the true optimum, even with not using any search constraints)

**To run**: edit the configuration section at the top of `main()`, then:

```bash
python -m src.algorithms.MH_algorithm
```

Key settings to edit in `main()`:

```python
n = 194                        # bit-string length to optimize
num_restarts = 15_000          # Ctrl+C at any time — partial results are saved
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

Results are saved to `results/mh_unbounded/n_XXXX_results.json`. Press Ctrl+C at any time — progress is saved automatically.

### 5. Greedy Search (`src/algorithms/greedy_search.py`)

Template-based construction exploiting the discovered linear scaling pattern:

- First 8 zero-blocks seem to be always `[1, 2, 1, 1, 1, 1, 1, 1]` for large n
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
- Asymptotic extrapolation with log₂ tail models
- De Bruijn embedding analysis with worked example

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

MH results with structural constraints (block ranges from `config/search_constraints.json`, with only loose density bounds) for n = 81–200. The constraints reduce the search space, allowing more restarts to focus on the most promising region. Same format as `mh_unbounded/` above.

### `results/large_n/`

Lower-bound computations for n = 3,000 up to 2 billion:

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

- **K_common = 8**: the first 8 separators are fixed across all optimal solutions
- **common_seps = [1, 2, 1, 1, 1, 1, 1, 1]**: the shared separator values
- **block_ranges**: `[min, max]` for each 1-block position
- **density bounds**: min/max 1-bit density

### `config/search_constraints.json`

Per-n structural constraints for n ≥ 81, used by the bounded MH algorithm. Includes tailored block ranges, density bounds, and separator constraints per value of n.

## Typical Workflows

### Compute exact values for a new n (≤ 150)

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

## Licence

This project is for research and educational purposes. The OEIS reference values are from [A112509](https://oeis.org/A112509) (Martin Fuller).
