"""Large-n lower-bound search for OEIS A112509.

This module is independent from the CBD proxy objective. It optimizes a
*certified lower bound* on the true distinct-substring-value count by counting
only a chosen subset of substrings exactly.

If S is any subset of all contiguous substrings of a bitstring x, then

    |{ value(t) : t in S }| <= a_true(x)

so the subset distinct-value count is always a valid lower bound.

The default subset is:
- all starts for short lengths (1..short_full_len)
- strided starts for longer lengths (short_full_len+1..max_len)

This scales to n=10_000+ while remaining mathematically valid.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json
import math
import os
import random
import time
from typing import Dict, List, Sequence, Tuple


def _load_seed_store(seed_path: str) -> Dict[str, str]:
    if not os.path.exists(seed_path):
        return {}
    with open(seed_path, "r", encoding="utf-8-sig") as f:
        data = json.load(f)
    return {str(k): str(v) for k, v in data.items()}


def _bitstring_to_runs(bits: str) -> List[Tuple[str, int]]:
    if not bits:
        return []
    runs: List[Tuple[str, int]] = []
    ch = bits[0]
    count = 1
    for nxt in bits[1:]:
        if nxt == ch:
            count += 1
        else:
            runs.append((ch, count))
            ch = nxt
            count = 1
    runs.append((ch, count))
    return runs


def _runs_to_bits(runs: Sequence[Tuple[str, int]]) -> str:
    return "".join(ch * count for ch, count in runs if count > 0)


def _extract_prefix_runs(bits: str, max_runs: int, prefix_ratio: float) -> List[Tuple[str, int]]:
    runs = _bitstring_to_runs(bits)
    if not runs:
        return []
    target_len = max(1, int(len(bits) * prefix_ratio))
    chosen: List[Tuple[str, int]] = []
    built = 0
    for ch, count in runs:
        if len(chosen) >= max_runs or built >= target_len:
            break
        chosen.append((ch, count))
        built += count
    return chosen


def _scaled_prefix(
    source_prefix_runs: Sequence[Tuple[str, int]],
    source_n: int,
    target_n: int,
    max_zero_run: int,
) -> List[Tuple[str, int]]:
    scale = target_n / max(1, source_n)
    scaled: List[Tuple[str, int]] = []
    for ch, count in source_prefix_runs:
        if ch == "1":
            new_count = max(1, int(round(count * scale)))
        else:
            # Keep separators short; this matches observed structure in high-n seeds.
            new_count = max(1, min(max_zero_run, count))
        scaled.append((ch, new_count))
    return scaled


def _build_tail(remaining: int, one_prob: float, rng: random.Random) -> str:
    if remaining <= 0:
        return ""

    bits: List[str] = []
    i = 0
    while i < remaining:
        if rng.random() < one_prob:
            bits.append("1")
            i += 1
            continue

        # Occasional short 0-bursts in tail.
        burst = 1 + int(rng.random() < 0.35) + int(rng.random() < 0.1)
        burst = min(burst, remaining - i)
        bits.extend("0" for _ in range(burst))
        i += burst

    # Ensure final bitstring starts with 1 when this is the whole string.
    if bits and bits[0] != "1":
        bits[0] = "1"

    return "".join(bits)


def build_structured_seed(
    n: int,
    reference_seed: str,
    rng: random.Random,
    prefix_runs: int = 28,
    prefix_ratio: float = 0.72,
    max_zero_run: int = 3,
    tail_one_prob: float = 0.90,
) -> str:
    """Construct a large-n candidate with structured prefix and stochastic tail."""
    source_n = len(reference_seed)
    prefix_template = _extract_prefix_runs(reference_seed, max_runs=prefix_runs, prefix_ratio=prefix_ratio)
    if not prefix_template:
        # Fallback: dense-random if no reference structure is available.
        bits = _build_tail(n, one_prob=tail_one_prob, rng=rng)
        if bits and bits[0] != "1":
            bits = "1" + bits[1:]
        return bits

    scaled_prefix_runs = _scaled_prefix(prefix_template, source_n=source_n, target_n=n, max_zero_run=max_zero_run)

    prefix_bits = _runs_to_bits(scaled_prefix_runs)
    if len(prefix_bits) > n:
        prefix_bits = prefix_bits[:n]
    remaining = n - len(prefix_bits)
    tail = _build_tail(remaining, one_prob=tail_one_prob, rng=rng)

    candidate = prefix_bits + tail
    if len(candidate) < n:
        candidate += "1" * (n - len(candidate))
    elif len(candidate) > n:
        candidate = candidate[:n]

    if candidate and candidate[0] != "1":
        candidate = "1" + candidate[1:]
    return candidate


def certified_lower_bound(
    bits: str,
    max_len: int = 62,
    short_full_len: int = 16,
    long_stride: int = 2,
) -> int:
    """Certified lower bound on true distinct substring values.

    Counts exact integer values for a subset of substrings only.
    Since this subset is included in the full substring set, this is always
    a valid lower bound on the true metric.
    """
    n = len(bits)
    if n == 0:
        return 0

    max_len = max(1, min(max_len, n))
    short_full_len = max(1, min(short_full_len, max_len))
    long_stride = max(1, long_stride)

    values = set()

    # Full starts for short lengths.
    for length in range(1, short_full_len + 1):
        limit = n - length + 1
        for start in range(limit):
            v = 0
            for idx in range(start, start + length):
                v = (v << 1) | (1 if bits[idx] == "1" else 0)
            values.add(v)

    # Strided starts for longer lengths.
    for length in range(short_full_len + 1, max_len + 1):
        limit = n - length + 1
        for start in range(0, limit, long_stride):
            v = 0
            for idx in range(start, start + length):
                v = (v << 1) | (1 if bits[idx] == "1" else 0)
            values.add(v)

    return len(values)


def _mutate_tail(
    bits: str,
    prefix_locked: int,
    rng: random.Random,
    flip_prob: float,
) -> str:
    n = len(bits)
    if prefix_locked >= n:
        return bits

    arr = list(bits)
    tail_len = n - prefix_locked

    op = rng.random()

    if op < 0.65:
        # Single-bit flip
        idx = prefix_locked + rng.randrange(tail_len)
        arr[idx] = "0" if arr[idx] == "1" else "1"
    elif op < 0.9:
        # Short burst flip
        start = prefix_locked + rng.randrange(tail_len)
        burst = 1 + int(rng.random() < 0.5) + int(rng.random() < 0.2)
        end = min(n, start + burst)
        for idx in range(start, end):
            arr[idx] = "0" if arr[idx] == "1" else "1"
    else:
        # Density-aware move: set one tail 0->1 and one 1->0 (if possible)
        zeros = [i for i in range(prefix_locked, n) if arr[i] == "0"]
        ones = [i for i in range(prefix_locked, n) if arr[i] == "1"]
        if zeros and ones:
            z = zeros[rng.randrange(len(zeros))]
            o = ones[rng.randrange(len(ones))]
            arr[z] = "1"
            arr[o] = "0"
        else:
            idx = prefix_locked + rng.randrange(tail_len)
            arr[idx] = "0" if arr[idx] == "1" else "1"

    # Mild bias back toward high 1-density in tail.
    if rng.random() < flip_prob:
        idx = prefix_locked + rng.randrange(tail_len)
        arr[idx] = "1"

    if arr[0] != "1":
        arr[0] = "1"
    return "".join(arr)


@dataclass
class LBSearchConfig:
    n: int
    restarts: int = 12
    steps_per_restart: int = 1200
    max_len: int = 62
    short_full_len: int = 16
    long_stride: int = 2
    prefix_runs: int = 28
    prefix_ratio: float = 0.72
    max_zero_run: int = 3
    tail_one_prob: float = 0.90
    mutate_bias_one_prob: float = 0.08
    anneal_temp0: float = 0.8
    anneal_cooling: float = 0.998
    random_seed: int = 0


@dataclass
class LBSearchResult:
    n: int
    best_lower_bound: int
    best_bits: str
    elapsed_seconds: float
    config: Dict[str, object]


def _choose_reference_seed(seed_store: Dict[str, str], n: int) -> str:
    if not seed_store:
        return ""
    available = sorted((int(k), v) for k, v in seed_store.items() if str(k).isdigit())
    if not available:
        return ""
    # Prefer nearest by n, tie-break larger n for richer prefix structure.
    nearest = min(available, key=lambda kv: (abs(kv[0] - n), -kv[0]))
    return nearest[1]


def run_large_n_lower_bound_search(
    config: LBSearchConfig,
    seed_store_path: str,
    progress_every: int = 200,
) -> LBSearchResult:
    """Run restart-based local search for a certified lower bound."""
    rng = random.Random(config.random_seed)
    seed_store = _load_seed_store(seed_store_path)
    reference = _choose_reference_seed(seed_store, config.n)

    t0 = time.time()
    global_best_score = -1
    global_best_bits = ""

    for restart in range(config.restarts):
        seed = build_structured_seed(
            n=config.n,
            reference_seed=reference,
            rng=rng,
            prefix_runs=config.prefix_runs,
            prefix_ratio=config.prefix_ratio,
            max_zero_run=config.max_zero_run,
            tail_one_prob=config.tail_one_prob,
        )

        prefix_locked = min(len(seed), int(config.n * config.prefix_ratio))

        current = seed
        current_score = certified_lower_bound(
            current,
            max_len=config.max_len,
            short_full_len=config.short_full_len,
            long_stride=config.long_stride,
        )

        best_local = current
        best_local_score = current_score

        temp = config.anneal_temp0

        for step in range(1, config.steps_per_restart + 1):
            proposal = _mutate_tail(
                current,
                prefix_locked=prefix_locked,
                rng=rng,
                flip_prob=config.mutate_bias_one_prob,
            )
            proposal_score = certified_lower_bound(
                proposal,
                max_len=config.max_len,
                short_full_len=config.short_full_len,
                long_stride=config.long_stride,
            )

            delta = proposal_score - current_score
            accept = delta >= 0
            if not accept and temp > 1e-12:
                threshold = math.exp(delta / max(temp, 1e-12))
                accept = rng.random() < threshold

            if accept:
                current = proposal
                current_score = proposal_score

            if current_score > best_local_score:
                best_local_score = current_score
                best_local = current

            temp *= config.anneal_cooling

            if progress_every > 0 and step % progress_every == 0:
                print(
                    f"[restart {restart + 1}/{config.restarts}] step={step} "
                    f"local_best={best_local_score} global_best={global_best_score}",
                    flush=True,
                )

        if best_local_score > global_best_score:
            global_best_score = best_local_score
            global_best_bits = best_local

        print(
            f"[restart {restart + 1}/{config.restarts}] done "
            f"best_local={best_local_score} global_best={global_best_score}",
            flush=True,
        )

    elapsed = time.time() - t0
    return LBSearchResult(
        n=config.n,
        best_lower_bound=global_best_score,
        best_bits=global_best_bits,
        elapsed_seconds=elapsed,
        config=asdict(config),
    )


def save_lb_result(result: LBSearchResult, output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    payload = {
        "n": result.n,
        "best_lower_bound": result.best_lower_bound,
        "best_bits": result.best_bits,
        "elapsed_seconds": result.elapsed_seconds,
        "ones_density": (result.best_bits.count("1") / result.n) if result.n else 0.0,
        "config": result.config,
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
