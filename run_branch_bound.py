"""
run_branch_bound.py — Entry point for exact branch-and-bound computation of A112509(n).

Proves the exact value of a(n) by exhaustive parallel DFS with SAM-based
upper-bound pruning.  When run with --collect-all (the default), also finds
every optimal n-bit string.

Usage
-----
    # Solve a single n with default settings:
    python run_branch_bound.py --n 112

    # Override worker count and lookahead:
    python run_branch_bound.py --n 100 --workers 8 --lookahead 3

    # Dry run (no result file written):
    python run_branch_bound.py --n 20 --no-save

    # Multi-machine sharding:
    python run_branch_bound.py --n 115 --shard-mode coordinator
    python run_branch_bound.py --n 115 --shard-mode worker      # on each machine
    python run_branch_bound.py --n 115 --shard-mode merge

Results are saved to results/branch_bound_exact/n_NNNN_results.json.
"""

import argparse
import json
import os
import re
import time

from data.reference.known_values import KNOWN_VALUES
from src.algorithms.branch_bound_search import (
    branch_bound_a112509,
    branch_bound_a112509_parallel,
    enumerate_shard_tasks,
    merge_shard_results,
    run_shard_worker_loop,
)


DEFAULT_COLLECT_ALL = True


def build_auto_config(
    n: int,
    collect_all: bool = DEFAULT_COLLECT_ALL,
    preset: str = "stable",
) -> dict[str, int | bool | str]:
    """Return stable defaults intended for hands-off exact runs.

    The point of this runner is to avoid per-n retuning.  Keep the defaults
    conservative and let branch_bound_a112509_parallel auto-select the
    parameters it already knows how to choose internally.
    """
    cfg: dict[str, int | bool | str] = {
        "n": n,
        "collect_all": collect_all,
        "lookahead": 2 if collect_all else 1,
        "num_workers": 0,
        "split_depth": 0,
        "hard_split_gap": 5,
        "hard_split_extra": 7,
        "dfs_lookahead": 0,
        "dfs_refine_lookahead": 0,
        "dfs_refine_margin": 2,
        "probe_nodes": 0,
        "max_nodes_per_task": 0,
        "preset": preset,
    }

    if collect_all and preset == "full-enum-fast":
        cfg["split_depth"] = 0
        cfg["hard_split_gap"] = 6
        cfg["hard_split_extra"] = 9
        cfg["probe_nodes"] = 2500 if n >= 80 else 0
        cfg["max_nodes_per_task"] = max(200_000, 2500 * n) if n >= 80 else 0
        cfg["dfs_lookahead"] = 3 if n >= 90 else 0
        cfg["dfs_refine_lookahead"] = 4 if n >= 100 else 0
        cfg["dfs_refine_margin"] = 2

        if n >= 100:
            cfg["hard_split_gap"] = 5
            cfg["hard_split_extra"] = 7
            cfg["probe_nodes"] = 0
            cfg["max_nodes_per_task"] = 0
            cfg["dfs_lookahead"] = 0
            cfg["dfs_refine_lookahead"] = 0
            cfg["dfs_refine_margin"] = 2

    return cfg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exact branch-and-bound solver for A112509(n).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_branch_bound.py --n 112
  python run_branch_bound.py --n 100 --workers 8 --lookahead 3
  python run_branch_bound.py --n 115 --shard-mode coordinator
  python run_branch_bound.py --n 115 --shard-mode worker
  python run_branch_bound.py --n 115 --shard-mode merge
""",
    )
    parser.add_argument("--n", type=int, required=True,
                        help="bit-string length to solve")
    parser.add_argument(
        "--preset",
        choices=("stable", "full-enum-fast"),
        default="stable",
        help=(
            "runner tuning preset. 'stable' keeps conservative defaults; "
            "'full-enum-fast' enables stronger full-enumeration balancing/pruning"
        ),
    )
    parser.add_argument("--workers", type=int, default=None,
                        help="override worker count (default: all CPUs minus one)")
    parser.add_argument("--lookahead", type=int, default=None,
                        help="override lookahead depth (default: 2 with collect-all)")
    parser.add_argument("--split-depth", type=int, default=None,
                        help="override BFS split depth (0 = auto)")
    parser.add_argument("--hard-split-gap", type=int, default=None,
                        help="hard-subtree gap threshold for extra re-splitting")
    parser.add_argument("--hard-split-extra", type=int, default=None,
                        help="extra BFS levels for hard-subtree re-splitting")
    parser.add_argument("--dfs-lookahead", type=int, default=None,
                        help="DFS worker lookahead (0 = use main lookahead)")
    parser.add_argument("--dfs-refine-lookahead", type=int, default=None,
                        help="DFS adaptive refinement lookahead (0 = auto, <0 = disable)")
    parser.add_argument("--dfs-refine-margin", type=int, default=None,
                        help="DFS refinement trigger margin (<0 = disable)")
    parser.add_argument("--probe-nodes", type=int, default=None,
                        help="shallow probe budget per subtree (0 = off)")
    parser.add_argument("--chunk-nodes", type=int, default=None,
                        help="DFS chunk budget for cooperative chunking (0 = off)")
    parser.add_argument("--no-collect-all", action="store_true",
                        help="find only a(n), not all optimal strings (faster)")
    parser.add_argument("--no-save", action="store_true",
                        help="run search but do not write result file")
    # ── Sharding options ──────────────────────────────────────────────────
    parser.add_argument(
        "--shard-mode",
        choices=("coordinator", "worker", "merge"),
        default=None,
        help=(
            "multi-machine sharding mode. "
            "'coordinator': enumerate task files; "
            "'worker': claim and process tasks; "
            "'merge': combine done-task files into final result"
        ),
    )
    parser.add_argument(
        "--shard-dir",
        default=None,
        help="base directory for shard work files (default: results/shard_work/)",
    )
    return parser.parse_args()


def compute_k_common(strs):
    if not strs:
        return None

    def sep_list(s):
        return [len(m.group()) for m in re.finditer(r"0+", s)]

    sep_lists = [sep_list(s) for s in strs]
    k = 0
    for seps in zip(*sep_lists):
        if len(set(seps)) == 1:
            k += 1
        else:
            break
    return k


def compute_common_seps(strs, k_common):
    if not strs or k_common is None or k_common == 0:
        return []

    def sep_list(s):
        return [len(m.group()) for m in re.finditer(r"0+", s)]

    return sep_list(strs[0])[:k_common]


def _save_result(
    n: int,
    a: int,
    optimals: list,
    stats: dict,
    config: dict,
    known: dict,
    elapsed: float = 0.0,
) -> None:
    out_dir = os.path.join(os.path.dirname(__file__), "results", "branch_bound_exact")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"n_{n:04d}_results.json")
    k_common = compute_k_common(optimals)
    common_seps = compute_common_seps(optimals, k_common)
    cpu_count = os.cpu_count() or 1
    resolved_workers = (
        int(config["num_workers"]) if int(config["num_workers"]) > 0 else max(1, cpu_count - 1)
    )
    entry = {
        "a(n)": a,
        "num_optimal": len(optimals),
        "optimal_strings": list(optimals),
        "K_common": k_common,
        "common_seps": common_seps,
        "meta": {
            "n": n,
            "preset": str(config.get("preset", "stable")),
            "lookahead": int(config.get("lookahead", 2)),
            "collect_all": bool(config.get("collect_all", True)),
            "requested_workers": int(config.get("num_workers", 0)),
            "resolved_workers": resolved_workers,
            "split_depth": int(config.get("split_depth", 0)),
            "hard_split_gap": int(config.get("hard_split_gap", 5)),
            "hard_split_extra": int(config.get("hard_split_extra", 7)),
            "dfs_lookahead": int(config.get("dfs_lookahead", 0)),
            "dfs_refine_lookahead": int(config.get("dfs_refine_lookahead", 0)),
            "dfs_refine_margin": int(config.get("dfs_refine_margin", 2)),
            "probe_nodes": int(config.get("probe_nodes", 0)),
            "max_nodes_per_task": int(config.get("max_nodes_per_task", 0)),
            "nodes_expanded": stats.get("nodes_expanded", 0),
            "nodes_pruned": stats.get("nodes_pruned", 0),
            "num_subtrees": stats.get("num_subtrees", stats.get("num_tasks", 0)),
            "elapsed_s": round(elapsed, 2),
        },
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(entry, f, indent=2)
    print(f"\nSaved to {out_path}")


def main() -> None:
    args = parse_args()
    known = {i + 1: v for i, v in enumerate(KNOWN_VALUES)}
    collect_all = not args.no_collect_all
    config = build_auto_config(args.n, collect_all=collect_all, preset=args.preset)
    if args.workers is not None:
        config["num_workers"] = args.workers
    if args.lookahead is not None:
        config["lookahead"] = args.lookahead
    if args.split_depth is not None:
        config["split_depth"] = args.split_depth
    if args.hard_split_gap is not None:
        config["hard_split_gap"] = args.hard_split_gap
    if args.hard_split_extra is not None:
        config["hard_split_extra"] = args.hard_split_extra
    if args.dfs_lookahead is not None:
        config["dfs_lookahead"] = args.dfs_lookahead
    if args.dfs_refine_lookahead is not None:
        config["dfs_refine_lookahead"] = args.dfs_refine_lookahead
    if args.dfs_refine_margin is not None:
        config["dfs_refine_margin"] = args.dfs_refine_margin
    if args.probe_nodes is not None:
        config["probe_nodes"] = args.probe_nodes
    if args.chunk_nodes is not None:
        config["max_nodes_per_task"] = args.chunk_nodes

    n = int(config["n"])

    shard_dir = args.shard_dir or os.path.join(
        os.path.dirname(__file__), "results", "shard_work"
    )

    # ── Shard: coordinator ────────────────────────────────────────────────
    if getattr(args, "shard_mode", None) == "coordinator":
        tasks = enumerate_shard_tasks(
            n,
            known_sequence=known,
            incumbent_hint=0,
            lookahead=int(config["lookahead"]),
            split_depth=int(config["split_depth"]),
            hard_split_gap=int(config["hard_split_gap"]),
            hard_split_extra=int(config["hard_split_extra"]),
            dfs_lookahead=int(config["dfs_lookahead"]),
            dfs_refine_lookahead=int(config["dfs_refine_lookahead"]),
            dfs_refine_margin=int(config["dfs_refine_margin"]),
            verbose=True,
        )
        shard_n_dir = os.path.join(shard_dir, f"n_{n:04d}")
        tasks_dir = os.path.join(shard_n_dir, "tasks")
        os.makedirs(tasks_dir, exist_ok=True)
        tasks.sort(key=lambda t: int(t.get("ub", 0)), reverse=True)
        for idx, task in enumerate(tasks):
            task["id"] = idx
            fname = os.path.join(tasks_dir, f"task_{idx:06d}.json")
            with open(fname, "w", encoding="utf-8") as f:
                json.dump(task, f)
        init_inc = int(tasks[0]["incumbent"]) if tasks else int(known.get(n, 0) or 0)
        inc_path = os.path.join(shard_n_dir, "incumbent.json")
        with open(inc_path, "w", encoding="utf-8") as f:
            json.dump({"incumbent": init_inc}, f)
        print(f"Wrote {len(tasks)} task files to {tasks_dir}")
        print(f"Initialized shard incumbent at {inc_path}: {init_inc}")
        print("Start workers with: python run_branch_bound.py --n <n> --shard-mode worker --shard-dir <dir>")
        return

    # ── Shard: worker ─────────────────────────────────────────────────────
    if getattr(args, "shard_mode", None) == "worker":
        chunk_nodes = int(config["max_nodes_per_task"]) if int(config["max_nodes_per_task"]) > 0 else 0
        n_processed = run_shard_worker_loop(
            shard_dir=shard_dir,
            n=n,
            chunk_nodes=chunk_nodes,
            verbose=True,
        )
        print(f"Worker done: processed {n_processed} tasks.")
        return

    # ── Shard: merge ──────────────────────────────────────────────────────
    if getattr(args, "shard_mode", None) == "merge":
        a, optimals, stats = merge_shard_results(shard_dir=shard_dir, n=n, verbose=True)
        expected = known.get(n)
        print(f"\na({n}) = {a}")
        if expected:
            print(f"correct = {a == expected}  (expected {expected})")
        print(f"#optimals = {len(optimals)}")
        if not args.no_save:
            _save_result(n, a, optimals, stats, config, known)
        return

    # ── Standard single-machine run ───────────────────────────────────────
    expected = known.get(n)
    cpu_count = os.cpu_count() or 1
    resolved_workers = int(config["num_workers"]) if int(config["num_workers"]) > 0 else max(1, cpu_count - 1)
    chunk_text = config["max_nodes_per_task"] if int(config["max_nodes_per_task"]) > 0 else "off"
    print(f"n={n}, expected a({n})={expected}, CPUs={cpu_count}")
    print(
        "config: "
        f"preset={config['preset']}, "
        f"collect_all={config['collect_all']}, "
        f"lookahead=L{config['lookahead']}, "
        f"workers={resolved_workers}, "
        f"split_depth={config['split_depth']}, "
        f"hard_split=+{config['hard_split_extra']}@{config['hard_split_gap']}, "
        f"dfs_L={config['dfs_lookahead']}, "
        f"dfs_refine_L={config['dfs_refine_lookahead']}±{config['dfs_refine_margin']}, "
        f"probe_nodes={config['probe_nodes']}, "
        f"chunk_nodes={chunk_text}"
    )

    # Warm up Numba in the main process before launching workers.
    branch_bound_a112509(10, known_sequence=known, lookahead=int(config["lookahead"]), verbose=False)

    t0 = time.perf_counter()
    a, optimals, stats = branch_bound_a112509_parallel(
        n,
        known_sequence=known,
        lookahead=int(config["lookahead"]),
        split_depth=int(config["split_depth"]),
        num_workers=int(config["num_workers"]),
        verbose=True,
        collect_all=bool(config["collect_all"]),
        hard_split_gap=int(config["hard_split_gap"]),
        hard_split_extra=int(config["hard_split_extra"]),
        dfs_lookahead=int(config["dfs_lookahead"]),
        dfs_refine_lookahead=int(config["dfs_refine_lookahead"]),
        dfs_refine_margin=int(config["dfs_refine_margin"]),
        probe_nodes=int(config["probe_nodes"]),
        max_nodes_per_task=int(config["max_nodes_per_task"]),
    )
    elapsed = time.perf_counter() - t0

    print(f"\n=== n={n} DONE ===")
    print(f"  a({n})     = {a}")
    correct_str = (str(a == expected) + f"  (expected {expected})") if expected else "unknown"
    print(f"  correct   = {correct_str}")
    print(f"  #optimals = {len(optimals)}")
    print(f"  nodes     = {stats['nodes_expanded']:,}")
    print(f"  pruned    = {stats['nodes_pruned']:,}")
    print(f"  subtrees  = {stats['num_subtrees']}")
    print(f"  time      = {elapsed:.1f}s  ({elapsed/60:.2f} min)")
    if optimals:
        print(f"\nFirst 5 optimal strings (out of {len(optimals)}):")
        for s in optimals[:5]:
            print(f"  {s}")

    if args.no_save:
        print("\nSkipping result write (--no-save).")
        return

    _save_result(n, a, optimals, stats, config, known, elapsed=elapsed)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        raise SystemExit(0)
