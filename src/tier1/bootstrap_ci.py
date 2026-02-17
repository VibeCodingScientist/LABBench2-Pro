"""Compute BCa bootstrap confidence intervals + pairwise model comparisons.

Usage: python -m src.tier1.bootstrap_ci --category LitQA2
"""

import argparse
import asyncio
import csv
import itertools
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import bootstrap

from src.db import fetch


def compute_ci(correct_vector: list[bool], confidence: float = 0.95) -> tuple[float, float, float]:
    """Compute BCa bootstrap CI for accuracy."""
    arr = np.array(correct_vector, dtype=float)
    if len(arr) < 2:
        return float(arr.mean()), 0.0, 1.0
    result = bootstrap(
        (arr,), np.mean, n_resamples=10000,
        confidence_level=confidence, method="BCa",
    )
    return float(arr.mean()), float(result.confidence_interval.low), float(result.confidence_interval.high)


def pairwise_bootstrap_test(vec_a: list[bool], vec_b: list[bool], n_resamples: int = 10000) -> float:
    """Bootstrap test for difference in accuracy between two models.

    Returns p-value (proportion of bootstrap samples where difference flips sign).
    """
    a = np.array(vec_a, dtype=float)
    b = np.array(vec_b, dtype=float)
    observed_diff = a.mean() - b.mean()
    if observed_diff == 0:
        return 1.0

    # Paired bootstrap on shared items (same tasks for both models)
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    count_flipped = 0
    rng = np.random.default_rng(42)
    for _ in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        boot_diff = a[idx].mean() - b[idx].mean()
        if (observed_diff > 0 and boot_diff <= 0) or (observed_diff < 0 and boot_diff >= 0):
            count_flipped += 1
    return count_flipped / n_resamples


async def main():
    parser = argparse.ArgumentParser(description="Bootstrap CIs + pairwise model comparisons")
    parser.add_argument("--category", required=True)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--output-csv", default=None)
    parser.add_argument("--output-plot", default=None)
    args = parser.parse_args()

    rows = await fetch(
        """SELECT m.model_name, e.task_id, e.correct
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           JOIN tasks t ON e.task_id = t.id
           WHERE t.category = $1
           ORDER BY m.model_name, e.task_id""",
        args.category,
    )

    if not rows:
        print(f"No eval results found for category '{args.category}'", file=sys.stderr)
        sys.exit(1)

    # Group by model
    by_model: dict[str, list[bool]] = {}
    # Also track per-task results for paired comparisons
    by_model_task: dict[str, dict[str, bool]] = {}
    for row in rows:
        by_model.setdefault(row["model_name"], []).append(row["correct"])
        by_model_task.setdefault(row["model_name"], {})[row["task_id"]] = row["correct"]

    # --- Per-model CIs ---
    print(f"=== {args.category} — {args.confidence*100:.0f}% BCa Bootstrap CIs ===\n")
    results = []
    for model, corrects in sorted(by_model.items()):
        mean, ci_low, ci_high = compute_ci(corrects, args.confidence)
        results.append({"model": model, "n": len(corrects), "accuracy": mean, "ci_low": ci_low, "ci_high": ci_high})
        print(f"  {model}: {mean:.3f} [{ci_low:.3f}, {ci_high:.3f}] (n={len(corrects)})")

    # --- Pairwise comparisons ---
    models = sorted(by_model.keys())
    if len(models) >= 2:
        print(f"\n=== Pairwise Comparisons ===\n")
        print(f"  {'Model A':<25} {'Model B':<25} {'Diff':>7} {'p-value':>8} {'Significant':>12}")
        print(f"  {'-'*80}")

        pairs = []
        for a, b in itertools.combinations(models, 2):
            # Get shared tasks
            shared_tasks = sorted(set(by_model_task[a].keys()) & set(by_model_task[b].keys()))
            if len(shared_tasks) < 5:
                continue
            vec_a = [by_model_task[a][t] for t in shared_tasks]
            vec_b = [by_model_task[b][t] for t in shared_tasks]
            diff = np.mean(vec_a) - np.mean(vec_b)
            p_val = pairwise_bootstrap_test(vec_a, vec_b)

            # Bonferroni correction
            n_comparisons = len(models) * (len(models) - 1) // 2
            sig = p_val < (0.05 / n_comparisons)
            sig_str = "YES" if sig else "no (overlap)"

            pairs.append({"model_a": a, "model_b": b, "diff": diff, "p_value": p_val, "significant": sig})
            print(f"  {a:<25} {b:<25} {diff:>+7.3f} {p_val:>8.4f} {sig_str:>12}")

        # Summary
        overlapping = [p for p in pairs if not p["significant"]]
        if overlapping:
            print(f"\n  WARNING: {len(overlapping)} model pair(s) have overlapping CIs:")
            for p in overlapping:
                print(f"    {p['model_a']} vs {p['model_b']} (diff={p['diff']:+.3f}, p={p['p_value']:.4f})")

    # CSV output
    if args.output_csv:
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["model", "n", "accuracy", "ci_low", "ci_high"])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nCSV saved to {args.output_csv}")

    # Plot
    if args.output_plot:
        models_sorted = sorted(results, key=lambda r: r["accuracy"])
        names = [r["model"] for r in models_sorted]
        accs = [r["accuracy"] for r in models_sorted]
        errs_low = [r["accuracy"] - r["ci_low"] for r in models_sorted]
        errs_high = [r["ci_high"] - r["accuracy"] for r in models_sorted]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(names, accs, xerr=[errs_low, errs_high], capsize=5, color="steelblue")
        ax.set_xlabel("Accuracy")
        ax.set_title(f"LABBench2 — {args.category} ({args.confidence*100:.0f}% BCa CI)")
        ax.set_xlim(0, 1)
        plt.tight_layout()
        fig.savefig(args.output_plot, dpi=150)
        print(f"Plot saved to {args.output_plot}")


if __name__ == "__main__":
    asyncio.run(main())
