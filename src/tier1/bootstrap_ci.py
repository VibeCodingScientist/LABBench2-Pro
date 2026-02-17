"""Compute BCa bootstrap confidence intervals from eval results.

Usage: python -m src.tier1.bootstrap_ci --category LitQA3
"""

import argparse
import asyncio
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import bootstrap

from src.db import fetch


def compute_ci(correct_vector: list[bool], confidence: float = 0.95) -> tuple[float, float, float]:
    """Compute BCa bootstrap CI for accuracy."""
    arr = np.array(correct_vector, dtype=float)
    if len(arr) < 2:
        return arr.mean(), 0.0, 1.0
    result = bootstrap(
        (arr,), np.mean, n_resamples=10000,
        confidence_level=confidence, method="BCa",
    )
    return float(arr.mean()), float(result.confidence_interval.low), float(result.confidence_interval.high)


async def main():
    parser = argparse.ArgumentParser(description="Bootstrap CIs for eval results")
    parser.add_argument("--category", required=True)
    parser.add_argument("--output-csv", default=None)
    parser.add_argument("--output-plot", default=None)
    args = parser.parse_args()

    rows = await fetch(
        """SELECT m.model_name, e.correct
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           JOIN tasks t ON e.task_id = t.id
           WHERE t.category = $1""",
        args.category,
    )

    if not rows:
        print(f"No eval results found for category '{args.category}'", file=sys.stderr)
        sys.exit(1)

    # Group by model
    by_model: dict[str, list[bool]] = {}
    for row in rows:
        by_model.setdefault(row["model_name"], []).append(row["correct"])

    # Compute CIs
    results = []
    for model, corrects in sorted(by_model.items()):
        mean, ci_low, ci_high = compute_ci(corrects)
        results.append({"model": model, "n": len(corrects), "accuracy": mean, "ci_low": ci_low, "ci_high": ci_high})
        print(f"{model}: {mean:.3f} [{ci_low:.3f}, {ci_high:.3f}] (n={len(corrects)})")

    # CSV output
    if args.output_csv:
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["model", "n", "accuracy", "ci_low", "ci_high"])
            writer.writeheader()
            writer.writerows(results)
        print(f"CSV saved to {args.output_csv}")

    # Plot
    if args.output_plot:
        models = [r["model"] for r in results]
        accs = [r["accuracy"] for r in results]
        errs_low = [r["accuracy"] - r["ci_low"] for r in results]
        errs_high = [r["ci_high"] - r["accuracy"] for r in results]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(models, accs, xerr=[errs_low, errs_high], capsize=5, color="steelblue")
        ax.set_xlabel("Accuracy")
        ax.set_title(f"LABBench2 — {args.category} (95% BCa CI)")
        ax.set_xlim(0, 1)
        plt.tight_layout()
        fig.savefig(args.output_plot, dpi=150)
        print(f"Plot saved to {args.output_plot}")


if __name__ == "__main__":
    asyncio.run(main())
