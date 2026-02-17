"""Item Response Theory analysis — fit 2PL model, flag low-discrimination items.

Usage: python -m src.tier1.irt_analysis
"""

import argparse
import asyncio
import csv
import sys

import numpy as np

from src.db import fetch


async def build_response_matrix() -> tuple[list[str], list[str], np.ndarray]:
    """Build task x model response matrix from eval_runs.

    Returns: (task_ids, model_names, matrix[n_tasks, n_models])
    """
    rows = await fetch(
        """SELECT t.id as task_id, m.model_name, e.correct
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           JOIN tasks t ON e.task_id = t.id
           ORDER BY t.id, m.model_name"""
    )
    if not rows:
        print("No eval results found.", file=sys.stderr)
        sys.exit(1)

    # Collect unique task_ids and models
    task_set: dict[str, int] = {}
    model_set: dict[str, int] = {}
    for row in rows:
        if row["task_id"] not in task_set:
            task_set[row["task_id"]] = len(task_set)
        if row["model_name"] not in model_set:
            model_set[row["model_name"]] = len(model_set)

    task_ids = list(task_set.keys())
    model_names = list(model_set.keys())

    matrix = np.full((len(task_ids), len(model_names)), np.nan)
    for row in rows:
        ti = task_set[row["task_id"]]
        mi = model_set[row["model_name"]]
        matrix[ti, mi] = 1.0 if row["correct"] else 0.0

    return task_ids, model_names, matrix


async def main():
    parser = argparse.ArgumentParser(description="IRT analysis of eval results")
    parser.add_argument("--discrimination-threshold", type=float, default=0.3)
    parser.add_argument("--output-csv", default=None)
    args = parser.parse_args()

    task_ids, model_names, matrix = await build_response_matrix()
    print(f"Response matrix: {matrix.shape[0]} tasks x {matrix.shape[1]} models")

    if matrix.shape[1] < 2:
        print("Need at least 2 models for IRT analysis.", file=sys.stderr)
        sys.exit(1)

    try:
        from py_irt.models import TwoParamLog

        model = TwoParamLog(matrix.tolist())
        model.fit()
        difficulties = model.item_difficulties
        discriminations = model.item_discriminations
    except ImportError:
        print("py-irt not installed, falling back to classical item analysis.", file=sys.stderr)
        # Classical fallback: difficulty = 1 - p(correct), discrimination = point-biserial
        difficulties = 1.0 - np.nanmean(matrix, axis=1)
        total_scores = np.nansum(matrix, axis=0)
        discriminations = np.array([
            np.corrcoef(matrix[i, ~np.isnan(matrix[i])], total_scores[~np.isnan(matrix[i])])[0, 1]
            if np.sum(~np.isnan(matrix[i])) > 1 else 0.0
            for i in range(matrix.shape[0])
        ])

    # Output results
    low_disc = []
    results = []
    for i, task_id in enumerate(task_ids):
        disc = float(discriminations[i]) if not np.isnan(discriminations[i]) else 0.0
        diff = float(difficulties[i]) if not np.isnan(difficulties[i]) else 0.5
        flag = disc < args.discrimination_threshold
        results.append({"task_id": task_id, "difficulty": diff, "discrimination": disc, "low_signal": flag})
        if flag:
            low_disc.append(task_id)

    print(f"\nItem parameters computed for {len(task_ids)} tasks")
    print(f"Low discrimination (<{args.discrimination_threshold}): {len(low_disc)} items")
    if low_disc[:10]:
        print(f"  Examples: {low_disc[:10]}")

    pruned = [r["task_id"] for r in results if not r["low_signal"]]
    print(f"Recommended pruned set: {len(pruned)} items")

    if args.output_csv:
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["task_id", "difficulty", "discrimination", "low_signal"])
            writer.writeheader()
            writer.writerows(results)
        print(f"CSV saved to {args.output_csv}")


if __name__ == "__main__":
    asyncio.run(main())
