"""Item Response Theory analysis — fit 2PL model, information curves, flag low-discrimination items.

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


def classical_item_analysis(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Classical test theory fallback: difficulty (p-value) + point-biserial discrimination."""
    n_items, n_models = matrix.shape

    # Difficulty = proportion correct (inverted: higher = harder)
    difficulties = 1.0 - np.nanmean(matrix, axis=1)

    # Discrimination: point-biserial correlation between item and total score
    # Total score per model (excluding the item being analyzed)
    discriminations = np.zeros(n_items)
    for i in range(n_items):
        valid_mask = ~np.isnan(matrix[i])
        if valid_mask.sum() < 2:
            discriminations[i] = 0.0
            continue

        item_scores = matrix[i, valid_mask]
        # Total score for each model EXCLUDING this item
        other_items = np.delete(matrix, i, axis=0)
        total_scores = np.nansum(other_items[:, valid_mask], axis=0)

        if np.std(item_scores) == 0 or np.std(total_scores) == 0:
            discriminations[i] = 0.0
            continue

        discriminations[i] = np.corrcoef(item_scores, total_scores)[0, 1]

    return difficulties, discriminations


def compute_information(difficulties: np.ndarray, discriminations: np.ndarray, theta_range: np.ndarray = None) -> np.ndarray:
    """Compute test information function across ability levels.

    Uses 2PL IRT information formula: I(theta) = sum_i( a_i^2 * P_i * (1 - P_i) )
    where P_i = 1 / (1 + exp(-a_i * (theta - b_i)))
    """
    if theta_range is None:
        theta_range = np.linspace(-3, 3, 61)

    info = np.zeros_like(theta_range)
    for b, a in zip(difficulties, discriminations):
        if np.isnan(a) or np.isnan(b) or a <= 0:
            continue
        for j, theta in enumerate(theta_range):
            p = 1.0 / (1.0 + np.exp(-a * (theta - b)))
            info[j] += a ** 2 * p * (1 - p)

    return info


async def main():
    parser = argparse.ArgumentParser(description="IRT analysis of eval results")
    parser.add_argument("--discrimination-threshold", type=float, default=0.3)
    parser.add_argument("--output-csv", default=None)
    parser.add_argument("--output-info-plot", default=None)
    args = parser.parse_args()

    task_ids, model_names, matrix = await build_response_matrix()
    print(f"Response matrix: {matrix.shape[0]} tasks x {matrix.shape[1]} models")

    if matrix.shape[1] < 2:
        print("Need at least 2 models for IRT analysis.", file=sys.stderr)
        sys.exit(1)

    # Try py-irt first, fall back to classical
    try:
        from py_irt.models import TwoParamLog
        model = TwoParamLog(matrix.tolist())
        model.fit()
        difficulties = np.array(model.item_difficulties)
        discriminations = np.array(model.item_discriminations)
        print("Fitted 2PL IRT model via py-irt")
    except (ImportError, Exception) as e:
        print(f"py-irt unavailable ({e}), using classical item analysis")
        difficulties, discriminations = classical_item_analysis(matrix)

    # Output results
    low_disc = []
    too_easy = []
    too_hard = []
    results = []

    for i, task_id in enumerate(task_ids):
        disc = float(discriminations[i]) if not np.isnan(discriminations[i]) else 0.0
        diff = float(difficulties[i]) if not np.isnan(difficulties[i]) else 0.5
        p_correct = 1.0 - diff

        flag_disc = disc < args.discrimination_threshold
        flag_easy = p_correct > 0.95  # All models get it right
        flag_hard = p_correct < 0.05  # All models get it wrong

        results.append({
            "task_id": task_id, "difficulty": round(diff, 3),
            "discrimination": round(disc, 3), "p_correct": round(p_correct, 3),
            "low_discrimination": flag_disc, "too_easy": flag_easy, "too_hard": flag_hard,
        })
        if flag_disc:
            low_disc.append(task_id)
        if flag_easy:
            too_easy.append(task_id)
        if flag_hard:
            too_hard.append(task_id)

    print(f"\n=== Item Analysis Summary ===")
    print(f"Total items: {len(task_ids)}")
    print(f"Low discrimination (<{args.discrimination_threshold}): {len(low_disc)}")
    print(f"Too easy (p>0.95): {len(too_easy)}")
    print(f"Too hard (p<0.05): {len(too_hard)}")

    pruned = [r["task_id"] for r in results if not r["low_discrimination"] and not r["too_easy"] and not r["too_hard"]]
    print(f"Recommended pruned set: {len(pruned)} items (removed {len(task_ids) - len(pruned)})")

    # --- Test Information Function ---
    theta_range = np.linspace(-3, 3, 61)
    info = compute_information(difficulties, discriminations, theta_range)
    peak_theta = theta_range[np.argmax(info)]
    print(f"\nTest information peaks at theta={peak_theta:.2f} (ability level where benchmark is most discriminative)")

    if info.max() > 0:
        # Find where information drops below 50% of peak
        half_info = info.max() / 2
        above = theta_range[info >= half_info]
        if len(above) >= 2:
            print(f"Effective range: theta=[{above[0]:.2f}, {above[-1]:.2f}] (>50% peak information)")

    if args.output_csv:
        with open(args.output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "task_id", "difficulty", "discrimination", "p_correct",
                "low_discrimination", "too_easy", "too_hard",
            ])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nCSV saved to {args.output_csv}")

    if args.output_info_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(theta_range, info, color="steelblue", linewidth=2)
        ax.axvline(x=peak_theta, color="red", linestyle="--", alpha=0.5, label=f"Peak at θ={peak_theta:.2f}")
        ax.set_xlabel("Ability (θ)")
        ax.set_ylabel("Test Information")
        ax.set_title("Test Information Function")
        ax.legend()
        plt.tight_layout()
        fig.savefig(args.output_info_plot, dpi=150)
        print(f"Information plot saved to {args.output_info_plot}")


if __name__ == "__main__":
    asyncio.run(main())
