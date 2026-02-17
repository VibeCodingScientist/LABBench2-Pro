"""Token/cost accounting + cost-accuracy Pareto frontier.

Usage: python -m src.tier3.cost_tracker
"""

import asyncio

import numpy as np

from src.db import fetch


def compute_pareto_frontier(points: list[dict]) -> list[str]:
    """Find Pareto-optimal models (max accuracy, min cost).

    Returns list of model names on the frontier.
    """
    if not points:
        return []

    # Sort by cost ascending
    sorted_pts = sorted(points, key=lambda p: p["cost_per_task"])
    frontier = []
    max_acc = -1.0

    for pt in sorted_pts:
        if pt["accuracy"] > max_acc:
            frontier.append(pt["model"])
            max_acc = pt["accuracy"]

    return frontier


async def main():
    # --- Eval runs cost summary ---
    eval_rows = await fetch(
        """SELECT m.model_name,
                  COUNT(*) as runs,
                  SUM(e.tokens_in) as total_tokens_in,
                  SUM(e.tokens_out) as total_tokens_out,
                  SUM(e.cost_usd) as total_cost,
                  AVG(e.latency_ms) as avg_latency,
                  COUNT(*) FILTER (WHERE e.correct) as correct_count
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           GROUP BY m.model_name
           ORDER BY total_cost DESC"""
    )

    print("=== Eval Runs Cost Summary ===")
    print(f"{'Model':<25} {'Runs':>6} {'Tokens In':>12} {'Tokens Out':>12} {'Cost ($)':>10} {'Avg Latency':>12} {'Accuracy':>9}")
    print("-" * 90)

    grand_total = 0.0
    pareto_points = []

    for row in eval_rows:
        cost = row["total_cost"] or 0
        grand_total += cost
        acc = row["correct_count"] / row["runs"] if row["runs"] > 0 else 0
        cost_per_task = cost / row["runs"] if row["runs"] > 0 else 0
        print(
            f"{row['model_name']:<25} {row['runs']:>6} "
            f"{row['total_tokens_in'] or 0:>12,} {row['total_tokens_out'] or 0:>12,} "
            f"{cost:>10.4f} {row['avg_latency'] or 0:>10.0f}ms {acc:>8.1%}"
        )
        pareto_points.append({
            "model": row["model_name"],
            "accuracy": acc,
            "cost_per_task": cost_per_task,
            "total_cost": cost,
        })

    # --- Pareto frontier ---
    if len(pareto_points) >= 2:
        frontier = compute_pareto_frontier(pareto_points)
        print(f"\n=== Cost-Accuracy Pareto Frontier ===")
        print(f"Pareto-optimal models (best accuracy at each cost level):")
        for model in frontier:
            pt = next(p for p in pareto_points if p["model"] == model)
            print(f"  {model}: accuracy={pt['accuracy']:.1%}, cost/task=${pt['cost_per_task']:.4f}")

        dominated = [p["model"] for p in pareto_points if p["model"] not in frontier]
        if dominated:
            print(f"Dominated (never optimal): {', '.join(dominated)}")

    # --- Cost by category ---
    cat_rows = await fetch(
        """SELECT t.category,
                  COUNT(*) as runs,
                  SUM(e.cost_usd) as total_cost,
                  COUNT(*) FILTER (WHERE e.correct) as correct_count
           FROM eval_runs e
           JOIN tasks t ON e.task_id = t.id
           GROUP BY t.category
           ORDER BY total_cost DESC"""
    )

    if cat_rows:
        print(f"\n=== Cost by Category ===")
        print(f"{'Category':<25} {'Runs':>6} {'Cost ($)':>10} {'$/correct':>10}")
        print("-" * 55)
        for row in cat_rows:
            cost_per_correct = (row["total_cost"] or 0) / row["correct_count"] if row["correct_count"] else float("inf")
            print(f"{row['category']:<25} {row['runs']:>6} {row['total_cost'] or 0:>10.4f} {cost_per_correct:>10.4f}")

    # --- Chain runs summary ---
    chain_rows = await fetch(
        """SELECT c.chain_id, m.model_name,
                  COUNT(*) as steps,
                  COUNT(*) FILTER (WHERE c.correct) as correct_steps
           FROM chain_runs c
           JOIN models m ON c.model_id = m.id
           GROUP BY c.chain_id, m.model_name
           ORDER BY c.chain_id"""
    )

    if chain_rows:
        print(f"\n=== Chain Runs Summary ===")
        print(f"{'Chain':<30} {'Model':<25} {'Steps':>6} {'Correct':>8} {'End2End':>8}")
        print("-" * 80)
        for row in chain_rows:
            e2e = "YES" if row["correct_steps"] == row["steps"] else "no"
            print(
                f"{row['chain_id']:<30} {row['model_name']:<25} "
                f"{row['steps']:>6} {row['correct_steps']:>8} {e2e:>8}"
            )

    print(f"\n{'Grand Total Cost:':<25} ${grand_total:>10.4f}")

    if not eval_rows and not chain_rows:
        print("\nNo run data found. Run some evals first!")


if __name__ == "__main__":
    asyncio.run(main())
