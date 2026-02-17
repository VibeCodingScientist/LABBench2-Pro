"""Token and cost accounting across all runs.

Usage: python -m src.tier3.cost_tracker
"""

import asyncio

from src.db import fetch


async def main():
    # Eval runs cost summary
    eval_rows = await fetch(
        """SELECT m.model_name,
                  COUNT(*) as runs,
                  SUM(e.tokens_in) as total_tokens_in,
                  SUM(e.tokens_out) as total_tokens_out,
                  SUM(e.cost_usd) as total_cost,
                  AVG(e.latency_ms) as avg_latency
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           GROUP BY m.model_name
           ORDER BY total_cost DESC"""
    )

    print("=== Eval Runs Cost Summary ===")
    print(f"{'Model':<25} {'Runs':>6} {'Tokens In':>12} {'Tokens Out':>12} {'Cost ($)':>10} {'Avg Latency':>12}")
    print("-" * 80)

    grand_total = 0.0
    for row in eval_rows:
        cost = row["total_cost"] or 0
        grand_total += cost
        print(
            f"{row['model_name']:<25} {row['runs']:>6} "
            f"{row['total_tokens_in'] or 0:>12,} {row['total_tokens_out'] or 0:>12,} "
            f"{cost:>10.4f} {row['avg_latency'] or 0:>10.0f}ms"
        )

    # Cost by category
    cat_rows = await fetch(
        """SELECT t.category,
                  COUNT(*) as runs,
                  SUM(e.cost_usd) as total_cost
           FROM eval_runs e
           JOIN tasks t ON e.task_id = t.id
           GROUP BY t.category
           ORDER BY total_cost DESC"""
    )

    if cat_rows:
        print(f"\n=== Cost by Category ===")
        print(f"{'Category':<25} {'Runs':>6} {'Cost ($)':>10}")
        print("-" * 45)
        for row in cat_rows:
            print(f"{row['category']:<25} {row['runs']:>6} {row['total_cost'] or 0:>10.4f}")

    # Chain runs summary
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
        print(f"{'Chain':<30} {'Model':<25} {'Steps':>6} {'Correct':>8}")
        print("-" * 75)
        for row in chain_rows:
            print(
                f"{row['chain_id']:<30} {row['model_name']:<25} "
                f"{row['steps']:>6} {row['correct_steps']:>8}"
            )

    print(f"\n{'Grand Total Cost:':<25} ${grand_total:>10.4f}")

    if not eval_rows and not chain_rows:
        print("\nNo run data found. Run some evals first!")


if __name__ == "__main__":
    asyncio.run(main())
