"""Feedback loop simulation — re-run chains with correctness signal.

Usage: python -m src.tier3.feedback_sim --model claude-opus-4.5 --chain lit_to_primer_01
"""

import argparse
import asyncio
import json
import sys

from src.config import MODEL_REGISTRY
from src.db import execute, fetchrow, fetchval
from src.models import call_model
from src.tier1.grading import grade
from src.tier3.chains import get_chain, list_chain_ids


async def ensure_model(model_key: str) -> int:
    from src.tier3.run_chains import ensure_model as _em
    return await _em(model_key)


async def get_task(task_id: str) -> dict | None:
    from src.tier3.run_chains import get_task as _gt
    return await _gt(task_id)


async def run_chain_with_feedback(model_key: str, chain_id: str):
    """Re-run a chain, but feed back correctness after each step."""
    chain = get_chain(chain_id)
    if not chain:
        print(f"Chain '{chain_id}' not found.", file=sys.stderr)
        sys.exit(1)

    model_id = await ensure_model(model_key)
    print(f"Running chain '{chain_id}' WITH feedback")
    print(f"Model: {model_key}")
    print()

    previous_response = None
    previous_correct = None
    results = []

    for i, step in enumerate(chain["steps"]):
        task_id = step["task_id"]
        input_type = step["input_type"]

        task = await get_task(task_id)
        if not task:
            print(f"  Step {i+1}: Task '{task_id}' not found, skipping")
            results.append({"step": i + 1, "task_id": task_id, "error": "not found"})
            continue

        # Build prompt with feedback
        prompt = task.get("question", "")
        if input_type == "previous_response" and previous_response:
            feedback_note = ""
            if previous_correct is not None:
                if previous_correct:
                    feedback_note = "\n[Note: Your previous answer was CORRECT.]\n"
                else:
                    feedback_note = "\n[Note: Your previous answer was INCORRECT. Please reconsider your approach.]\n"

            prompt = (
                f"Previous context:\n{previous_response}\n"
                f"{feedback_note}\n---\n\n{prompt}"
            )

        system = "You are an expert biologist working through a multi-step research task."

        try:
            result = await call_model(model_key, prompt, system=system, cache=False)
        except Exception as e:
            print(f"  Step {i+1}: ERROR — {e}")
            results.append({"step": i + 1, "task_id": task_id, "error": str(e)})
            previous_response = None
            previous_correct = None
            continue

        grading = await grade(task, result["response"])

        # Save as a separate chain run with feedback prefix
        await execute(
            """INSERT INTO chain_runs (chain_id, model_id, step_num, task_id,
               input_from, response, correct)
               VALUES ($1, $2, $3, $4, $5, $6, $7)""",
            f"{chain_id}_feedback", model_id, i + 1, task_id,
            previous_response[:500] if previous_response else None,
            result["response"], grading["correct"],
        )

        status = "CORRECT" if grading["correct"] else "WRONG"
        print(f"  Step {i+1} [{status}]: {task_id}")

        previous_response = result["response"]
        previous_correct = grading["correct"]
        results.append({"step": i + 1, "task_id": task_id, "correct": grading["correct"]})

    # Compare with no-feedback run
    no_fb_rows = await fetchrow(
        """SELECT COUNT(*) FILTER (WHERE correct) as correct_count,
                  COUNT(*) as total
           FROM chain_runs
           WHERE chain_id = $1 AND model_id = $2""",
        chain_id, model_id,
    )

    completed = [r for r in results if "correct" in r]
    fb_correct = sum(1 for r in completed if r["correct"])

    print(f"\n--- Feedback Comparison ---")
    print(f"With feedback: {fb_correct}/{len(completed)} correct")
    if no_fb_rows and no_fb_rows["total"] > 0:
        print(f"Without feedback: {no_fb_rows['correct_count']}/{no_fb_rows['total']} correct")
        improvement = fb_correct / len(completed) - no_fb_rows["correct_count"] / no_fb_rows["total"]
        print(f"Improvement: {improvement:+.3f}")
    else:
        print("No baseline run found — run without feedback first.")


async def main():
    parser = argparse.ArgumentParser(description="Run chains with feedback")
    parser.add_argument("--model", required=True, choices=list(MODEL_REGISTRY.keys()))
    parser.add_argument("--chain", required=True)
    args = parser.parse_args()

    await run_chain_with_feedback(args.model, args.chain)


if __name__ == "__main__":
    asyncio.run(main())
