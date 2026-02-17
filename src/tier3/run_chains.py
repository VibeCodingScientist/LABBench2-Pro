"""Execute compositional chains, measure error propagation.

Usage: python -m src.tier3.run_chains --model claude-opus-4.5 --chain lit_to_primer_01
"""

import argparse
import asyncio
import json
import sys

from src.config import MODEL_REGISTRY
from src.db import execute, fetch, fetchrow, fetchval
from src.models import call_model
from src.tier1.grading import grade
from src.tier3.chains import get_chain, list_chain_ids


async def ensure_model(model_key: str) -> int:
    info = MODEL_REGISTRY[model_key]
    row = await fetchrow(
        "SELECT id FROM models WHERE provider=$1 AND model_name=$2 AND tools=$3",
        info["provider"], model_key, False,
    )
    if row:
        return row["id"]
    return await fetchval(
        "INSERT INTO models (provider, model_name, tools) VALUES ($1, $2, $3) RETURNING id",
        info["provider"], model_key, False,
    )


async def get_task(task_id: str) -> dict | None:
    row = await fetchrow("SELECT * FROM tasks WHERE id=$1", task_id)
    if not row:
        return None
    meta = json.loads(row["meta"]) if row["meta"] else {}
    return {
        "id": row["id"],
        "category": row["category"],
        "question": meta.get("question", f"Task {row['id']}"),
        "ideal": meta.get("ideal", ""),
        "verification": meta.get("verification", "llm-judge"),
        **meta,
    }


async def run_chain(model_key: str, chain_id: str):
    chain = get_chain(chain_id)
    if not chain:
        available = list_chain_ids()
        print(f"Chain '{chain_id}' not found. Available: {available}", file=sys.stderr)
        sys.exit(1)

    model_id = await ensure_model(model_key)
    print(f"Running chain '{chain_id}': {chain['description']}")
    print(f"Model: {model_key} (id={model_id})")
    print(f"Steps: {len(chain['steps'])}")
    print()

    # Resume support: find already-completed steps for this chain+model
    existing = await fetch(
        """SELECT step_num, response, correct FROM chain_runs
           WHERE chain_id=$1 AND model_id=$2 ORDER BY step_num""",
        chain_id, model_id,
    )
    completed_steps = {row["step_num"]: row for row in existing}
    if completed_steps:
        print(f"Resuming — {len(completed_steps)} steps already completed")

    previous_response = None
    results = []

    for i, step in enumerate(chain["steps"]):
        step_num = i + 1
        task_id = step["task_id"]
        input_type = step["input_type"]

        # Skip already-completed steps, but carry forward the response
        if step_num in completed_steps:
            prev = completed_steps[step_num]
            previous_response = prev["response"]
            status = "CORRECT" if prev["correct"] else "WRONG"
            print(f"  Step {step_num} [{status}]: {task_id} (cached)")
            results.append({"step": step_num, "task_id": task_id, "correct": prev["correct"], "skipped": True})
            continue

        task = await get_task(task_id)
        if not task:
            print(f"  Step {step_num}: Task '{task_id}' not found in DB, skipping")
            results.append({"step": step_num, "task_id": task_id, "error": "task not found"})
            continue

        # Build prompt
        prompt = task.get("question", "")
        if input_type == "previous_response" and previous_response:
            prompt = f"Previous context:\n{previous_response}\n\n---\n\n{prompt}"

        system = "You are an expert biologist working through a multi-step research task."

        try:
            result = await call_model(model_key, prompt, system=system, cache=False)
        except Exception as e:
            print(f"  Step {i+1}: ERROR — {e}")
            results.append({"step": i + 1, "task_id": task_id, "error": str(e)})
            previous_response = None
            continue

        grading = await grade(task, result["response"])

        # Save to DB
        await execute(
            """INSERT INTO chain_runs (chain_id, model_id, step_num, task_id,
               input_from, response, correct)
               VALUES ($1, $2, $3, $4, $5, $6, $7)""",
            chain_id, model_id, i + 1, task_id,
            previous_response[:500] if previous_response else None,
            result["response"], grading["correct"],
        )

        status = "CORRECT" if grading["correct"] else "WRONG"
        print(f"  Step {i+1} [{status}]: {task_id} ({result['latency_ms']}ms)")

        previous_response = result["response"]
        results.append({
            "step": i + 1, "task_id": task_id,
            "correct": grading["correct"], "latency_ms": result["latency_ms"],
        })

    # Summary
    completed = [r for r in results if "correct" in r]
    correct = sum(1 for r in completed if r["correct"])
    chain_correct = all(r["correct"] for r in completed) if completed else False

    print(f"\n--- Chain Summary ---")
    print(f"Steps completed: {len(completed)}/{len(chain['steps'])}")
    print(f"Steps correct: {correct}/{len(completed)}")
    print(f"End-to-end correct: {chain_correct}")


async def main():
    parser = argparse.ArgumentParser(description="Run compositional chains")
    parser.add_argument("--model", required=True, choices=list(MODEL_REGISTRY.keys()))
    parser.add_argument("--chain", required=True)
    args = parser.parse_args()

    await run_chain(args.model, args.chain)


if __name__ == "__main__":
    asyncio.run(main())
