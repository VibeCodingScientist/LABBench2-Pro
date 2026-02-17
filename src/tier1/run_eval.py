"""Run models against LABBench2 tasks. Save results to DB.

Usage: python -m src.tier1.run_eval --model claude-opus-4.5 --category LitQA3 --concurrency 5
"""

import argparse
import asyncio
import json
import sys

from datasets import load_dataset

from src.config import HF_TOKEN, MODEL_REGISTRY
from src.db import execute, fetch, fetchrow, fetchval
from src.models import call_model
from src.tier1.grading import grade


async def ensure_model(model_key: str) -> int:
    """Insert model into DB if not exists, return model_id."""
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


async def upsert_task(task: dict) -> None:
    """Insert task into DB (skip if exists)."""
    await execute(
        """INSERT INTO tasks (id, category, variant, source, meta)
           VALUES ($1, $2, $3, $4, $5)
           ON CONFLICT (id) DO NOTHING""",
        task["id"], task["category"], task.get("variant"),
        task.get("source", "labbench2"), json.dumps(task.get("meta", {})),
    )


async def already_evaluated(model_id: int, task_id: str) -> bool:
    """Check if this model+task combo already has a result in the DB."""
    row = await fetchrow(
        "SELECT id FROM eval_runs WHERE model_id=$1 AND task_id=$2 LIMIT 1",
        model_id, task_id,
    )
    return row is not None


async def run_single(model_key: str, model_id: int, task: dict, sem: asyncio.Semaphore, force: bool = False) -> dict:
    """Run a single eval: call model, grade, insert result. Skips if already done."""
    async with sem:
        # Resume support: skip tasks already evaluated for this model
        if not force and await already_evaluated(model_id, task["id"]):
            row = await fetchrow(
                "SELECT correct FROM eval_runs WHERE model_id=$1 AND task_id=$2 LIMIT 1",
                model_id, task["id"],
            )
            print(f"  [SKIP] {task['id']} (already evaluated)")
            return {"task_id": task["id"], "correct": row["correct"], "skipped": True}

        prompt = task.get("question", task.get("prompt", ""))
        system = task.get("system", "You are an expert biologist answering scientific questions.")

        try:
            result = await call_model(model_key, prompt, system=system)
        except Exception as e:
            print(f"  ERROR on {task['id']}: {e}", file=sys.stderr)
            return {"task_id": task["id"], "error": str(e)}

        grading = await grade(task, result["response"])

        await execute(
            """INSERT INTO eval_runs (model_id, task_id, response, correct, score, grader,
               tokens_in, tokens_out, cost_usd, latency_ms, meta)
               VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)""",
            model_id, task["id"], result["response"],
            grading["correct"], grading["score"], grading["grader"],
            result["tokens_in"], result["tokens_out"], result["cost_usd"], result["latency_ms"],
            json.dumps({"explanation": grading.get("explanation", "")}),
        )

        status = "CORRECT" if grading["correct"] else "WRONG"
        print(f"  [{status}] {task['id']} ({result['latency_ms']}ms, ${result['cost_usd']:.4f})")
        return {"task_id": task["id"], "correct": grading["correct"]}


async def load_tasks(category: str) -> list[dict]:
    """Load tasks from HuggingFace datasets."""
    token = HF_TOKEN if HF_TOKEN else None
    ds = load_dataset("futurehouse/labbench2", category, token=token, split="test")
    tasks = []
    for i, item in enumerate(ds):
        task = {
            "id": item.get("id", f"{category}_{i:04d}"),
            "category": category,
            "question": item.get("question", item.get("prompt", "")),
            "ideal": item.get("ideal", item.get("answer", "")),
            "verification": item.get("verification", "llm-judge"),
            "source": "labbench2",
            "meta": {k: v for k, v in item.items() if k not in ("id", "question", "prompt", "ideal", "answer", "verification")},
        }
        tasks.append(task)
    return tasks


async def main():
    parser = argparse.ArgumentParser(description="Run LABBench2 eval")
    parser.add_argument("--model", required=True, choices=list(MODEL_REGISTRY.keys()))
    parser.add_argument("--category", required=True)
    parser.add_argument("--concurrency", type=int, default=5)
    parser.add_argument("--limit", type=int, default=None, help="Limit number of tasks (for testing)")
    parser.add_argument("--force", action="store_true", help="Re-run tasks even if already evaluated")
    args = parser.parse_args()

    print(f"Loading tasks for {args.category}...")
    tasks = await load_tasks(args.category)
    if args.limit:
        tasks = tasks[:args.limit]
    print(f"Loaded {len(tasks)} tasks")

    model_id = await ensure_model(args.model)
    print(f"Model: {args.model} (id={model_id})")

    # Upsert all tasks
    for task in tasks:
        await upsert_task(task)

    # Run evals with concurrency limit
    sem = asyncio.Semaphore(args.concurrency)
    print(f"Running eval (concurrency={args.concurrency})...")
    results = await asyncio.gather(*[run_single(args.model, model_id, t, sem, force=args.force) for t in tasks])

    # Summary
    correct = sum(1 for r in results if r.get("correct"))
    errors = sum(1 for r in results if "error" in r)
    skipped = sum(1 for r in results if r.get("skipped"))
    total = len(results) - errors
    print(f"\nResults: {correct}/{total} correct ({correct/total*100:.1f}%)" if total > 0 else "\nNo results.")
    if skipped:
        print(f"Skipped (already done): {skipped}")
    if errors:
        print(f"Errors: {errors}")


if __name__ == "__main__":
    asyncio.run(main())
