"""Run models against LABBench tasks. Save results to DB.

Supports both:
  - futurehouse/labbench2 (gated, requires HF_TOKEN + access approval)
  - futurehouse/lab-bench (open, CC-BY-SA-4.0 fallback)

Also runs Tier 2 generated tasks from tasks/ directory.

Usage:
  python -m src.tier1.run_eval --model claude-opus-4.6 --category LitQA2 --concurrency 5
  python -m src.tier1.run_eval --model claude-opus-4.6 --tasks-dir tasks/stats_reasoning
"""

import argparse
import asyncio
import glob
import json
import os
import sys

from src.config import HF_TOKEN, MODEL_REGISTRY
from src.db import execute, fetch, fetchrow, fetchval
from src.models import call_model
from src.tier1.grading import grade

# LABBench2 categories (gated dataset)
LABBENCH2_CATEGORIES = [
    "LitQA3", "FigQA2", "SeqQA2", "SuppQA2", "ProtocolQA2",
    "DbQA2", "TableQA2", "CloningScenarios2",
]

# Original LAB-Bench categories (open dataset)
LABBENCH_CATEGORIES = [
    "LitQA2", "FigQA", "SeqQA", "SuppQA", "ProtocolQA",
    "DbQA", "TableQA", "CloningScenarios",
]


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

        prompt = task.get("question", "")
        system = "You are an expert biologist answering scientific questions. Give a concise answer."

        # For multiple-choice: append choices to prompt
        if task.get("choices"):
            choices_str = "\n".join(f"  {c}" for c in task["choices"])
            prompt = f"{prompt}\n\nChoices:\n{choices_str}\n\nRespond with only the letter of the correct answer."

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


def load_tasks_from_hf(category: str) -> list[dict]:
    """Load tasks from HuggingFace. Tries labbench2 first, falls back to lab-bench."""
    from datasets import load_dataset

    token = HF_TOKEN if HF_TOKEN else None

    # Try labbench2 (gated) first
    if category in LABBENCH2_CATEGORIES and token:
        try:
            ds = load_dataset("futurehouse/labbench2", category, token=token, split="test")
            print(f"  Loaded from futurehouse/labbench2")
            return _parse_hf_dataset(ds, category, source="labbench2")
        except Exception as e:
            print(f"  labbench2 failed ({e}), trying lab-bench...")

    # Fall back to original lab-bench (open)
    fallback_cat = category
    # Map labbench2 category names to lab-bench names
    category_map = {
        "LitQA3": "LitQA2", "FigQA2": "FigQA", "SeqQA2": "SeqQA",
        "SuppQA2": "SuppQA", "ProtocolQA2": "ProtocolQA",
        "DbQA2": "DbQA", "TableQA2": "TableQA",
        "CloningScenarios2": "CloningScenarios",
    }
    if category in category_map:
        fallback_cat = category_map[category]

    try:
        ds = load_dataset("futurehouse/lab-bench", fallback_cat, token=token, split="test")
        print(f"  Loaded from futurehouse/lab-bench ({fallback_cat})")
        return _parse_hf_dataset(ds, category, source="lab-bench")
    except Exception as e:
        print(f"  lab-bench also failed ({e})")
        return []


def _parse_hf_dataset(ds, category: str, source: str) -> list[dict]:
    """Parse HuggingFace dataset into our task format."""
    tasks = []
    for i, item in enumerate(ds):
        task_id = item.get("id", f"{category}_{i:04d}")

        # LAB-Bench format: question + choices (multiple choice, answer is a letter)
        question = item.get("question", item.get("prompt", ""))
        choices = item.get("choices", [])
        ideal = item.get("ideal", item.get("answer", ""))

        # Determine grading type
        if choices:
            verification = "programmatic"
            verification_fn = "multiple_choice"
        else:
            verification = item.get("verification", "llm-judge")
            verification_fn = "exact_match"

        # Collect remaining fields as metadata
        skip_keys = {"id", "question", "prompt", "choices", "ideal", "answer", "verification", "figures"}
        meta = {k: v for k, v in item.items() if k not in skip_keys and not hasattr(v, 'size')}
        meta["question"] = question
        meta["ideal"] = ideal
        if choices:
            meta["choices"] = choices

        task = {
            "id": task_id,
            "category": category,
            "question": question,
            "choices": choices if choices else None,
            "ideal": ideal,
            "verification": verification,
            "verification_fn": verification_fn,
            "source": source,
            "meta": meta,
        }
        tasks.append(task)
    return tasks


def load_tasks_from_dir(tasks_dir: str) -> list[dict]:
    """Load generated tasks from a local JSON directory."""
    tasks = []
    for filepath in sorted(glob.glob(os.path.join(tasks_dir, "*.json"))):
        with open(filepath) as f:
            task = json.load(f)
            tasks.append(task)
    return tasks


async def run_eval(model_key: str, tasks: list[dict], concurrency: int = 5, force: bool = False) -> list[dict]:
    """Run eval on a list of tasks. Returns results list."""
    model_id = await ensure_model(model_key)
    print(f"Model: {model_key} (id={model_id})")

    # Upsert all tasks
    for task in tasks:
        await upsert_task(task)

    sem = asyncio.Semaphore(concurrency)
    print(f"Running eval (concurrency={concurrency})...")
    results = await asyncio.gather(*[run_single(model_key, model_id, t, sem, force=force) for t in tasks])

    correct = sum(1 for r in results if r.get("correct"))
    errors = sum(1 for r in results if "error" in r)
    skipped = sum(1 for r in results if r.get("skipped"))
    total = len(results) - errors
    if total > 0:
        print(f"\nResults: {correct}/{total} correct ({correct/total*100:.1f}%)")
    else:
        print("\nNo results.")
    if skipped:
        print(f"Skipped (already done): {skipped}")
    if errors:
        print(f"Errors: {errors}")

    return results


async def main():
    parser = argparse.ArgumentParser(description="Run LABBench eval")
    parser.add_argument("--model", required=True, choices=list(MODEL_REGISTRY.keys()))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--category", help="HuggingFace dataset category")
    group.add_argument("--tasks-dir", help="Local directory with task JSON files")
    parser.add_argument("--concurrency", type=int, default=5)
    parser.add_argument("--limit", type=int, default=None, help="Limit number of tasks (for testing)")
    parser.add_argument("--force", action="store_true", help="Re-run tasks even if already evaluated")
    args = parser.parse_args()

    if args.category:
        print(f"Loading tasks for {args.category} from HuggingFace...")
        tasks = load_tasks_from_hf(args.category)
    else:
        print(f"Loading tasks from {args.tasks_dir}...")
        tasks = load_tasks_from_dir(args.tasks_dir)

    if not tasks:
        print("No tasks loaded.", file=sys.stderr)
        sys.exit(1)

    if args.limit:
        tasks = tasks[:args.limit]
    print(f"Loaded {len(tasks)} tasks")

    await run_eval(args.model, tasks, args.concurrency, args.force)


if __name__ == "__main__":
    asyncio.run(main())
