"""Contamination probes — cloze, reverse, temporal split.

Usage: python -m src.tier1.contamination --model claude-opus-4.5
"""

import argparse
import asyncio
import difflib
import sys

import numpy as np
from scipy.stats import chi2_contingency

from src.db import fetch
from src.models import call_model


async def cloze_probe(model_key: str, tasks: list[dict], sample_size: int = 50) -> dict:
    """Cloze completion: can the model complete a truncated benchmark question?"""
    import random
    sample = random.sample(tasks, min(sample_size, len(tasks)))
    matches = 0

    for task in sample:
        question = task.get("question", "")
        if len(question) < 20:
            continue
        half = len(question) // 2
        first_half = question[:half]
        second_half = question[half:]

        prompt = f"Complete this scientific benchmark question: '{first_half}...'"
        result = await call_model(model_key, prompt, cache=False)
        completion = result["response"]

        similarity = difflib.SequenceMatcher(None, second_half.lower(), completion.lower()).ratio()
        if similarity > 0.8:
            matches += 1

    rate = matches / len(sample) if sample else 0
    print(f"  Cloze probe: {matches}/{len(sample)} matched (rate={rate:.3f})")
    return {"probe": "cloze", "matches": matches, "total": len(sample), "rate": rate}


async def reverse_probe(model_key: str, tasks: list[dict], sample_size: int = 50) -> dict:
    """Reverse: can the model guess the question from the answer?"""
    import random
    sample = random.sample(tasks, min(sample_size, len(tasks)))
    matches = 0

    for task in sample:
        ideal = task.get("ideal", "")
        question = task.get("question", "")
        if not ideal or not question:
            continue

        prompt = f"What question would produce this answer in a biology benchmark: '{ideal}'"
        result = await call_model(model_key, prompt, cache=False)
        generated_q = result["response"]

        similarity = difflib.SequenceMatcher(None, question.lower(), generated_q.lower()).ratio()
        if similarity > 0.6:
            matches += 1

    rate = matches / len(sample) if sample else 0
    print(f"  Reverse probe: {matches}/{len(sample)} matched (rate={rate:.3f})")
    return {"probe": "reverse", "matches": matches, "total": len(sample), "rate": rate}


async def temporal_probe(model_key: str, tasks: list[dict]) -> dict:
    """Compare accuracy on pre- vs post-cutoff tasks (by source paper date)."""
    # Split tasks by presence of date metadata
    pre_cutoff = []
    post_cutoff = []
    for task in tasks:
        meta = task.get("meta", {})
        if isinstance(meta, str):
            import json
            try:
                meta = json.loads(meta)
            except (json.JSONDecodeError, TypeError):
                meta = {}
        year = meta.get("year", meta.get("pub_year"))
        if year and int(year) < 2024:
            pre_cutoff.append(task)
        elif year and int(year) >= 2024:
            post_cutoff.append(task)

    if not pre_cutoff or not post_cutoff:
        print("  Temporal probe: insufficient date metadata for split")
        return {"probe": "temporal", "chi2": None, "p_value": None, "note": "insufficient metadata"}

    # Get accuracy from DB for these tasks
    pre_ids = [t["task_id"] for t in pre_cutoff]
    post_ids = [t["task_id"] for t in post_cutoff]

    pre_correct = sum(1 for t in pre_cutoff if t.get("correct"))
    pre_wrong = len(pre_cutoff) - pre_correct
    post_correct = sum(1 for t in post_cutoff if t.get("correct"))
    post_wrong = len(post_cutoff) - post_correct

    table = np.array([[pre_correct, pre_wrong], [post_correct, post_wrong]])

    if table.min() == 0:
        print("  Temporal probe: zero cell in contingency table")
        return {"probe": "temporal", "chi2": None, "p_value": None, "note": "zero cell"}

    chi2, p_value, dof, expected = chi2_contingency(table)
    print(f"  Temporal probe: chi2={chi2:.3f}, p={p_value:.4f}")
    print(f"    Pre-cutoff: {pre_correct}/{len(pre_cutoff)} = {pre_correct/len(pre_cutoff):.3f}")
    print(f"    Post-cutoff: {post_correct}/{len(post_cutoff)} = {post_correct/len(post_cutoff):.3f}")
    return {"probe": "temporal", "chi2": float(chi2), "p_value": float(p_value)}


async def main():
    parser = argparse.ArgumentParser(description="Contamination probes")
    parser.add_argument("--model", required=True)
    parser.add_argument("--sample-size", type=int, default=50)
    args = parser.parse_args()

    # Load tasks with their eval results for this model
    rows = await fetch(
        """SELECT t.id as task_id, t.category, t.meta,
                  e.correct, e.response
           FROM eval_runs e
           JOIN tasks t ON e.task_id = t.id
           JOIN models m ON e.model_id = m.id
           WHERE m.model_name = $1""",
        args.model,
    )

    if not rows:
        print(f"No eval results for model '{args.model}'", file=sys.stderr)
        sys.exit(1)

    tasks = [dict(row) for row in rows]
    # Parse meta from JSON string if needed
    import json
    for t in tasks:
        if isinstance(t.get("meta"), str):
            try:
                t["meta"] = json.loads(t["meta"])
            except (json.JSONDecodeError, TypeError):
                t["meta"] = {}
        t["question"] = t.get("meta", {}).get("question", "")
        t["ideal"] = t.get("meta", {}).get("ideal", "")

    print(f"Running contamination probes for {args.model} ({len(tasks)} tasks)")

    # Run all three probes
    print("\n1. Cloze completion probe:")
    cloze = await cloze_probe(args.model, tasks, args.sample_size)

    print("\n2. Reverse question probe:")
    reverse = await reverse_probe(args.model, tasks, args.sample_size)

    print("\n3. Temporal split probe:")
    temporal = await temporal_probe(args.model, tasks)

    print("\n--- Summary ---")
    print(f"Cloze match rate: {cloze['rate']:.3f}")
    print(f"Reverse match rate: {reverse['rate']:.3f}")
    if temporal.get("p_value") is not None:
        sig = "SIGNIFICANT" if temporal["p_value"] < 0.05 else "not significant"
        print(f"Temporal chi2: {temporal['chi2']:.3f} (p={temporal['p_value']:.4f}, {sig})")


if __name__ == "__main__":
    asyncio.run(main())
