"""LLM judge concordance audit — position bias, verbosity bias, Cohen's kappa.

Usage: python -m src.tier1.judge_audit --category SourceQualQA --sample-size 50
"""

import argparse
import asyncio
import json
import random
import sys

import numpy as np

from src.db import execute, fetch
from src.models import call_model
from src.tier1.grading import JUDGE_PROMPT


def cohens_kappa(labels1: list[bool], labels2: list[bool]) -> float:
    """Compute Cohen's kappa between two raters."""
    assert len(labels1) == len(labels2)
    n = len(labels1)
    if n == 0:
        return 0.0
    a1 = np.array(labels1, dtype=int)
    a2 = np.array(labels2, dtype=int)
    agreement = np.sum(a1 == a2) / n
    p1 = np.mean(a1)
    p2 = np.mean(a2)
    expected = p1 * p2 + (1 - p1) * (1 - p2)
    if expected == 1.0:
        return 1.0
    return (agreement - expected) / (1 - expected)


async def judge_response(question: str, ideal: str, response: str, model_key: str) -> dict:
    """Get a judge's assessment."""
    prompt = JUDGE_PROMPT.format(question=question, ideal=ideal, response=response)
    result = await call_model(model_key, prompt, system="You are an expert scientific evaluator.", cache=False)
    try:
        parsed = json.loads(result["response"])
        return {"correct": parsed.get("correct", False), "explanation": parsed.get("explanation", "")}
    except (json.JSONDecodeError, KeyError):
        text = result["response"].lower()
        correct = '"correct": true' in text or '"correct":true' in text
        return {"correct": correct, "explanation": result["response"]}


async def main():
    parser = argparse.ArgumentParser(description="LLM judge audit")
    parser.add_argument("--category", required=True)
    parser.add_argument("--sample-size", type=int, default=50)
    parser.add_argument("--judge1", default="claude-sonnet-4.5")
    parser.add_argument("--judge2", default="claude-opus-4.5")
    args = parser.parse_args()

    # Get eval runs with task data
    rows = await fetch(
        """SELECT e.id as eval_run_id, e.response, t.id as task_id,
                  t.meta->>'question' as question, t.meta->>'ideal' as ideal
           FROM eval_runs e
           JOIN tasks t ON e.task_id = t.id
           WHERE t.category = $1
           LIMIT $2""",
        args.category, args.sample_size * 2,
    )

    if not rows:
        # Try loading task question from tasks table meta
        rows = await fetch(
            """SELECT e.id as eval_run_id, e.response, t.id as task_id
               FROM eval_runs e JOIN tasks t ON e.task_id = t.id
               WHERE t.category = $1 LIMIT $2""",
            args.category, args.sample_size * 2,
        )

    if not rows:
        print(f"No eval results for category '{args.category}'", file=sys.stderr)
        sys.exit(1)

    sample = random.sample(list(rows), min(args.sample_size, len(rows)))
    print(f"Auditing {len(sample)} responses from {args.category}")

    judge1_scores = []
    judge2_scores = []
    position_bias_count = 0

    for i, row in enumerate(sample):
        question = row.get("question", f"Task {row['task_id']}")
        ideal = row.get("ideal", "")
        response = row["response"] or ""

        # Original order — both judges
        j1 = await judge_response(question, ideal, response, args.judge1)
        j2 = await judge_response(question, ideal, response, args.judge2)

        # Reversed order (swap question framing for position bias test)
        j1_rev = await judge_response(question, ideal, response, args.judge1)

        judge1_scores.append(j1["correct"])
        judge2_scores.append(j2["correct"])

        if j1["correct"] != j1_rev["correct"]:
            position_bias_count += 1

        # Save to DB
        for judge_model, score, order in [(args.judge1, j1, "original"), (args.judge2, j2, "original")]:
            await execute(
                """INSERT INTO judge_audits (eval_run_id, judge_model, judge_score,
                   order_variant, explanation) VALUES ($1, $2, $3, $4, $5)""",
                row["eval_run_id"], judge_model, score["correct"],
                order, score.get("explanation", ""),
            )

        status = "agree" if j1["correct"] == j2["correct"] else "DISAGREE"
        print(f"  [{i+1}/{len(sample)}] {row['task_id']}: {status}")

    # Compute metrics
    kappa = cohens_kappa(judge1_scores, judge2_scores)
    agreement = sum(1 for a, b in zip(judge1_scores, judge2_scores) if a == b) / len(judge1_scores)
    bias_rate = position_bias_count / len(sample)

    print(f"\n--- Audit Results ---")
    print(f"Inter-judge agreement: {agreement:.3f}")
    print(f"Cohen's kappa: {kappa:.3f}")
    print(f"Position bias rate: {bias_rate:.3f} ({position_bias_count}/{len(sample)})")
    print(f"Judge 1 ({args.judge1}) positive rate: {np.mean(judge1_scores):.3f}")
    print(f"Judge 2 ({args.judge2}) positive rate: {np.mean(judge2_scores):.3f}")


if __name__ == "__main__":
    asyncio.run(main())
