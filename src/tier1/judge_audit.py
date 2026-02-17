"""LLM judge concordance audit — position bias, verbosity bias, Cohen's kappa.

Usage: python -m src.tier1.judge_audit --category LitQA2 --sample-size 50
"""

import argparse
import asyncio
import json
import random
import sys

import numpy as np

from src.db import execute, fetch
from src.models import call_model

JUDGE_PROMPT = """You are grading a scientific benchmark response.

Question: {question}
Reference answer: {ideal}
Model response: {response}

Is the model response correct? It does not need to match word-for-word,
but must contain the key factual content of the reference answer.

Respond with JSON only: {{"correct": true, "explanation": "brief reason"}} or {{"correct": false, "explanation": "brief reason"}}"""

# Position bias test: swap reference and response positions
JUDGE_PROMPT_REVERSED = """You are grading a scientific benchmark response.

Question: {question}
Model response: {response}
Reference answer: {ideal}

Is the model response correct? It does not need to match word-for-word,
but must contain the key factual content of the reference answer.

Respond with JSON only: {{"correct": true, "explanation": "brief reason"}} or {{"correct": false, "explanation": "brief reason"}}"""

# Verbosity bias test: instruct the model to be brief
VERBOSITY_PROMPT = """Rewrite this response more concisely, keeping only the core factual claim in 1-2 sentences:

{response}

Concise version:"""


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


def parse_judge_response(raw: str) -> dict:
    """Parse a judge response, handling JSON and non-JSON."""
    try:
        parsed = json.loads(raw)
        return {"correct": parsed.get("correct", False), "explanation": parsed.get("explanation", "")}
    except (json.JSONDecodeError, KeyError):
        text = raw.lower()
        correct = '"correct": true' in text or '"correct":true' in text or "is correct" in text
        return {"correct": correct, "explanation": raw}


async def judge_response(question: str, ideal: str, response: str, model_key: str, reversed_order: bool = False) -> dict:
    """Get a judge's assessment. reversed_order=True swaps reference/response position."""
    template = JUDGE_PROMPT_REVERSED if reversed_order else JUDGE_PROMPT
    prompt = template.format(question=question, ideal=ideal, response=response)
    result = await call_model(model_key, prompt, system="You are an expert scientific evaluator.", cache=False)
    return parse_judge_response(result["response"])


async def shorten_response(response: str, model_key: str = "claude-sonnet-4.5") -> str:
    """Create a shortened version of a response for verbosity bias testing."""
    prompt = VERBOSITY_PROMPT.format(response=response)
    result = await call_model(model_key, prompt, cache=False)
    return result["response"]


async def main():
    parser = argparse.ArgumentParser(description="LLM judge audit")
    parser.add_argument("--category", required=True)
    parser.add_argument("--sample-size", type=int, default=50)
    parser.add_argument("--judge1", default="claude-sonnet-4.5")
    parser.add_argument("--judge2", default="claude-opus-4.6")
    args = parser.parse_args()

    # Get eval runs with task data (question/ideal stored in task meta)
    rows = await fetch(
        """SELECT e.id as eval_run_id, e.response, t.id as task_id, t.meta
           FROM eval_runs e
           JOIN tasks t ON e.task_id = t.id
           WHERE t.category = $1 AND e.response IS NOT NULL""",
        args.category,
    )

    if not rows:
        print(f"No eval results for category '{args.category}'", file=sys.stderr)
        sys.exit(1)

    # Parse meta to get question/ideal
    parsed_rows = []
    for row in rows:
        meta = json.loads(row["meta"]) if row["meta"] else {}
        question = meta.get("question", "")
        ideal = meta.get("ideal", "")
        if question and ideal and row["response"]:
            parsed_rows.append({
                "eval_run_id": row["eval_run_id"],
                "task_id": row["task_id"],
                "question": question,
                "ideal": ideal,
                "response": row["response"],
            })

    if not parsed_rows:
        print(f"No tasks with valid question/ideal/response found.", file=sys.stderr)
        sys.exit(1)

    sample = random.sample(parsed_rows, min(args.sample_size, len(parsed_rows)))
    print(f"Auditing {len(sample)} responses from {args.category}")
    print(f"Judges: {args.judge1} vs {args.judge2}")

    judge1_scores = []
    judge2_scores = []
    position_bias_count = 0
    verbosity_bias_diffs = []

    for i, row in enumerate(sample):
        question = row["question"]
        ideal = row["ideal"]
        response = row["response"]

        # --- Inter-judge agreement: both judges, original order ---
        j1 = await judge_response(question, ideal, response, args.judge1)
        j2 = await judge_response(question, ideal, response, args.judge2)
        judge1_scores.append(j1["correct"])
        judge2_scores.append(j2["correct"])

        # --- Position bias: same judge, reversed order ---
        j1_rev = await judge_response(question, ideal, response, args.judge1, reversed_order=True)
        if j1["correct"] != j1_rev["correct"]:
            position_bias_count += 1

        # --- Verbosity bias: shorten the response, re-judge ---
        if i < min(20, len(sample)):  # Test on subset to save cost
            try:
                short_response = await shorten_response(response)
                j1_short = await judge_response(question, ideal, short_response, args.judge1)
                verbosity_bias_diffs.append(int(j1["correct"]) - int(j1_short["correct"]))
            except Exception:
                pass

        # Save to DB
        for judge_model, score, order in [
            (args.judge1, j1, "original"),
            (args.judge2, j2, "original"),
            (args.judge1, j1_rev, "reversed"),
        ]:
            await execute(
                """INSERT INTO judge_audits (eval_run_id, judge_model, judge_score,
                   order_variant, length_variant, explanation) VALUES ($1, $2, $3, $4, $5, $6)""",
                row["eval_run_id"], judge_model, score["correct"],
                order, "original" if order == "original" else None,
                score.get("explanation", ""),
            )

        agree = "agree" if j1["correct"] == j2["correct"] else "DISAGREE"
        flip = " POSITION-FLIP" if j1["correct"] != j1_rev["correct"] else ""
        print(f"  [{i+1}/{len(sample)}] {row['task_id']}: {agree}{flip}")

    # Compute metrics
    kappa = cohens_kappa(judge1_scores, judge2_scores)
    agreement = sum(1 for a, b in zip(judge1_scores, judge2_scores) if a == b) / len(judge1_scores)
    position_bias_rate = position_bias_count / len(sample)

    print(f"\n--- Audit Results ---")
    print(f"Inter-judge agreement: {agreement:.3f}")
    print(f"Cohen's kappa: {kappa:.3f}")
    print(f"Position bias rate: {position_bias_rate:.3f} ({position_bias_count}/{len(sample)})")
    print(f"Judge 1 ({args.judge1}) positive rate: {np.mean(judge1_scores):.3f}")
    print(f"Judge 2 ({args.judge2}) positive rate: {np.mean(judge2_scores):.3f}")

    if verbosity_bias_diffs:
        mean_diff = np.mean(verbosity_bias_diffs)
        print(f"Verbosity bias (original - shortened): {mean_diff:+.3f} "
              f"(>0 = verbose answers favored, n={len(verbosity_bias_diffs)})")


if __name__ == "__main__":
    asyncio.run(main())
