"""Re-grade hypothesis generation tasks with a strict rubric.

The original LLM-judge rubric produced 100% accuracy because it checked whether
the response "contained the key factual content" of a rubric description, not
whether the hypotheses actually met the rubric criteria. This script applies a
rigorous per-criterion evaluation.

Usage:
  python -m src.tier1.regrade_hypothesis
"""

import asyncio
import json
import re
import sys

from src.db import execute, fetch, fetchrow
from src.models import call_model

STRICT_JUDGE_PROMPT = """You are a senior scientific reviewer evaluating AI-generated hypotheses.

A model was given this research abstract and asked to propose 3 testable hypotheses:

--- ABSTRACT ---
{abstract}
--- END ABSTRACT ---

--- MODEL RESPONSE ---
{response}
--- END RESPONSE ---

Evaluate the response on EACH of these 4 criteria INDEPENDENTLY. Be strict — most responses should NOT get full marks.

## Criterion 1: Falsifiable
Each hypothesis must describe a CONCRETE experiment that could disprove it.
- FAIL if hypotheses say "this could be investigated" or "future studies should examine" without specifying an experimental approach
- FAIL if the experiment described is vague (e.g., "use CRISPR" without specifying what gene, what cell type, what readout)
- PASS only if at least 2 of 3 hypotheses include a specific experimental design with method, subject, and measurable outcome

## Criterion 2: Connected
Hypotheses must logically follow from specific findings in the abstract.
- FAIL if hypotheses are generic to the field but not specifically motivated by THIS abstract's findings
- FAIL if the logical chain from finding to hypothesis requires unstated assumptions
- PASS only if at least 2 of 3 hypotheses explicitly reference a specific finding from the abstract and explain the logical connection

## Criterion 3: Novel
Hypotheses must go BEYOND restating or inverting the abstract's conclusions.
- FAIL if any hypothesis simply restates what the abstract already says or suggests
- FAIL if hypotheses propose what the abstract's own "future work" section suggests
- FAIL if hypotheses are obvious next steps that any scientist in the field would immediately suggest
- PASS only if at least 2 of 3 hypotheses propose genuinely non-obvious extensions

## Criterion 4: Specific
Hypotheses must name concrete biological entities.
- FAIL if hypotheses refer to vague categories ("certain genes", "relevant pathways", "key proteins")
- FAIL if specificity comes only from entities already named in the abstract (no new entities)
- PASS only if at least 2 of 3 hypotheses name specific genes, proteins, pathways, cell types, or mechanisms NOT already in the abstract

---

Score each criterion as PASS or FAIL. A response needs 3/4 PASS to be marked correct.

Respond with JSON only:
{{"criterion_1_falsifiable": true/false, "criterion_1_reason": "...", "criterion_2_connected": true/false, "criterion_2_reason": "...", "criterion_3_novel": true/false, "criterion_3_reason": "...", "criterion_4_specific": true/false, "criterion_4_reason": "...", "total_pass": <0-4>, "correct": true/false}}"""


async def regrade():
    # Get all hypothesis eval runs with their responses and task metadata
    rows = await fetch("""
        SELECT e.id as eval_id, e.task_id, e.response, e.correct as old_correct,
               e.score as old_score, t.meta as task_meta
        FROM eval_runs e
        JOIN tasks t ON t.id = e.task_id
        WHERE t.category = 'HypothesisGeneration'
        ORDER BY e.task_id
    """)

    print(f"Re-grading {len(rows)} hypothesis responses with strict rubric")
    print(f"Original: {sum(1 for r in rows if r['old_correct'])}/{len(rows)} correct")
    print()

    new_correct = 0
    new_wrong = 0
    criteria_pass_counts = {1: 0, 2: 0, 3: 0, 4: 0}
    results = []

    sem = asyncio.Semaphore(3)

    async def grade_one(row):
        nonlocal new_correct, new_wrong
        async with sem:
            meta = json.loads(row["task_meta"]) if isinstance(row["task_meta"], str) else row["task_meta"]

            # Extract the abstract from the question
            question = meta.get("question", "")
            # The abstract is in the question text after "Abstract:"
            abstract = question
            if "Abstract:" in question:
                abstract = question.split("Abstract:", 1)[1].strip()
            elif "abstract" in question.lower():
                abstract = question

            prompt = STRICT_JUDGE_PROMPT.format(
                abstract=abstract[:3000],  # Truncate very long abstracts
                response=row["response"][:4000],  # Truncate very long responses
            )

            try:
                result = await call_model(
                    "claude-sonnet-4.5", prompt,
                    system="You are a strict scientific reviewer. Most AI-generated hypotheses are too vague or obvious to pass rigorous review. Respond with JSON only, no markdown.",
                    cache=False,
                )
                resp_text = result["response"].strip()
                # Strip markdown code fences if present
                if resp_text.startswith("```"):
                    resp_text = re.sub(r"^```(?:json)?\s*", "", resp_text)
                    resp_text = re.sub(r"\s*```$", "", resp_text)
                parsed = json.loads(resp_text)
            except json.JSONDecodeError as exc:
                # Try to extract JSON from mixed text
                import re as re_mod
                match = re_mod.search(r'\{[^{}]*"correct"[^{}]*\}', result["response"], re_mod.DOTALL)
                if match:
                    try:
                        parsed = json.loads(match.group())
                    except json.JSONDecodeError:
                        print(f"  ERROR {row['task_id']}: could not parse JSON")
                        return
                else:
                    print(f"  ERROR {row['task_id']}: no JSON found in response")
                    return
            except Exception as exc:
                print(f"  ERROR {row['task_id']}: {exc}")
                return

            is_correct = parsed.get("correct", False)
            total_pass = parsed.get("total_pass", 0)

            for i in range(1, 5):
                key = f"criterion_{i}_{'falsifiable' if i==1 else 'connected' if i==2 else 'novel' if i==3 else 'specific'}"
                if parsed.get(key, False):
                    criteria_pass_counts[i] += 1

            # Update the eval run with new grade
            await execute("""
                UPDATE eval_runs
                SET correct = $1, score = $2, grader = 'llm-judge-strict-rubric',
                    meta = jsonb_set(
                        COALESCE(meta, '{}'::jsonb),
                        '{strict_regrade}',
                        $3::jsonb
                    )
                WHERE id = $4
            """,
                is_correct,
                total_pass / 4.0,
                json.dumps({
                    "criteria": {
                        "falsifiable": parsed.get("criterion_1_falsifiable", False),
                        "connected": parsed.get("criterion_2_connected", False),
                        "novel": parsed.get("criterion_3_novel", False),
                        "specific": parsed.get("criterion_4_specific", False),
                    },
                    "reasons": {
                        "falsifiable": parsed.get("criterion_1_reason", ""),
                        "connected": parsed.get("criterion_2_reason", ""),
                        "novel": parsed.get("criterion_3_reason", ""),
                        "specific": parsed.get("criterion_4_reason", ""),
                    },
                    "total_pass": total_pass,
                    "old_correct": row["old_correct"],
                }),
                row["eval_id"],
            )

            status = "CORRECT" if is_correct else "WRONG"
            if is_correct:
                new_correct += 1
            else:
                new_wrong += 1

            print(f"  [{status}] {row['task_id']}: {total_pass}/4 criteria passed "
                  f"(F:{'Y' if parsed.get('criterion_1_falsifiable') else 'N'} "
                  f"C:{'Y' if parsed.get('criterion_2_connected') else 'N'} "
                  f"N:{'Y' if parsed.get('criterion_3_novel') else 'N'} "
                  f"S:{'Y' if parsed.get('criterion_4_specific') else 'N'})")

            results.append({
                "task_id": row["task_id"],
                "old_correct": row["old_correct"],
                "new_correct": is_correct,
                "total_pass": total_pass,
                "cost": result.get("cost_usd", 0),
            })

    await asyncio.gather(*[grade_one(row) for row in rows])

    total = new_correct + new_wrong
    print(f"\n{'='*60}")
    print(f"STRICT RE-GRADING RESULTS")
    print(f"{'='*60}")
    print(f"Original: 100/100 correct (100.0%) — lenient rubric")
    print(f"Strict:   {new_correct}/{total} correct ({new_correct/total*100:.1f}%) — strict rubric")
    print(f"Drop:     {100 - new_correct} responses failed strict criteria")
    print(f"\nPer-criterion pass rates:")
    print(f"  1. Falsifiable:  {criteria_pass_counts[1]}/{total} ({criteria_pass_counts[1]/total*100:.1f}%)")
    print(f"  2. Connected:    {criteria_pass_counts[2]}/{total} ({criteria_pass_counts[2]/total*100:.1f}%)")
    print(f"  3. Novel:        {criteria_pass_counts[3]}/{total} ({criteria_pass_counts[3]/total*100:.1f}%)")
    print(f"  4. Specific:     {criteria_pass_counts[4]}/{total} ({criteria_pass_counts[4]/total*100:.1f}%)")

    total_cost = sum(r.get("cost", 0) for r in results)
    print(f"\nRe-grading cost: ${total_cost:.2f}")


if __name__ == "__main__":
    asyncio.run(regrade())
