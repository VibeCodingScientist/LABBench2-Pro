"""Auto-generate compositional chains from tasks in the database.

Composes multi-step chains by picking tasks from different categories
that form logical research workflows.

Usage: python -m src.tier3.gen_chains
"""

import asyncio
import json
import os
import random

from src.db import fetch

# Chain templates: each defines a category sequence that forms a logical workflow
CHAIN_TEMPLATES = [
    {
        "name": "lit_to_stats",
        "description": "Literature comprehension -> statistical reasoning",
        "categories": ["LitQA2", "StatisticalReasoning"],
    },
    {
        "name": "struct_to_hypothesis",
        "description": "Structure analysis -> hypothesis generation",
        "categories": ["StructureAnalysis", "HypothesisGeneration"],
    },
    {
        "name": "stats_pipeline",
        "description": "Test selection -> p-value computation -> significance count",
        "categories": ["StatisticalReasoning", "StatisticalReasoning", "StatisticalReasoning"],
    },
    {
        "name": "lit_to_struct_to_hyp",
        "description": "Literature -> structure analysis -> hypothesis",
        "categories": ["LitQA2", "StructureAnalysis", "HypothesisGeneration"],
    },
    {
        "name": "calibration_chain",
        "description": "Calibration -> literature -> calibration (uncertainty reasoning)",
        "categories": ["Calibration", "LitQA2", "Calibration"],
    },
    {
        "name": "full_pipeline",
        "description": "Literature -> stats -> structure -> hypothesis",
        "categories": ["LitQA2", "StatisticalReasoning", "StructureAnalysis", "HypothesisGeneration"],
    },
]


async def get_tasks_by_category() -> dict[str, list[str]]:
    """Get all task IDs grouped by category from the database."""
    rows = await fetch("SELECT id, category FROM tasks ORDER BY category, id")
    by_cat: dict[str, list[str]] = {}
    for row in rows:
        by_cat.setdefault(row["category"], []).append(row["id"])
    return by_cat


async def generate_chains(chains_per_template: int = 5) -> list[dict]:
    """Generate chain definitions from tasks in the DB."""
    tasks_by_cat = await get_tasks_by_category()

    if not tasks_by_cat:
        print("No tasks in database. Run task generation and eval first.")
        return []

    print(f"Tasks in DB by category:")
    for cat, ids in sorted(tasks_by_cat.items()):
        print(f"  {cat}: {len(ids)} tasks")

    chains = []
    chain_num = 0

    for template in CHAIN_TEMPLATES:
        # Check all categories in this template have tasks
        available = all(cat in tasks_by_cat for cat in template["categories"])
        if not available:
            missing = [c for c in template["categories"] if c not in tasks_by_cat]
            print(f"  Skipping '{template['name']}' — missing categories: {missing}")
            continue

        for i in range(chains_per_template):
            chain_num += 1
            steps = []
            used_ids = set()

            for j, cat in enumerate(template["categories"]):
                # Pick a random task from this category (avoid repeats within a chain)
                candidates = [t for t in tasks_by_cat[cat] if t not in used_ids]
                if not candidates:
                    candidates = tasks_by_cat[cat]  # allow repeats if exhausted

                task_id = random.choice(candidates)
                used_ids.add(task_id)

                steps.append({
                    "task_id": task_id,
                    "input_type": "none" if j == 0 else "previous_response",
                })

            chains.append({
                "chain_id": f"{template['name']}_{chain_num:03d}",
                "description": template["description"],
                "steps": steps,
            })

    return chains


async def main():
    random.seed(42)  # Reproducible chain composition

    chains = await generate_chains(chains_per_template=5)

    if not chains:
        print("No chains generated.")
        return

    output_path = os.path.join(os.path.dirname(__file__), "../../tasks/chains/chain_definitions.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(chains, f, indent=2)

    print(f"\nGenerated {len(chains)} chains -> {output_path}")
    for chain in chains:
        step_ids = [s["task_id"] for s in chain["steps"]]
        print(f"  {chain['chain_id']}: {' -> '.join(step_ids)}")


if __name__ == "__main__":
    asyncio.run(main())
