"""Generate calibration/uncertainty tasks — modified LABBench2 questions with info stripped.

Usage: python -m src.tier2.gen_calibration --output-dir tasks/calibration --count 100
"""

import argparse
import json
import os

from datasets import load_dataset

HF_TOKEN = os.environ.get("HF_TOKEN", "")


def strip_critical_info(question: str) -> tuple[str, str]:
    """Remove critical information from a question and note what was removed.

    Returns (modified_question, removed_info_description).
    """
    # Simple heuristic: remove the last sentence of context before the question mark
    lines = question.strip().split("\n")

    if len(lines) < 2:
        # Single line: truncate middle portion
        words = question.split()
        if len(words) > 10:
            cut_start = len(words) // 3
            cut_end = 2 * len(words) // 3
            removed = " ".join(words[cut_start:cut_end])
            modified = " ".join(words[:cut_start]) + " [information redacted] " + " ".join(words[cut_end:])
            return modified, removed
        return question, ""

    # Multi-line: remove a paragraph from the middle
    if len(lines) >= 3:
        mid = len(lines) // 2
        removed = lines[mid]
        modified_lines = lines[:mid] + ["[Critical experimental details have been redacted]"] + lines[mid+1:]
        return "\n".join(modified_lines), removed
    else:
        removed = lines[0]
        modified_lines = ["[Key context has been redacted from this question]"] + lines[1:]
        return "\n".join(modified_lines), removed


def main():
    parser = argparse.ArgumentParser(description="Generate calibration tasks")
    parser.add_argument("--output-dir", default="tasks/calibration")
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--source-category", default="LitQA3")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    token = HF_TOKEN if HF_TOKEN else None
    try:
        ds = load_dataset("futurehouse/labbench2", args.source_category, token=token, split="test")
    except Exception as e:
        print(f"Could not load dataset: {e}")
        print("Generating placeholder calibration tasks instead...")
        ds = None

    tasks = []
    count = min(args.count, len(ds) if ds else args.count)

    for i in range(count):
        if ds and i < len(ds):
            item = ds[i]
            orig_question = item.get("question", item.get("prompt", f"Question {i}"))
            orig_answer = item.get("ideal", item.get("answer", ""))
            orig_id = item.get("id", f"{args.source_category}_{i:04d}")
        else:
            # Placeholder if dataset unavailable
            orig_question = f"Placeholder question {i} about a biological mechanism."
            orig_answer = f"Answer {i}"
            orig_id = f"placeholder_{i:04d}"

        modified_q, removed = strip_critical_info(orig_question)

        task = {
            "id": f"pro_cal_{i+1:04d}",
            "category": "Calibration",
            "question": modified_q + "\n\nBased on the information provided, answer the question. "
                        "If there is insufficient information to determine the answer, "
                        "state 'Insufficient information' and explain what is missing.",
            "ideal": "The model should recognize that critical information has been removed and express appropriate uncertainty. Correct responses include stating 'insufficient information', refusing to commit to a definitive answer, or explaining what additional data would be needed. Incorrect responses confidently provide a specific answer despite missing context.",
            "verification": "llm-judge",
            "source": "pro",
            "meta": {
                "template": "calibration",
                "original_task_id": orig_id,
                "removed_info": removed,
                "original_answer": orig_answer,
                "difficulty": "medium",
            },
        }
        tasks.append(task)

        filepath = os.path.join(args.output_dir, f"{task['id']}.json")
        with open(filepath, "w") as f:
            json.dump(task, f, indent=2)

    print(f"Generated {len(tasks)} calibration tasks in {args.output_dir}")


if __name__ == "__main__":
    main()
