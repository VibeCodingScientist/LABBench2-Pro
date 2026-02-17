"""Validate generated tasks — schema check, ground truth verification.

Usage: python -m src.tier2.validate_tasks
"""

import argparse
import json
import os
import sys

REQUIRED_FIELDS = {"id", "category", "question", "ideal", "verification", "source"}
VALID_VERIFICATIONS = {"programmatic", "llm-judge"}
VALID_SOURCES = {"pro", "labbench2"}

TASK_DIRS = [
    "tasks/stats_reasoning",
    "tasks/structures",
    "tasks/calibration",
    "tasks/hypothesis",
]


def validate_task(filepath: str) -> list[str]:
    """Validate a single task file. Returns list of errors."""
    errors = []

    try:
        with open(filepath) as f:
            task = json.load(f)
    except json.JSONDecodeError as e:
        return [f"Invalid JSON: {e}"]

    if not isinstance(task, dict):
        return ["Task must be a JSON object"]

    # Required fields
    missing = REQUIRED_FIELDS - set(task.keys())
    if missing:
        errors.append(f"Missing fields: {missing}")

    # Field types
    if "id" in task and not isinstance(task["id"], str):
        errors.append(f"id must be string, got {type(task['id']).__name__}")
    if "question" in task and not isinstance(task["question"], str):
        errors.append(f"question must be string, got {type(task['question']).__name__}")
    if "ideal" in task and not isinstance(task["ideal"], str):
        errors.append(f"ideal must be string, got {type(task['ideal']).__name__}")

    # Verification type
    if "verification" in task and task["verification"] not in VALID_VERIFICATIONS:
        errors.append(f"Invalid verification: {task['verification']}")

    # Source
    if "source" in task and task["source"] not in VALID_SOURCES:
        errors.append(f"Invalid source: {task['source']}")

    # Non-empty question and ideal
    if task.get("question", "").strip() == "":
        errors.append("Empty question")
    if task.get("ideal", "").strip() == "":
        errors.append("Empty ideal answer")

    # Meta should be a dict if present
    if "meta" in task and not isinstance(task["meta"], dict):
        errors.append(f"meta must be dict, got {type(task['meta']).__name__}")

    return errors


def main():
    parser = argparse.ArgumentParser(description="Validate generated tasks")
    parser.add_argument("--dirs", nargs="*", default=TASK_DIRS)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    total = 0
    passed = 0
    failed = 0
    all_errors = []

    for task_dir in args.dirs:
        if not os.path.isdir(task_dir):
            print(f"  SKIP {task_dir} (not found)")
            continue

        dir_total = 0
        dir_passed = 0

        for filename in sorted(os.listdir(task_dir)):
            if not filename.endswith(".json"):
                continue

            filepath = os.path.join(task_dir, filename)
            errors = validate_task(filepath)
            total += 1
            dir_total += 1

            if errors:
                failed += 1
                all_errors.append((filepath, errors))
                if args.verbose:
                    print(f"  FAIL {filepath}: {errors}")
            else:
                passed += 1
                dir_passed += 1

        if dir_total > 0:
            print(f"  {task_dir}: {dir_passed}/{dir_total} passed")

    print(f"\nTotal: {passed}/{total} passed, {failed} failed")

    if all_errors and not args.verbose:
        print(f"\nFirst 5 failures:")
        for filepath, errors in all_errors[:5]:
            print(f"  {filepath}: {errors}")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
