"""Grading module — programmatic and LLM-judge grading."""

import json
import re

from src.models import call_model

JUDGE_PROMPT = """You are grading a scientific benchmark response.

Question: {question}
Reference answer: {ideal}
Model response: {response}

Is the model response correct? It does not need to match word-for-word,
but must contain the key factual content of the reference answer.

Respond with JSON only: {{"correct": true, "explanation": "brief reason"}} or {{"correct": false, "explanation": "brief reason"}}"""


def grade_programmatic(task: dict, response: str) -> dict:
    """Grade by exact match, numeric tolerance, multiple choice, or custom verifier."""
    ideal = task.get("ideal", "")
    verification = task.get("verification_fn", "exact_match")

    if verification == "multiple_choice":
        return _grade_multiple_choice(ideal, response)
    elif verification == "numeric_tolerance":
        return _grade_numeric(ideal, response, tolerance=task.get("meta", {}).get("tolerance", 0.01))
    elif verification == "integer_match":
        return _grade_integer(ideal, response)
    else:
        return _grade_exact(ideal, response)


def _grade_multiple_choice(ideal: str, response: str) -> dict:
    """Grade multiple choice: extract the letter from the response and compare."""
    ideal_letter = ideal.strip().upper()[:1]
    # Extract the first standalone letter from the response
    response_clean = response.strip()
    # Try: starts with a letter, or "The answer is X", or just a single letter
    match = re.search(r'\b([A-Z])\b', response_clean.upper())
    if match:
        response_letter = match.group(1)
    else:
        response_letter = response_clean.upper()[:1]
    correct = ideal_letter == response_letter
    return {"correct": correct, "score": 1.0 if correct else 0.0, "grader": "programmatic_mc"}


def _grade_exact(ideal: str, response: str) -> dict:
    ideal_clean = ideal.strip().lower()
    response_clean = response.strip().lower()
    # Check if ideal answer is contained in response
    correct = ideal_clean in response_clean
    return {"correct": correct, "score": 1.0 if correct else 0.0, "grader": "programmatic_exact"}


def _grade_numeric(ideal: str, response: str, tolerance: float = 0.01) -> dict:
    try:
        ideal_val = float(re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", ideal).group())
        response_val = float(re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", response).group())
        correct = abs(ideal_val - response_val) <= tolerance
        return {"correct": correct, "score": 1.0 if correct else 0.0, "grader": "programmatic_numeric"}
    except (AttributeError, ValueError):
        return {"correct": False, "score": 0.0, "grader": "programmatic_numeric"}


def _grade_integer(ideal: str, response: str) -> dict:
    try:
        ideal_int = int(re.search(r"\d+", ideal).group())
        response_int = int(re.search(r"\d+", response).group())
        correct = ideal_int == response_int
        return {"correct": correct, "score": 1.0 if correct else 0.0, "grader": "programmatic_integer"}
    except (AttributeError, ValueError):
        return {"correct": False, "score": 0.0, "grader": "programmatic_integer"}


async def grade_llm_judge(task: dict, response: str, model_key: str = "claude-sonnet-4.5") -> dict:
    """Grade using an LLM judge."""
    prompt = JUDGE_PROMPT.format(
        question=task.get("question", ""),
        ideal=task.get("ideal", ""),
        response=response,
    )
    result = await call_model(model_key, prompt, system="You are an expert scientific evaluator.", cache=False)
    try:
        parsed = json.loads(result["response"])
        return {
            "correct": parsed.get("correct", False),
            "score": 1.0 if parsed.get("correct") else 0.0,
            "grader": f"llm-judge-{model_key}",
            "explanation": parsed.get("explanation", ""),
        }
    except (json.JSONDecodeError, KeyError):
        # Try to extract from non-JSON response
        text = result["response"].lower()
        correct = '"correct": true' in text or '"correct":true' in text
        return {
            "correct": correct,
            "score": 1.0 if correct else 0.0,
            "grader": f"llm-judge-{model_key}",
            "explanation": result["response"],
        }


async def grade(task: dict, response: str) -> dict:
    """Dispatch to the right grader based on task verification type."""
    verification = task.get("verification", "llm-judge")
    if verification == "programmatic":
        return grade_programmatic(task, response)
    else:
        return await grade_llm_judge(task, response)
