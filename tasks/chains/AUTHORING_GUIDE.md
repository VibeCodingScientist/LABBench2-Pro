# Chain Authoring Guide

## How Chains Work

Each chain is a multi-step task where the model's answer to step N is prepended
as context to step N+1. This tests **compositional reasoning** — can the model
carry information forward and build on it?

## What Makes a Good Chain

1. **Each step genuinely depends on the previous.** If you can answer step 3
   without step 2's output, the chain is too loose.

2. **Wrong answers should cascade.** If the model misidentifies a protein in
   step 1, it should get step 2 wrong too. This measures error propagation.

3. **Steps should cross categories.** Lit → Stats → Hypothesis is better than
   Lit → Lit → Lit because it tests transfer across reasoning types.

## How to Fill In Templates

Each task file has `[FILL: description]` markers. Replace them with real content:

```json
{
  "question": "[FILL: A question about a specific paper finding]",
  "ideal": "[FILL: The correct answer]"
}
```

becomes:

```json
{
  "question": "A 2023 Nature study found that TRIM28 ...",
  "ideal": "TRIM28 acts as a transcriptional repressor by ..."
}
```

## Verification Types

- `"programmatic"` + `"verification_fn": "exact_match"` — for factual answers (numbers, names)
- `"programmatic"` + `"verification_fn": "numeric_tolerance"` — for numerical answers (p-values, concentrations)
- `"programmatic"` + `"verification_fn": "integer_match"` — for counts
- `"programmatic"` + `"verification_fn": "multiple_choice"` — for A/B/C/D answers
- `"llm-judge"` — for open-ended answers (hypotheses, explanations, experimental designs)

## Testing a Chain

After filling in tasks:

```bash
python -m src.tier1.run_eval --model claude-opus-4.6 --tasks-dir tasks/chains/tasks --concurrency 3
python -m src.tier3.run_chains --model claude-opus-4.6 --chain <chain_id>
```
