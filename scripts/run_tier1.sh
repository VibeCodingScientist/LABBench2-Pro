#!/usr/bin/env bash
set -euo pipefail

# Convenience script: run all Tier 1 analysis for a given model and category.
# Usage: ./scripts/run_tier1.sh --model claude-opus-4.5 --category LitQA3

MODEL="${1:?Usage: $0 --model MODEL --category CATEGORY}"
shift
CATEGORY="${1:?Usage: $0 --model MODEL --category CATEGORY}"
shift

echo "=== Tier 1: ${MODEL} on ${CATEGORY} ==="

echo ""
echo "--- Step 1: Run eval ---"
python -m src.tier1.run_eval --model "$MODEL" --category "$CATEGORY" --concurrency 5

echo ""
echo "--- Step 2: Bootstrap CIs ---"
python -m src.tier1.bootstrap_ci --category "$CATEGORY" --output-csv "results_${CATEGORY}_ci.csv"

echo ""
echo "--- Step 3: IRT analysis ---"
python -m src.tier1.irt_analysis --output-csv "results_irt.csv"

echo ""
echo "--- Step 4: Judge audit ---"
python -m src.tier1.judge_audit --category "$CATEGORY" --sample-size 20

echo ""
echo "--- Step 5: Contamination probes ---"
python -m src.tier1.contamination --model "$MODEL" --sample-size 20

echo ""
echo "=== Tier 1 complete ==="
