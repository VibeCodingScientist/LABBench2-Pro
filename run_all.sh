#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# LABBench2-Pro — Fully Automated Pipeline
#
# One command: ./run_all.sh --model claude-opus-4.6
#
# Prerequisites:
#   1. Docker running
#   2. .env file with ANTHROPIC_API_KEY (and HF_TOKEN for LABBench2)
#   3. pip install -e .
###############################################################################

MODEL="${MODEL:-claude-opus-4.6}"
CONCURRENCY="${CONCURRENCY:-5}"
SKIP_DOCKER="${SKIP_DOCKER:-false}"

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --model) MODEL="$2"; shift 2 ;;
        --concurrency) CONCURRENCY="$2"; shift 2 ;;
        --skip-docker) SKIP_DOCKER=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

log() { echo ""; echo "========== $1 =========="; echo ""; }
timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

echo "LABBench2-Pro Pipeline"
echo "Model: $MODEL | Concurrency: $CONCURRENCY"
echo "Started: $(timestamp)"
echo ""

###############################################################################
# Phase 0: Infrastructure
###############################################################################
log "Phase 0: Infrastructure"

if [ "$SKIP_DOCKER" = false ]; then
    echo "Starting Postgres + Redis..."
    docker compose up -d
    echo "Waiting for services to be healthy..."
    sleep 5
    # Wait for Postgres
    for i in $(seq 1 30); do
        if docker compose exec -T postgres pg_isready -U dev -d labbench2pro > /dev/null 2>&1; then
            echo "Postgres ready."
            break
        fi
        sleep 1
    done
fi

# Schema is auto-applied via docker-entrypoint-initdb.d, but re-apply if tables missing
echo "Ensuring schema..."
PGPASSWORD=dev psql -h localhost -U dev -d labbench2pro -f db/schema.sql 2>/dev/null || true

###############################################################################
# Phase 1: Generate Tier 2 Tasks (local, no API cost)
###############################################################################
log "Phase 1: Generate Tier 2 Tasks"

echo "[1/4] Statistical reasoning tasks..."
python -m src.tier2.gen_stats_tasks --output-dir tasks/stats_reasoning --count 200

echo "[2/4] Structure analysis tasks (PDB download + gel images)..."
python -m src.tier2.gen_structure --output-dir tasks/structures --pdb-count 50 --gel-count 60

echo "[3/4] Hypothesis generation tasks (PubMed fetch)..."
python -m src.tier2.gen_hypothesis --output-dir tasks/hypothesis --count 100

echo "[4/4] Calibration tasks..."
python -m src.tier2.gen_calibration --output-dir tasks/calibration --count 100

echo ""
echo "Validating all generated tasks..."
python -m src.tier2.validate_tasks

###############################################################################
# Phase 2: Run Tier 1 Evals (LABBench from HuggingFace)
###############################################################################
log "Phase 2: Tier 1 — LABBench Evals"

# Try all available categories
CATEGORIES="LitQA2 FigQA SeqQA SuppQA ProtocolQA DbQA TableQA CloningScenarios"
for CAT in $CATEGORIES; do
    echo ""
    echo "--- $CAT ---"
    python -m src.tier1.run_eval --model "$MODEL" --category "$CAT" --concurrency "$CONCURRENCY" || \
        echo "  Warning: $CAT failed or unavailable, continuing..."
done

###############################################################################
# Phase 3: Run Tier 2 Generated Tasks Through Eval
###############################################################################
log "Phase 3: Tier 2 — Generated Task Evals"

TASK_DIRS="tasks/stats_reasoning tasks/structures tasks/calibration tasks/hypothesis"
for DIR in $TASK_DIRS; do
    if [ -d "$DIR" ] && ls "$DIR"/*.json 1>/dev/null 2>&1; then
        echo ""
        echo "--- $DIR ---"
        python -m src.tier1.run_eval --model "$MODEL" --tasks-dir "$DIR" --concurrency "$CONCURRENCY"
    else
        echo "  Skipping $DIR (no tasks)"
    fi
done

###############################################################################
# Phase 4: Auto-generate and Run Chains
###############################################################################
log "Phase 4: Tier 3 — Compositional Chains"

echo "Auto-generating chain definitions from DB tasks..."
python -m src.tier3.gen_chains

echo "Running chains..."
CHAIN_FILE="tasks/chains/chain_definitions.json"
if [ -f "$CHAIN_FILE" ]; then
    # Extract chain IDs and run each
    CHAIN_IDS=$(python -c "import json; chains=json.load(open('$CHAIN_FILE')); [print(c['chain_id']) for c in chains]")
    for CHAIN_ID in $CHAIN_IDS; do
        echo ""
        echo "--- Chain: $CHAIN_ID ---"
        python -m src.tier3.run_chains --model "$MODEL" --chain "$CHAIN_ID" || \
            echo "  Warning: chain $CHAIN_ID failed, continuing..."
    done
fi

###############################################################################
# Phase 5: Analysis
###############################################################################
log "Phase 5: Analysis"

echo "[1/4] Bootstrap CIs..."
for CAT in $CATEGORIES; do
    python -m src.tier1.bootstrap_ci --category "$CAT" 2>/dev/null || true
done

echo ""
echo "[2/4] IRT analysis..."
python -m src.tier1.irt_analysis 2>/dev/null || echo "  (needs 2+ models for IRT)"

echo ""
echo "[3/4] Judge audit (sample of 20)..."
python -m src.tier1.judge_audit --category LitQA2 --sample-size 20 2>/dev/null || \
    echo "  (needs eval results first)"

echo ""
echo "[4/4] Contamination probes..."
python -m src.tier1.contamination --model "$MODEL" --sample-size 20 2>/dev/null || \
    echo "  (needs eval results first)"

###############################################################################
# Phase 6: Cost Summary
###############################################################################
log "Phase 6: Cost Summary"
python -m src.tier3.cost_tracker

###############################################################################
# Done
###############################################################################
log "PIPELINE COMPLETE"
echo "Finished: $(timestamp)"
echo ""
echo "Results are in PostgreSQL. Query with:"
echo "  PGPASSWORD=dev psql -h localhost -U dev -d labbench2pro"
echo ""
echo "Or start the API:"
echo "  uvicorn src.api:app --host 0.0.0.0 --port 8000"
