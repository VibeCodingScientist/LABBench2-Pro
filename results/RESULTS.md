# LABBench2-Pro Results

## Run Details

| Parameter | Value |
|---|---|
| Model | Claude Opus 4.6 (claude-opus-4-6) |
| Provider | Anthropic |
| Date | 2026-02-18 |
| Runtime | 5h 5min (07:29 - 12:34 UTC) |
| Total cost | $102.81 |
| Total eval runs | 2,454 |
| Total tokens (in) | 766,478 |
| Total tokens (out) | 1,217,518 |
| Infrastructure | PostgreSQL 16 + Redis 7, DigitalOcean 4vCPU/8GB |
| Concurrency | 3 parallel API calls |

## Tier 1: LABBench2 Categories

Accuracy with BCa 95% bootstrap confidence intervals (10,000 resamples):

| Category | n | Correct | Accuracy | 95% CI | Cost |
|---|---|---|---|---|---|
| CloningScenarios | 33 | 13 | **39.4%** | [24.2%, 57.6%] | $5.25 |
| LitQA2 | 186 | 58 | **31.2%** | [24.7%, 38.2%] | $3.72 |
| SeqQA | 585 | 104 | **17.8%** | [14.9%, 21.0%] | $38.62 |
| ProtocolQA | 100 | 15 | **15.0%** | [9.0%, 23.0%] | $3.16 |
| SuppQA | 81 | 9 | **11.1%** | [4.9%, 19.8%] | $1.15 |
| FigQA | 181 | 19 | **10.5%** | [6.6%, 15.5%] | $2.84 |
| DbQA | 492 | 23 | **4.7%** | [3.0%, 6.7%] | $13.20 |
| TableQA | — | — | skipped | — | — |

TableQA was skipped (image loading issue with the HuggingFace dataset).

## Tier 2: New Task Categories

| Category | n | Correct | Accuracy | Verification | Cost |
|---|---|---|---|---|---|
| Calibration | 100 | 100 | **100.0%** | LLM-judge | $1.11 |
| Hypothesis Generation | 100 | 100 | **100.0%** | LLM-judge | $9.95 |
| Structure Analysis | 303 | 138 | **45.5%** | Programmatic | $2.30 |
| Statistical Reasoning | 200 | 46 | **23.0%** | Programmatic | $11.22 |

**Note:** Calibration and Hypothesis Generation are LLM-judge graded. The 100% accuracy should be interpreted with caution—see Judge Audit below.

## Tier 3: Compositional Chains

### Summary

| Metric | Value |
|---|---|
| Chains evaluated | 30 (29 completed, 1 partial) |
| Total steps | 93 |
| Step-level accuracy | 83/93 = **89.2%** |
| End-to-end accuracy | 17/29 = **58.6%** |
| Error propagation gap | **30.6 pp** |

### Per-Chain Results

| Chain | Template | Steps | Correct | E2E |
|---|---|---|---|---|
| chain01 | paper_to_experiment | 4 | 3 | No |
| chain02 | paper_to_experiment | 4 | 4 | **Yes** |
| chain03 | paper_to_experiment | 3 | 2 | No |
| chain04 | critical_appraisal | 3 | 3 | **Yes** |
| chain05 | critical_appraisal | 3 | 1 | No |
| chain06 | critical_appraisal | 3 | 3 | **Yes** |
| chain07 | genetics_to_therapy | 3 | 2 | No |
| chain08 | genetics_to_therapy | 3 | 1 | No |
| chain09 | genetics_to_therapy | 3 | 2 | No |
| chain10 | structure_to_drug | 3 | 3 | **Yes** |
| chain11 | structure_to_drug | 4 | 4 | **Yes** |
| chain12 | structure_to_drug | 4 | 4 | **Yes** |
| chain13 | stats_pipeline | 3 | 3 | **Yes** |
| chain14 | stats_pipeline | 3 | 3 | **Yes** |
| chain15 | stats_pipeline | 3 | 1 | No |
| chain16 | protocol_troubleshoot | 3 | 2 | No |
| chain17 | protocol_troubleshoot | 3 | 3 | **Yes** |
| chain18 | protocol_troubleshoot | 2 | 2 | **Yes** |
| chain19 | paradox_resolution | 3 | 3 | **Yes** |
| chain20 | paradox_resolution | 3 | 3 | **Yes** |
| chain21 | paradox_resolution | 3 | 2 | No |
| chain22 | sequence_to_function | 3 | 3 | **Yes** |
| chain23 | sequence_to_function | 2 | 2 | **Yes** |
| chain24 | sequence_to_function | 4 | 4 | **Yes** |
| chain25 | data_to_mechanism | 3 | 2 | No |
| chain26 | data_to_mechanism | 3 | 3 | **Yes** |
| chain27 | data_to_mechanism | 3 | 3 | **Yes** |
| chain28 | evidence_synthesis | 3 | 2 | No |
| chain29 | evidence_synthesis | 3 | 2 | No |
| chain30 | evidence_synthesis | 3 | 3 | **Yes** |

### By Template Type

| Template | Chains | E2E Correct | Rate |
|---|---|---|---|
| Evidence Synthesis | 3 | 1 | 33% |
| Paper to Experiment | 3 | 1 | 33% |
| Genetics to Therapy | 3 | 0 | 0% |
| Protocol Troubleshoot | 3 | 2 | 67% |
| Structure to Drug | 3 | 3 | 100% |
| Stats Pipeline | 3 | 2 | 67% |
| Critical Appraisal | 3 | 2 | 67% |
| Paradox Resolution | 3 | 2 | 67% |
| Sequence to Function | 3 | 3 | 100% |
| Data to Mechanism | 3 | 2 | 67% |

## Methodological Audit

### Judge Audit (n=20 LitQA2 responses)

| Metric | Value |
|---|---|
| Inter-judge agreement | 90.0% |
| Cohen's kappa | 0.765 (substantial) |
| Position bias rate | 10.0% (2/20 changed when order swapped) |
| Verbosity bias | +5.0% (verbose answers scored higher) |
| Judge 1 (Sonnet 4.5) positive rate | 35.0% |
| Judge 2 (Opus 4.6) positive rate | 25.0% |

### Contamination Probes

| Probe | Sample | Match Rate |
|---|---|---|
| Cloze completion | 20 | **0.0%** |
| Reverse question reconstruction | 20 | **0.0%** |
| Temporal split | — | Insufficient date metadata |

**Interpretation:** No evidence of benchmark memorization detected. The model could not complete truncated questions or reconstruct questions from answers alone.

### IRT Analysis

Requires 2+ models for meaningful item parameter estimation. Single-model run collected; awaiting additional model evaluations.

## Cost Analysis

### By Category (sorted by total cost)

| Category | Runs | Total Cost | $/Correct | Avg Latency |
|---|---|---|---|---|
| SeqQA | 585 | $38.62 | $0.37 | 16,282ms |
| DbQA | 492 | $13.20 | $0.57 | 8,585ms |
| Statistical Reasoning | 200 | $11.22 | $0.24 | 11,591ms |
| Chain Tasks | 93 | $10.30 | $0.12 | 33,476ms |
| Hypothesis Gen. | 100 | $9.95 | $0.10 | 27,966ms |
| Cloning Scenarios | 33 | $5.25 | $0.40 | 24,384ms |
| LitQA2 | 186 | $3.72 | $0.06 | 7,491ms |
| ProtocolQA | 100 | $3.16 | $0.21 | 11,443ms |
| FigQA | 181 | $2.84 | $0.15 | 6,217ms |
| Structure Analysis | 303 | $2.30 | $0.02 | 3,387ms |
| SuppQA | 81 | $1.15 | $0.13 | 5,857ms |
| Calibration | 100 | $1.11 | $0.01 | 4,826ms |

### Token Breakdown

| Metric | Value |
|---|---|
| Total input tokens | 766,478 |
| Total output tokens | 1,217,518 |
| Input cost | $11.50 |
| Output cost | $91.31 |
| **Total** | **$102.81** |

Output tokens dominate cost (88.8%), reflecting the model's verbose scientific responses at Opus pricing ($75/M output).

## Reproducing Results

```bash
# Clone
git clone https://github.com/VibeCodingScientist/LABBench2-Pro.git
cd LABBench2-Pro

# Setup
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY
pip install -e .
docker compose up -d

# Run full pipeline
./run_all.sh --model claude-opus-4.6 --concurrency 3

# Generate figures and tables from results
python results/generate_all.py
```

## Raw Data Files

| File | Description | Size |
|---|---|---|
| `results/raw/eval_runs.csv` | All 2,454 eval runs with model responses, grades, costs, latency | Full model traces |
| `results/raw/chain_runs.csv` | All 93 chain step executions with responses | Full chain traces |
| `results/raw/judge_audits.csv` | 60 judge audit records (20 items x 3 judge variants) | Judge comparison data |
| `results/raw/tasks.csv` | All 2,522 tasks with metadata | Task definitions |
| `results/raw/pipeline.log` | Complete pipeline execution log | Execution trace |
| `results/summary.json` | Machine-readable summary of all results | Structured data |

## Figures

| Figure | Description |
|---|---|
| `fig1_accuracy_overview` | Bar chart with 95% CIs across all tiers |
| `fig2_chain_error_propagation` | Per-step and cumulative accuracy; step vs end-to-end comparison |
| `fig3_cost_accuracy` | Cost-accuracy scatter plot by category |
| `fig4_judge_audit` | Agreement rates and Cohen's kappa visualization |
| `fig5_latency` | Latency distribution box plots by category |
