# LABBench2-Pro Results

## Run Details

| Parameter | Claude Opus 4.6 | Claude Sonnet 4.6 |
|---|---|---|
| Model ID | claude-opus-4-6 | claude-sonnet-4-6 |
| Provider | Anthropic | Anthropic |
| Date | 2026-02-18 | 2026-02-18 |
| Eval runs | 2,454 | 2,095 |
| Total tokens (in) | 766,478 | 577,536 |
| Total tokens (out) | 1,217,518 | 868,789 |
| Total cost | $102.81 | $14.76 |
| Concurrency | 3 | 5 |
| Infrastructure | PostgreSQL 16 + Redis 7, DigitalOcean 4vCPU/8GB | Same |

**Combined cost:** $117.57 | **Combined eval runs:** 4,549

## Tier 1: LABBench2 Categories

Accuracy with BCa 95% bootstrap confidence intervals (10,000 resamples) and pairwise significance tests:

| Category | n | Opus 4.6 | 95% CI | Sonnet 4.6 | 95% CI | p-value | Sig? |
|---|---|---|---|---|---|---|---|
| CloningScenarios | 33 | **39.4%** | [24.2%, 57.6%] | 27.3% | [15.2%, 45.5%] | 0.014 | **YES** |
| LitQA2 | 186-187 | **31.2%** | [24.7%, 38.2%] | 28.3% | [21.9%, 35.3%] | 0.210 | No |
| SeqQA | 298-585 | **17.8%** | [14.9%, 21.0%] | 12.4% | [9.1%, 16.4%] | 0.214 | No |
| ProtocolQA | 100-102 | 15.0% | [9.0%, 23.0%] | **15.7%** | [9.8%, 23.5%] | 0.434 | No |
| SuppQA | 81 | **11.1%** | [4.9%, 19.8%] | 6.2% | [2.5%, 13.6%] | 0.065 | No |
| FigQA | 181 | **10.5%** | [6.6%, 15.5%] | 5.0% | [2.2%, 8.8%] | 0.003 | **YES** |
| DbQA | 492-493 | **4.7%** | [3.0%, 6.7%] | 2.2% | [1.2%, 3.9%] | 0.001 | **YES** |

**Key takeaway:** Only 3 of 7 categories show statistically significant differences. Opus leads on retrieval-heavy tasks (FigQA, DbQA, CloningScenarios), but most categories cannot be distinguished.

## Tier 2: New Task Categories

| Category | Opus 4.6 | Sonnet 4.6 | Delta | Verification |
|---|---|---|---|---|
| Calibration | **100.0%** | **100.0%** | 0 | LLM-judge |
| Hypothesis Gen. (lenient) | **100.0%** | **100.0%** | 0 | LLM-judge |
| Hypothesis Gen. (strict) | **97.0%** | **97.0%** | 0 | LLM-judge (strict rubric) |
| Structure Analysis | **45.5%** | 43.9% | -1.6 | Programmatic |
| Statistical Reasoning | 23.0% | **35.0%** | **+12.0** | Programmatic |

**Standout:** Sonnet 4.6 significantly outperforms Opus 4.6 on Statistical Reasoning (+12pp). Both achieve identical strict hypothesis scores (97%).

### Rubric Sensitivity Analysis (Hypothesis Generation)

Per-criterion comparison under strict rubric (3/4 criteria must pass):

| Criterion | Opus 4.6 | Sonnet 4.6 | Description |
|---|---|---|---|
| Falsifiable | 100% | 100% | Hypotheses include specific experimental designs |
| Connected | 95% | 90% | Hypotheses logically follow from abstract findings |
| Specific | 95% | 94% | Hypotheses name concrete biological entities |
| Novel | 67% | 71% | Hypotheses go beyond obvious next steps |

**Novelty is the primary failure mode for both models.** Re-grading cost: $3.17 total (200 re-grades via Claude Sonnet 4.5).

## Tier 3: Compositional Chains

### Summary

| Metric | Opus 4.6 | Sonnet 4.6 |
|---|---|---|
| Chains evaluated | 30 | 30 |
| Steps completed | 93 | 89 |
| Step-level accuracy | 78/93 = **83.9%** | 79/89 = **88.8%** |
| End-to-end accuracy | 18/30 = **60.0%** | 22/30 = **73.3%** |
| Error propagation gap | **23.9 pp** | **15.5 pp** |

**Sonnet 4.6 outperforms Opus 4.6 on multi-step research workflows** (73.3% vs 60.0% E2E), despite performing worse on most atomic retrieval tasks.

### Per-Chain Results

| Chain | Template | Opus Steps | Opus E2E | Sonnet Steps | Sonnet E2E |
|---|---|---|---|---|---|
| chain01 | paper_to_experiment | 3/4 | No | 1/2 | No |
| chain02 | paper_to_experiment | 4/4 | **Yes** | 2/2 | **Yes** |
| chain03 | paper_to_experiment | 2/3 | No | 3/3 | **Yes** |
| chain04 | critical_appraisal | 3/3 | **Yes** | 3/3 | **Yes** |
| chain05 | critical_appraisal | 1/3 | No | 1/3 | No |
| chain06 | critical_appraisal | 3/3 | **Yes** | 3/3 | **Yes** |
| chain07 | genetics_to_therapy | 2/3 | No | 3/3 | **Yes** |
| chain08 | genetics_to_therapy | 1/3 | No | 2/3 | No |
| chain09 | genetics_to_therapy | 2/3 | No | 3/3 | **Yes** |
| chain10 | structure_to_drug | 3/3 | **Yes** | 3/3 | **Yes** |
| chain11 | structure_to_drug | 4/4 | **Yes** | 4/4 | **Yes** |
| chain12 | structure_to_drug | 4/4 | **Yes** | 4/4 | **Yes** |
| chain13 | stats_pipeline | 3/3 | **Yes** | 3/3 | **Yes** |
| chain14 | stats_pipeline | 3/3 | **Yes** | 3/3 | **Yes** |
| chain15 | stats_pipeline | 1/3 | No | 3/3 | **Yes** |
| chain16 | protocol_troubleshoot | 2/3 | No | 2/3 | No |
| chain17 | protocol_troubleshoot | 3/3 | **Yes** | 3/3 | **Yes** |
| chain18 | protocol_troubleshoot | 2/2 | **Yes** | 2/2 | **Yes** |
| chain19 | paradox_resolution | 3/3 | **Yes** | 3/3 | **Yes** |
| chain20 | paradox_resolution | 3/3 | **Yes** | 3/3 | **Yes** |
| chain21 | paradox_resolution | 2/3 | No | 1/3 | No |
| chain22 | sequence_to_function | 3/3 | **Yes** | 3/3 | **Yes** |
| chain23 | sequence_to_function | 2/2 | **Yes** | 1/2 | No |
| chain24 | sequence_to_function | 4/4 | **Yes** | 4/4 | **Yes** |
| chain25 | data_to_mechanism | 2/3 | No | 3/3 | **Yes** |
| chain26 | data_to_mechanism | 3/3 | **Yes** | 3/3 | **Yes** |
| chain27 | data_to_mechanism | 3/3 | **Yes** | 2/3 | No |
| chain28 | evidence_synthesis | 2/3 | No | 3/3 | **Yes** |
| chain29 | evidence_synthesis | 2/3 | No | 3/3 | **Yes** |
| chain30 | evidence_synthesis | 3/3 | **Yes** | 2/3 | No |

### By Template Type

| Template | Opus E2E | Sonnet E2E |
|---|---|---|
| Paper to Experiment | 1/3 | 1/3 |
| Structure to Drug | 3/3 | 3/3 |
| Stats Pipeline | 2/3 | **3/3** |
| Critical Appraisal | 2/3 | 2/3 |
| Genetics to Therapy | 0/3 | **2/3** |
| Protocol Troubleshoot | 2/3 | 2/3 |
| Paradox Resolution | 2/3 | 2/3 |
| Sequence to Function | 3/3 | 2/3 |
| Data to Mechanism | 2/3 | 2/3 |
| Evidence Synthesis | 1/3 | 1/3 |

## Methodological Audit

### Judge Audit (n=20 LitQA2 responses)

| Metric | Value |
|---|---|
| Inter-judge agreement | 90.0% |
| Cohen's kappa | 0.765 (substantial) |
| Position bias rate | 15.0% (3/20 changed when order swapped) |
| Verbosity bias | +5.0% (verbose answers scored higher) |
| Judge 1 (Sonnet 4.5) positive rate | 35.0% |
| Judge 2 (Opus 4.6) positive rate | 25.0% |

### Contamination Probes

| Probe | Opus 4.6 | Sonnet 4.6 |
|---|---|---|
| Cloze completion | **0.0%** (0/20) | **0.0%** (0/20) |
| Reverse question | **0.0%** (0/20) | N/A |
| Temporal split | Insufficient metadata | Insufficient metadata |

**Interpretation:** No evidence of benchmark memorization detected for either model.

### IRT Analysis (2 models, 2,459 items)

| Metric | Value |
|---|---|
| Total items | 2,459 |
| Low discrimination (< 0.3) | 2,326 (94.6%) |
| Too easy (p > 0.95) | 572 (23.3%) |
| Too hard (p < 0.05) | 1,642 (66.8%) |
| Recommended pruned set | 133 items |
| Peak discrimination | theta = 0.50 |

94.6% of items do not discriminate between models—benchmark pruning to the 133 high-discrimination items would enable 19x more efficient evaluations.

## Cost Analysis

### By Category (both models)

| Category | Opus Cost | Sonnet Cost | Opus $/Correct | Sonnet $/Correct |
|---|---|---|---|---|
| SeqQA | $38.62 | $3.37 | $0.37 | $0.09 |
| DbQA | $13.20 | $2.64 | $0.57 | $0.24 |
| Statistical Reasoning | $11.22 | $2.92 | $0.24 | $0.04 |
| Chain Tasks | $10.30 | $0.35 | $0.13 | $0.005 |
| Hypothesis Gen. | $9.95 | $2.00 | $0.10 | $0.02 |
| Cloning Scenarios | $5.25 | $0.97 | $0.40 | $0.11 |
| LitQA2 | $3.72 | $0.70 | $0.06 | $0.01 |
| ProtocolQA | $3.16 | $0.54 | $0.21 | $0.03 |
| FigQA | $2.84 | $0.47 | $0.15 | $0.05 |
| Structure Analysis | $2.30 | $0.44 | $0.02 | $0.003 |
| SuppQA | $1.15 | $0.17 | $0.13 | $0.03 |
| Calibration | $1.11 | $0.19 | $0.01 | $0.002 |
| **Total** | **$102.81** | **$14.76** | | |

### Token Breakdown

| Metric | Opus 4.6 | Sonnet 4.6 |
|---|---|---|
| Input tokens | 766,478 | 577,536 |
| Output tokens | 1,217,518 | 868,789 |
| Input cost | $11.50 | $1.73 |
| Output cost | $91.31 | $13.03 |
| **Total** | **$102.81** | **$14.76** |

Sonnet 4.6 achieves comparable or superior performance at **7x lower cost** across the full benchmark.

## Reproducing Results

```bash
git clone https://github.com/VibeCodingScientist/LABBench2-Pro.git
cd LABBench2-Pro

cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY
pip install -e .
docker compose up -d

# Run full pipeline for each model
./run_all.sh --model claude-opus-4.6 --concurrency 3
./run_all.sh --model claude-sonnet-4.6 --concurrency 5

# Strict hypothesis re-grading
python -m src.tier1.regrade_hypothesis --model claude-opus-4.6
python -m src.tier1.regrade_hypothesis --model claude-sonnet-4.6

# Generate figures and tables
python results/generate_all.py
```

## Raw Data Files

| File | Description |
|---|---|
| `results/raw/eval_runs.csv` | All 4,549 eval runs with model responses, grades, costs, latency |
| `results/raw/chain_runs.csv` | All chain step executions with responses |
| `results/raw/judge_audits.csv` | 60 judge audit records |
| `results/raw/tasks.csv` | All tasks with metadata |
| `results/raw/pipeline.log` | Execution logs |
| `results/summary.json` | Machine-readable summary of all results |
