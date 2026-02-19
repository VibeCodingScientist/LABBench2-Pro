# LABBench2-Pro Results

## Run Details

| Parameter | Claude Opus 4.6 | Claude Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Provider | Anthropic | Anthropic | OpenAI | Google |
| Model ID | claude-opus-4-6 | claude-sonnet-4-6 | gpt-5.2 | gemini-2.5-pro |
| Date | 2026-02-18 | 2026-02-18 | 2026-02-19 | 2026-02-19 |
| Eval runs | 2,454 | 2,095 | 2,521 | 2,521 |
| Total cost | $102.81 | $14.76 | ~$8.50 | ~$6.50 |
| Concurrency | 3 | 5 | 5 | 5 |
| Infrastructure | PostgreSQL 16 + Redis 7, DigitalOcean 4vCPU/8GB | Same | Same | Same |

**Combined:** 9,591 eval runs across 4 models, 3 providers | **Total cost:** ~$132.57

## Tier 1: LABBench2 Categories

Accuracy with BCa 95% bootstrap confidence intervals (10,000 resamples):

| Category | n | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|---|
| CloningScenarios | 33 | **39.4%** | 27.3% | 15.2% | 27.3% | Opus |
| LitQA2 | 186-199 | **31.2%** | 28.3% | 30.2% | 23.1% | Opus |
| SeqQA | 298-600 | **17.8%** | 12.4% | 10.8% | 15.5% | Opus |
| ProtocolQA | 100-108 | 15.0% | 15.7% | **18.5%** | 13.9% | GPT-5.2 |
| SuppQA | 81-82 | 11.1% | 6.2% | **12.2%** | 7.3% | GPT-5.2 |
| FigQA | 181 | **10.5%** | 5.0% | 6.1% | 8.8% | Opus |
| DbQA | 492-520 | **4.7%** | 2.2% | 1.7% | 2.1% | Opus |

**Key findings:**
- Opus 4.6 leads in 5 of 7 Tier 1 categories (retrieval-heavy tasks)
- GPT-5.2 wins on ProtocolQA and SuppQA (protocol reasoning and supplementary data)
- Gemini 2.5 Pro is competitive on SeqQA (15.5%) but struggles on LitQA2 (23.1%)
- No model achieves >40% on any Tier 1 category except CloningScenarios

## Tier 2: New Task Categories

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|
| Calibration | **100.0%** | **100.0%** | **100.0%** | 94.0% | Tie (Anthropic + OpenAI) |
| Hypothesis Gen. | 97.0% | 97.0% | **100.0%** | 95.0% | GPT-5.2 |
| Structure Analysis | 45.5% | 43.9% | **52.8%** | 51.2% | GPT-5.2 |
| Statistical Reasoning | 23.0% | 35.0% | 19.5% | **67.5%** | Gemini |

**Standout results:**
- **Gemini 2.5 Pro dominates Statistical Reasoning** at 67.5% — nearly 2x the next best model (Sonnet 4.6 at 35.0%). This is the largest single-category advantage in the entire benchmark.
- **GPT-5.2 leads on Structure Analysis** (52.8%) — both non-Anthropic models outperform both Anthropic models on structural biology tasks.
- **Calibration ceiling effect** — 3 of 4 models score 100%, suggesting these tasks are too easy and need refinement. Gemini's 94% is the only discriminative result.
- **GPT-5.2 achieves perfect hypothesis scores** — 100% vs 97% for Anthropic models, though this is a near-ceiling category.

## Tier 3: Compositional Chains

### Summary

| Metric | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Chains evaluated | 30 | 30 | 30 | 30 |
| Steps completed | 93 | 89 | 96 | 96 |
| Step-level accuracy | 78/93 = **83.9%** | 79/89 = **88.8%** | 78/96 = **81.2%** | 69/96 = **71.9%** |
| End-to-end accuracy | 18/30 = **60.0%** | 22/30 = **73.3%** | 11/30 = **36.7%** | 11/30 = **36.7%** |
| Error propagation gap | 23.9 pp | **15.5 pp** | 44.5 pp | 35.2 pp |

**Anthropic models dramatically outperform on multi-step research workflows:**
- Sonnet 4.6 achieves 73.3% E2E — double the rate of GPT-5.2 and Gemini (both 36.7%)
- Even Opus 4.6 at 60.0% significantly outperforms non-Anthropic models
- GPT-5.2 has the largest error propagation gap (44.5pp) — high step accuracy does not translate to chain completion
- Gemini has the lowest step-level accuracy (71.9%) AND high error propagation

### Per-Chain Results

| Chain | Template | Opus | Sonnet | GPT-5.2 | Gemini |
|---|---|---|---|---|---|
| chain01 | paper_to_experiment | 3/4 No | 1/2 No | 1/4 No | 0/4 No |
| chain02 | paper_to_experiment | 4/4 **Yes** | 2/2 **Yes** | 3/4 No | 4/4 **Yes** |
| chain03 | paper_to_experiment | 2/3 No | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** |
| chain04 | critical_appraisal | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 3/3 **Yes** |
| chain05 | critical_appraisal | 1/3 No | 1/3 No | 1/3 No | 3/3 **Yes** |
| chain06 | critical_appraisal | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 1/3 No |
| chain07 | genetics_to_therapy | 2/3 No | 3/3 **Yes** | 2/3 No | 3/3 **Yes** |
| chain08 | genetics_to_therapy | 1/3 No | 2/3 No | 0/3 No | 1/2 No |
| chain09 | genetics_to_therapy | 2/3 No | 3/3 **Yes** | 2/3 No | 3/3 **Yes** |
| chain10 | structure_to_drug | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** |
| chain11 | structure_to_drug | 4/4 **Yes** | 4/4 **Yes** | 4/4 **Yes** | 3/4 No |
| chain12 | structure_to_drug | 4/4 **Yes** | 4/4 **Yes** | 4/4 **Yes** | 4/4 **Yes** |
| chain13 | stats_pipeline | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 3/3 **Yes** |
| chain14 | stats_pipeline | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 2/3 No |
| chain15 | stats_pipeline | 1/3 No | 3/3 **Yes** | 3/3 **Yes** | 2/3 No |
| chain16 | protocol_troubleshoot | 2/3 No | 2/3 No | 1/3 No | 1/3 No |
| chain17 | protocol_troubleshoot | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** |
| chain18 | protocol_troubleshoot | 2/2 **Yes** | 2/2 **Yes** | 2/4 No | 3/4 No |
| chain19 | paradox_resolution | 3/3 **Yes** | 3/3 **Yes** | 1/3 No | 2/3 No |
| chain20 | paradox_resolution | 3/3 **Yes** | 3/3 **Yes** | 3/3 **Yes** | 2/3 No |
| chain21 | paradox_resolution | 2/3 No | 1/3 No | 2/3 No | 1/3 No |
| chain22 | sequence_to_function | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 2/3 No |
| chain23 | sequence_to_function | 2/2 **Yes** | 1/2 No | 3/3 **Yes** | 2/3 No |
| chain24 | sequence_to_function | 4/4 **Yes** | 4/4 **Yes** | 3/4 No | 3/4 No |
| chain25 | data_to_mechanism | 2/3 No | 3/3 **Yes** | 3/3 **Yes** | 1/3 No |
| chain26 | data_to_mechanism | 3/3 **Yes** | 3/3 **Yes** | 2/3 No | 3/3 **Yes** |
| chain27 | data_to_mechanism | 3/3 **Yes** | 2/3 No | 3/3 **Yes** | 2/3 No |
| chain28 | evidence_synthesis | 2/3 No | 3/3 **Yes** | 2/3 No | 2/3 No |
| chain29 | evidence_synthesis | 2/3 No | 3/3 **Yes** | 3/3 **Yes** | 2/3 No |
| chain30 | evidence_synthesis | 3/3 **Yes** | 2/3 No | 2/3 No | 2/3 No |

### By Template Type (E2E counts)

| Template | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Paper to Experiment | 1/3 | 1/3 | 1/3 | 1/3 |
| Critical Appraisal | 2/3 | 2/3 | 0/3 | 2/3 |
| Genetics to Therapy | 0/3 | **2/3** | 0/3 | **2/3** |
| Structure to Drug | **3/3** | **3/3** | **3/3** | 2/3 |
| Stats Pipeline | 2/3 | **3/3** | 1/3 | 1/3 |
| Protocol Troubleshoot | 2/3 | 2/3 | 1/3 | 1/3 |
| Paradox Resolution | 2/3 | 2/3 | 1/3 | 0/3 |
| Sequence to Function | **3/3** | 2/3 | 1/3 | 0/3 |
| Data to Mechanism | 2/3 | 2/3 | 2/3 | 1/3 |
| Evidence Synthesis | 1/3 | **2/3** | 1/3 | 0/3 |

**Template-level findings:**
- Structure-to-Drug is the "easy" chain template: 3 models achieve 3/3
- Evidence Synthesis and Paradox Resolution are hard across the board
- Sonnet 4.6 is the only model to achieve 3/3 on Stats Pipeline and 2/3 on Evidence Synthesis
- GPT-5.2 and Gemini struggle on sequence analysis and paradox resolution chains

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

**Interpretation:** No evidence of benchmark memorization detected. GPT-5.2 and Gemini contamination probes not run (requires Anthropic-specific judge setup).

### IRT Analysis (4 models, 2,522 items)

| Metric | Value |
|---|---|
| Total items | 2,522 |
| Low discrimination (< 0.3) | 2,219 (88.0%) |
| Too easy (p > 0.95) | 402 (15.9%) |
| Too hard (p < 0.05) | 1,481 (58.7%) |
| Recommended pruned set | 303 items |
| Peak discrimination | theta = 0.50 |
| Effective range | theta = [-1.90, 2.90] |

With 4 models, 88.0% of items still do not discriminate between models. The recommended pruned set of 303 items would enable 8x more efficient evaluations. The effective range expanded from [-1.20, 2.20] (2-model) to [-1.90, 2.90] (4-model), indicating that the additional models stretch the ability spectrum.

## Cost Analysis

### By Category (all models)

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| SeqQA | $38.62 | $3.37 | ~$2.50 | ~$2.00 |
| DbQA | $13.20 | $2.64 | ~$1.50 | ~$1.20 |
| Statistical Reasoning | $11.22 | $2.92 | ~$1.00 | ~$0.80 |
| Hypothesis Gen. | $9.95 | $2.00 | ~$1.50 | ~$1.20 |
| Cloning Scenarios | $5.25 | $0.97 | ~$0.30 | ~$0.25 |
| LitQA2 | $3.72 | $0.70 | ~$0.50 | ~$0.40 |
| ProtocolQA | $3.16 | $0.54 | ~$0.30 | ~$0.25 |
| FigQA | $2.84 | $0.47 | ~$0.25 | ~$0.20 |
| Structure Analysis | $2.30 | $0.44 | ~$0.30 | ~$0.25 |
| SuppQA | $1.15 | $0.17 | ~$0.12 | ~$0.10 |
| Calibration | $1.11 | $0.19 | ~$0.10 | ~$0.08 |
| Chain Tasks | $10.30 | $0.35 | ~$0.30 | ~$0.25 |
| **Total** | **$102.81** | **$14.76** | **~$8.50** | **~$6.50** |

### Cost Efficiency

| Model | Total Cost | Avg Tier 1 Accuracy | Cost per Tier 1 % Point |
|---|---|---|---|
| Gemini 2.5 Pro | ~$6.50 | 14.0% | $0.46/pp |
| GPT-5.2 | ~$8.50 | 13.6% | $0.63/pp |
| Claude Sonnet 4.6 | $14.76 | 13.9% | $1.06/pp |
| Claude Opus 4.6 | $102.81 | 18.5% | $5.56/pp |

Opus 4.6 achieves the highest Tier 1 accuracy but at 12x the cost of Gemini. For Tier 2 tasks, Gemini's Statistical Reasoning dominance makes it the most cost-effective model overall.

## Reproducing Results

```bash
git clone https://github.com/VibeCodingScientist/LABBench2-Pro.git
cd LABBench2-Pro

cp .env.example .env
# Edit .env with your API keys:
#   ANTHROPIC_API_KEY (required for judge + Anthropic models)
#   OPENAI_API_KEY (for GPT-5.2)
#   GOOGLE_API_KEY (for Gemini 2.5 Pro)
pip install -e .
docker compose up -d

# Run full pipeline for each model
./run_all.sh --model claude-opus-4.6 --concurrency 3
./run_all.sh --model claude-sonnet-4.6 --concurrency 5
./run_all.sh --model gpt-5.2 --concurrency 5
./run_all.sh --model gemini-2.5-pro --concurrency 5

# Generate figures and tables
python results/generate_all.py
```

## Raw Data Files

| File | Description |
|---|---|
| `results/raw/eval_runs.csv` | All 9,591 eval runs with model responses, grades, costs, latency |
| `results/raw/chain_runs.csv` | All chain step executions with responses (4 models x 30 chains) |
| `results/raw/judge_audits.csv` | 60 judge audit records |
| `results/raw/tasks.csv` | All tasks with metadata |
| `results/raw/pipeline.log` | Execution logs |
| `results/summary.json` | Machine-readable summary of all results |
