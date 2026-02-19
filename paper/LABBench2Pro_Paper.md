# LABBench2-Pro: A Methodological Audit and Extension of Scientific LLM Evaluation

## Abstract

We present LABBench2-Pro, a methodological audit and extension of FutureHouse's LABBench2 benchmark for evaluating large language models on biology research tasks. Our analysis identifies five categories of gaps in current scientific LLM evaluation: scoring reliability, statistical reporting, contamination risk, coverage gaps, and atomic-only evaluation. We address these through a three-tier framework: (1) a methodological audit of existing LABBench2 categories with bootstrap confidence intervals, IRT analysis, judge reliability testing, and contamination probes; (2) 703 new programmatically-generated tasks spanning statistical reasoning, structural biology, uncertainty calibration, and hypothesis generation; and (3) 30 compositional chains testing multi-step scientific reasoning where errors compound. Evaluating four frontier models from three providers---Claude Opus 4.6, Claude Sonnet 4.6 (Anthropic), GPT-5.2 (OpenAI), and Gemini 2.5 Pro (Google)---across 9,591 total evaluation runs ($132.57 total cost), we find: (a) no single model dominates across all tiers---Opus leads 5 of 7 LABBench categories, GPT-5.2 wins protocol and supplementary reasoning, Gemini dominates statistical reasoning at 67.5% (nearly 2x the next best), and GPT-5.2 leads structural biology; (b) Anthropic models dramatically outperform on compositional chains (Sonnet 73.3%, Opus 60.0% end-to-end accuracy vs. GPT-5.2 and Gemini both at 36.7%), demonstrating that atomic task performance does not predict multi-step reasoning ability; (c) cost varies 16x across providers ($102.81 for Opus vs. ~$6.50 for Gemini) with the most expensive model not uniformly the best; (d) LLM-as-judge grading shows substantial inter-judge agreement (Cohen's kappa = 0.765) with measurable position bias (15%) and verbosity bias (+5%); (e) rubric sensitivity analysis reveals that both Anthropic models achieve 97% on hypothesis generation under strict grading, with novelty as the primary failure mode (67--71% pass rate); and (f) zero contamination signal across cloze probes for Anthropic models. IRT analysis across 2,522 items with 4 models identifies 303 high-discrimination items from a pool of 2,219 low-discrimination items (88%), enabling 8x more efficient future evaluations. All code, data, model traces, and analysis are publicly available.

## 1. Introduction

Benchmarks for evaluating LLMs on scientific tasks are proliferating rapidly. FutureHouse's LABBench2 established a comprehensive evaluation framework covering literature comprehension (LitQA2), figure interpretation (FigQA), sequence analysis (SeqQA), supplementary material reasoning (SuppQA), protocol understanding (ProtocolQA), database queries (DbQA), table interpretation (TableQA), and molecular cloning scenarios (CloningScenarios). However, benchmarks themselves require benchmarking.

Our gap analysis identified five categories of methodological issues in current scientific LLM evaluation:

1. **Scoring reliability.** LLM-as-judge grading introduces unquantified noise. Position bias (does answer order affect grading?), verbosity bias (are longer answers scored more favorably?), and inter-judge disagreement are never measured.

2. **Statistical reporting.** Results are reported as point estimates without confidence intervals. With small per-category sample sizes, two models can appear different when their confidence intervals overlap entirely.

3. **Contamination risk.** No probes test whether models have memorized benchmark questions from training data.

4. **Coverage gaps.** LABBench2 tests retrieval and comprehension but not statistical reasoning, uncertainty calibration, hypothesis generation, or structural biology interpretation.

5. **Atomic-only evaluation.** Every task is independent. Real research requires multi-step reasoning where errors compound---reading a paper, interpreting its figures, choosing the right statistical test, designing a follow-up experiment.

LABBench2-Pro addresses each gap systematically. Critically, by evaluating four models from three independent providers, we move beyond single-provider comparisons to reveal cross-provider capability patterns that no two-model study can capture.

## 2. Methods

### 2.1 Tier 1: Methodological Audit

We ran four models against all available LABBench2 categories from HuggingFace, then applied four analyses:

**Bootstrap confidence intervals.** BCa (bias-corrected and accelerated) 95% confidence intervals computed via 10,000 bootstrap resamples on accuracy within each category, with pairwise significance tests (Bonferroni-corrected).

**Item Response Theory.** Classical item analysis with discrimination indices and difficulty parameters across 2,522 items and 4 models, used to identify low-discrimination items and compute a recommended pruned item set.

**Judge audit.** We sampled 20 LitQA2 responses and graded each with two LLM judges (Claude Sonnet 4.5 and Claude Opus 4.6). We measured inter-judge agreement (Cohen's kappa), position bias (swapping reference/response order), and verbosity bias (comparing grades on original vs. shortened responses).

**Contamination probes.** Two probes per Anthropic model: (1) cloze completion---can the model complete a truncated benchmark question?; (2) reverse reconstruction---can it guess the question from only the answer? Contamination probes were run for Anthropic models only; cross-provider contamination testing was not performed.

### 2.2 Tier 2: New Task Categories (703 tasks)

All tasks are programmatically generated with deterministic or rubric-graded ground truth:

| Category | Count | Verification | Source |
|---|---|---|---|
| Statistical Reasoning | 200 | Programmatic (scipy) | Synthetic gene expression data |
| Structure Analysis | 303 | Programmatic + LLM-judge | Real PDB structures (BioPython) + synthetic gel images |
| Uncertainty Calibration | 100 | LLM-judge | Placeholder questions with deliberately insufficient information |
| Hypothesis Generation | 100 | LLM-judge (rubric) | Real PubMed abstracts (NCBI Entrez) |

**Statistical Reasoning** tasks present gene expression data and ask the model to select appropriate statistical tests, compute p-values, or count significant genes. Ground truth is computed via scipy.stats.

**Structure Analysis** tasks query real PDB protein structures (downloaded via BioPython) for properties such as chain count, residue count, secondary structure content, and resolution. Gel image tasks present synthetic gel electrophoresis images with known band patterns.

**Uncertainty Calibration** tasks present questions with critical information deliberately removed. The correct response is to identify the information gap rather than guess.

**Hypothesis Generation** tasks provide real PubMed abstracts and ask the model to generate testable hypotheses. Graded by LLM-judge with structured rubric; re-graded with a strict per-criterion rubric (Section 3.6).

### 2.3 Tier 3: Compositional Chains (30 chains)

Multi-step research workflows where each step's output is prepended to the next step's prompt. If the model answers step 1 incorrectly, the error propagates. We designed 30 chains across 10 template types:

| Template | Steps | Chains | Example Topics |
|---|---|---|---|
| Paper to Experiment | 4 | 3 | SHP2 allosteric inhibition, JAK2 V617F, PCSK9 |
| Structure to Drug | 4 | 3 | EGFR T790M, KRAS G12C, SARS-CoV-2 Mpro |
| Stats Pipeline | 3 | 3 | TNBC RNA-seq, T2D GWAS, scRNA-seq TILs |
| Critical Appraisal | 3 | 3 | IDH1 glioma, Lecanemab, Microbiome-melanoma |
| Genetics to Therapy | 3 | 3 | PINK1 Parkinson's, CFTR DF508, SCN1A Dravet |
| Protocol Troubleshoot | 3 | 3 | KRAS-BRAF co-IP, ChIP-seq, CRISPR base editing |
| Paradox Resolution | 3 | 3 | ZEB1 EMT, PD-1 hyperprogression, Exercise immunosuppression |
| Sequence to Function | 3 | 3 | Psychrophilic LDH, AmpC beta-lactamase, Cas effector |
| Data to Mechanism | 3 | 3 | GAPDH confound, Imatinib resistance, Venetoclax synergy |
| Evidence Synthesis | 3 | 3 | ctDNA lung cancer, FLT3 AML, BRAF melanoma |

All chains use verified data: 95/95 data points checked against PDB, UniProt, ChEMBL, ClinVar, ClinicalTrials.gov, and Open Targets.

### 2.4 Experimental Setup

| Parameter | Claude Opus 4.6 | Claude Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Model ID | claude-opus-4-6 | claude-sonnet-4-6 | gpt-5.2 | gemini-2.5-pro |
| Provider | Anthropic | Anthropic | OpenAI | Google |
| Input cost | $15.00/M tokens | $3.00/M tokens | $1.75/M tokens | $1.25/M tokens |
| Output cost | $75.00/M tokens | $15.00/M tokens | $14.00/M tokens | $10.00/M tokens |
| Eval runs | 2,454 | 2,095 | 2,521 | 2,521 |
| Total cost | $102.81 | $14.76 | ~$8.50 | ~$6.50 |
| Concurrency | 3 | 5 | 5 | 5 |

- **Infrastructure:** PostgreSQL 16 + Redis 7 on DigitalOcean (4 vCPU, 8 GB RAM)
- **Combined cost:** ~$132.57 (4 models, 3 providers)
- **Total eval runs:** 9,591
- **Resume safety:** All evaluations check the database before calling the API; interrupted runs resume without re-evaluating completed tasks

## 3. Results

### 3.1 Tier 1: LABBench Category Performance

| Category | n | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|---|
| CloningScenarios | 33 | **39.4%** | 27.3% | 15.2% | 27.3% | Opus |
| LitQA2 | 186--199 | **31.2%** | 28.3% | 30.2% | 23.1% | Opus |
| SeqQA | 298--600 | **17.8%** | 12.4% | 10.8% | 15.5% | Opus |
| ProtocolQA | 100--108 | 15.0% | 15.7% | **18.5%** | 13.9% | GPT-5.2 |
| SuppQA | 81--82 | 11.1% | 6.2% | **12.2%** | 7.3% | GPT-5.2 |
| FigQA | 181 | **10.5%** | 5.0% | 6.1% | 8.8% | Opus |
| DbQA | 492--520 | **4.7%** | 2.2% | 1.7% | 2.1% | Opus |

Key observations:

- **Opus 4.6 leads 5 of 7 categories** but does not dominate universally. GPT-5.2 wins ProtocolQA and SuppQA---tasks requiring protocol reasoning and supplementary data extraction.

- **No model achieves >40%** on any Tier 1 category except CloningScenarios, underscoring the difficulty of these tasks for all frontier models.

- **Gemini 2.5 Pro is competitive on SeqQA** (15.5%, second to Opus) but struggles on LitQA2 (23.1%, the lowest score).

- **Opus advantages cluster in retrieval-heavy tasks.** The categories where Opus leads most clearly (FigQA, DbQA, CloningScenarios) all require extracting specific information from structured inputs.

### 3.2 Tier 2: New Task Categories

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|
| Calibration | **100.0%** | **100.0%** | **100.0%** | 94.0% | Tie (3 models) |
| Hypothesis Generation (strict) | 97.0% | 97.0% | **100.0%** | 95.0% | GPT-5.2 |
| Structure Analysis | 45.5% | 43.9% | **52.8%** | 51.2% | GPT-5.2 |
| **Statistical Reasoning** | 23.0% | 35.0% | 19.5% | **67.5%** | **Gemini** |

**Statistical Reasoning is the standout result.** Gemini 2.5 Pro achieves 67.5%---nearly 2x the next best model (Sonnet 4.6 at 35.0%). This is the largest single-category advantage in the entire benchmark and reveals a clear domain specialization: Gemini's training or architecture confers a substantial edge on tasks requiring test selection, p-value computation, and multiple testing correction.

**Structure Analysis favors non-Anthropic models.** Both GPT-5.2 (52.8%) and Gemini (51.2%) outperform both Anthropic models (Opus 45.5%, Sonnet 43.9%) on structural biology tasks, suggesting that protein structure reasoning is differentially developed across providers.

**Calibration** achieved 100% on three of four models, suggesting frontier models reliably identify when insufficient information is available. Gemini's 94% is the only discriminative result.

**Hypothesis Generation** reached a near-ceiling: GPT-5.2 achieves perfect strict scores (100%), while Anthropic models match at 97% and Gemini at 95%. See Section 3.6 for per-criterion breakdown.

### 3.3 Tier 3: Compositional Chains

| Metric | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Chains evaluated | 30 | 30 | 30 | 30 |
| Steps completed | 93 | 89 | 96 | 96 |
| Step-level accuracy | 78/93 = 83.9% | 79/89 = **88.8%** | 78/96 = 81.2% | 69/96 = 71.9% |
| End-to-end accuracy | 18/30 = 60.0% | 22/30 = **73.3%** | 11/30 = 36.7% | 11/30 = 36.7% |
| Error propagation gap | 23.9 pp | **15.5 pp** | 44.5 pp | 35.2 pp |

**Anthropic models dramatically outperform on compositional chains.** Sonnet 4.6 achieves 73.3% end-to-end accuracy---double the rate of GPT-5.2 and Gemini (both 36.7%). Even Opus 4.6 at 60.0% significantly exceeds non-Anthropic models. This is the most striking cross-provider finding: a clear Anthropic advantage on multi-step research workflows that does not track with Tier 1 or Tier 2 performance.

**GPT-5.2 has competitive step-level accuracy (81.2%) but the worst error propagation gap (44.5 pp).** This means individual steps are often correct, but errors at critical junctures cascade more severely. Atomic task performance does not predict compositional ability.

Chain-level results by template (E2E correct out of 3):

| Template | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Paper to Experiment (4 steps) | 1/3 | 1/3 | 1/3 | 1/3 |
| Structure to Drug (4 steps) | 3/3 | 3/3 | 3/3 | 2/3 |
| Stats Pipeline | 2/3 | **3/3** | 1/3 | 1/3 |
| Critical Appraisal | 2/3 | 2/3 | 0/3 | 2/3 |
| Genetics to Therapy | 0/3 | **2/3** | 0/3 | **2/3** |
| Protocol Troubleshoot | 2/3 | 2/3 | 1/3 | 1/3 |
| Paradox Resolution | 2/3 | 2/3 | 1/3 | 0/3 |
| Sequence to Function | **3/3** | 2/3 | 1/3 | 0/3 |
| Data to Mechanism | 2/3 | 2/3 | 2/3 | 1/3 |
| Evidence Synthesis | 1/3 | **2/3** | 1/3 | 0/3 |

**Template-level observations:**

- **Structure-to-Drug** is the "easy" template: three models achieve 3/3, with Gemini at 2/3.
- **Evidence Synthesis** and **Paradox Resolution** are the hardest templates across the board.
- **Sonnet 4.6** is the only model to achieve 3/3 on Stats Pipeline and 2/3 on Evidence Synthesis.
- **Gemini** and **GPT-5.2** struggle particularly on Sequence to Function and Paradox Resolution chains, achieving 0/3 and 1/3 respectively.

### 3.4 Judge Audit

From 20 sampled LitQA2 responses graded by two LLM judges:

| Metric | Value |
|---|---|
| Inter-judge agreement | 90.0% |
| Cohen's kappa | 0.765 (substantial) |
| Position bias rate | 15.0% (3/20) |
| Verbosity bias | +5.0% (verbose answers favored) |
| Sonnet 4.5 positive rate | 35.0% |
| Opus 4.6 positive rate | 25.0% |

Inter-judge agreement is substantial (kappa = 0.765) but not perfect. Position bias at 15% and verbosity bias at +5% are measurable and should be reported alongside any LLM-judged results.

### 3.5 Contamination Probes

| Probe | Opus 4.6 | Sonnet 4.6 |
|---|---|---|
| Cloze completion | 0/20 = 0.0% | 0/20 = 0.0% |
| Reverse question | 0/20 = 0.0% | N/A |
| Temporal split | Insufficient metadata | Insufficient metadata |

No contamination signal was detected for Anthropic models. Neither model could complete truncated questions or reconstruct questions from answers, suggesting the benchmark content has not been memorized from training data. Contamination probes were not run for GPT-5.2 or Gemini 2.5 Pro (requires Anthropic-specific judge setup).

### 3.6 Rubric Sensitivity Analysis

Anthropic models achieved 100% under the lenient rubric and 97% under the strict per-criterion rubric. GPT-5.2 achieved 100% under strict grading; Gemini 95%. Per-criterion comparison (Anthropic models):

| Criterion | Opus 4.6 | Sonnet 4.6 | Description |
|---|---|---|---|
| Falsifiable | 100% | 100% | Hypotheses include specific experimental designs |
| Connected | 95% | 90% | Hypotheses logically follow from abstract findings |
| Specific | 95% | 94% | Hypotheses name concrete biological entities |
| Novel | 67% | 71% | Hypotheses go beyond obvious next steps |

Both Anthropic models show the same pattern: **novelty is the primary failure mode.** Falsifiability is universally passed, while novelty hovers around 67--71%. The models reliably generate well-structured, specific, connected hypotheses but frequently propose obvious extensions rather than genuinely creative follow-ups.

### 3.7 IRT Analysis

Classical item analysis across 2,522 items and 4 models:

| Metric | Value |
|---|---|
| Total items | 2,522 |
| Low discrimination (< 0.3) | 2,219 (88.0%) |
| Too easy (p > 0.95) | 402 (15.9%) |
| Too hard (p < 0.05) | 1,481 (58.7%) |
| Recommended pruned set | 303 items |
| Peak discrimination | theta = 0.50 |
| Effective range | theta = [-1.90, 2.90] |

With four models (up from two), 88.0% of items still do not discriminate between models. However, the effective range expanded from [-1.20, 2.20] (2-model analysis) to [-1.90, 2.90], indicating that the additional models stretch the ability spectrum. The recommended 303-item pruned set retains only items where model responses diverge, enabling 8x more efficient future evaluations.

### 3.8 Cost Analysis

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| SeqQA | $38.62 | $3.37 | ~$2.50 | ~$2.00 |
| DbQA | $13.20 | $2.64 | ~$1.50 | ~$1.20 |
| Statistical Reasoning | $11.22 | $2.92 | ~$1.00 | ~$0.80 |
| Chain Tasks | $10.30 | $0.35 | ~$0.30 | ~$0.25 |
| Hypothesis Gen. | $9.95 | $2.00 | ~$1.50 | ~$1.20 |
| Cloning Scenarios | $5.25 | $0.97 | ~$0.30 | ~$0.25 |
| LitQA2 | $3.72 | $0.70 | ~$0.50 | ~$0.40 |
| Protocol QA | $3.16 | $0.54 | ~$0.30 | ~$0.25 |
| FigQA | $2.84 | $0.47 | ~$0.25 | ~$0.20 |
| Structure Analysis | $2.30 | $0.44 | ~$0.30 | ~$0.25 |
| SuppQA | $1.15 | $0.17 | ~$0.12 | ~$0.10 |
| Calibration | $1.11 | $0.19 | ~$0.10 | ~$0.08 |
| **Total** | **$102.81** | **$14.76** | **~$8.50** | **~$6.50** |

**Cost varies 16x across providers.** Opus 4.6 at $102.81 is 16x more expensive than Gemini 2.5 Pro at ~$6.50 for running the identical benchmark suite. Cost efficiency by model:

| Model | Total Cost | Avg Tier 1 Accuracy | Cost per Tier 1 % Point |
|---|---|---|---|
| Gemini 2.5 Pro | ~$6.50 | 14.0% | $0.46/pp |
| GPT-5.2 | ~$8.50 | 13.6% | $0.63/pp |
| Claude Sonnet 4.6 | $14.76 | 13.9% | $1.06/pp |
| Claude Opus 4.6 | $102.81 | 18.5% | $5.56/pp |

Opus achieves the highest Tier 1 accuracy but at 12x the cost per percentage point of Gemini. For Tier 2 tasks, Gemini's statistical reasoning dominance at 67.5% makes it the most cost-effective model for that category. For compositional chains, Sonnet at $14.76 total delivers the best end-to-end accuracy at a fraction of Opus's cost.

## 4. Discussion

### 4.1 No Single Model Dominates

The most important finding from the four-model comparison is that **no single model is best across all evaluation tiers.** Opus 4.6 leads Tier 1 retrieval-heavy tasks (5/7 categories). GPT-5.2 leads protocol reasoning and structural biology. Gemini 2.5 Pro dominates statistical reasoning at 67.5%---nearly 2x the next best model. And Sonnet 4.6 achieves the highest compositional chain accuracy at 73.3%. This pattern of domain specialization has practical implications: research teams should match model selection to their specific scientific task profile rather than defaulting to any single "best" model.

### 4.2 Atomic Performance Does Not Predict Compositional Ability

The gap between step-level and end-to-end accuracy persists for all four models, ranging from 15.5 pp (Sonnet) to 44.5 pp (GPT-5.2). The most striking finding is that **GPT-5.2 achieves competitive step-level accuracy (81.2%) but has the worst error propagation gap (44.5 pp)**, resulting in only 36.7% end-to-end accuracy. Conversely, Sonnet 4.6 has the highest step-level accuracy (88.8%) and the smallest propagation gap (15.5 pp), achieving 73.3% end-to-end.

This suggests that multi-step scientific reasoning draws on capabilities---consistency, coherent context integration, calibrated confidence at decision points---that are not well-predicted by single-question benchmarks. Compositional evaluation captures something that atomic evaluation misses, and the difference is not simply about accuracy levels: it is about how gracefully errors compound.

### 4.3 Domain Specialization Across Providers

The cross-provider comparison reveals clear domain-specific strengths:

- **Gemini: statistical reasoning.** At 67.5%, Gemini's advantage on test selection, p-value computation, and multiple testing correction is the largest single-category gap in the benchmark. This likely reflects differences in training data or architecture tuned for quantitative reasoning.
- **GPT-5.2: structural biology.** At 52.8% on Structure Analysis, GPT-5.2 leads both Anthropic models by 7+ percentage points. Gemini is close behind at 51.2%.
- **Opus: retrieval and cloning.** Opus leads on tasks requiring extraction from structured inputs (FigQA, DbQA, CloningScenarios).
- **Sonnet: compositional chains.** Sonnet's chain superiority (73.3% E2E) is consistent with its statistical reasoning advantage and suggests stronger error suppression during multi-step workflows.

These patterns would be invisible in a single-provider comparison.

### 4.4 Cost-Performance Tradeoffs

Cost varies 16x across models ($6.50 to $102.81 for the full benchmark). **The most expensive model is not the best for most task types.** Opus's cost premium is justified only for retrieval-heavy Tier 1 tasks where it shows significant advantages. For statistical reasoning, Gemini delivers 3x higher accuracy at 16x lower cost. For compositional chains, Sonnet delivers the highest end-to-end accuracy at $14.76---7x less than Opus. For most research teams, a mixed-model strategy matching model to task type would optimize both accuracy and cost.

### 4.5 Judge Reliability Is Quantifiable

LLM-as-judge grading is increasingly common but rarely audited. Our results show it is reliable (kappa = 0.765) but not unbiased. Position and verbosity biases, while small, are systematic and should be measured and reported.

The rubric sensitivity analysis on Hypothesis Generation demonstrates a subtler problem: even with perfect inter-judge agreement, the rubric itself can be the source of inflation. A generic "does the response address the rubric" prompt produces 100% accuracy for both Anthropic models; a strict per-criterion decomposition reveals a 67--71% novelty pass rate hiding underneath. This suggests that **benchmarks using LLM-as-judge should report rubric sensitivity alongside inter-judge agreement**---the rubric matters as much as the judge.

### 4.6 Limitations

- **Judge audit sample size.** Twenty items provide directional signal but not precise estimates of bias rates.
- **Temporal probe.** Insufficient date metadata in the dataset prevented the temporal contamination split.
- **Contamination probes limited to Anthropic.** Cloze and reverse probes were only run for Claude models. Cross-provider contamination testing would require adapting the probe methodology.
- **LLM-judge ceiling.** Calibration remains at 100% under LLM-judge grading for three of four models and may benefit from a rubric sensitivity analysis similar to Hypothesis Generation.
- **IRT with 4 models.** While improved from 2-model analysis, classical item analysis with only 4 respondents still produces noisy discrimination estimates; 8+ models would yield more reliable item parameters.
- **Cost approximations.** GPT-5.2 and Gemini cost figures are approximate (token-level tracking used estimated counts rather than exact API-reported values).

## 5. Conclusion

LABBench2-Pro demonstrates that methodological rigor and cross-provider evaluation are both essential for scientific LLM benchmarking. Evaluating four frontier models from three providers across 9,591 runs reveals a landscape of domain specialization: no single model dominates all tiers, atomic performance does not predict compositional ability, and cost varies 16x with the most expensive model not uniformly the best.

The cross-provider findings are the most actionable. Gemini 2.5 Pro's statistical reasoning dominance (67.5%, nearly 2x the next best) suggests significant architectural or training differences for quantitative tasks. Anthropic models' compositional chain superiority (Sonnet 73.3%, Opus 60.0% vs. both competitors at 36.7%) reveals capabilities---consistency, error suppression, coherent context integration---that are invisible to atomic benchmarks. GPT-5.2's structural biology lead and protocol reasoning advantage add further nuance.

For the research community, the practical implication is clear: benchmark reports based on a single model or provider are insufficient. The methodological tools we introduce---bootstrap CIs, IRT-based pruning, judge audits, contamination probes, compositional chains, and cost-accuracy analysis---provide a template for more rigorous evaluation. IRT analysis identifies that 88% of items do not discriminate between models, enabling 8x more efficient evaluations through pruning. We release all code, 9,591 evaluation traces, 30 compositional chains, and analysis scripts to support reproducible scientific LLM evaluation.

## Data Availability

- **Code:** https://github.com/VibeCodingScientist/LABBench2-Pro
- **Raw results:** `results/raw/` (9,591 eval traces across 4 models, chain execution traces, judge audits)
- **Figures:** `results/figures/` (publication-ready figures)
- **Tasks:** `tasks/` (703 generated tasks + 30 compositional chains with 95/95 verified data points)
- **Analysis:** `results/generate_all.py` (reproducible figure and table generation)
- **Machine-readable summary:** `results/summary.json`

## Acknowledgments

Built on FutureHouse's LABBench2 benchmark. Not affiliated with FutureHouse.
