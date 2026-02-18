# LABBench2-Pro: A Methodological Audit and Extension of Scientific LLM Evaluation

## Abstract

We present LABBench2-Pro, a methodological audit and extension of FutureHouse's LABBench2 benchmark for evaluating large language models on biology research tasks. Our analysis identifies five categories of gaps in current scientific LLM evaluation: scoring reliability, statistical reporting, contamination risk, coverage gaps, and atomic-only evaluation. We address these through a three-tier framework: (1) a methodological audit of existing LABBench2 categories with bootstrap confidence intervals, IRT analysis, judge reliability testing, and contamination probes; (2) 703 new programmatically-generated tasks spanning statistical reasoning, structural biology, uncertainty calibration, and hypothesis generation; and (3) 30 compositional chains (93 steps) testing multi-step scientific reasoning where errors compound. Evaluating Claude Opus 4.6 across 2,454 tasks ($102.81 total cost), we find: (a) performance varies dramatically across categories (4.7%--100%), with confidence intervals revealing that many apparent differences are not statistically significant; (b) the model achieves 89.2% step-level accuracy on compositional chains but only 58.6% end-to-end, quantifying error propagation; (c) LLM-as-judge grading shows substantial inter-judge agreement (Cohen's kappa = 0.765) with measurable position bias (10%) and verbosity bias (+5%); (d) zero contamination signal across cloze and reverse probes; and (e) rubric sensitivity analysis reveals that a lenient judge prompt inflates Hypothesis Generation from 97% (strict per-criterion rubric) to 100% (generic rubric), with novelty as the primary failure mode (67% pass rate). All code, data, model traces, and analysis are publicly available.

## 1. Introduction

Benchmarks for evaluating LLMs on scientific tasks are proliferating rapidly. FutureHouse's LABBench2 established a comprehensive evaluation framework covering literature comprehension (LitQA2), figure interpretation (FigQA), sequence analysis (SeqQA), supplementary material reasoning (SuppQA), protocol understanding (ProtocolQA), database queries (DbQA), table interpretation (TableQA), and molecular cloning scenarios (CloningScenarios). However, benchmarks themselves require benchmarking.

Our gap analysis identified five categories of methodological issues in current scientific LLM evaluation:

1. **Scoring reliability.** LLM-as-judge grading introduces unquantified noise. Position bias (does answer order affect grading?), verbosity bias (are longer answers scored more favorably?), and inter-judge disagreement are never measured.

2. **Statistical reporting.** Results are reported as point estimates without confidence intervals. With small per-category sample sizes, two models can appear different when their confidence intervals overlap entirely.

3. **Contamination risk.** No probes test whether models have memorized benchmark questions from training data.

4. **Coverage gaps.** LABBench2 tests retrieval and comprehension but not statistical reasoning, uncertainty calibration, hypothesis generation, or structural biology interpretation.

5. **Atomic-only evaluation.** Every task is independent. Real research requires multi-step reasoning where errors compound---reading a paper, interpreting its figures, choosing the right statistical test, designing a follow-up experiment.

LABBench2-Pro addresses each gap systematically.

## 2. Methods

### 2.1 Tier 1: Methodological Audit

We ran Claude Opus 4.6 against all available LABBench2 categories from HuggingFace, then applied four analyses:

**Bootstrap confidence intervals.** BCa (bias-corrected and accelerated) 95% confidence intervals computed via 10,000 bootstrap resamples on accuracy within each category.

**Item Response Theory.** Two-parameter logistic (2PL) IRT models to estimate item difficulty and discrimination parameters, identify low-discrimination items, and compute test information functions. (Requires 2+ models; single-model results reported here.)

**Judge audit.** We sampled 20 LitQA2 responses and graded each with two LLM judges (Claude Sonnet 4.5 and Claude Opus 4.6). We measured inter-judge agreement (Cohen's kappa), position bias (swapping reference/response order), and verbosity bias (comparing grades on original vs. shortened responses).

**Contamination probes.** Three probes: (1) cloze completion---can the model complete a truncated benchmark question?; (2) reverse reconstruction---can it guess the question from only the answer?; (3) temporal split---chi-squared test on accuracy for pre- vs. post-training-cutoff items (insufficient date metadata in current dataset).

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

**Hypothesis Generation** tasks provide real PubMed abstracts and ask the model to generate testable hypotheses. Graded by LLM-judge with structured rubric.

### 2.3 Tier 3: Compositional Chains (30 chains, 93 steps)

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

- **Model:** Claude Opus 4.6 (claude-opus-4-6), Anthropic API
- **Cost:** $15.00/M input tokens, $75.00/M output tokens
- **Concurrency:** 3 parallel API calls
- **Infrastructure:** PostgreSQL 16 + Redis 7 on DigitalOcean (4 vCPU, 8 GB RAM)
- **Runtime:** 5 hours 5 minutes
- **Total cost:** $102.81 (766,478 input tokens, 1,217,518 output tokens)
- **Resume safety:** All evaluations check the database before calling the API; interrupted runs resume without re-evaluating completed tasks

## 3. Results

### 3.1 Tier 1: LABBench Category Performance

| Category | n | Correct | Accuracy | 95% CI |
|---|---|---|---|---|
| CloningScenarios | 33 | 13 | 39.4% | [24.2%, 57.6%] |
| LitQA2 | 186 | 58 | 31.2% | [24.7%, 38.2%] |
| SeqQA | 585 | 104 | 17.8% | [14.9%, 21.0%] |
| ProtocolQA | 100 | 15 | 15.0% | [9.0%, 23.0%] |
| SuppQA | 81 | 9 | 11.1% | [4.9%, 19.8%] |
| FigQA | 181 | 19 | 10.5% | [6.6%, 15.5%] |
| DbQA | 492 | 23 | 4.7% | [3.0%, 6.7%] |

Key observations:

- **Confidence intervals matter.** SuppQA (11.1% [4.9%, 19.8%]) and FigQA (10.5% [6.6%, 15.5%]) have overlapping CIs---these categories cannot be meaningfully distinguished without more data. CloningScenarios, with only 33 items, has a CI spanning 33 percentage points.

- **DbQA is nearly at chance.** At 4.7% accuracy, database query tasks represent the largest performance gap, likely because these require structured query generation from natural language descriptions.

- **LitQA2 leads retrieval tasks** at 31.2%, consistent with LLMs' strength in literature comprehension.

### 3.2 Tier 2: New Task Categories

| Category | n | Correct | Accuracy |
|---|---|---|---|
| Calibration | 100 | 100 | 100.0% |
| Hypothesis Generation (lenient) | 100 | 100 | 100.0% |
| Hypothesis Generation (strict) | 100 | 97 | 97.0% |
| Structure Analysis | 303 | 138 | 45.5% |
| Statistical Reasoning | 200 | 46 | 23.0% |

**Calibration** achieved 100% accuracy, suggesting frontier models reliably identify when insufficient information is available.

**Hypothesis Generation** achieved 100% accuracy under the original lenient rubric. Re-grading with a strict per-criterion rubric (Section 3.6) reduced this to 97%, with per-criterion analysis revealing that novelty is the primary failure mode.

**Structure Analysis** at 45.5% shows moderate capability on protein structure questions and gel image interpretation. This is the strongest performance on any programmatically-verified category.

**Statistical Reasoning** at 23.0% reveals a substantial gap in quantitative scientific reasoning---test selection, p-value computation, and multiple testing correction remain challenging.

### 3.3 Tier 3: Compositional Chains

Across 30 chains (93 total steps):

- **Step-level accuracy:** 83/93 = 89.2%
- **End-to-end accuracy:** 17/29 = 58.6% (one chain had a step loading error)
- **Error propagation gap:** 30.6 percentage points between step and end-to-end accuracy

This is the central finding: **a model that gets individual steps right 89% of the time only completes entire research workflows 59% of the time.** For a 3-step chain at 89% per-step accuracy, the expected end-to-end rate is 70.5% (0.89^3), and for a 4-step chain, 62.7% (0.89^4). The observed 58.6% is slightly below the independence assumption, suggesting correlated errors.

Chain-level results by template:

| Template | Chains | E2E Correct | E2E Rate |
|---|---|---|---|
| Paper to Experiment (4 steps) | 3 | 1 | 33% |
| Structure to Drug (4 steps) | 3 | 2 | 67% |
| Stats Pipeline | 3 | 2 | 67% |
| Critical Appraisal | 3 | 2 | 67% |
| Genetics to Therapy | 3 | 1 | 33% |
| Protocol Troubleshoot | 3 | 1 | 33% |
| Paradox Resolution | 3 | 2 | 67% |
| Sequence to Function | 3 | 1 | 33% |
| Data to Mechanism | 3 | 2 | 67% |
| Evidence Synthesis | 3 | 3 | 100% |

The 4-step templates (Paper to Experiment, Structure to Drug) show lower end-to-end rates than 3-step templates, consistent with compounding error.

### 3.4 Judge Audit

From 20 sampled LitQA2 responses graded by two LLM judges:

| Metric | Value |
|---|---|
| Inter-judge agreement | 90.0% |
| Cohen's kappa | 0.765 (substantial) |
| Position bias rate | 10.0% (2/20) |
| Verbosity bias | +5.0% (verbose answers favored) |
| Sonnet positive rate | 35.0% |
| Opus positive rate | 25.0% |

Inter-judge agreement is substantial (kappa = 0.765) but not perfect. Position bias at 10% and verbosity bias at +5% are measurable and should be reported alongside any LLM-judged results.

### 3.5 Contamination Probes

| Probe | Match Rate |
|---|---|
| Cloze completion | 0/20 = 0.0% |
| Reverse question | 0/20 = 0.0% |
| Temporal split | Insufficient date metadata |

No contamination signal was detected. The model could not complete truncated questions or reconstruct questions from answers, suggesting the benchmark content has not been memorized from training data.

### 3.6 Rubric Sensitivity Analysis

The 100% accuracy on Hypothesis Generation under the original LLM-judge rubric prompted a rubric sensitivity analysis. The original grading prompt asked whether the response "contained the key factual content" of a rubric description---a low bar that any coherent response could clear. We designed a strict per-criterion rubric evaluating four dimensions independently:

| Criterion | Pass Rate | Description |
|---|---|---|
| Falsifiable | 100/100 (100%) | Hypotheses include specific experimental designs |
| Connected | 95/100 (95%) | Hypotheses logically follow from abstract findings |
| Specific | 95/100 (95%) | Hypotheses name concrete biological entities |
| Novel | 67/100 (67%) | Hypotheses go beyond obvious next steps |

A response requires 3/4 criteria to pass. Under strict grading, accuracy dropped from 100% to **97%** (97/100). The per-criterion analysis reveals that **novelty is the primary failure mode**: while the model reliably generates falsifiable, connected, and specific hypotheses, only 67% of responses propose genuinely non-obvious extensions of the source material. The remaining 33% restate the abstract's implications or suggest experiments that any domain scientist would immediately propose.

This result validates the paper's thesis about scoring reliability: a 3-percentage-point drop may seem small, but the per-criterion decomposition reveals that the lenient rubric was masking a substantial novelty deficit. The strict rubric cost $1.59 for 100 re-grades (Claude Sonnet 4.5 as judge).

### 3.7 Cost Analysis

Total evaluation cost: **$102.81** across 2,454 runs.

| Category | Cost | $/Correct Answer |
|---|---|---|
| SeqQA | $38.62 | $0.37 |
| DbQA | $13.20 | $0.57 |
| Statistical Reasoning | $11.22 | $0.24 |
| Chain Tasks | $10.30 | $0.12 |
| Hypothesis Gen. | $9.95 | $0.10 |
| Cloning Scenarios | $5.25 | $0.40 |
| LitQA2 | $3.72 | $0.06 |
| Protocol QA | $3.16 | $0.21 |
| FigQA | $2.84 | $0.15 |
| Structure Analysis | $2.30 | $0.02 |
| SuppQA | $1.15 | $0.13 |
| Calibration | $1.11 | $0.01 |

SeqQA dominates cost (37.6% of total) due to long sequence contexts. Structure Analysis offers the best cost-efficiency at $0.02 per correct answer.

## 4. Discussion

### 4.1 The Confidence Interval Problem

Point estimates without confidence intervals are the most widespread issue in LLM benchmarking. Our results demonstrate why: SuppQA (11.1%) and FigQA (10.5%) appear nearly identical, but their overlapping CIs ([4.9%, 19.8%] and [6.6%, 15.5%]) make ranking impossible. CloningScenarios at 39.4% has a CI spanning [24.2%, 57.6%]---a 33-point range that renders the point estimate nearly meaningless. Every benchmark should report confidence intervals.

### 4.2 The Compositionality Gap

The 30.6 percentage-point gap between step-level (89.2%) and end-to-end (58.6%) accuracy is the most practically significant finding. Real scientific research is inherently compositional---a literature review informs experimental design, which determines the statistical test, which shapes the interpretation. A model that answers individual questions well but fails to chain reasoning together is fundamentally limited as a research tool.

### 4.3 Judge Reliability Is Quantifiable

LLM-as-judge grading is increasingly common but rarely audited. Our results show it is reliable (kappa = 0.765) but not unbiased. Position and verbosity biases, while small, are systematic and should be measured and reported.

The rubric sensitivity analysis on Hypothesis Generation demonstrates a subtler problem: even with perfect inter-judge agreement, the rubric itself can be the source of inflation. A generic "does the response address the rubric" prompt produces 100% accuracy; a strict per-criterion decomposition reveals a 67% novelty pass rate hiding underneath. This suggests that **benchmarks using LLM-as-judge should report rubric sensitivity alongside inter-judge agreement**---the rubric matters as much as the judge.

### 4.4 Limitations

- **Single model.** We evaluate only Claude Opus 4.6. Multi-model comparison is needed for IRT analysis and relative performance ranking.
- **Judge audit sample size.** Twenty items provide directional signal but not precise estimates of bias rates.
- **Temporal probe.** Insufficient date metadata in the dataset prevented the temporal contamination split.
- **LLM-judge ceiling.** Calibration remains at 100% under LLM-judge grading and may benefit from a similar rubric sensitivity analysis. Hypothesis Generation dropped to 97% under strict grading, but the novelty criterion (67%) warrants further investigation with human expert judges.

## 5. Conclusion

LABBench2-Pro demonstrates that methodological rigor is not optional in LLM benchmarking. Confidence intervals reveal that many reported performance differences are not statistically significant. Compositional evaluation exposes a 31-point gap that atomic tasks cannot capture. Judge audits quantify biases that are otherwise invisible. Rubric sensitivity analysis shows that grading rubric design can mask substantial capability gaps---a 100% headline accuracy concealed a 67% novelty pass rate. Contamination probes provide a necessary baseline of trust. We release all code, tasks, model traces, and analysis to support reproducible scientific LLM evaluation.

## Data Availability

- **Code:** https://github.com/VibeCodingScientist/LABBench2-Pro
- **Raw results:** `results/raw/` (2,454 eval traces, 93 chain execution traces, 60 judge audits)
- **Figures:** `results/figures/` (5 publication-ready figures)
- **Tasks:** `tasks/` (703 generated tasks + 30 compositional chains with 95/95 verified data points)
- **Analysis:** `results/generate_all.py` (reproducible figure and table generation)

## Acknowledgments

Built on FutureHouse's LABBench2 benchmark. Not affiliated with FutureHouse.
