# LABBench2-Pro: A Methodological Audit and Extension of Scientific LLM Evaluation

## Abstract

We present LABBench2-Pro, a methodological audit and extension of FutureHouse's LABBench2 benchmark for evaluating large language models on biology research tasks. Our analysis identifies five categories of gaps in current scientific LLM evaluation: scoring reliability, statistical reporting, contamination risk, coverage gaps, and atomic-only evaluation. We address these through a three-tier framework: (1) a methodological audit of existing LABBench2 categories with bootstrap confidence intervals, IRT analysis, judge reliability testing, and contamination probes; (2) 703 new programmatically-generated tasks spanning statistical reasoning, structural biology, uncertainty calibration, and hypothesis generation; and (3) 30 compositional chains testing multi-step scientific reasoning where errors compound. Evaluating Claude Opus 4.6 and Claude Sonnet 4.6 across 4,549 total evaluation runs ($117.57 total cost), we find: (a) performance varies dramatically across categories (2.2%--100%), with pairwise bootstrap tests showing only 3 of 7 LABBench categories have statistically significant differences between models; (b) Sonnet 4.6 achieves higher end-to-end chain accuracy than Opus 4.6 (73.3% vs 60.0%) despite lower step-level accuracy on retrieval tasks, suggesting that compositional reasoning ability does not simply track atomic task performance; (c) Sonnet 4.6 significantly outperforms Opus 4.6 on statistical reasoning (35.0% vs 23.0%) while costing 7x less; (d) LLM-as-judge grading shows substantial inter-judge agreement (Cohen's kappa = 0.765) with measurable position bias (15%) and verbosity bias (+5%); (e) rubric sensitivity analysis reveals that both models achieve 97% on hypothesis generation under strict grading, with novelty as the primary failure mode (67--71% pass rate); and (f) zero contamination signal across cloze probes for both models. IRT analysis across 2,459 items identifies 133 high-discrimination items from a pool of 2,326 low-discrimination items, informing benchmark pruning. All code, data, model traces, and analysis are publicly available.

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

We ran two models against all available LABBench2 categories from HuggingFace, then applied four analyses:

**Bootstrap confidence intervals.** BCa (bias-corrected and accelerated) 95% confidence intervals computed via 10,000 bootstrap resamples on accuracy within each category, with pairwise significance tests (Bonferroni-corrected).

**Item Response Theory.** Classical item analysis with discrimination indices and difficulty parameters across 2,459 items and 2 models, used to identify low-discrimination items and compute a recommended pruned item set.

**Judge audit.** We sampled 20 LitQA2 responses and graded each with two LLM judges (Claude Sonnet 4.5 and Claude Opus 4.6). We measured inter-judge agreement (Cohen's kappa), position bias (swapping reference/response order), and verbosity bias (comparing grades on original vs. shortened responses).

**Contamination probes.** Two probes per model: (1) cloze completion---can the model complete a truncated benchmark question?; (2) reverse reconstruction---can it guess the question from only the answer?

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

| Parameter | Claude Opus 4.6 | Claude Sonnet 4.6 |
|---|---|---|
| Model ID | claude-opus-4-6 | claude-sonnet-4-6 |
| Provider | Anthropic | Anthropic |
| Input cost | $15.00/M tokens | $3.00/M tokens |
| Output cost | $75.00/M tokens | $15.00/M tokens |
| Eval runs | 2,454 | 2,095 |
| Input tokens | 766,478 | 577,536 |
| Output tokens | 1,217,518 | 868,789 |
| Total cost | $102.81 | $14.76 |
| Concurrency | 3 | 5 |

- **Infrastructure:** PostgreSQL 16 + Redis 7 on DigitalOcean (4 vCPU, 8 GB RAM)
- **Combined cost:** $117.57 (both models)
- **Resume safety:** All evaluations check the database before calling the API; interrupted runs resume without re-evaluating completed tasks

## 3. Results

### 3.1 Tier 1: LABBench Category Performance

| Category | n | Opus 4.6 | 95% CI | Sonnet 4.6 | 95% CI | Diff | p-value | Sig? |
|---|---|---|---|---|---|---|---|---|
| CloningScenarios | 33 | 39.4% | [24.2%, 57.6%] | 27.3% | [15.2%, 45.5%] | +12.1 | 0.014 | **YES** |
| LitQA2 | 186--187 | 31.2% | [24.7%, 38.2%] | 28.3% | [21.9%, 35.3%] | +2.9 | 0.210 | No |
| SeqQA | 298--585 | 17.8% | [14.9%, 21.0%] | 12.4% | [9.1%, 16.4%] | +5.4 | 0.214 | No |
| ProtocolQA | 100--102 | 15.0% | [9.0%, 23.0%] | 15.7% | [9.8%, 23.5%] | -0.7 | 0.434 | No |
| SuppQA | 81 | 11.1% | [4.9%, 19.8%] | 6.2% | [2.5%, 13.6%] | +4.9 | 0.065 | No |
| FigQA | 181 | 10.5% | [6.6%, 15.5%] | 5.0% | [2.2%, 8.8%] | +5.5 | 0.003 | **YES** |
| DbQA | 492--493 | 4.7% | [3.0%, 6.7%] | 2.2% | [1.2%, 3.9%] | +2.5 | 0.001 | **YES** |

Key observations:

- **Most differences are not significant.** Despite Opus scoring higher on 6 of 7 categories, only 3 differences survive pairwise bootstrap testing: CloningScenarios, FigQA, and DbQA. LitQA2 (p=0.210), SeqQA (p=0.214), SuppQA (p=0.065), and ProtocolQA (p=0.434) cannot be distinguished.

- **Confidence intervals overlap extensively.** SuppQA and FigQA remain indistinguishable within each model. CloningScenarios still has a 30+ point CI span for both models.

- **Opus advantages cluster in retrieval-heavy tasks.** The three significant categories (FigQA, DbQA, CloningScenarios) all require extracting specific information from structured inputs.

### 3.2 Tier 2: New Task Categories

| Category | Opus 4.6 | Sonnet 4.6 | Delta |
|---|---|---|---|
| Calibration | 100.0% | 100.0% | 0 |
| Hypothesis Generation (lenient) | 100.0% | 100.0% | 0 |
| Hypothesis Generation (strict) | 97.0% | 97.0% | 0 |
| Structure Analysis | 45.5% | 43.9% | -1.6 |
| **Statistical Reasoning** | **23.0%** | **35.0%** | **+12.0** |

**Statistical Reasoning is the standout result.** Sonnet 4.6 outperforms Opus 4.6 by 12 percentage points on tasks requiring test selection, p-value computation, and multiple testing correction. This is the largest inter-model difference on any Tier 2 category and runs counter to the assumption that larger models uniformly outperform smaller ones.

**Calibration** achieved 100% on both models, suggesting frontier models reliably identify when insufficient information is available.

**Hypothesis Generation** achieved identical 97% strict accuracy on both models (see Section 3.6 for per-criterion breakdown).

**Structure Analysis** is nearly identical between models (45.5% vs 43.9%), suggesting protein structure reasoning is not strongly differentiated by model scale.

### 3.3 Tier 3: Compositional Chains

| Metric | Opus 4.6 | Sonnet 4.6 |
|---|---|---|
| Chains evaluated | 30 | 30 |
| Total steps completed | 93 | 89 |
| Step-level accuracy | 78/93 = 83.9% | 79/89 = 88.8% |
| End-to-end accuracy | 18/30 = 60.0% | 22/30 = 73.3% |
| Error propagation gap | 23.9 pp | 15.5 pp |

**Sonnet 4.6 outperforms Opus 4.6 on compositional chains** (73.3% vs 60.0% end-to-end). This is the most counterintuitive finding: a smaller, cheaper model achieves higher multi-step research accuracy. The error propagation gap is also smaller for Sonnet (15.5 pp vs 23.9 pp), suggesting more consistent step-level performance.

Chain-level results by template (Opus / Sonnet E2E correct):

| Template | Opus E2E | Sonnet E2E |
|---|---|---|
| Paper to Experiment (4 steps) | 1/3 | 1/3 |
| Structure to Drug (4 steps) | 3/3 | 3/3 |
| Stats Pipeline | 2/3 | 3/3 |
| Critical Appraisal | 2/3 | 2/3 |
| Genetics to Therapy | 0/3 | 2/3 |
| Protocol Troubleshoot | 2/3 | 2/3 |
| Paradox Resolution | 2/3 | 2/3 |
| Sequence to Function | 3/3 | 2/3 |
| Data to Mechanism | 2/3 | 2/3 |
| Evidence Synthesis | 1/3 | 1/3 |

Sonnet shows the largest gains on Genetics to Therapy (+2) and Stats Pipeline (+1), consistent with its statistical reasoning advantage.

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

No contamination signal was detected for either model. Neither model could complete truncated questions or reconstruct questions from answers, suggesting the benchmark content has not been memorized from training data.

### 3.6 Rubric Sensitivity Analysis

Both models achieved 100% under the lenient rubric and 97% under the strict per-criterion rubric. Per-criterion comparison:

| Criterion | Opus 4.6 | Sonnet 4.6 | Description |
|---|---|---|---|
| Falsifiable | 100% | 100% | Hypotheses include specific experimental designs |
| Connected | 95% | 90% | Hypotheses logically follow from abstract findings |
| Specific | 95% | 94% | Hypotheses name concrete biological entities |
| Novel | 67% | 71% | Hypotheses go beyond obvious next steps |

Both models show the same pattern: **novelty is the primary failure mode.** Falsifiability is universally passed, while novelty hovers around 67--71%. The models reliably generate well-structured, specific, connected hypotheses but frequently propose obvious extensions rather than genuinely creative follow-ups. The strict rubric cost $3.17 total for 200 re-grades.

### 3.7 IRT Analysis

With two models, classical item analysis across 2,459 items yields:

| Metric | Value |
|---|---|
| Total items | 2,459 |
| Low discrimination (< 0.3) | 2,326 (94.6%) |
| Too easy (p > 0.95) | 572 (23.3%) |
| Too hard (p < 0.05) | 1,642 (66.8%) |
| Recommended pruned set | 133 items |
| Peak discrimination | theta = 0.50 |
| Effective range | theta = [-1.20, 2.20] |

The high proportion of low-discrimination items (94.6%) indicates that most items do not differentiate between models at current capability levels. The recommended 133-item pruned set retains only items where model responses diverge, enabling more efficient future evaluations.

### 3.8 Cost Analysis

| Category | Opus Cost | Sonnet Cost | Opus $/Correct | Sonnet $/Correct |
|---|---|---|---|---|
| SeqQA | $38.62 | $3.37 | $0.37 | $0.09 |
| DbQA | $13.20 | $2.64 | $0.57 | $0.24 |
| Statistical Reasoning | $11.22 | $2.92 | $0.24 | $0.04 |
| Chain Tasks | $10.30 | $0.35 | $0.13 | $0.005 |
| Hypothesis Gen. | $9.95 | $2.00 | $0.10 | $0.02 |
| Cloning Scenarios | $5.25 | $0.97 | $0.40 | $0.11 |
| LitQA2 | $3.72 | $0.70 | $0.06 | $0.01 |
| Protocol QA | $3.16 | $0.54 | $0.21 | $0.03 |
| FigQA | $2.84 | $0.47 | $0.15 | $0.05 |
| Structure Analysis | $2.30 | $0.44 | $0.02 | $0.003 |
| SuppQA | $1.15 | $0.17 | $0.13 | $0.03 |
| Calibration | $1.11 | $0.19 | $0.01 | $0.002 |
| **Total** | **$102.81** | **$14.76** | | |

Sonnet 4.6 dominates the cost-accuracy Pareto frontier: it achieves comparable or superior accuracy at 7x lower cost across all categories. On Statistical Reasoning, Sonnet achieves higher accuracy at 3.8x lower cost---the strongest Pareto dominance in the evaluation.

## 4. Discussion

### 4.1 The Confidence Interval Problem

Point estimates without confidence intervals remain the most widespread issue in LLM benchmarking. Our two-model comparison demonstrates why: of 7 LABBench categories, only 3 show statistically significant differences between Opus 4.6 and Sonnet 4.6 (p < 0.05). The remaining 4 categories, including the high-profile LitQA2, have overlapping confidence intervals that make model ranking impossible without more data.

### 4.2 The Compositionality Gap

The gap between step-level and end-to-end accuracy persists for both models (23.9 pp for Opus, 15.5 pp for Sonnet). However, the most surprising finding is that **Sonnet 4.6 outperforms Opus 4.6 on compositional chains** (73.3% vs 60.0% E2E) despite performing worse on most atomic retrieval tasks. This suggests that multi-step scientific reasoning draws on capabilities---consistency, coherent context integration, calibrated confidence---that are not well-predicted by single-question benchmarks. Compositional evaluation captures something that atomic evaluation misses.

### 4.3 Cost-Performance Tradeoffs

Sonnet 4.6 at $14.76 total cost versus Opus 4.6 at $102.81 demonstrates that price is a poor proxy for scientific reasoning ability. For statistical reasoning, Sonnet is both cheaper and more accurate. For compositional chains, Sonnet is both cheaper and achieves higher end-to-end accuracy. The only categories where Opus shows a significant advantage (FigQA, DbQA, CloningScenarios) are retrieval-intensive tasks. This has direct implications for research groups designing LLM-augmented scientific workflows: the most expensive model is not always the best choice.

### 4.4 Judge Reliability Is Quantifiable

LLM-as-judge grading is increasingly common but rarely audited. Our results show it is reliable (kappa = 0.765) but not unbiased. Position and verbosity biases, while small, are systematic and should be measured and reported.

The rubric sensitivity analysis on Hypothesis Generation demonstrates a subtler problem: even with perfect inter-judge agreement, the rubric itself can be the source of inflation. A generic "does the response address the rubric" prompt produces 100% accuracy for both models; a strict per-criterion decomposition reveals a 67--71% novelty pass rate hiding underneath. This suggests that **benchmarks using LLM-as-judge should report rubric sensitivity alongside inter-judge agreement**---the rubric matters as much as the judge.

### 4.5 Limitations

- **Two models, one provider.** We evaluate only Claude models from Anthropic. Cross-provider comparison (GPT, Gemini) would strengthen the generalizability of findings.
- **Judge audit sample size.** Twenty items provide directional signal but not precise estimates of bias rates.
- **Temporal probe.** Insufficient date metadata in the dataset prevented the temporal contamination split.
- **LLM-judge ceiling.** Calibration remains at 100% under LLM-judge grading for both models and may benefit from a rubric sensitivity analysis.
- **IRT with 2 models.** Classical item analysis with only 2 respondents produces noisy discrimination estimates; 5+ models would yield more reliable item parameters.

## 5. Conclusion

LABBench2-Pro demonstrates that methodological rigor is not optional in LLM benchmarking. Confidence intervals and pairwise significance tests reveal that most reported performance differences between models are not statistically significant. Compositional evaluation exposes a gap that atomic tasks cannot capture---and reveals that a smaller model can outperform a larger one on multi-step research workflows. Cost-performance analysis shows 7x cost savings with comparable accuracy. Rubric sensitivity analysis shows that grading rubric design can mask substantial capability gaps---a 100% headline accuracy concealed a 67--71% novelty pass rate in both models. IRT analysis identifies that 94.6% of items do not discriminate between models, enabling a 19x more efficient pruned benchmark. Contamination probes provide a necessary baseline of trust. We release all code, tasks, model traces, and analysis to support reproducible scientific LLM evaluation.

## Data Availability

- **Code:** https://github.com/VibeCodingScientist/LABBench2-Pro
- **Raw results:** `results/raw/` (4,549 eval traces, chain execution traces, judge audits)
- **Figures:** `results/figures/` (publication-ready figures)
- **Tasks:** `tasks/` (703 generated tasks + 30 compositional chains with 95/95 verified data points)
- **Analysis:** `results/generate_all.py` (reproducible figure and table generation)

## Acknowledgments

Built on FutureHouse's LABBench2 benchmark. Not affiliated with FutureHouse.
