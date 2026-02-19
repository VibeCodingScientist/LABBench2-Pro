# LABBench2-Pro

A methodological audit and extension of [LABBench2](https://huggingface.co/datasets/futurehouse/labbench2) — FutureHouse's benchmark for evaluating LLMs on biology research tasks.

LABBench2 established the standard for measuring how well frontier models handle real scientific work: literature comprehension, figure interpretation, sequence analysis, protocol reasoning. But benchmarks themselves need benchmarking. LABBench2-Pro identifies and addresses specific methodological gaps in the current evaluation framework.

## Motivation

Our [gap analysis](labbench2_pro_dev_plan.md) identified five categories of issues:

1. **Scoring reliability** — LLM-as-judge grading introduces unquantified noise. Position bias (does answer order matter?), verbosity bias (are longer answers favored?), and inter-judge disagreement are never measured.

2. **Statistical reporting** — Results are reported as point estimates without confidence intervals. With small per-category sample sizes, two models can appear different when their CIs overlap entirely.

3. **Contamination risk** — No probes test whether models have memorized benchmark questions from training data.

4. **Coverage gaps** — LABBench2 tests retrieval and comprehension but not statistical reasoning, uncertainty calibration, hypothesis generation, or structural biology interpretation.

5. **Atomic-only evaluation** — Every task is independent. Real research requires multi-step reasoning where errors compound — reading a paper, interpreting its figures, choosing the right statistical test, designing a follow-up experiment.

## What LABBench2-Pro Does

### Tier 1: Methodological Audit of LABBench2

Run frontier models against the existing LABBench2 categories, then apply rigorous statistical analysis that the original benchmark lacks:

- **Bootstrap CIs** — BCa 95% confidence intervals on accuracy with pairwise significance tests (Bonferroni-corrected)
- **IRT Analysis** — 2-parameter logistic Item Response Theory to identify low-discrimination items, compute test information functions, and recommend a pruned item set
- **Judge Audit** — Run two LLM judges on the same responses, measure Cohen's kappa, test position bias (swap reference/response order) and verbosity bias (shorten responses, re-judge)
- **Contamination Probes** — Cloze completion (can the model finish a truncated question?), reverse reconstruction (can it guess the question from the answer?), temporal split (chi-squared test on pre- vs post-cutoff accuracy)

### Tier 2: New Task Categories (703 tasks)

Programmatically generated tasks with deterministic or rubric-graded ground truth:

| Category | Count | Verification | Source |
|---|---|---|---|
| Statistical Reasoning | 200 | Programmatic (scipy) | Synthetic gene expression data |
| Structure Analysis | 303 | Programmatic + LLM-judge | Real PDB structures (BioPython) + synthetic gel images |
| Uncertainty Calibration | 100 | LLM-judge | LABBench2 questions with critical info stripped |
| Hypothesis Generation | 100 | LLM-judge (rubric) | Real PubMed abstracts (NCBI Entrez) |

### Tier 3: Compositional Chains (30 chains, 96 steps)

Multi-step research workflows where each step depends on the previous answer. If the model gets step 1 wrong, it cascades. **30 hand-authored chains** across 10 template types, with 3 variants each covering different biomedical domains:

| Template | Workflow | Steps | Example Topics |
|---|---|---|---|
| Paper to Experiment | Paper finding → data interpretation → stats test → hypothesis | 4 | SHP2, JAK2 V617F, PCSK9 |
| Structure to Drug | Protein structure → binding mechanism → SAR prediction → validation | 4 | EGFR, KRAS G12C, SARS-CoV-2 Mpro |
| Stats Pipeline | Test selection → multiple testing correction → pathway interpretation | 3 | TNBC RNA-seq, T2D GWAS, scRNA-seq TILs |
| Critical Appraisal | Evaluate weak evidence → integrate conflicting data → design definitive experiment | 3 | IDH1 glioma, Lecanemab, Microbiome |
| Genetics to Therapy | Genetic finding → structural impact → therapeutic strategy | 3 | PINK1 Parkinson's, CFTR CF, SCN1A Dravet |
| Protocol Troubleshoot | Diagnose error → interpret fix → quantitative follow-up | 3 | KRAS-BRAF co-IP, ChIP-seq, CRISPR base editing |
| Paradox Resolution | Explain paradox → discriminating experiment → synthesize conclusion | 3 | ZEB1 EMT, PD-1 hyperprogression, Exercise immunosuppression |
| Sequence to Function | Identify protein → predict adaptations → design validation | 3 | Psychrophilic LDH, AmpC β-lactamase, Novel Cas effector |
| Data to Mechanism | Interpret ambiguous data → update with evidence → correct prior analysis | 3 | GAPDH confound, Imatinib resistance, Venetoclax synergy |
| Evidence Synthesis | Compare conflicting studies → meta-analysis → clinical recommendation | 3 | ctDNA lung cancer, FLT3 AML, BRAF melanoma |

All chains use real data verified against 6 databases: PDB, UniProt, ChEMBL, ClinVar, ClinicalTrials.gov, and Open Targets. See [`tasks/chains/LABBench2Pro_AllExamples.md`](tasks/chains/LABBench2Pro_AllExamples.md) for the complete chain content and [`tasks/chains/VERIFICATION_REPORT.md`](tasks/chains/VERIFICATION_REPORT.md) for the data verification audit (95/95 data points verified correct).

## Architecture

```
labbench2-pro/
├── src/
│   ├── config.py              # Model registry, API keys, cost table
│   ├── models.py              # Unified model caller (Anthropic/OpenAI/Google)
│   ├── db.py                  # Thin asyncpg wrapper
│   ├── cache.py               # Redis response cache
│   ├── api.py                 # FastAPI REST endpoints
│   │
│   ├── tier1/                 # Methodological audit
│   │   ├── run_eval.py        # Run models against tasks (HF or local)
│   │   ├── grading.py         # Programmatic + LLM-judge grading
│   │   ├── bootstrap_ci.py    # BCa CIs + pairwise tests
│   │   ├── irt_analysis.py    # 2PL IRT + test information
│   │   ├── judge_audit.py     # Inter-judge agreement + bias tests
│   │   └── contamination.py   # Cloze, reverse, temporal probes
│   │
│   ├── tier2/                 # Task generation
│   │   ├── gen_stats_tasks.py # Statistical reasoning (scipy ground truth)
│   │   ├── gen_structure.py   # PDB parsing + gel images
│   │   ├── gen_calibration.py # Uncertainty calibration
│   │   ├── gen_hypothesis.py  # PubMed-based hypothesis tasks
│   │   └── validate_tasks.py  # Schema validation
│   │
│   └── tier3/                 # Compositional chains
│       ├── run_chains.py      # Execute chains, measure error propagation
│       ├── feedback_sim.py    # Re-run with correctness signal
│       ├── gen_chains.py      # Auto-generate chains (scaffolding)
│       └── cost_tracker.py    # Pareto frontier + cost analysis
│
├── tasks/chains/              # 30 chains (96 steps), verified data
│   ├── tasks/                 # Individual step JSON files
│   ├── chain_definitions.json # Chain wiring (step order, dependencies)
│   ├── LABBench2Pro_AllExamples.md  # Full chain content (readable)
│   └── VERIFICATION_REPORT.md # Data verification audit
├── db/schema.sql              # PostgreSQL schema (5 tables)
├── docker-compose.yml         # Postgres 16 + Redis 7
└── run_all.sh                 # Single-command pipeline
```

No LangChain, no orchestration frameworks. Direct SDK calls, raw SQL, standalone scripts.

## Quick Start

### Prerequisites

- Docker (for Postgres + Redis)
- Python 3.11+
- API keys: `ANTHROPIC_API_KEY` (required for judge), plus any of `OPENAI_API_KEY`, `GOOGLE_API_KEY`

### Setup

```bash
# Clone
git clone https://github.com/VibeCodingScientist/LABBench2-Pro.git
cd LABBench2-Pro

# Environment
cp .env.example .env
# Edit .env — add your API keys

# Install
pip install -e .

# Start services
docker compose up -d
```

### Run the Full Pipeline

```bash
# Run each model (or just the ones you have API keys for)
./run_all.sh --model claude-opus-4.6
./run_all.sh --model gpt-5.2
./run_all.sh --model gemini-2.5-pro
```

This runs all 6 phases automatically:
1. Start Postgres + Redis, apply schema
2. Generate Tier 2 tasks (stats, structures, calibration, hypothesis)
3. Run Tier 1 evals against LABBench2 categories from HuggingFace
4. Run Tier 2 generated tasks through eval
5. Run compositional chains
6. Analysis: bootstrap CIs, IRT, judge audit, contamination probes, cost summary

Results are stored in PostgreSQL. Resume-safe — if interrupted, re-running skips already-completed tasks.

### Run Individual Components

```bash
# Single category eval
python -m src.tier1.run_eval --model claude-opus-4.6 --category LitQA2 --concurrency 5

# Bootstrap CIs for a category
python -m src.tier1.bootstrap_ci --category LitQA2

# IRT analysis across all results
python -m src.tier1.irt_analysis

# Judge audit (sample 50 responses)
python -m src.tier1.judge_audit --category LitQA2 --sample-size 50

# Contamination probes
python -m src.tier1.contamination --model claude-opus-4.6

# Generate tasks
python -m src.tier2.gen_stats_tasks --output-dir tasks/stats_reasoning --count 200
python -m src.tier2.gen_structure --output-dir tasks/structures --pdb-count 50 --gel-count 60

# Run a specific chain
python -m src.tier3.run_chains --model claude-opus-4.6 --chain chain01

# Cost summary
python -m src.tier3.cost_tracker
```

### Query Results

```bash
# Direct SQL
PGPASSWORD=dev psql -h localhost -p 5433 -U dev -d labbench2pro

# Or start the API
uvicorn src.api:app --host 0.0.0.0 --port 8000
# GET /results/eval?category=LitQA2
# GET /results/ci?category=LitQA2
# GET /status
```

## Models Tested

| Model | Provider | Input $/1M | Output $/1M | Status |
|---|---|---|---|---|
| claude-opus-4.6 | Anthropic | $15.00 | $75.00 | Complete |
| claude-sonnet-4.6 | Anthropic | $3.00 | $15.00 | Complete |
| gpt-5.2 | OpenAI | $1.75 | $14.00 | Complete |
| gemini-2.5-pro | Google | $1.25 | $10.00 | Complete |

## Results (4 Models, 3 Providers — Feb 2026)

Cross-provider comparison: 9,591 eval runs, ~$132.57 combined cost. See [`results/RESULTS.md`](results/RESULTS.md) for complete analysis.

### Tier 1: LABBench Categories

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|
| CloningScenarios | **39.4%** | 27.3% | 15.2% | 27.3% | Opus |
| LitQA2 | **31.2%** | 28.3% | 30.2% | 23.1% | Opus |
| SeqQA | **17.8%** | 12.4% | 10.8% | 15.5% | Opus |
| ProtocolQA | 15.0% | 15.7% | **18.5%** | 13.9% | GPT-5.2 |
| SuppQA | 11.1% | 6.2% | **12.2%** | 7.3% | GPT-5.2 |
| FigQA | **10.5%** | 5.0% | 6.1% | 8.8% | Opus |
| DbQA | **4.7%** | 2.2% | 1.7% | 2.1% | Opus |

### Tier 2: New Tasks

| Category | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro | Best |
|---|---|---|---|---|---|
| Calibration | **100%** | **100%** | **100%** | 94.0% | Tie |
| Hypothesis Gen. | 97.0% | 97.0% | **100%** | 95.0% | GPT-5.2 |
| Structure Analysis | 45.5% | 43.9% | **52.8%** | 51.2% | GPT-5.2 |
| Statistical Reasoning | 23.0% | 35.0% | 19.5% | **67.5%** | Gemini |

### Tier 3: Compositional Chains

| Metric | Opus 4.6 | Sonnet 4.6 | GPT-5.2 | Gemini 2.5 Pro |
|---|---|---|---|---|
| Step-level accuracy | 83.9% | **88.8%** | 81.2% | 71.9% |
| End-to-end accuracy | 60.0% | **73.3%** | 36.7% | 36.7% |
| Error propagation gap | 23.9 pp | **15.5 pp** | 44.5 pp | 35.2 pp |
| Total cost | $102.81 | $14.76 | ~$8.50 | **~$6.50** |

### Key Findings

1. **No single model dominates.** Opus leads Tier 1 retrieval (5/7 categories), GPT-5.2 leads structural biology, Gemini crushes statistical reasoning (67.5%), and Sonnet wins compositional chains.
2. **Anthropic models dominate multi-step workflows.** Sonnet 4.6 (73.3% E2E) and Opus 4.6 (60.0%) dramatically outperform GPT-5.2 and Gemini (both 36.7%) on compositional chains.
3. **Gemini's statistical reasoning is a standout.** 67.5% — nearly 2x the next best model (Sonnet at 35.0%). The largest single-category advantage in the entire benchmark.
4. **Atomic task performance does not predict compositional ability.** GPT-5.2 has competitive step-level accuracy (81.2%) but the worst error propagation gap (44.5pp).
5. **Cost varies 16x across providers.** Opus costs $102.81 vs Gemini's ~$6.50 for the full benchmark. For most research teams, non-Opus models offer better value.
6. **88% of benchmark items don't discriminate between models.** IRT analysis recommends pruning to 303 high-discrimination items for 8x more efficient evaluations.
7. **Zero contamination detected.** 0% cloze match rates for Anthropic models.

### Methodological Audit

| Metric | Value |
|---|---|
| Inter-judge agreement | 90.0% |
| Cohen's kappa | 0.765 |
| Position bias | 15.0% |
| Verbosity bias | +5.0% |
| IRT items analyzed | 2,522 |
| High-discrimination items | 303 (12.0%) |
| Cloze contamination | 0.0% |

## Design Principles

- **Simplest thing that works.** No abstraction until it's needed twice.
- **Every script runnable standalone.** No hidden dependencies between modules.
- **Provider-agnostic.** Unified `call_model()` interface works with Anthropic, OpenAI, and Google.
- **Resume-safe.** Every eval checks the DB before calling the API. Interrupted runs pick up where they left off.
- **Raw SQL, no ORM.** Five tables, four indexes. Schema changes = drop and recreate (results are reproducible).
- **Cost-aware.** Every API call is tracked. Cost-accuracy Pareto frontier identifies dominated models.

## Artifacts

| Directory | Contents |
|---|---|
| `results/raw/` | All 9,591 eval traces, chain traces, 60 judge audits (CSV) |
| `results/figures/` | Publication-ready figures (PDF + PNG) |
| `results/tables/` | LaTeX tables + supplementary chain traces |
| `results/summary.json` | Machine-readable results summary |
| `paper/` | Draft manuscript |
| `tasks/chains/` | 30 chains with verification report |

## Contributing

All 30 compositional chains are authored and verified. Contributions welcome for:
- Additional chain variants (new biomedical topics using the 10 existing templates)
- New chain templates (novel multi-step reasoning patterns)
- Running the benchmark against additional models
- Running with tool-augmented models to measure skill-based improvement

## Citation

```bibtex
@misc{labbench2pro2026,
  title={LABBench2-Pro: A Methodological Audit and Extension of Scientific LLM Evaluation},
  author={Weidener, Lukas},
  year={2026},
  url={https://github.com/VibeCodingScientist/LABBench2-Pro}
}
```

## License

Research use. Not affiliated with FutureHouse.
