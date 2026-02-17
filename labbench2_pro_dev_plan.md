# LABBench2-Pro — Technical Development Plan

**Audience:** Developer implementing the evaluation infrastructure.
**Principle:** Simplest thing that works. No abstraction until the second time you need it. Every script should be runnable standalone.

---

## Stack

| Layer | Choice | Why |
|-------|--------|-----|
| Language | Python 3.11+ | Everything lives here |
| DB | PostgreSQL 16 | Single source of truth for all results |
| Cache / Queue | Redis 7 | Cache API responses, queue eval jobs |
| API | FastAPI | Thin REST layer to trigger runs + query results |
| Model SDKs | `openai`, `anthropic`, `google-genai` | Direct SDK calls, no wrapper frameworks |
| Data | HuggingFace `datasets` | Load LABBench2 tasks |
| Stats | `scipy`, `numpy`, `pandas` | Bootstrap CIs, basic stats |
| IRT | `py-irt` | Item Response Theory models |
| Bio tools | `biopython` | PDB parsing, sequence work |
| Plotting | `matplotlib`, `seaborn` | Figures for the paper |

No LangChain, no LlamaIndex, no orchestration frameworks. Direct SDK calls wrapped in thin async functions.

---

## Project Structure

```
labbench2-pro/
├── README.md
├── pyproject.toml              # single dependency file, use uv or pip
├── .env.example                # API keys template
├── docker-compose.yml          # postgres + redis
│
├── db/
│   └── schema.sql              # one file, all tables
│
├── src/
│   ├── config.py               # env vars, model registry
│   ├── db.py                   # thin postgres wrapper (asyncpg)
│   ├── cache.py                # redis get/set helpers
│   ├── models.py               # unified model caller
│   ├── api.py                  # FastAPI app
│   │
│   ├── tier1/
│   │   ├── run_eval.py         # run models against LABBench2, save to DB
│   │   ├── bootstrap_ci.py     # compute CIs from DB results
│   │   ├── irt_analysis.py     # fit IRT model
│   │   ├── judge_audit.py      # LLM-judge concordance
│   │   └── contamination.py    # contamination probes
│   │
│   ├── tier2/
│   │   ├── gen_stats_tasks.py  # generate statistical reasoning tasks
│   │   ├── gen_structure.py    # generate PDB / gel image tasks
│   │   ├── gen_calibration.py  # generate uncertainty tasks
│   │   ├── gen_hypothesis.py   # generate hypothesis tasks
│   │   └── validate_tasks.py   # verify generated tasks have correct answers
│   │
│   └── tier3/
│       ├── chains.py           # define compositional task chains
│       ├── run_chains.py       # execute chains, measure decay
│       ├── feedback_sim.py     # feedback loop simulation
│       └── cost_tracker.py     # token/cost accounting
│
├── tasks/                      # generated task data (Tier 2 outputs)
│   ├── stats_reasoning/
│   ├── structures/
│   ├── calibration/
│   └── hypothesis/
│
├── notebooks/                  # analysis + paper figures
│   ├── 01_bootstrap_results.ipynb
│   ├── 02_irt_results.ipynb
│   ├── 03_judge_audit.ipynb
│   ├── 04_contamination.ipynb
│   ├── 05_chain_decay.ipynb
│   └── 06_cost_frontier.ipynb
│
└── scripts/
    ├── setup_db.sh             # createdb + run schema.sql
    └── run_tier1.sh            # convenience: run all Tier 1 in sequence
```

---

## Database Schema

One schema. No migrations framework needed — just `schema.sql`. If the schema changes, drop and recreate (results are reproducible from reruns).

```sql
-- schema.sql

CREATE TABLE IF NOT EXISTS models (
    id          SERIAL PRIMARY KEY,
    provider    TEXT NOT NULL,           -- 'openai', 'anthropic', 'google'
    model_name  TEXT NOT NULL,           -- 'gpt-5.2-pro', 'claude-opus-4.5', etc.
    tools       BOOLEAN DEFAULT FALSE,  -- with or without tool augmentation
    UNIQUE(provider, model_name, tools)
);

CREATE TABLE IF NOT EXISTS tasks (
    id          TEXT PRIMARY KEY,        -- original LABBench2 task ID
    category    TEXT NOT NULL,           -- 'LitQA3', 'FigQA2', 'SeqQA2', etc.
    variant     TEXT,                    -- 'img', 'pdf', 'retrieve', 'inject', 'file'
    source      TEXT DEFAULT 'labbench2', -- 'labbench2' or 'pro' for new tasks
    meta        JSONB                   -- any extra metadata (source paper DOI, etc.)
);

CREATE TABLE IF NOT EXISTS eval_runs (
    id          SERIAL PRIMARY KEY,
    model_id    INTEGER REFERENCES models(id),
    task_id     TEXT REFERENCES tasks(id),
    response    TEXT,                    -- raw model response
    correct     BOOLEAN,                -- graded result
    score       FLOAT,                  -- partial credit if applicable
    grader      TEXT,                   -- 'programmatic', 'llm-judge-claude', etc.
    tokens_in   INTEGER,
    tokens_out  INTEGER,
    cost_usd    FLOAT,
    latency_ms  INTEGER,
    run_at      TIMESTAMPTZ DEFAULT NOW(),
    meta        JSONB                   -- extra (e.g., judge explanation)
);

CREATE TABLE IF NOT EXISTS chain_runs (
    id          SERIAL PRIMARY KEY,
    chain_id    TEXT NOT NULL,           -- e.g., 'chain_litqa_to_primer_01'
    model_id    INTEGER REFERENCES models(id),
    step_num    INTEGER NOT NULL,
    task_id     TEXT REFERENCES tasks(id),
    input_from  TEXT,                   -- previous step's response (or NULL for step 1)
    response    TEXT,
    correct     BOOLEAN,
    run_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE IF NOT EXISTS judge_audits (
    id           SERIAL PRIMARY KEY,
    eval_run_id  INTEGER REFERENCES eval_runs(id),
    judge_model  TEXT NOT NULL,          -- which model judged
    judge_score  BOOLEAN,
    order_variant TEXT,                  -- 'original', 'reversed' (for position bias)
    length_variant TEXT,                 -- 'original', 'shortened', 'lengthened'
    explanation  TEXT,
    run_at       TIMESTAMPTZ DEFAULT NOW()
);

-- Indexes for common queries
CREATE INDEX idx_eval_model_task ON eval_runs(model_id, task_id);
CREATE INDEX idx_eval_task_correct ON eval_runs(task_id, correct);
CREATE INDEX idx_chain_model ON chain_runs(chain_id, model_id);
CREATE INDEX idx_tasks_category ON tasks(category);
```

---

## Redis Usage

Keep it simple. Two use cases only:

**1. API response cache** — avoid re-calling a model for the same prompt.

```
Key:    cache:{model_name}:{sha256(prompt)}
Value:  JSON string of response
TTL:    30 days
```

**2. Job queue** — for batch eval runs.

```
Key:    queue:eval_jobs
Type:   Redis list (LPUSH / BRPOP)
Item:   JSON {"model_id": 1, "task_ids": ["t001", "t002", ...]}
```

No Celery. No RQ. Just a simple worker loop:

```python
# worker.py (pseudocode)
async def worker():
    while True:
        job = await redis.brpop("queue:eval_jobs")
        job = json.loads(job)
        await run_eval_batch(job["model_id"], job["task_ids"])
```

---

## Core Module: Unified Model Caller

One function to call any model. No abstraction beyond this.

```python
# src/models.py

import os, hashlib, json
from openai import AsyncOpenAI
from anthropic import AsyncAnthropic
from google import genai

openai_client = AsyncOpenAI()
anthropic_client = AsyncAnthropic()
google_client = genai.Client()

MODEL_REGISTRY = {
    "gpt-5.2-pro":      {"provider": "openai",    "model": "gpt-5.2-pro"},
    "claude-opus-4.5":   {"provider": "anthropic", "model": "claude-opus-4-5-20250929"},
    "gemini-3-pro":      {"provider": "google",    "model": "gemini-3.0-pro"},
    # add more as needed
}

# Cost per 1M tokens (input/output) — update as pricing changes
COST_TABLE = {
    "gpt-5.2-pro":     {"input": 2.50, "output": 10.00},
    "claude-opus-4.5":  {"input": 15.00, "output": 75.00},
    "gemini-3-pro":     {"input": 1.25, "output": 5.00},
}


async def call_model(
    model_key: str,
    prompt: str,
    system: str = "",
    images: list[bytes] | None = None,
    max_tokens: int = 2048,
    cache: Redis | None = None,
) -> dict:
    """
    Returns: {
        "response": str,
        "tokens_in": int,
        "tokens_out": int,
        "cost_usd": float,
        "latency_ms": int,
    }
    """
    # Check cache first
    if cache:
        cache_key = f"cache:{model_key}:{hashlib.sha256(prompt.encode()).hexdigest()}"
        cached = await cache.get(cache_key)
        if cached:
            return json.loads(cached)

    info = MODEL_REGISTRY[model_key]
    start = time.monotonic()

    if info["provider"] == "openai":
        resp = await openai_client.chat.completions.create(
            model=info["model"],
            messages=_build_openai_messages(system, prompt, images),
            max_tokens=max_tokens,
        )
        text = resp.choices[0].message.content
        t_in, t_out = resp.usage.prompt_tokens, resp.usage.completion_tokens

    elif info["provider"] == "anthropic":
        resp = await anthropic_client.messages.create(
            model=info["model"],
            system=system,
            messages=_build_anthropic_messages(prompt, images),
            max_tokens=max_tokens,
        )
        text = resp.content[0].text
        t_in, t_out = resp.usage.input_tokens, resp.usage.output_tokens

    elif info["provider"] == "google":
        resp = await google_client.aio.models.generate_content(
            model=info["model"],
            contents=_build_google_messages(system, prompt, images),
        )
        text = resp.text
        t_in = resp.usage_metadata.prompt_token_count
        t_out = resp.usage_metadata.candidates_token_count

    elapsed = int((time.monotonic() - start) * 1000)
    costs = COST_TABLE[model_key]
    cost = (t_in * costs["input"] + t_out * costs["output"]) / 1_000_000

    result = {
        "response": text,
        "tokens_in": t_in,
        "tokens_out": t_out,
        "cost_usd": round(cost, 6),
        "latency_ms": elapsed,
    }

    if cache:
        await cache.set(cache_key, json.dumps(result), ex=86400 * 30)

    return result
```

The `_build_*_messages` helpers format the prompt + images per provider. Keep them in the same file. ~20 lines each.

---

## Tier 1 — Implementation Details

### 1.1 `run_eval.py` — Run models against LABBench2

```
Usage: python -m src.tier1.run_eval --model gpt-5.2-pro --category LitQA3 --concurrency 8
```

Steps:
1. Load tasks from HuggingFace (`datasets.load_dataset("futurehouse/labbench2", "LitQA3")`)
2. Insert tasks into `tasks` table (idempotent — skip if exists)
3. For each task, call `call_model()` with the task prompt
4. Grade the response (see grading section below)
5. Insert result into `eval_runs`

Run all models × all categories. This is the foundation everything else builds on.

**Concurrency:** Use `asyncio.Semaphore(concurrency)` to limit parallel API calls. No need for anything fancier.

### 1.2 `bootstrap_ci.py` — Confidence Intervals

```
Usage: python -m src.tier1.bootstrap_ci --category LitQA3
```

Steps:
1. Query `eval_runs` for all models on the given category
2. For each model, get the binary vector of `correct` values
3. Compute 95% BCa bootstrap CI:

```python
from scipy.stats import bootstrap

def compute_ci(correct_vector: list[bool]) -> tuple[float, float, float]:
    arr = np.array(correct_vector, dtype=float)
    result = bootstrap(
        (arr,), np.mean, n_resamples=10000,
        confidence_level=0.95, method='BCa'
    )
    return arr.mean(), result.confidence_interval.low, result.confidence_interval.high
```

4. Output: CSV table + matplotlib figure with error bars.

### 1.3 `irt_analysis.py` — Item Response Theory

```
Usage: python -m src.tier1.irt_analysis
```

Steps:
1. Query all `eval_runs`, build a response matrix: rows = tasks, cols = models, values = 0/1
2. Fit 2PL IRT model:

```python
from py_irt.models import TwoParamLog

model = TwoParamLog(response_matrix)
model.fit()
difficulties = model.item_difficulties
discriminations = model.item_discriminations
```

3. Flag items with discrimination < 0.3 as "low signal"
4. Output: item parameter distributions, recommended pruned set

### 1.4 `judge_audit.py` — LLM Judge Concordance

```
Usage: python -m src.tier1.judge_audit --category SourceQualQA --sample-size 50
```

Steps:
1. Sample N tasks from the category
2. For each task, get the existing model response from `eval_runs`
3. Send to two different LLM judges (e.g., Claude and GPT) with the grading rubric
4. For position bias: also send the response with answer options reversed
5. For verbosity bias: send a shortened paraphrase of the same response
6. Insert all results into `judge_audits`
7. Compute: Cohen's κ between judges, position bias rate, verbosity correlation

### 1.5 `contamination.py` — Contamination Probes

```
Usage: python -m src.tier1.contamination --model gpt-5.2-pro
```

Three probes, implemented as simple prompt templates:

**Cloze:** `"Complete this scientific benchmark question: '{first_half_of_question}...'"`
— Score: does the completion match the actual second half? (fuzzy string match, >0.8 similarity)

**Reverse:** `"What question would produce this answer in a biology benchmark: '{ideal_answer}'"`
— Score: does the generated question match the real question? (semantic similarity via embedding cosine)

**Temporal:** Split tasks by source paper date vs. model training cutoff. Compare accuracy on pre-cutoff vs. post-cutoff items using a chi-squared test.

---

## Tier 2 — Task Generation

Each generator is a standalone script that outputs JSON files to `tasks/`.

### Common task format

```json
{
    "id": "pro_stats_001",
    "category": "StatisticalReasoning",
    "question": "Given the following gene expression matrix...",
    "ideal": "Wilcoxon rank-sum test; top genes: BRCA1 (p=0.003), ...",
    "verification": "programmatic",
    "verification_fn": "verify_stats_001",
    "source": "pro",
    "meta": {"geo_accession": "GSE12345", "difficulty": "medium"}
}
```

Each script follows the same pattern:
1. Fetch data from public source (GEO, PDB, PMC)
2. Generate question from template
3. Compute ground-truth answer programmatically
4. Write to `tasks/{category}/` as JSON
5. Run `validate_tasks.py` to verify all ground truths are correct

### 2.1 `gen_stats_tasks.py`

Data source: GEO (via `GEOparse` library) or pre-downloaded CSVs.

Template types:
- "Which statistical test is appropriate for this design?" (categorical answer)
- "Compute the p-value for differential expression of gene X" (numerical, tolerance ±0.01)
- "How many genes are significant at FDR < 0.05?" (exact integer)

Ground truth: pre-computed via `scipy.stats` (t-test, Mann-Whitney, Fisher's exact).

Target: 200 tasks. Budget: ~2 weeks.

### 2.2 `gen_structure.py`

**PDB tasks:** Download 50 structures from RCSB. Use BioPython to extract facts (number of chains, active site residues, secondary structure fractions). Generate questions. Answer is deterministic from the PDB file.

**Gel images:** Use matplotlib to programmatically draw gel images with known band patterns. Questions: "What is the approximate size of band 3?" "Is there evidence of degradation?" Answers are deterministic from the generation parameters.

Target: 150 tasks. Budget: ~2 weeks.

### 2.3 `gen_calibration.py`

Take 100 existing LABBench2 questions. For each, create a modified version where one critical piece of information is removed from the context. The correct answer becomes "insufficient information to determine." Score: binary — did the model refuse to answer or express appropriate uncertainty?

Target: 100 tasks. Budget: ~1 week.

### 2.4 `gen_hypothesis.py`

Data source: PubMed Central Open Access (via `biopython.Entrez`).

Steps:
1. Fetch recent papers with "unexpected finding" or "future direction" in text
2. Extract the key finding paragraph
3. Format question: "Based on this finding, propose 3 testable hypotheses"
4. Ground truth: scored by LLM judge against a rubric (falsifiable? connected? novel?)

Rubric must be validated on a 50-task pilot with manual expert review before scaling.

Target: 100 tasks. Budget: ~2 weeks.

---

## Tier 3 — Compositional Chains

### 3.1 Chain Definitions

Define chains in a simple JSON file:

```json
// tasks/chains/chain_definitions.json
[
    {
        "chain_id": "lit_to_primer_01",
        "description": "Retrieve paper → read figure → design primers → check protocol",
        "steps": [
            {"task_id": "litqa3_042", "input_type": "none"},
            {"task_id": "figqa2_017", "input_type": "previous_response"},
            {"task_id": "seqqa2_108", "input_type": "previous_response"},
            {"task_id": "protqa2_055", "input_type": "previous_response"}
        ]
    }
]
```

Design 25–30 chains manually. This is a research design task, not an engineering task — Lukas will specify the chains, developer implements the runner.

### 3.2 `run_chains.py`

```
Usage: python -m src.tier3.run_chains --model claude-opus-4.5 --chain lit_to_primer_01
```

For each step:
1. If step 1: use the original task prompt
2. If step N>1: prepend the previous step's response to the task prompt
3. Call model, grade, save to `chain_runs`
4. If the model got step N wrong, still continue the chain (to measure error propagation)

### 3.3 Analysis

Query `chain_runs` to compute:
- Per-step accuracy within chains
- End-to-end chain accuracy (all steps correct)
- Theoretical independent baseline: product of per-step accuracies from Tier 1
- Error propagation: how often does a wrong step N cause wrong step N+1?

---

## Grading

Two grading paths. No third option.

**Programmatic** (SeqQA2, CloningQA, StatsReasoning, Structure tasks, Calibration):
- Exact match, numerical tolerance, or custom verifier function
- Deterministic, fast, no API cost

**LLM Judge** (LitQA3, FigQA2-retrieve, PatentQA, TrialQA, SourceQualQA, SuppQA2, ProtocolQA2, Hypothesis):
- Single prompt: provide the question, ideal answer, and model response
- Ask judge to output JSON: `{"correct": true/false, "explanation": "..."}`
- Use Claude Sonnet as default judge (fast, cheap, good enough)
- Cache judge results in Redis

```python
JUDGE_PROMPT = """You are grading a scientific benchmark response.

Question: {question}
Reference answer: {ideal}
Model response: {response}

Is the model response correct? It does not need to match word-for-word,
but must contain the key factual content of the reference answer.

Respond with JSON only: {{"correct": true/false, "explanation": "brief reason"}}"""
```

---

## API Endpoints

Minimal. The API exists to trigger runs and query results — not to be a product.

```
POST /runs/eval          — start an eval run (model + category)
POST /runs/chain         — start a chain run (model + chain_id)
GET  /results/eval       — query eval results (filter by model, category)
GET  /results/chains     — query chain results
GET  /results/ci         — get bootstrap CIs for a category
GET  /status             — job queue depth, running jobs
```

All return JSON. No auth needed (this is a research tool, not production).

---

## Order of Implementation

| Step | What | Depends on | Est. days |
|------|------|------------|-----------|
| 1 | `docker-compose.yml`, `schema.sql`, `setup_db.sh` | Nothing | 1 |
| 2 | `config.py`, `db.py`, `cache.py` | Step 1 | 1 |
| 3 | `models.py` — unified model caller | Step 2 | 2 |
| 4 | `run_eval.py` — run models on LABBench2 | Step 3 | 3 |
| 5 | Grading functions (programmatic + LLM judge) | Step 3 | 2 |
| 6 | `bootstrap_ci.py` | Step 4 | 1 |
| 7 | `irt_analysis.py` | Step 4 | 2 |
| 8 | `judge_audit.py` | Step 5 | 2 |
| 9 | `contamination.py` | Step 3 | 2 |
| 10 | `api.py` — FastAPI endpoints | Steps 4-9 | 2 |
| 11 | Tier 2 task generators (parallel work) | Step 2 | 8 |
| 12 | `validate_tasks.py` | Step 11 | 1 |
| 13 | Run Tier 2 tasks through eval pipeline | Steps 4, 11 | 2 |
| 14 | Chain definitions (Lukas provides) | None | — |
| 15 | `run_chains.py` + `feedback_sim.py` | Steps 4, 14 | 3 |
| 16 | `cost_tracker.py` + analysis notebooks | All above | 3 |

**Total estimate: ~5 weeks** for a single developer working full-time.

**Critical path:** Steps 1–4 must be sequential. After step 4, Tier 2 generation (step 11) can run in parallel with Tier 1 analysis (steps 6–9).

---

## Environment & Deployment

```yaml
# docker-compose.yml
services:
  postgres:
    image: postgres:16
    environment:
      POSTGRES_DB: labbench2pro
      POSTGRES_USER: dev
      POSTGRES_PASSWORD: dev
    ports: ["5432:5432"]
    volumes: ["pgdata:/var/lib/postgresql/data"]

  redis:
    image: redis:7-alpine
    ports: ["6379:6379"]

volumes:
  pgdata:
```

```bash
# .env.example
OPENAI_API_KEY=sk-...
ANTHROPIC_API_KEY=sk-ant-...
GOOGLE_API_KEY=...
HF_TOKEN=hf_...
DATABASE_URL=postgresql://dev:dev@localhost:5432/labbench2pro
REDIS_URL=redis://localhost:6379/0
```

No Kubernetes. No cloud deployment. Everything runs locally or on a single VM. If compute budget requires it, spin up a bigger VM — don't add infrastructure complexity.
