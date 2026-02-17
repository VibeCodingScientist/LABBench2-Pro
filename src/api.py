"""FastAPI app — trigger runs and query results.

Usage: uvicorn src.api:app --reload
"""

import json

from fastapi import BackgroundTasks, FastAPI, Query

from src.cache import get_redis
from src.db import fetch, fetchrow

app = FastAPI(title="LABBench2-Pro", version="0.1.0")


@app.post("/runs/eval")
async def start_eval_run(
    model: str,
    category: str,
    concurrency: int = 5,
    background_tasks: BackgroundTasks = None,
):
    """Queue an eval run via Redis."""
    r = await get_redis()
    if r is None:
        return {"error": "Redis unavailable — run eval directly via CLI"}

    job = json.dumps({"type": "eval", "model": model, "category": category, "concurrency": concurrency})
    await r.lpush("queue:eval_jobs", job)
    return {"status": "queued", "model": model, "category": category}


@app.post("/runs/chain")
async def start_chain_run(model: str, chain_id: str):
    """Queue a chain run via Redis."""
    r = await get_redis()
    if r is None:
        return {"error": "Redis unavailable — run chain directly via CLI"}

    job = json.dumps({"type": "chain", "model": model, "chain_id": chain_id})
    await r.lpush("queue:eval_jobs", job)
    return {"status": "queued", "model": model, "chain_id": chain_id}


@app.get("/results/eval")
async def get_eval_results(
    model: str | None = None,
    category: str | None = None,
    limit: int = Query(default=100, le=1000),
):
    """Query eval results, optionally filtered by model and category."""
    query = """
        SELECT e.id, m.model_name, t.id as task_id, t.category,
               e.correct, e.score, e.grader, e.cost_usd, e.latency_ms, e.run_at
        FROM eval_runs e
        JOIN models m ON e.model_id = m.id
        JOIN tasks t ON e.task_id = t.id
        WHERE 1=1
    """
    args = []
    i = 1

    if model:
        query += f" AND m.model_name = ${i}"
        args.append(model)
        i += 1
    if category:
        query += f" AND t.category = ${i}"
        args.append(category)
        i += 1

    query += f" ORDER BY e.run_at DESC LIMIT ${i}"
    args.append(limit)

    rows = await fetch(query, *args)
    return [dict(row) for row in rows]


@app.get("/results/chains")
async def get_chain_results(chain_id: str | None = None, model: str | None = None):
    """Query chain run results."""
    query = """
        SELECT c.chain_id, m.model_name, c.step_num, c.task_id,
               c.correct, c.run_at
        FROM chain_runs c
        JOIN models m ON c.model_id = m.id
        WHERE 1=1
    """
    args = []
    i = 1

    if chain_id:
        query += f" AND c.chain_id = ${i}"
        args.append(chain_id)
        i += 1
    if model:
        query += f" AND m.model_name = ${i}"
        args.append(model)
        i += 1

    query += " ORDER BY c.chain_id, c.step_num"

    rows = await fetch(query, *args)
    return [dict(row) for row in rows]


@app.get("/results/ci")
async def get_bootstrap_ci(category: str):
    """Get accuracy + bootstrap CIs (pre-computed, or compute on the fly)."""
    rows = await fetch(
        """SELECT m.model_name, e.correct
           FROM eval_runs e
           JOIN models m ON e.model_id = m.id
           JOIN tasks t ON e.task_id = t.id
           WHERE t.category = $1""",
        category,
    )

    if not rows:
        return {"error": f"No results for category '{category}'"}

    import numpy as np
    from scipy.stats import bootstrap

    by_model: dict[str, list[bool]] = {}
    for row in rows:
        by_model.setdefault(row["model_name"], []).append(row["correct"])

    results = []
    for model, corrects in sorted(by_model.items()):
        arr = np.array(corrects, dtype=float)
        if len(arr) < 2:
            results.append({"model": model, "accuracy": float(arr.mean()), "ci_low": 0.0, "ci_high": 1.0, "n": len(arr)})
            continue
        ci = bootstrap((arr,), np.mean, n_resamples=10000, confidence_level=0.95, method="BCa")
        results.append({
            "model": model,
            "accuracy": float(arr.mean()),
            "ci_low": float(ci.confidence_interval.low),
            "ci_high": float(ci.confidence_interval.high),
            "n": len(arr),
        })

    return results


@app.get("/status")
async def get_status():
    """Job queue depth and system status."""
    r = await get_redis()
    queue_depth = 0
    redis_status = "disconnected"

    if r:
        try:
            queue_depth = await r.llen("queue:eval_jobs")
            redis_status = "connected"
        except Exception:
            pass

    # DB check
    db_status = "disconnected"
    try:
        row = await fetchrow("SELECT COUNT(*) as n FROM eval_runs")
        eval_count = row["n"] if row else 0
        db_status = "connected"
    except Exception:
        eval_count = 0

    return {
        "redis": redis_status,
        "database": db_status,
        "queue_depth": queue_depth,
        "total_eval_runs": eval_count,
    }
