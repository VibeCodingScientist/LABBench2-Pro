"""Simple Redis worker — BRPOP loop for eval/chain jobs.

Usage: python -m src.worker
"""

import asyncio
import json
import sys

from src.cache import get_redis


async def process_eval_job(job: dict):
    """Run an eval batch."""
    from src.tier1.run_eval import load_tasks, ensure_model, upsert_task, run_single

    model_key = job["model"]
    category = job["category"]
    concurrency = job.get("concurrency", 5)

    print(f"[worker] Starting eval: model={model_key}, category={category}")
    tasks = await load_tasks(category)
    model_id = await ensure_model(model_key)

    for task in tasks:
        await upsert_task(task)

    sem = asyncio.Semaphore(concurrency)
    results = await asyncio.gather(*[run_single(model_key, model_id, t, sem) for t in tasks])

    correct = sum(1 for r in results if r.get("correct"))
    total = len(results)
    print(f"[worker] Eval complete: {correct}/{total}")


async def process_chain_job(job: dict):
    """Run a chain."""
    from src.tier3.run_chains import run_chain

    model_key = job["model"]
    chain_id = job["chain_id"]

    print(f"[worker] Starting chain: model={model_key}, chain={chain_id}")
    await run_chain(model_key, chain_id)
    print(f"[worker] Chain complete.")


async def worker_loop():
    """Main worker loop — blocks on Redis queue."""
    r = await get_redis()
    if r is None:
        print("Redis unavailable. Cannot start worker.", file=sys.stderr)
        sys.exit(1)

    print("[worker] Listening for jobs on queue:eval_jobs...")

    while True:
        try:
            result = await r.brpop("queue:eval_jobs", timeout=30)
            if result is None:
                continue

            _, raw = result
            job = json.loads(raw)
            job_type = job.get("type", "eval")

            if job_type == "eval":
                await process_eval_job(job)
            elif job_type == "chain":
                await process_chain_job(job)
            else:
                print(f"[worker] Unknown job type: {job_type}")

        except KeyboardInterrupt:
            print("\n[worker] Shutting down.")
            break
        except Exception as e:
            print(f"[worker] Error processing job: {e}", file=sys.stderr)
            await asyncio.sleep(1)


if __name__ == "__main__":
    asyncio.run(worker_loop())
