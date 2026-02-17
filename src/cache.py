"""Redis cache — graceful degradation if unavailable."""

import json
import redis.asyncio as aioredis
from src.config import REDIS_URL

_redis: aioredis.Redis | None = None


async def get_redis() -> aioredis.Redis | None:
    global _redis
    if _redis is None:
        try:
            _redis = aioredis.from_url(REDIS_URL, decode_responses=True)
            await _redis.ping()
        except Exception:
            _redis = None
    return _redis


async def cache_get(key: str) -> dict | None:
    r = await get_redis()
    if r is None:
        return None
    try:
        val = await r.get(key)
        return json.loads(val) if val else None
    except Exception:
        return None


async def cache_set(key: str, value: dict, ttl: int = 30 * 86400) -> None:
    r = await get_redis()
    if r is None:
        return
    try:
        await r.set(key, json.dumps(value), ex=ttl)
    except Exception:
        pass
