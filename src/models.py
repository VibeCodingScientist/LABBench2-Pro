"""Unified model caller. One function, all providers."""

import asyncio
import base64
import hashlib
import time

import httpx
from anthropic import AsyncAnthropic

from src.config import (
    ANTHROPIC_API_KEY,
    BIOXYZ_API_KEY,
    COST_TABLE,
    GOOGLE_API_KEY,
    MODEL_REGISTRY,
    OPENAI_API_KEY,
)
from src.cache import cache_get, cache_set

anthropic_client = AsyncAnthropic(api_key=ANTHROPIC_API_KEY)


def _build_anthropic_messages(prompt: str, images: list[bytes] | None = None) -> list[dict]:
    content = []
    if images:
        for img in images:
            content.append({
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": "image/png",
                    "data": base64.b64encode(img).decode(),
                },
            })
    content.append({"type": "text", "text": prompt})
    return [{"role": "user", "content": content}]


def _build_openai_messages(system: str, prompt: str, images: list[bytes] | None = None) -> list[dict]:
    msgs = []
    if system:
        msgs.append({"role": "system", "content": system})
    content = []
    if images:
        for img in images:
            content.append({
                "type": "image_url",
                "image_url": {"url": f"data:image/png;base64,{base64.b64encode(img).decode()}"},
            })
    content.append({"type": "text", "text": prompt})
    msgs.append({"role": "user", "content": content})
    return msgs


def _build_google_messages(system: str, prompt: str, images: list[bytes] | None = None) -> list:
    parts = []
    if system:
        parts.append(system + "\n\n")
    if images:
        for img in images:
            parts.append({"inline_data": {"mime_type": "image/png", "data": base64.b64encode(img).decode()}})
    parts.append(prompt)
    return parts


_BIOXYZ_BASE = "https://api.ai.bio.xyz"
_BIOXYZ_POLL_INTERVAL = 15   # seconds between status checks
_BIOXYZ_TIMEOUT = 45 * 60   # 45 minutes max


async def _call_bioxyz(prompt: str, system: str) -> tuple[str, int, int]:
    """Call bio.xyz BioAgent deep-research API (async polling)."""
    if not BIOXYZ_API_KEY:
        raise NotImplementedError("bio.xyz provider not configured — set BIOXYZ_API_KEY")

    message = f"{system}\n\n{prompt}" if system else prompt
    headers = {"Authorization": f"Bearer {BIOXYZ_API_KEY}"}

    async with httpx.AsyncClient(timeout=60) as client:
        # Start research job
        resp = await client.post(
            f"{_BIOXYZ_BASE}/deep-research/start",
            data={"message": message, "researchMode": "steering"},
            headers=headers,
        )
        resp.raise_for_status()
        conversation_id = resp.json()["conversationId"]

        # Poll until completed
        deadline = time.monotonic() + _BIOXYZ_TIMEOUT
        while True:
            await asyncio.sleep(_BIOXYZ_POLL_INTERVAL)
            if time.monotonic() > deadline:
                raise TimeoutError(
                    f"BioAgent job {conversation_id} timed out after {_BIOXYZ_TIMEOUT}s"
                )
            poll = await client.get(
                f"{_BIOXYZ_BASE}/deep-research/{conversation_id}",
                headers=headers,
            )
            poll.raise_for_status()
            data = poll.json()
            if data.get("status") == "completed":
                messages = data.get("messages", [])
                if messages:
                    return messages[-1].get("content", ""), 0, 0
                return "", 0, 0

    # unreachable, but satisfies type checkers
    return "", 0, 0


async def call_model(
    model_key: str,
    prompt: str,
    system: str = "",
    images: list[bytes] | None = None,
    max_tokens: int = 2048,
    cache: bool = True,
) -> dict:
    """Call a model and return response + usage metadata.

    Returns: {"response", "tokens_in", "tokens_out", "cost_usd", "latency_ms"}
    """
    if model_key not in MODEL_REGISTRY:
        raise ValueError(f"Unknown model: {model_key}. Available: {list(MODEL_REGISTRY)}")

    # Check cache
    if cache:
        cache_key = f"cache:{model_key}:{hashlib.sha256(prompt.encode()).hexdigest()}"
        cached = await cache_get(cache_key)
        if cached:
            return cached

    info = MODEL_REGISTRY[model_key]
    start = time.monotonic()

    if info["provider"] == "anthropic":
        resp = await anthropic_client.messages.create(
            model=info["model"],
            system=system if system else [],
            messages=_build_anthropic_messages(prompt, images),
            max_tokens=max_tokens,
        )
        text = resp.content[0].text if resp.content else ""
        t_in, t_out = resp.usage.input_tokens, resp.usage.output_tokens

    elif info["provider"] == "openai":
        if not OPENAI_API_KEY:
            raise NotImplementedError("OpenAI provider not configured — set OPENAI_API_KEY")
        from openai import AsyncOpenAI
        client = AsyncOpenAI(api_key=OPENAI_API_KEY)
        # GPT-5.x models require max_completion_tokens; older models use max_tokens
        use_new_param = info["model"].startswith("gpt-5") or info["model"].startswith("o")
        token_kwarg = {"max_completion_tokens": max_tokens} if use_new_param else {"max_tokens": max_tokens}
        resp = await client.chat.completions.create(
            model=info["model"],
            messages=_build_openai_messages(system, prompt, images),
            **token_kwarg,
        )
        text = resp.choices[0].message.content or ""
        t_in, t_out = resp.usage.prompt_tokens, resp.usage.completion_tokens

    elif info["provider"] == "google":
        if not GOOGLE_API_KEY:
            raise NotImplementedError("Google provider not configured — set GOOGLE_API_KEY")
        from google import genai
        client = genai.Client(api_key=GOOGLE_API_KEY)
        resp = await client.aio.models.generate_content(
            model=info["model"],
            contents=_build_google_messages(system, prompt, images),
        )
        text = resp.text or ""
        t_in = resp.usage_metadata.prompt_token_count or 0
        t_out = resp.usage_metadata.candidates_token_count or 0

    elif info["provider"] == "bioxyz":
        text, t_in, t_out = await _call_bioxyz(prompt, system)

    else:
        raise ValueError(f"Unknown provider: {info['provider']}")

    elapsed = int((time.monotonic() - start) * 1000)
    costs = COST_TABLE.get(model_key, {"input": 0, "output": 0})
    cost = (t_in * costs["input"] + t_out * costs["output"]) / 1_000_000

    result = {
        "response": text,
        "tokens_in": t_in,
        "tokens_out": t_out,
        "cost_usd": round(cost, 6),
        "latency_ms": elapsed,
    }

    if cache:
        await cache_set(cache_key, result)

    return result
