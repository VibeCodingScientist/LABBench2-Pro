"""Unified model caller. One function, all providers."""

import base64
import hashlib
import time

from anthropic import AsyncAnthropic

from src.config import (
    ANTHROPIC_API_KEY,
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
        text = resp.text
        t_in = resp.usage_metadata.prompt_token_count
        t_out = resp.usage_metadata.candidates_token_count

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
