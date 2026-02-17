"""Configuration — env vars, model registry, cost table."""

import os
from dotenv import load_dotenv

load_dotenv()

# --- Required ---
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
if not ANTHROPIC_API_KEY:
    raise RuntimeError("ANTHROPIC_API_KEY is required. Set it in .env or environment.")

# --- Optional (stubbed providers) ---
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")
GOOGLE_API_KEY = os.environ.get("GOOGLE_API_KEY", "")
HF_TOKEN = os.environ.get("HF_TOKEN", "")
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

# --- Infrastructure ---
DATABASE_URL = os.environ.get("DATABASE_URL", "postgresql://dev:dev@localhost:5432/labbench2pro")
REDIS_URL = os.environ.get("REDIS_URL", "redis://localhost:6379/0")

# --- Model registry ---
MODEL_REGISTRY = {
    "claude-opus-4.6": {"provider": "anthropic", "model": "claude-opus-4-6"},
    "claude-opus-4.5": {"provider": "anthropic", "model": "claude-opus-4-5-20250929"},
    "claude-sonnet-4.5": {"provider": "anthropic", "model": "claude-sonnet-4-5-20250929"},
    "gpt-5.2-pro": {"provider": "openai", "model": "gpt-5.2-pro"},
    "gemini-3-pro": {"provider": "google", "model": "gemini-3.0-pro"},
}

# Cost per 1M tokens (input / output)
COST_TABLE = {
    "claude-opus-4.6": {"input": 15.00, "output": 75.00},
    "claude-opus-4.5": {"input": 15.00, "output": 75.00},
    "claude-sonnet-4.5": {"input": 3.00, "output": 15.00},
    "gpt-5.2-pro": {"input": 2.50, "output": 10.00},
    "gemini-3-pro": {"input": 1.25, "output": 5.00},
}
