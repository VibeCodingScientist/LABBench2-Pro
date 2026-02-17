"""Load and query chain definitions."""

import json
import os

CHAIN_FILE = os.path.join(os.path.dirname(__file__), "../../tasks/chains/chain_definitions.json")


def load_chains() -> list[dict]:
    with open(CHAIN_FILE) as f:
        return json.load(f)


def get_chain(chain_id: str) -> dict | None:
    for chain in load_chains():
        if chain["chain_id"] == chain_id:
            return chain
    return None


def list_chain_ids() -> list[str]:
    return [c["chain_id"] for c in load_chains()]
