"""Prompt registry for report summary generation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from importlib.resources.abc import Traversable
from typing import Dict, Iterable, Optional


DEFAULT_PROMPT_IDS = {
    "run": "run_overview",
    "collection": "collection_summary",
    "comparison": "comparison_summary",
    "leading_edge": "leading_edge_summary",
}


@dataclass(frozen=True)
class PromptVariant:
    name: str
    system_prompt: str
    instruction: str


@dataclass(frozen=True)
class PromptSpec:
    prompt_id: str
    version: str
    scope: str
    variants: Dict[str, PromptVariant]


@dataclass(frozen=True)
class PromptSelection:
    prompt_id: str
    version: str
    scope: str
    variant: str
    system_prompt: str
    instruction: str


def _prompt_dir():
    return resources.files(__package__) / "prompt_specs"


def _iter_prompt_files() -> Iterable[Traversable]:
    for entry in _prompt_dir().iterdir():
        if entry.is_file() and entry.name.endswith(".json"):
            yield entry


@lru_cache(maxsize=1)
def load_prompt_specs() -> Dict[str, PromptSpec]:
    registry: Dict[str, PromptSpec] = {}
    for path in _iter_prompt_files():
        raw = json.loads(path.read_text(encoding="utf-8"))
        variants = {
            str(name): PromptVariant(
                name=str(name),
                system_prompt=str(payload["system_prompt"]).strip(),
                instruction=str(payload["instruction"]).strip(),
            )
            for name, payload in (raw.get("variants") or {}).items()
        }
        spec = PromptSpec(
            prompt_id=str(raw["id"]).strip(),
            version=str(raw["version"]).strip(),
            scope=str(raw["scope"]).strip(),
            variants=variants,
        )
        registry[spec.prompt_id] = spec
    return registry


def list_prompt_specs() -> list[PromptSpec]:
    return sorted(load_prompt_specs().values(), key=lambda spec: (spec.scope, spec.prompt_id))


def resolve_prompt(
    *,
    scope: str,
    prompt_id: Optional[str] = None,
    variant: str = "default",
) -> PromptSelection:
    registry = load_prompt_specs()
    target_id = (prompt_id or DEFAULT_PROMPT_IDS.get(scope) or "").strip()
    if not target_id:
        raise KeyError(f"No default prompt id configured for scope {scope!r}")
    if target_id not in registry:
        raise KeyError(f"Unknown prompt id: {target_id}")
    spec = registry[target_id]
    if spec.scope != scope:
        raise KeyError(f"Prompt {target_id!r} has scope {spec.scope!r}, expected {scope!r}")
    target_variant = variant.strip() or "default"
    if target_variant not in spec.variants:
        available = ", ".join(sorted(spec.variants))
        raise KeyError(f"Prompt {target_id!r} does not define variant {target_variant!r}. Available: {available}")
    selected = spec.variants[target_variant]
    return PromptSelection(
        prompt_id=spec.prompt_id,
        version=spec.version,
        scope=spec.scope,
        variant=selected.name,
        system_prompt=selected.system_prompt,
        instruction=selected.instruction,
    )


__all__ = [
    "DEFAULT_PROMPT_IDS",
    "PromptSelection",
    "PromptSpec",
    "PromptVariant",
    "list_prompt_specs",
    "load_prompt_specs",
    "resolve_prompt",
]
