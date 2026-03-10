"""LLM-backed summarisation helpers for tackle2 reports."""

from __future__ import annotations

import configparser
import hashlib
import json
import logging
import math
import os
import socket
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Protocol, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import ollama
from ollama import Client, RequestError, ResponseError

from .prompts import resolve_prompt
from .summary import CollectionSummary, ComparisonSummary, GeneratedSummary, TopPathway

logger = logging.getLogger(__name__)

_STATE_DIR_ENV = "ISPEC_STATE_DIR"
_AGENT_CONF_ENV = "TACKLE_AGENT_CONF"
_DEFAULT_AGENT_CONF_FILENAME = "tackle-agent.conf"
_FALLBACK_AGENT_CONF_FILENAME = "ispec.conf"
_ALLOWED_AGENT_CONF_KEYS = {
    "agent_api",
    "agent_id",
    "api_key",
    "ispec_api_key",
    "bearer_token",
    "token",
    "tackle_agent_api",
    "tackle_agent_id",
    "tackle_agent_api_key",
    "tackle_agent_bearer_token",
    "tackle_agent_timeout_seconds",
    "tackle_agent_model_label",
    "timeout_seconds",
    "model_label",
}

_PROMPT_FLOAT_SIG_DIGITS = 5


class SummarisationError(RuntimeError):
    """Raised when the configured LLM fails in an unexpected way."""


class ReportSummarizer(Protocol):
    """Interface implemented by helpers capable of generating report prose."""

    def set_cache_dir(self, cache_dir: Path | None) -> None:
        ...

    def build_requests(
        self,
        collections: Sequence[CollectionSummary],
        *,
        include_run: bool = True,
        include_collections: bool = True,
        include_comparisons: bool = False,
        include_leading_edges: bool = False,
        user_notes: Sequence[str] = (),
    ) -> list[SummaryRequest]:
        ...

    def summarise_request(self, request: SummaryRequest) -> Optional[GeneratedSummary]:
        ...

    def summarise_run(self, collections: Sequence[CollectionSummary]) -> Optional[GeneratedSummary]:
        ...

    def summarise_collections(self, collections: Sequence[CollectionSummary]) -> Dict[str, GeneratedSummary]:
        ...

    def summarise_comparisons(self, collection: CollectionSummary) -> Dict[str, GeneratedSummary]:
        ...


@dataclass(frozen=True)
class SummaryRequest:
    """A single summary generation request with fully rendered prompt context."""

    scope: str
    key: str
    title: str
    target: Mapping[str, str]
    prompt_id: str
    prompt_version: str
    prompt_variant: str
    system_prompt: str
    instruction: str
    payload: Mapping[str, object]
    prompt_text: str
    user_notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class OllamaConfig:
    """Configuration for :class:`OllamaSummarizer`."""

    model: str = "llama3.2:3b"
    host: Optional[str] = None
    temperature: float = 0.2
    max_tokens: int = 512
    run_collection_limit: int = 6
    comparisons_per_collection: int = 3
    pathways_per_comparison: int = 3
    enable_collection_summaries: bool = True
    collection_summary_limit: int = 3
    enable_comparison_summaries: bool = False
    comparison_summary_limit: int = 0
    enable_leading_edge_summaries: bool = False
    leading_edge_summary_limit: int = 3
    leading_edge_gene_limit: int = 12
    prompt_variant: str = "default"
    run_prompt_id: str = "run_overview"
    collection_prompt_id: str = "collection_summary"
    comparison_prompt_id: str = "comparison_summary"
    leading_edge_prompt_id: str = "leading_edge_summary"


def _strip_wrapping_quotes(value: str) -> str:
    value = str(value).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1].strip()
    return value


def _agent_conf_paths() -> list[Path]:
    override = (os.environ.get(_AGENT_CONF_ENV) or "").strip()
    if override:
        return [Path(override).expanduser()]

    base_dir = Path(os.environ.get(_STATE_DIR_ENV, Path.home() / ".ispec")).expanduser()
    return [
        base_dir / _DEFAULT_AGENT_CONF_FILENAME,
        base_dir / _FALLBACK_AGENT_CONF_FILENAME,
    ]


def _load_ini_config(path: Path) -> dict[str, str]:
    try:
        raw = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return {}
    except OSError:
        return {}

    stripped_lines = [
        line
        for line in raw.splitlines()
        if line.strip() and not line.lstrip().startswith(("#", ";"))
    ]
    if not stripped_lines:
        return {}

    normalized = raw
    if not any(line.startswith("[") and "]" in line for line in stripped_lines[:5]):
        normalized = "[DEFAULT]\n" + raw

    parser = configparser.ConfigParser(interpolation=None)
    try:
        parser.read_string(normalized)
    except configparser.Error:
        return {}

    sections_by_lower = {section.lower(): section for section in parser.sections()}
    preferred_sections = [
        sections_by_lower.get("agent"),
        sections_by_lower.get("tackle"),
        sections_by_lower.get("ispec"),
    ]
    if "DEFAULT" in parser:
        preferred_sections.append("DEFAULT")

    result: dict[str, str] = {}
    for section in [item for item in preferred_sections if item]:
        for key, value in parser[section].items():
            if not isinstance(key, str) or not isinstance(value, str):
                continue
            if key not in _ALLOWED_AGENT_CONF_KEYS:
                continue
            cleaned = _strip_wrapping_quotes(value)
            if key not in result and cleaned:
                result[key] = cleaned
    return result


def _load_agent_config() -> dict[str, str]:
    merged: dict[str, str] = {}
    for path in _agent_conf_paths():
        merged.update(_load_ini_config(path))
    return merged


def _default_agent_id() -> str:
    user = os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"
    try:
        host = socket.gethostname() or "unknown"
    except Exception:
        host = "unknown"
    return f"{user}@{host}"


def _build_support_chat_url(agent_api: str) -> str:
    agent_api = str(agent_api).strip()
    if not agent_api:
        raise ValueError("agent_api is empty")
    if agent_api.endswith("/api/support/chat"):
        return agent_api
    return agent_api.rstrip("/") + "/api/support/chat"


def _summary_cache_path(cache_dir: Path, *, message: str) -> Path:
    digest = hashlib.sha256(message.encode("utf-8")).hexdigest()
    return cache_dir / f"{digest}.json"


def _coerce_summary_payload(text: str) -> Optional[Mapping[str, object]]:
    raw = str(text or "").strip()
    if not raw:
        return None

    if raw.startswith("```"):
        lines = raw.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip() == "```":
            lines = lines[:-1]
        raw = "\n".join(lines).strip()

    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return {"summary": raw, "bullets": []}

    if isinstance(payload, Mapping):
        return payload
    if isinstance(payload, list):
        bullets = [str(item).strip() for item in payload if str(item).strip()]
        if bullets:
            return {"summary": bullets[0], "bullets": bullets}
    return {"summary": raw, "bullets": []}


def _compact_prompt_value(value: object, *, digits: int = _PROMPT_FLOAT_SIG_DIGITS) -> object:
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return value
        return float(f"{value:.{digits}g}")
    if isinstance(value, Mapping):
        return {
            str(key): _compact_prompt_value(item, digits=digits)
            for key, item in value.items()
        }
    if isinstance(value, tuple):
        return tuple(_compact_prompt_value(item, digits=digits) for item in value)
    if isinstance(value, list):
        return [_compact_prompt_value(item, digits=digits) for item in value]
    return value


def _compact_prompt_payload(
    payload: Mapping[str, object],
    *,
    digits: int = _PROMPT_FLOAT_SIG_DIGITS,
) -> Dict[str, object]:
    compacted = _compact_prompt_value(payload, digits=digits)
    if isinstance(compacted, dict):
        return compacted
    return dict(compacted)


def _load_summary_cache(cache_dir: Path, *, message: str) -> Optional[Mapping[str, object]]:
    cache_path = _summary_cache_path(cache_dir, message=message)
    if not cache_path.exists() or not cache_path.is_file():
        return None
    try:
        obj = json.loads(cache_path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("AI summary cache read failed for %s: %s", cache_path, exc)
        return None

    if str(obj.get("prompt_message") or "") != message:
        return None
    payload = obj.get("summary_payload")
    if isinstance(payload, Mapping):
        return payload
    text = str(obj.get("response_text") or "").strip()
    return _coerce_summary_payload(text)


def _write_summary_cache(
    cache_dir: Path,
    *,
    message: str,
    response: Mapping[str, Any],
    summary_payload: Mapping[str, object],
    model_label: str,
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = _summary_cache_path(cache_dir, message=message)
    payload = {
        "schema_version": 1,
        "message_sha256": cache_path.stem,
        "prompt_message": message,
        "response": dict(response),
        "response_text": str(response.get("message") or ""),
        "summary_payload": dict(summary_payload),
        "model_label": model_label,
    }
    tmp_path = cache_dir / f".{cache_path.stem}.{uuid.uuid4().hex}.tmp"
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(cache_path)
    return cache_path


@dataclass(frozen=True)
class AgentApiConfig:
    """Configuration for :class:`AgentApiSummarizer`."""

    agent_api: str
    agent_id: str
    api_key: Optional[str] = None
    bearer_token: Optional[str] = None
    timeout_seconds: float = 120.0
    model: str = "llama3.1 tulu"
    run_collection_limit: int = 6
    comparisons_per_collection: int = 3
    pathways_per_comparison: int = 3
    enable_collection_summaries: bool = True
    collection_summary_limit: int = 3
    enable_comparison_summaries: bool = False
    comparison_summary_limit: int = 0
    force_refresh: bool = False
    enable_leading_edge_summaries: bool = False
    leading_edge_summary_limit: int = 3
    leading_edge_gene_limit: int = 12
    prompt_variant: str = "default"
    run_prompt_id: str = "run_overview"
    collection_prompt_id: str = "collection_summary"
    comparison_prompt_id: str = "comparison_summary"
    leading_edge_prompt_id: str = "leading_edge_summary"

    @classmethod
    def from_env(
        cls,
        *,
        agent_api: Optional[str] = None,
        agent_id: Optional[str] = None,
        model: Optional[str] = None,
        run_collection_limit: int = 6,
        comparisons_per_collection: int = 3,
        pathways_per_comparison: int = 3,
        enable_collection_summaries: bool = True,
        collection_summary_limit: int = 3,
        enable_comparison_summaries: bool = False,
        comparison_summary_limit: int = 0,
        force_refresh: bool = False,
        enable_leading_edge_summaries: bool = False,
        leading_edge_summary_limit: int = 3,
        leading_edge_gene_limit: int = 12,
        prompt_variant: str = "default",
        run_prompt_id: str = "run_overview",
        collection_prompt_id: str = "collection_summary",
        comparison_prompt_id: str = "comparison_summary",
        leading_edge_prompt_id: str = "leading_edge_summary",
    ) -> "AgentApiConfig":
        conf = _load_agent_config()

        resolved_agent_api = (
            agent_api
            or os.environ.get("TACKLE_AGENT_API")
            or conf.get("tackle_agent_api")
            or conf.get("agent_api")
            or ""
        ).strip()
        if not resolved_agent_api:
            raise SummarisationError(
                "AI summary requested but TACKLE_AGENT_API is not configured "
                "(set it in ~/.ispec/tackle-agent.conf or as an env var)."
            )

        timeout_raw = (
            os.environ.get("TACKLE_AGENT_TIMEOUT_SECONDS")
            or conf.get("tackle_agent_timeout_seconds")
            or conf.get("timeout_seconds")
            or "120.0"
        )
        try:
            timeout_seconds = float(timeout_raw)
        except ValueError:
            timeout_seconds = 120.0

        resolved_model = (
            model
            or os.environ.get("TACKLE_AGENT_MODEL_LABEL")
            or conf.get("tackle_agent_model_label")
            or conf.get("model_label")
            or "llama3.1 tulu"
        ).strip()

        return cls(
            agent_api=resolved_agent_api,
            agent_id=(
                agent_id
                or os.environ.get("TACKLE_AGENT_ID")
                or conf.get("tackle_agent_id")
                or conf.get("agent_id")
                or _default_agent_id()
            ),
            api_key=(
                os.environ.get("TACKLE_AGENT_API_KEY")
                or conf.get("tackle_agent_api_key")
                or conf.get("api_key")
                or conf.get("ispec_api_key")
                or None
            ),
            bearer_token=(
                os.environ.get("TACKLE_AGENT_BEARER_TOKEN")
                or conf.get("tackle_agent_bearer_token")
                or conf.get("bearer_token")
                or conf.get("token")
                or None
            ),
            timeout_seconds=timeout_seconds,
            model=resolved_model,
            run_collection_limit=run_collection_limit,
            comparisons_per_collection=comparisons_per_collection,
            pathways_per_comparison=pathways_per_comparison,
            enable_collection_summaries=enable_collection_summaries,
            collection_summary_limit=collection_summary_limit,
            enable_comparison_summaries=enable_comparison_summaries,
            comparison_summary_limit=comparison_summary_limit,
            force_refresh=force_refresh,
            enable_leading_edge_summaries=enable_leading_edge_summaries,
            leading_edge_summary_limit=leading_edge_summary_limit,
            leading_edge_gene_limit=leading_edge_gene_limit,
            prompt_variant=str(prompt_variant or "default"),
            run_prompt_id=str(run_prompt_id or "run_overview"),
            collection_prompt_id=str(collection_prompt_id or "collection_summary"),
            comparison_prompt_id=str(comparison_prompt_id or "comparison_summary"),
            leading_edge_prompt_id=str(leading_edge_prompt_id or "leading_edge_summary"),
        )


class OllamaSummarizer:
    """Generate narrative summaries using a locally hosted Ollama model."""

    _run_system_prompt = (
        "You are assisting with summarising gene set enrichment analysis (GSEA) results. "
        "Use only the provided data to produce concise, factual summaries. "
        "Highlight comparisons with the most significant pathways (padj < 0.05)."
    )

    _collection_system_prompt = (
        "You are writing short summaries for a single gene set enrichment analysis collection. "
        "Summaries must stay grounded in the provided JSON data and emphasise key pathways."
    )

    _comparison_system_prompt = (
        "You are summarising findings for a single gene set enrichment comparison. "
        "Base your statements strictly on the supplied JSON data."
    )

    def __init__(self, config: OllamaConfig | None = None) -> None:
        self._config = config or OllamaConfig()
        self._client: Any
        if self._config.host:
            self._client = Client(host=self._config.host)
        else:
            self._client = ollama

    # ------------------------------------------------------------------
    @property
    def config(self) -> OllamaConfig:
        return self._config

    def set_cache_dir(self, cache_dir: Path | None) -> None:
        _ = cache_dir

    # ------------------------------------------------------------------
    def summarise_run(self, collections: Sequence[CollectionSummary]) -> Optional[GeneratedSummary]:
        request = self.build_run_request(collections)
        if request is None:
            return None
        return self.summarise_request(request)

    # ------------------------------------------------------------------
    def summarise_collections(self, collections: Sequence[CollectionSummary]) -> Dict[str, GeneratedSummary]:
        if not self._config.enable_collection_summaries:
            return {}

        summaries: Dict[str, GeneratedSummary] = {}
        for request in self.build_collection_requests(collections):
            summary = self.summarise_request(request)
            if summary:
                target_key = request.target.get("collection_id") or request.key
                summaries[target_key] = summary
        return summaries

    # ------------------------------------------------------------------
    def summarise_comparisons(self, collection: CollectionSummary) -> Dict[str, GeneratedSummary]:
        if not self._config.enable_comparison_summaries:
            return {}

        summaries: Dict[str, GeneratedSummary] = {}
        for request in self.build_comparison_requests(collection):
            summary = self.summarise_request(request)
            if summary:
                target_key = request.target.get("comparison_id") or request.key
                summaries[target_key] = summary
        return summaries

    # ------------------------------------------------------------------
    def summarise_request(self, request: SummaryRequest) -> Optional[GeneratedSummary]:
        data = self._invoke_llm(request.prompt_text, system=request.system_prompt)
        return self._parse_summary(data)

    # ------------------------------------------------------------------
    def build_requests(
        self,
        collections: Sequence[CollectionSummary],
        *,
        include_run: bool = True,
        include_collections: bool = True,
        include_comparisons: bool = False,
        include_leading_edges: bool = False,
        user_notes: Sequence[str] = (),
    ) -> list[SummaryRequest]:
        requests: list[SummaryRequest] = []
        if include_run:
            request = self.build_run_request(collections, user_notes=user_notes)
            if request is not None:
                requests.append(request)
        if include_collections:
            requests.extend(self.build_collection_requests(collections, user_notes=user_notes))
        if include_comparisons:
            for collection in collections:
                requests.extend(self.build_comparison_requests(collection, user_notes=user_notes))
        if include_leading_edges:
            for collection in collections:
                requests.extend(self.build_leading_edge_requests(collection, user_notes=user_notes))
        return requests

    def build_run_request(
        self,
        collections: Sequence[CollectionSummary],
        *,
        user_notes: Sequence[str] = (),
    ) -> Optional[SummaryRequest]:
        selected = self._select_collections(collections, limit=self._config.run_collection_limit)
        if not selected:
            return None
        prompt = self._prompt_selection(scope="run")
        payload = _compact_prompt_payload(
            {
                "collections": [self._collection_payload(collection) for collection in selected],
                "total_collections": len(collections),
            }
        )
        prompt_text = self._render_prompt(prompt.instruction, payload, user_notes=user_notes)
        return SummaryRequest(
            scope="run",
            key="run",
            title="Overall report summary",
            target={"scope": "run"},
            prompt_id=prompt.prompt_id,
            prompt_version=prompt.version,
            prompt_variant=prompt.variant,
            system_prompt=prompt.system_prompt,
            instruction=prompt.instruction,
            payload=payload,
            prompt_text=prompt_text,
            user_notes=tuple(str(note).strip() for note in user_notes if str(note).strip()),
        )

    def build_collection_requests(
        self,
        collections: Sequence[CollectionSummary],
        *,
        user_notes: Sequence[str] = (),
    ) -> list[SummaryRequest]:
        selected = self._select_collections(collections, limit=self._config.collection_summary_limit)
        prompt = self._prompt_selection(scope="collection")
        requests: list[SummaryRequest] = []
        notes = tuple(str(note).strip() for note in user_notes if str(note).strip())
        for collection in selected:
            payload = _compact_prompt_payload(self._collection_payload(collection))
            prompt_text = self._render_prompt(prompt.instruction, payload, user_notes=notes)
            requests.append(
                SummaryRequest(
                    scope="collection",
                    key=collection.identifier,
                    title=collection.display_name,
                    target={
                        "scope": "collection",
                        "collection_id": collection.identifier,
                    },
                    prompt_id=prompt.prompt_id,
                    prompt_version=prompt.version,
                    prompt_variant=prompt.variant,
                    system_prompt=prompt.system_prompt,
                    instruction=prompt.instruction,
                    payload=payload,
                    prompt_text=prompt_text,
                    user_notes=notes,
                )
            )
        return requests

    def build_comparison_requests(
        self,
        collection: CollectionSummary,
        *,
        user_notes: Sequence[str] = (),
    ) -> list[SummaryRequest]:
        limit = self._config.comparison_summary_limit
        comparisons = collection.comparisons[:limit] if limit > 0 else collection.comparisons
        prompt = self._prompt_selection(scope="comparison")
        requests: list[SummaryRequest] = []
        notes = tuple(str(note).strip() for note in user_notes if str(note).strip())
        for comparison in comparisons:
            payload = _compact_prompt_payload(self._comparison_payload(collection, comparison))
            prompt_text = self._render_prompt(prompt.instruction, payload, user_notes=notes)
            requests.append(
                SummaryRequest(
                    scope="comparison",
                    key=f"{collection.identifier}:{comparison.identifier}",
                    title=f"{collection.display_name} / {comparison.display_name}",
                    target={
                        "scope": "comparison",
                        "collection_id": collection.identifier,
                        "comparison_id": comparison.identifier,
                    },
                    prompt_id=prompt.prompt_id,
                    prompt_version=prompt.version,
                    prompt_variant=prompt.variant,
                    system_prompt=prompt.system_prompt,
                    instruction=prompt.instruction,
                    payload=payload,
                    prompt_text=prompt_text,
                    user_notes=notes,
                )
            )
        return requests

    def build_leading_edge_requests(
        self,
        collection: CollectionSummary,
        *,
        user_notes: Sequence[str] = (),
    ) -> list[SummaryRequest]:
        if not self._config.enable_leading_edge_summaries:
            return []

        prompt = self._prompt_selection(scope="leading_edge")
        notes = tuple(str(note).strip() for note in user_notes if str(note).strip())
        requests: list[SummaryRequest] = []
        pathway_limit = self._config.leading_edge_summary_limit

        for comparison in collection.comparisons:
            pathways = comparison.top_pathways[:pathway_limit] if pathway_limit > 0 else comparison.top_pathways
            for pathway in pathways:
                payload = _compact_prompt_payload(self._leading_edge_payload(collection, comparison, pathway))
                if not payload.get("items"):
                    continue
                prompt_text = self._render_prompt(prompt.instruction, payload, user_notes=notes)
                requests.append(
                    SummaryRequest(
                        scope="leading_edge",
                        key=f"{collection.identifier}:{comparison.identifier}:{pathway.pathway}",
                        title=f"{collection.display_name} / {comparison.display_name} / {pathway.pathway}",
                        target={
                            "scope": "leading_edge",
                            "collection_id": collection.identifier,
                            "comparison_id": comparison.identifier,
                            "pathway_name": pathway.pathway,
                        },
                        prompt_id=prompt.prompt_id,
                        prompt_version=prompt.version,
                        prompt_variant=prompt.variant,
                        system_prompt=prompt.system_prompt,
                        instruction=prompt.instruction,
                        payload=payload,
                        prompt_text=prompt_text,
                        user_notes=notes,
                    )
                )
        return requests

    # ------------------------------------------------------------------
    def _select_collections(
        self,
        collections: Sequence[CollectionSummary],
        *,
        limit: int,
    ) -> Sequence[CollectionSummary]:
        if limit <= 0:
            return collections
        ordered = sorted(collections, key=self._collection_priority, reverse=True)
        return ordered[:limit]

    @staticmethod
    def _collection_priority(collection: CollectionSummary) -> tuple[int, int]:
        significant = sum(comp.significant_pathways for comp in collection.comparisons)
        total = sum(comp.total_pathways for comp in collection.comparisons)
        return significant, total

    def _render_prompt(
        self,
        instruction: str,
        payload: Mapping[str, object],
        *,
        user_notes: Sequence[str] = (),
    ) -> str:
        data = json.dumps(payload, indent=2, ensure_ascii=False)
        notes = [str(note).strip() for note in user_notes if str(note).strip()]
        if not notes:
            return f"{instruction}\n\nDATA:\n{data}"
        notes_text = "\n".join(f"- {note}" for note in notes)
        return f"{instruction}\n\nAdditional user notes:\n{notes_text}\n\nDATA:\n{data}"

    def _invoke_llm(self, prompt: str, *, system: str) -> Optional[Mapping[str, object]]:
        options: Dict[str, Any] = {}
        if self._config.temperature is not None:
            options["temperature"] = self._config.temperature
        if self._config.max_tokens > 0:
            options["num_predict"] = self._config.max_tokens
        try:
            response = self._client.generate(  # type: ignore[attr-defined]
                model=self._config.model,
                prompt=prompt,
                system=system,
                format="json",
                options=options or None,
            )
        except (RequestError, ResponseError, ConnectionError, TimeoutError, OSError) as exc:
            logger.warning("Failed to obtain summary from model %s: %s", self._config.model, exc)
            return None

        text = getattr(response, "response", None)
        if text is None and isinstance(response, Mapping):
            text = response.get("response")
        if not text:
            logger.info("Model %s returned an empty summary response.", self._config.model)
            return None

        try:
            return json.loads(text)
        except json.JSONDecodeError:
            logger.warning(
                "Model %s returned invalid JSON payload: %s",
                self._config.model,
                text,
            )
            return None

    def _parse_summary(self, payload: Optional[Mapping[str, object]]) -> Optional[GeneratedSummary]:
        if not payload:
            return None

        summary_text = str(payload.get("summary", "")).strip()
        bullets_raw = payload.get("bullets", [])

        bullets: list[str] = []
        if isinstance(bullets_raw, str):
            bullets = [line.strip(" •-\t") for line in bullets_raw.splitlines() if line.strip()]
        elif isinstance(bullets_raw, Iterable):
            for entry in bullets_raw:
                text = str(entry).strip()
                if text:
                    bullets.append(text.lstrip("•- ").strip())

        if not summary_text and bullets:
            summary_text = bullets[0]

        if not summary_text:
            return None

        return GeneratedSummary(text=summary_text, bullets=tuple(bullets), model=self._config.model)

    def _prompt_selection(self, *, scope: str):
        if scope == "run":
            prompt_id = self._config.run_prompt_id
        elif scope == "collection":
            prompt_id = self._config.collection_prompt_id
        elif scope == "comparison":
            prompt_id = self._config.comparison_prompt_id
        elif scope == "leading_edge":
            prompt_id = self._config.leading_edge_prompt_id
        else:
            raise KeyError(f"Unsupported prompt scope: {scope}")

        return resolve_prompt(
            scope=scope,
            prompt_id=prompt_id,
            variant=self._config.prompt_variant,
        )

    def _collection_payload(self, collection: CollectionSummary) -> Dict[str, object]:
        comparisons = self._comparison_entries(collection.comparisons)
        return {
            "collection": collection.display_name,
            "identifier": collection.identifier,
            "comparisons": comparisons,
            "comparisons_total": len(collection.comparisons),
            "significant_pathways": sum(c["significant_pathways"] for c in comparisons),
        }

    def _comparison_entries(self, comparisons: Sequence[ComparisonSummary]) -> list[Dict[str, object]]:
        limit = self._config.comparisons_per_collection
        if limit > 0:
            candidates = list(comparisons)[:limit]
        else:
            candidates = list(comparisons)

        entries: list[Dict[str, object]] = []
        for comparison in candidates:
            entries.append(
                {
                    "comparison": comparison.display_name,
                    "identifier": comparison.identifier,
                    "significant_pathways": comparison.significant_pathways,
                    "total_pathways": comparison.total_pathways,
                    "top_pathways": self._pathway_entries(comparison.top_pathways),
                }
            )
        return entries

    def _pathway_entries(self, pathways: Sequence[TopPathway]) -> list[Dict[str, object]]:
        limit = self._config.pathways_per_comparison
        if limit > 0:
            candidates = list(pathways)[:limit]
        else:
            candidates = list(pathways)

        entries: list[Dict[str, object]] = []
        for pathway in candidates:
            entries.append(
                {
                    "name": pathway.pathway,
                    "nes": pathway.nes,
                    "padj": pathway.padj,
                    "leading_genes": pathway.leading_genes,
                    "leading_gene_count": len(pathway.leading_gene_items),
                    "peak_rank_pct": pathway.peak_rank_pct,
                    "leading_edge_fraction": pathway.leading_edge_fraction,
                    "leading_edge_span_pct": pathway.leading_edge_span_pct,
                    "front_loaded_score": pathway.front_loaded_score,
                }
            )
        return entries

    def _comparison_payload(self, collection: CollectionSummary, comparison: ComparisonSummary) -> Dict[str, object]:
        return {
            "collection": {
                "name": collection.display_name,
                "identifier": collection.identifier,
            },
            "comparison": {
                "name": comparison.display_name,
                "identifier": comparison.identifier,
                "significant_pathways": comparison.significant_pathways,
                "total_pathways": comparison.total_pathways,
                "top_pathways": self._pathway_entries(comparison.top_pathways),
            },
        }

    def _leading_edge_payload(
        self,
        collection: CollectionSummary,
        comparison: ComparisonSummary,
        pathway: TopPathway,
    ) -> Dict[str, object]:
        return {
            "context": "leading_edge",
            "collection": collection.identifier,
            "collection_label": collection.display_name,
            "group": comparison.identifier,
            "group_label": comparison.display_name,
            "geneset": pathway.pathway,
            "nes": pathway.nes,
            "padj": pathway.padj,
            "leading_gene_count": len(pathway.leading_gene_items),
            "leading_genes_display": pathway.leading_genes,
            "peak_rank_pct": pathway.peak_rank_pct,
            "leading_edge_fraction": pathway.leading_edge_fraction,
            "leading_edge_span_pct": pathway.leading_edge_span_pct,
            "front_loaded_score": pathway.front_loaded_score,
            "items": self._leading_edge_items(pathway),
        }

    def _leading_edge_items(self, pathway: TopPathway) -> list[Dict[str, object]]:
        limit = self._config.leading_edge_gene_limit
        if limit > 0:
            genes = list(pathway.leading_gene_items)[:limit]
        else:
            genes = list(pathway.leading_gene_items)
        return [
            {
                "gene": gene,
                "rank": index + 1,
            }
            for index, gene in enumerate(genes)
        ]


class AgentApiSummarizer(OllamaSummarizer):
    """Generate report summaries through the same agent API used by tackle make-html."""

    def __init__(self, config: AgentApiConfig | None = None, *, cache_dir: Path | None = None) -> None:
        self._config = config or AgentApiConfig.from_env()
        self._cache_dir = Path(cache_dir).expanduser().resolve() if cache_dir else None
        self._session_prefix = f"tackle2-report:{uuid.uuid4().hex}"

    @property
    def config(self) -> AgentApiConfig:
        return self._config

    def set_cache_dir(self, cache_dir: Path | None) -> None:
        self._cache_dir = Path(cache_dir).expanduser().resolve() if cache_dir else None

    def _invoke_llm(self, prompt: str, *, system: str) -> Optional[Mapping[str, object]]:
        prompt_message = f"SYSTEM:\n{system}\n\nUSER:\n{prompt}"
        cache_dir = self._cache_dir
        if cache_dir is not None and not self._config.force_refresh:
            cached = _load_summary_cache(cache_dir, message=prompt_message)
            if cached is not None:
                return cached

        url = _build_support_chat_url(self._config.agent_api)
        session_digest = hashlib.sha256(prompt_message.encode("utf-8")).hexdigest()[:12]
        payload_obj: dict[str, Any] = {
            "sessionId": f"{self._session_prefix}:{session_digest}",
            "message": prompt_message,
            "meta": {
                "_queue_force_inline": True,
                "source": "tackle2_report",
                "model_label": self._config.model,
            },
        }
        payload = json.dumps(payload_obj).encode("utf-8")
        headers = {"Content-Type": "application/json"}
        if self._config.api_key:
            headers["X-API-Key"] = self._config.api_key
        if self._config.bearer_token:
            headers["Authorization"] = f"Bearer {self._config.bearer_token}"

        req = Request(url, data=payload, headers=headers, method="POST")
        try:
            with urlopen(req, timeout=self._config.timeout_seconds) as response:
                body = response.read().decode("utf-8", errors="replace")
        except (HTTPError, URLError, ConnectionError, TimeoutError, OSError, ValueError) as exc:
            logger.warning("Failed to obtain summary from agent API %s: %s", self._config.agent_api, exc)
            return None

        try:
            resp_obj = json.loads(body)
        except json.JSONDecodeError:
            resp_obj = {"message": body}

        if not isinstance(resp_obj, Mapping):
            resp_obj = {"message": body}

        text = str(resp_obj.get("message") or "").strip()
        if not text:
            logger.info("Agent API %s returned an empty summary response.", self._config.agent_api)
            return None

        summary_payload = _coerce_summary_payload(text)
        if summary_payload is None:
            return None

        if cache_dir is not None:
            _write_summary_cache(
                cache_dir,
                message=prompt_message,
                response=resp_obj,
                summary_payload=summary_payload,
                model_label=self._config.model,
            )

        return summary_payload


__all__ = [
    "AgentApiConfig",
    "AgentApiSummarizer",
    "OllamaConfig",
    "OllamaSummarizer",
    "ReportSummarizer",
    "SummaryRequest",
    "SummarisationError",
]
