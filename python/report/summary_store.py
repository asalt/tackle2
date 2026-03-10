"""Persist and reload report summaries as JSON artifacts."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

from .. import export_packager
from . import catalog, summary
from .llm import ReportSummarizer, SummaryRequest
from .summary import CollectionSummary, ComparisonSummary, GeneratedSummary

logger = logging.getLogger(__name__)

SCHEMA_VERSION = 1


@dataclass(frozen=True)
class SummaryGenerationResult:
    savedir: Path
    output_dir: Path
    manifest_path: Path
    generated: int
    reused: int
    preview_only: bool
    files: tuple[Path, ...]


def default_summary_dir(savedir: Path) -> Path:
    return Path(savedir).expanduser().resolve() / "report" / "ai"


def resolve_savedir_and_summary_dir(
    *,
    savedir: Optional[Path] = None,
    config: Optional[Path] = None,
    output: Optional[Path] = None,
) -> tuple[Path, Optional[Path], Path]:
    resolved_savedir, resolved_config, _ = export_packager._find_savedir(  # type: ignore[attr-defined]
        Path(savedir) if savedir else None,
        Path(config) if config else None,
    )
    summary_dir = Path(output).expanduser().resolve() if output else default_summary_dir(resolved_savedir)
    summary_dir.mkdir(parents=True, exist_ok=True)
    return resolved_savedir, resolved_config, summary_dir


def _slugify_token(token: Optional[str], fallback: str = "item") -> str:
    if not token:
        token = fallback
    token = token.strip().lower()
    token = re.sub(r"[^a-z0-9]+", "-", token)
    token = re.sub(r"-{2,}", "-", token).strip("-")
    return token or fallback


def _prompt_hash(request: SummaryRequest, *, backend_label: str, model_label: str) -> str:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "backend": backend_label,
        "model_label": model_label,
        "scope": request.scope,
        "key": request.key,
        "prompt_id": request.prompt_id,
        "prompt_version": request.prompt_version,
        "prompt_variant": request.prompt_variant,
        "system_prompt": request.system_prompt,
        "instruction": request.instruction,
        "payload": request.payload,
        "user_notes": list(request.user_notes),
        "prompt_text": request.prompt_text,
    }
    raw = json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _request_doc_path(root: Path, request: SummaryRequest) -> Path:
    if request.scope == "run":
        return root / "run.json"
    if request.scope == "collection":
        slug = _slugify_token(request.target.get("collection_id"), fallback="collection")
        return root / "collections" / f"{slug}.json"
    if request.scope == "comparison":
        collection_slug = _slugify_token(request.target.get("collection_id"), fallback="collection")
        comparison_slug = _slugify_token(request.target.get("comparison_id"), fallback="comparison")
        return root / "comparisons" / collection_slug / f"{comparison_slug}.json"
    if request.scope == "leading_edge":
        collection_slug = _slugify_token(request.target.get("collection_id"), fallback="collection")
        comparison_slug = _slugify_token(request.target.get("comparison_id"), fallback="comparison")
        pathway_slug = _slugify_token(request.target.get("pathway_name"), fallback="pathway")
        return root / "leading_edge" / collection_slug / comparison_slug / f"{pathway_slug}.json"
    slug = _slugify_token(request.key, fallback="summary")
    return root / f"{slug}.json"


def _load_existing_doc(path: Path) -> Optional[Mapping[str, Any]]:
    if not path.exists() or not path.is_file():
        return None
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("Failed to read stored summary %s: %s", path, exc)
        return None
    if isinstance(doc, Mapping):
        return doc
    return None


def _generated_summary_from_doc(doc: Mapping[str, Any]) -> Optional[GeneratedSummary]:
    block = doc.get("summary")
    if not isinstance(block, Mapping):
        return None
    text = str(block.get("text") or "").strip()
    if not text:
        return None
    bullets_raw = block.get("bullets") or []
    bullets = tuple(str(item).strip() for item in bullets_raw if str(item).strip())
    model = str(block.get("model") or "").strip() or None
    return GeneratedSummary(text=text, bullets=bullets, model=model)


def _build_summary_doc(
    request: SummaryRequest,
    *,
    generated_summary: Optional[GeneratedSummary],
    prompt_hash: str,
    backend_label: str,
    model_label: str,
    status: str,
) -> Dict[str, Any]:
    summary_block = None
    if generated_summary is not None:
        summary_block = {
            "text": generated_summary.text,
            "bullets": list(generated_summary.bullets),
            "model": generated_summary.model,
        }

    return {
        "schema_version": SCHEMA_VERSION,
        "scope": request.scope,
        "key": request.key,
        "title": request.title,
        "target": dict(request.target),
        "status": status,
        "backend": backend_label,
        "model_label": model_label,
        "prompt": {
            "id": request.prompt_id,
            "version": request.prompt_version,
            "variant": request.prompt_variant,
            "system_prompt": request.system_prompt,
            "instruction": request.instruction,
            "payload": request.payload,
            "user_notes": list(request.user_notes),
            "prompt_text": request.prompt_text,
            "prompt_sha256": prompt_hash,
        },
        "summary": summary_block,
    }


def _write_summary_doc(path: Path, doc: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.parent / f".{path.name}.tmp"
    tmp_path.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def generate_and_store_summaries(
    *,
    savedir: Optional[Path] = None,
    config: Optional[Path] = None,
    output: Optional[Path] = None,
    summarizer: ReportSummarizer,
    include_run: bool = True,
    include_collections: bool = True,
    include_comparisons: bool = False,
    include_leading_edges: bool = False,
    prompt_notes: Sequence[str] = (),
    preview_only: bool = False,
    force: bool = False,
) -> SummaryGenerationResult:
    resolved_savedir, resolved_config, summary_dir = resolve_savedir_and_summary_dir(
        savedir=savedir,
        config=config,
        output=output,
    )
    artefacts = catalog.scan_savedir(resolved_savedir)
    context = summary.build_context(resolved_savedir, artefacts, config_path=resolved_config)
    collections = context.get("collections", [])

    summarizer.set_cache_dir(summary_dir / ".cache")
    requests = summarizer.build_requests(
        collections,
        include_run=include_run,
        include_collections=include_collections,
        include_comparisons=include_comparisons,
        include_leading_edges=include_leading_edges,
        user_notes=prompt_notes,
    )

    model_label = getattr(getattr(summarizer, "config", None), "model", None) or "agent"
    backend_label = summarizer.__class__.__name__
    generated = 0
    reused = 0
    written_files: list[Path] = []
    manifest_entries: list[dict[str, Any]] = []

    for request in requests:
        destination = _request_doc_path(summary_dir, request)
        prompt_hash = _prompt_hash(request, backend_label=backend_label, model_label=str(model_label))
        existing = _load_existing_doc(destination)

        if (
            not force
            and existing is not None
            and str(existing.get("prompt", {}).get("prompt_sha256") or "") == prompt_hash
        ):
            reused += 1
            written_files.append(destination)
            manifest_entries.append(
                {
                    "scope": request.scope,
                    "key": request.key,
                    "prompt_id": request.prompt_id,
                    "prompt_version": request.prompt_version,
                    "prompt_variant": request.prompt_variant,
                    "path": destination.relative_to(summary_dir).as_posix(),
                    "status": "reused",
                }
            )
            continue

        generated_summary = None if preview_only else summarizer.summarise_request(request)
        status = "preview_only" if preview_only else ("generated" if generated_summary else "empty")
        doc = _build_summary_doc(
            request,
            generated_summary=generated_summary,
            prompt_hash=prompt_hash,
            backend_label=backend_label,
            model_label=str(model_label),
            status=status,
        )
        _write_summary_doc(destination, doc)
        generated += 1
        written_files.append(destination)
        manifest_entries.append(
            {
                "scope": request.scope,
                "key": request.key,
                "prompt_id": request.prompt_id,
                "prompt_version": request.prompt_version,
                "prompt_variant": request.prompt_variant,
                "path": destination.relative_to(summary_dir).as_posix(),
                "status": status,
            }
        )

    manifest_path = summary_dir / "manifest.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "savedir": str(resolved_savedir),
        "config_path": str(resolved_config) if resolved_config else None,
        "output_dir": str(summary_dir),
        "backend": backend_label,
        "model_label": str(model_label),
        "preview_only": bool(preview_only),
        "force": bool(force),
        "prompt_notes": [str(note) for note in prompt_notes if str(note).strip()],
        "generated": generated,
        "reused": reused,
        "entries": manifest_entries,
    }
    _write_summary_doc(manifest_path, manifest)

    return SummaryGenerationResult(
        savedir=resolved_savedir,
        output_dir=summary_dir,
        manifest_path=manifest_path,
        generated=generated,
        reused=reused,
        preview_only=bool(preview_only),
        files=tuple(written_files),
    )


def apply_stored_summaries(
    *,
    context: Dict[str, Any],
    collections: Sequence[CollectionSummary],
    summary_dir: Optional[Path],
) -> None:
    if summary_dir is None:
        return
    root = Path(summary_dir).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        return

    run_doc = _load_existing_doc(root / "run.json")
    if run_doc:
        generated = _generated_summary_from_doc(run_doc)
        if generated is not None:
            context["llm_summary"] = generated

    by_collection = {collection.identifier: collection for collection in collections}
    for path in (root / "collections").glob("*.json"):
        doc = _load_existing_doc(path)
        if not doc:
            continue
        identifier = str(doc.get("target", {}).get("collection_id") or "").strip()
        collection = by_collection.get(identifier)
        if collection is None:
            continue
        generated = _generated_summary_from_doc(doc)
        if generated is not None:
            collection.llm_summary = generated

    for collection in collections:
        comparison_dir = root / "comparisons" / _slugify_token(collection.identifier, fallback="collection")
        if not comparison_dir.exists():
            continue
        by_comparison = {comparison.identifier: comparison for comparison in collection.comparisons}
        for path in comparison_dir.glob("*.json"):
            doc = _load_existing_doc(path)
            if not doc:
                continue
            comparison_id = str(doc.get("target", {}).get("comparison_id") or "").strip()
            comparison = by_comparison.get(comparison_id)
            if comparison is None:
                continue
            generated = _generated_summary_from_doc(doc)
            if generated is not None:
                comparison.llm_summary = generated

    leading_edge_root = root / "leading_edge"
    if not leading_edge_root.exists():
        return

    for collection in collections:
        comparison_root = leading_edge_root / _slugify_token(collection.identifier, fallback="collection")
        if not comparison_root.exists():
            continue
        by_comparison = {comparison.identifier: comparison for comparison in collection.comparisons}
        for comparison in collection.comparisons:
            pathway_root = comparison_root / _slugify_token(comparison.identifier, fallback="comparison")
            if not pathway_root.exists():
                continue
            by_pathway = {pathway.pathway: pathway for pathway in comparison.top_pathways}
            for path in pathway_root.glob("*.json"):
                doc = _load_existing_doc(path)
                if not doc:
                    continue
                comparison_id = str(doc.get("target", {}).get("comparison_id") or "").strip()
                if comparison_id and comparison_id != comparison.identifier:
                    continue
                pathway_name = str(doc.get("target", {}).get("pathway_name") or "").strip()
                pathway = by_pathway.get(pathway_name)
                if pathway is None:
                    continue
                generated = _generated_summary_from_doc(doc)
                if generated is not None:
                    pathway.llm_summary = generated


__all__ = [
    "SummaryGenerationResult",
    "apply_stored_summaries",
    "default_summary_dir",
    "generate_and_store_summaries",
    "resolve_savedir_and_summary_dir",
]
