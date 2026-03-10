"""Build the data context consumed by the HTML report template."""

from __future__ import annotations

import datetime as dt
import re
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from .catalog import SavedirArtefacts

try:  # Python >=3.11
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover - fallback for older runtimes
    import tomli as tomllib  # type: ignore[no-redef]


MAX_TOP_PATHWAYS = 12
LEADING_GENE_DISPLAY_LIMIT = 8
_NON_COLLECTION_DIRS = {
    "cache",
    "collections",
    "gsea_tables",
    "ranks",
    "report",
    "run_cache",
    "static",
    "tmp",
}


@dataclass
class GeneratedSummary:
    """Container for LLM-authored narrative snippets."""

    text: str
    bullets: Tuple[str, ...] = field(default_factory=tuple)
    model: Optional[str] = None


@dataclass
class TopPathway:
    pathway: str
    nes: Optional[float]
    padj: Optional[float]
    leading_genes: str
    peak_rank_pct: Optional[float] = None
    leading_edge_fraction: Optional[float] = None
    leading_edge_span_pct: Optional[float] = None
    front_loaded_score: Optional[float] = None
    leading_gene_items: Tuple[str, ...] = field(default_factory=tuple)
    llm_summary: Optional[GeneratedSummary] = None


@dataclass
class ComparisonSummary:
    identifier: str
    display_name: str
    table_path: Path
    relative_path: Path
    total_pathways: int
    significant_pathways: int
    top_pathways: List[TopPathway] = field(default_factory=list)
    llm_summary: Optional[GeneratedSummary] = None


@dataclass
class CollectionSummary:
    identifier: str
    display_name: str
    comparisons: List[ComparisonSummary] = field(default_factory=list)
    llm_summary: Optional[GeneratedSummary] = None


def _sanitize_collection_token(value: object) -> str:
    token = str(value or "").strip()
    if not token:
        return ""
    token = unicodedata.normalize("NFKD", token)
    token = token.encode("ascii", "ignore").decode("ascii")
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", token)
    token = re.sub(r"_+", "_", token).strip("_")
    return token


def _unique_tokens(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for raw in values:
        token = _sanitize_collection_token(raw)
        if not token or token in seen:
            continue
        seen.add(token)
        ordered.append(token)
    return ordered


def _config_collection_tokens(config_path: Optional[Path]) -> list[str]:
    if config_path is None:
        return []

    try:
        with Path(config_path).expanduser().resolve().open("rb") as fh:
            config_data = tomllib.load(fh)
    except (FileNotFoundError, OSError, tomllib.TOMLDecodeError):
        return []

    params = config_data.get("params") or {}
    genesets = params.get("genesets") or []
    if not isinstance(genesets, list):
        return []

    tokens: list[str] = []
    for entry in genesets:
        if not isinstance(entry, dict):
            continue
        collection_name = _sanitize_collection_token(entry.get("collection_name"))
        if collection_name:
            tokens.append(collection_name)
            continue

        category = str(entry.get("category") or "").strip()
        if not category:
            continue
        subcategory = str(entry.get("subcategory") or "").strip()
        tokens.append(_sanitize_collection_token(f"{category}_{subcategory}"))

    return _unique_tokens(tokens)


def _savedir_collection_tokens(savedir: Path, artefacts: SavedirArtefacts) -> list[str]:
    tokens: list[str] = []

    try:
        root_entries = sorted(savedir.iterdir(), key=lambda path: path.name.lower())
    except OSError:
        root_entries = []

    for entry in root_entries:
        if not entry.is_dir():
            continue
        if entry.name.startswith(".") or entry.name in _NON_COLLECTION_DIRS:
            continue
        tokens.append(entry.name)

    for file_path in list(artefacts.plot_files) + list(artefacts.log_files):
        try:
            rel = file_path.relative_to(savedir)
        except ValueError:
            continue
        if not rel.parts:
            continue
        root = rel.parts[0]
        if root.startswith(".") or root in _NON_COLLECTION_DIRS:
            continue
        tokens.append(root)

    for table_path in artefacts.gsea_tables:
        stem = table_path.stem
        if "__" in stem:
            tokens.append(stem.split("__", 1)[0])
        if stem.endswith("_all"):
            tokens.append(stem[:-4])

    return _unique_tokens(tokens)


def _known_collection_tokens(
    savedir: Path,
    artefacts: SavedirArtefacts,
    config_path: Optional[Path],
) -> tuple[str, ...]:
    tokens = _savedir_collection_tokens(savedir, artefacts)
    tokens.extend(_config_collection_tokens(config_path))
    deduped = _unique_tokens(tokens)
    deduped.sort(key=lambda token: (-len(token), token.lower()))
    return tuple(deduped)


def _split_table_name(
    table_path: Path,
    collection_tokens: Iterable[str] = (),
) -> tuple[str, str]:
    stem = table_path.stem
    if "__" in stem:
        collection, comparison = stem.split("__", 1)
        return collection, comparison

    for collection in collection_tokens:
        if stem == collection:
            return collection, "comparison"
        if stem.startswith(f"{collection}__"):
            comparison = stem[len(collection) + 2 :]
            return collection, comparison or "comparison"
        if stem.startswith(f"{collection}_"):
            comparison = stem[len(collection) + 1 :]
            return collection, comparison or "comparison"

    if "_" in stem:
        collection, comparison = stem.split("_", 1)
        return collection, comparison

    collection, comparison = stem, "comparison"
    return collection, comparison


def _reportable_table_entries(
    table_paths: Iterable[Path],
    collection_tokens: Iterable[str],
) -> list[tuple[Path, str, str]]:
    entries = [
        (table_path, *_split_table_name(table_path, collection_tokens))
        for table_path in table_paths
    ]
    collections_with_real_comparisons = {
        collection_token
        for _, collection_token, comparison_token in entries
        if comparison_token.lower() != "all"
    }
    return [
        (table_path, collection_token, comparison_token)
        for table_path, collection_token, comparison_token in entries
        if comparison_token.lower() != "all" or collection_token not in collections_with_real_comparisons
    ]


def _format_display(token: str) -> str:
    return token.replace("_", " ").strip() or token


def _display_path(path: Path) -> str:
    resolved = Path(path).expanduser().resolve()
    parts = resolved.parts
    lowered = [part.lower() for part in parts]

    for anchor in ("results", "gsea", "report"):
        if anchor in lowered:
            index = lowered.index(anchor)
            return "/".join(parts[index:])

    cwd = Path.cwd().resolve()
    try:
        relative = resolved.relative_to(cwd)
    except ValueError:
        relative = None
    if relative is not None:
        relative_text = relative.as_posix()
        if relative_text and relative_text != ".":
            return relative_text

    tail = [part for part in parts[-4:] if part not in (resolved.anchor, "/", "")]
    return "/".join(tail) or resolved.name


def _parse_leading_genes(value: object) -> Tuple[str, ...]:
    if not isinstance(value, str):
        return ()
    parts = [token.strip() for token in re.split(r"[\/,;]\s*", value) if token.strip()]
    return tuple(parts)


def _optional_float(value: object) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _summarise_table(
    table_path: Path,
    savedir: Path,
    collection_tokens: Iterable[str] = (),
) -> ComparisonSummary:
    df = pd.read_csv(table_path, sep="\t")
    total = len(df.index)
    padj_series = df.get("padj")

    significant = 0
    if padj_series is not None:
        significant = int(pd.to_numeric(padj_series, errors="coerce").lt(0.05).sum())

    padj_sorted = pd.to_numeric(df.get("padj"), errors="coerce")
    nes_sorted = pd.to_numeric(df.get("NES"), errors="coerce")
    df = df.assign(_padj_sorted=padj_sorted, _nes_sorted=nes_sorted)

    top_rows = df.sort_values(
        by=["_padj_sorted", "_nes_sorted"],
        ascending=[True, False],
        na_position="last",
    ).head(MAX_TOP_PATHWAYS)

    pathways: List[TopPathway] = []
    for _, row in top_rows.iterrows():
        leading = row.get("leadingEdge_genesymbol")
        if pd.isna(leading) or not isinstance(leading, str):
            leading = row.get("leadingEdge")
        genes = _parse_leading_genes(leading)
        leading_display = ", ".join(genes[:LEADING_GENE_DISPLAY_LIMIT])
        if len(genes) > LEADING_GENE_DISPLAY_LIMIT:
            leading_display += " …"

        nes_value = row.get("NES")
        padj_value = row.get("padj")

        if pd.isna(nes_value):
            nes_value = None
        else:
            nes_value = float(nes_value)

        if pd.isna(padj_value):
            padj_value = None
        else:
            padj_value = float(padj_value)

        pathways.append(
            TopPathway(
                pathway=str(row.get("pathway", "")),
                nes=nes_value,
                padj=padj_value,
                leading_genes=leading_display,
                peak_rank_pct=_optional_float(row.get("peak_rank_pct")),
                leading_edge_fraction=_optional_float(row.get("leading_edge_fraction")),
                leading_edge_span_pct=_optional_float(row.get("leading_edge_span_pct")),
                front_loaded_score=_optional_float(row.get("front_loaded_score")),
                leading_gene_items=genes,
            )
        )

    collection_token, comparison_token = _split_table_name(table_path, collection_tokens)

    return ComparisonSummary(
        identifier=comparison_token,
        display_name=_format_display(comparison_token),
        table_path=table_path,
        relative_path=table_path.relative_to(savedir),
        total_pathways=total,
        significant_pathways=significant,
        top_pathways=pathways,
    )


def build_context(
    savedir: Path,
    artefacts: SavedirArtefacts,
    config_path: Optional[Path] = None,
) -> Dict:
    """Assemble a JSON-serialisable context for the report template."""

    savedir = Path(savedir)
    collections: Dict[str, CollectionSummary] = {}
    collection_tokens = _known_collection_tokens(savedir, artefacts, config_path)
    table_entries = _reportable_table_entries(artefacts.gsea_tables, collection_tokens)

    for table_path, collection_token, _comparison_token in table_entries:
        comparison_summary = _summarise_table(table_path, artefacts.savedir, collection_tokens)
        collection_summary = collections.setdefault(
            collection_token,
            CollectionSummary(
                identifier=collection_token,
                display_name=_format_display(collection_token),
            ),
        )
        collection_summary.comparisons.append(comparison_summary)

    for collection in collections.values():
        collection.comparisons.sort(key=lambda cs: cs.display_name.lower())

    global_significant = sum(c.significant_pathways for coll in collections.values() for c in coll.comparisons)
    total_pathways = sum(c.total_pathways for coll in collections.values() for c in coll.comparisons)
    total_comparisons = sum(len(coll.comparisons) for coll in collections.values())

    plots = [
        {
            "path": path.relative_to(savedir),
            "name": path.relative_to(savedir).as_posix(),
        }
        for path in artefacts.plot_files
    ][:24]

    logs = [
        {
            "path": path.relative_to(savedir),
            "name": path.name,
        }
        for path in artefacts.log_files
    ]

    now_utc = dt.datetime.now(dt.timezone.utc)

    context = {
        "generated_at": now_utc.strftime("%Y-%m-%d %H:%M UTC"),
        "savedir_name": savedir.name,
        "savedir_path": _display_path(savedir),
        "stats": {
            "collections": len(collections),
            "comparisons": total_comparisons,
            "tables": len(table_entries),
            "pathways": total_pathways,
            "significant_pathways": global_significant,
        },
        "collections": list(collections.values()),
        "plots": plots,
        "logs": logs,
        "has_tables": len(table_entries) > 0,
        "config_path": _display_path(Path(config_path)) if config_path else None,
        "llm_summary": None,
    }

    return context


__all__ = ["build_context", "CollectionSummary", "ComparisonSummary", "TopPathway", "GeneratedSummary"]
