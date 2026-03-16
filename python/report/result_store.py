"""Result-store abstractions for report generation.

The HTML report still works from TSV files, but can now prefer SQLite when a
database is present and looks usable.
"""

from __future__ import annotations

import logging
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import pandas as pd

from .catalog import SavedirArtefacts

try:  # Python >=3.11
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


logger = logging.getLogger(__name__)

DEFAULT_DB_NAME = "gsea_results.sqlite"
MAX_TOP_PATHWAYS = 12
LEADING_GENE_DISPLAY_LIMIT = 8
NON_COLLECTION_DIRS = {
    "cache",
    "collections",
    "gsea_tables",
    "ranks",
    "report",
    "run_cache",
    "static",
    "tmp",
}


@dataclass(frozen=True)
class TopPathwayRecord:
    pathway: str
    nes: Optional[float]
    padj: Optional[float]
    leading_genes: str
    peak_rank_pct: Optional[float] = None
    leading_edge_fraction: Optional[float] = None
    leading_edge_span_pct: Optional[float] = None
    front_loaded_score: Optional[float] = None
    leading_gene_items: tuple[str, ...] = ()


@dataclass(frozen=True)
class ComparisonRecord:
    collection_id: str
    comparison_id: str
    table_path: Optional[Path]
    relative_path: Optional[Path]
    total_pathways: int
    significant_pathways: int
    top_pathways: tuple[TopPathwayRecord, ...] = ()


def _sanitize_collection_token(value: object) -> str:
    token = str(value or "").strip()
    return token


def _load_config_data(config_path: Optional[Path]) -> dict:
    if config_path is None:
        return {}

    try:
        with Path(config_path).expanduser().resolve().open("rb") as fh:
            config_data = tomllib.load(fh)
    except (FileNotFoundError, OSError, tomllib.TOMLDecodeError):
        return {}
    return config_data if isinstance(config_data, dict) else {}


def _config_db_path(savedir: Path, config_path: Optional[Path]) -> Optional[Path]:
    if config_path is None:
        return None

    config_data = _load_config_data(config_path)
    params = config_data.get("params") or {}
    db_cfg = params.get("db") or {}
    raw_path = str(db_cfg.get("path") or "").strip()
    if not raw_path or raw_path == "savedir":
        return savedir / DEFAULT_DB_NAME

    candidate = Path(raw_path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()

    return (Path(config_path).expanduser().resolve().parent / candidate).resolve()


def _config_collection_tokens(config_path: Optional[Path]) -> tuple[str, ...]:
    config_data = _load_config_data(config_path)
    params = config_data.get("params") or {}
    genesets = params.get("genesets") or []
    if not isinstance(genesets, list):
        return ()

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
    return tuple(sorted({token for token in tokens if token}, key=lambda token: (-len(token), token.lower())))


def _sorted_collection_tokens(tokens: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({token for token in tokens if token}, key=lambda token: (-len(token), token.lower())))


def resolve_db_path(savedir: Path, config_path: Optional[Path] = None) -> Path:
    db_path = _config_db_path(savedir, config_path)
    if db_path is not None:
        return db_path
    return savedir / DEFAULT_DB_NAME


def _split_table_name(stem: str, collection_tokens: Sequence[str] = ()) -> tuple[str, str]:
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
    return stem, "comparison"


def _known_collection_tokens(savedir: Path, artefacts: SavedirArtefacts, config_path: Optional[Path]) -> tuple[str, ...]:
    tokens: list[str] = []
    for entry in sorted(savedir.iterdir(), key=lambda path: path.name.lower()):
        if entry.is_dir() and not entry.name.startswith(".") and entry.name not in NON_COLLECTION_DIRS:
            tokens.append(entry.name)
    for table_path in artefacts.gsea_tables:
        stem = table_path.stem
        if "__" in stem:
            tokens.append(stem.split("__", 1)[0])
        if stem.endswith("_all"):
            tokens.append(stem[:-4])
    tokens.extend(_config_collection_tokens(config_path))
    return _sorted_collection_tokens(tokens)


def _parse_leading_genes(value: object) -> tuple[str, ...]:
    if not isinstance(value, str):
        return ()
    parts = [token.strip() for token in value.replace(";", "/").replace(",", "/").split("/") if token.strip()]
    return tuple(parts)


def _optional_float(value: object) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sorted_top_rows(df: pd.DataFrame) -> pd.DataFrame:
    padj_sorted = pd.to_numeric(df.get("padj"), errors="coerce")
    significant = int(padj_sorted.lt(0.05).sum()) if "padj" in df.columns else 0

    nes_sorted = pd.to_numeric(df.get("NES"), errors="coerce")
    if nes_sorted.isna().all():
        nes_sorted = pd.to_numeric(df.get("nes"), errors="coerce")

    ranked_df = df.assign(_padj_sorted=padj_sorted, _nes_sorted=nes_sorted)
    top_rows = ranked_df.sort_values(
        by=["_padj_sorted", "_nes_sorted"],
        ascending=[True, False],
        na_position="last",
    ).head(MAX_TOP_PATHWAYS)
    return top_rows.assign(_significant_count=significant)


def _top_pathway_record_from_row(row: pd.Series) -> TopPathwayRecord:
    leading = row.get("leadingEdge_genesymbol")
    if pd.isna(leading) or not isinstance(leading, str):
        leading = row.get("leadingEdge")
    genes = _parse_leading_genes(leading)
    leading_display = ", ".join(genes[:LEADING_GENE_DISPLAY_LIMIT])
    if len(genes) > LEADING_GENE_DISPLAY_LIMIT:
        leading_display += " …"

    nes_value = row.get("NES") if "NES" in row else row.get("nes")
    return TopPathwayRecord(
        pathway=str(row.get("pathway", "")),
        nes=_optional_float(nes_value),
        padj=_optional_float(row.get("padj")),
        leading_genes=leading_display,
        peak_rank_pct=_optional_float(row.get("peak_rank_pct")),
        leading_edge_fraction=_optional_float(row.get("leading_edge_fraction")),
        leading_edge_span_pct=_optional_float(row.get("leading_edge_span_pct")),
        front_loaded_score=_optional_float(row.get("front_loaded_score")),
        leading_gene_items=genes,
    )


def _comparison_record_from_dataframe(
    *,
    df: pd.DataFrame,
    collection_id: str,
    comparison_id: str,
    table_path: Optional[Path],
    relative_path: Optional[Path],
) -> ComparisonRecord:
    total = len(df.index)
    top_rows = _sorted_top_rows(df)
    significant = int(top_rows["_significant_count"].iloc[0]) if not top_rows.empty else 0
    top_pathways = [
        _top_pathway_record_from_row(row)
        for _, row in top_rows.iterrows()
    ]

    return ComparisonRecord(
        collection_id=collection_id,
        comparison_id=comparison_id,
        table_path=table_path,
        relative_path=relative_path,
        total_pathways=total,
        significant_pathways=significant,
        top_pathways=tuple(top_pathways),
    )


def _table_entries(
    *,
    table_paths: Sequence[Path],
    collection_tokens: Sequence[str],
) -> list[tuple[Path, str, str]]:
    return [
        (table_path, *_split_table_name(table_path.stem, collection_tokens))
        for table_path in table_paths
    ]


def _reportable_table_entries(
    entries: Sequence[tuple[Path, str, str]],
) -> list[tuple[Path, str, str]]:
    collections_with_real_comparisons = {
        collection_id
        for _, collection_id, comparison_id in entries
        if comparison_id.lower() != "all"
    }
    return [
        (table_path, collection_id, comparison_id)
        for table_path, collection_id, comparison_id in entries
        if comparison_id.lower() != "all" or collection_id not in collections_with_real_comparisons
    ]


def _read_tsv_comparison_record(
    *,
    savedir: Path,
    table_path: Path,
    collection_id: str,
    comparison_id: str,
) -> ComparisonRecord:
    df = pd.read_csv(table_path, sep="\t")
    return _comparison_record_from_dataframe(
        df=df,
        collection_id=collection_id,
        comparison_id=comparison_id,
        table_path=table_path,
        relative_path=table_path.relative_to(savedir),
    )


def _sqlite_table_names(con: sqlite3.Connection) -> set[str]:
    tables = pd.read_sql_query(
        sql="SELECT name FROM sqlite_master WHERE type='table' AND name IN ('CollectionResults', 'Collections', 'RankObjects')",
        con=con,
    )
    return set(tables["name"])


def _sqlite_results_query() -> str:
    return """
        SELECT
            collections.name AS collection_id,
            rank_objects.name AS comparison_id,
            COALESCE(collection_results.pathway, pathways.name) AS pathway,
            collection_results.pval AS pval,
            collection_results.padj AS padj,
            collection_results.log2err AS log2err,
            collection_results.es AS es,
            collection_results.nes AS NES,
            collection_results.size AS size,
            COALESCE(collection_results.leadingEdge, leading_edge_members.bio_ids) AS leadingEdge,
            COALESCE(collection_results.leadingEdge_genesymbol, leading_edge_members.gene_symbols) AS leadingEdge_genesymbol,
            collection_results.mainpathway AS mainpathway,
            collection_results.peak_rank_pct AS peak_rank_pct,
            collection_results.leading_edge_fraction AS leading_edge_fraction,
            collection_results.leading_edge_span_pct AS leading_edge_span_pct,
            collection_results.front_loaded_score AS front_loaded_score
        FROM CollectionResults AS collection_results
        JOIN Collections AS collections
          ON collections.collection_id = collection_results.collection_id
        JOIN RankObjects AS rank_objects
          ON rank_objects.rankobj_id = collection_results.rankobj_id
        LEFT JOIN Pathways AS pathways
          ON pathways.pathway_id = collection_results.pathway_id
        LEFT JOIN (
            SELECT
                ordered.result_id AS result_id,
                group_concat(ordered.bio_idstr, '/') AS bio_ids,
                group_concat(ordered.gene_symbol, '/') AS gene_symbols
            FROM (
                SELECT
                    result_id,
                    member_ordinal,
                    bio_idstr,
                    gene_symbol
                FROM LeadingEdgeMembers
                ORDER BY result_id, member_ordinal
            ) AS ordered
            GROUP BY ordered.result_id
        ) AS leading_edge_members
          ON leading_edge_members.result_id = collection_results.result_id
        ORDER BY collections.name, rank_objects.name
    """


def _load_sqlite_results_frame(db_path: Path) -> pd.DataFrame:
    with sqlite3.connect(str(db_path)) as con:
        return pd.read_sql_query(sql=_sqlite_results_query(), con=con)


class FileResultStore:
    def load_comparisons(
        self,
        *,
        savedir: Path,
        artefacts: SavedirArtefacts,
        config_path: Optional[Path] = None,
    ) -> list[ComparisonRecord]:
        collection_tokens = _known_collection_tokens(savedir, artefacts, config_path)
        entries = _table_entries(
            table_paths=artefacts.gsea_tables,
            collection_tokens=collection_tokens,
        )
        reportable_entries = _reportable_table_entries(entries)

        records: list[ComparisonRecord] = []
        for table_path, collection_id, comparison_id in reportable_entries:
            records.append(
                _read_tsv_comparison_record(
                    savedir=savedir,
                    table_path=table_path,
                    collection_id=collection_id,
                    comparison_id=comparison_id,
                )
            )
        return records


class SqliteResultStore:
    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = Path(db_path).expanduser().resolve() if db_path else None

    def resolved_db_path(self, savedir: Path, config_path: Optional[Path] = None) -> Path:
        if self.db_path is not None:
            return self.db_path
        return resolve_db_path(savedir=savedir, config_path=config_path)

    def is_available(self, savedir: Path, config_path: Optional[Path] = None) -> bool:
        db_path = self.resolved_db_path(savedir, config_path=config_path)
        if not db_path.exists() or not db_path.is_file():
            return False
        try:
            with sqlite3.connect(str(db_path)) as con:
                table_names = _sqlite_table_names(con=con)
        except sqlite3.DatabaseError:
            return False
        return {"CollectionResults", "Collections", "RankObjects"}.issubset(table_names)

    def load_comparisons(
        self,
        *,
        savedir: Path,
        artefacts: SavedirArtefacts,
        config_path: Optional[Path] = None,
    ) -> list[ComparisonRecord]:
        _ = artefacts
        db_path = self.resolved_db_path(savedir=savedir, config_path=config_path)
        if not self.is_available(savedir=savedir, config_path=config_path):
            raise FileNotFoundError(f"SQLite results DB not available: {db_path}")

        df = _load_sqlite_results_frame(db_path=db_path)

        if df.empty:
            return []

        records: list[ComparisonRecord] = []
        for (collection_id, comparison_id), group_df in df.groupby(["collection_id", "comparison_id"], sort=True):
            records.append(
                _comparison_record_from_dataframe(
                    df=group_df.reset_index(drop=True),
                    collection_id=_sanitize_collection_token(collection_id),
                    comparison_id=str(comparison_id),
                    table_path=None,
                    relative_path=None,
                )
            )
        return records


class HybridResultStore:
    def __init__(
        self,
        *,
        sqlite_store: Optional[SqliteResultStore] = None,
        file_store: Optional[FileResultStore] = None,
    ):
        self.sqlite_store = sqlite_store or SqliteResultStore()
        self.file_store = file_store or FileResultStore()

    def load_comparisons(
        self,
        *,
        savedir: Path,
        artefacts: SavedirArtefacts,
        config_path: Optional[Path] = None,
    ) -> list[ComparisonRecord]:
        if self.sqlite_store.is_available(savedir, config_path=config_path):
            try:
                return self.sqlite_store.load_comparisons(
                    savedir=savedir,
                    artefacts=artefacts,
                    config_path=config_path,
                )
            except Exception as exc:
                logger.warning("Falling back to TSV results after SQLite read failure: %s", exc)
        return self.file_store.load_comparisons(
            savedir=savedir,
            artefacts=artefacts,
            config_path=config_path,
        )


__all__ = [
    "ComparisonRecord",
    "FileResultStore",
    "HybridResultStore",
    "SqliteResultStore",
    "TopPathwayRecord",
    "resolve_db_path",
]
