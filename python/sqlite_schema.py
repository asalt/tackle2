"""Canonical SQLite schema definitions for tackle2 GSEA outputs.

This module keeps the SQLite schema formalized in one place while still
allowing the R runtime to bootstrap databases from committed SQL files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from sqlalchemy import (
    CheckConstraint,
    Column,
    Float,
    ForeignKey,
    Index,
    Integer,
    MetaData,
    Table,
    Text,
    text,
)
from sqlalchemy.dialects import sqlite
from sqlalchemy.schema import CreateIndex, CreateTable


SCHEMA_NAME = "bcm_gsea"
LATEST_SCHEMA_VERSION = 2


metadata = MetaData()


schema_info = Table(
    "SchemaInfo",
    metadata,
    Column("singleton_id", Integer, primary_key=True),
    Column("schema_name", Text, nullable=False, server_default=text(f"'{SCHEMA_NAME}'")),
    Column("schema_version", Integer, nullable=False),
    Column("updated_at", Text, nullable=False, server_default=text("CURRENT_TIMESTAMP")),
    CheckConstraint("singleton_id = 1", name="ck_schema_info_singleton"),
)


schema_migrations = Table(
    "SchemaMigrations",
    metadata,
    Column("version", Integer, primary_key=True),
    Column("name", Text, nullable=False),
    Column("applied_at", Text, nullable=False, server_default=text("CURRENT_TIMESTAMP")),
)


run_metadata = Table(
    "RunMetadata",
    metadata,
    Column("run_id", Integer, primary_key=True),
    Column("savedir_name", Text),
    Column("savedir_path", Text),
    Column("config_path", Text),
    Column("config_sha256", Text),
    Column("species", Text, server_default=text("'Homo sapiens'")),
    Column("ranks_from", Text),
    Column("created_at", Text, nullable=False, server_default=text("CURRENT_TIMESTAMP")),
    Column("updated_at", Text, nullable=False, server_default=text("CURRENT_TIMESTAMP")),
    CheckConstraint("run_id = 1", name="ck_run_metadata_singleton"),
)


rank_objects = Table(
    "RankObjects",
    metadata,
    Column("rankobj_id", Integer, primary_key=True),
    Column("name", Text, nullable=False),
    Column("species", Text, server_default=text("'Homo sapiens'")),
    sqlite_autoincrement=True,
)


ranks = Table(
    "Ranks",
    metadata,
    Column("rank_id", Integer, primary_key=True),
    Column("rankobj_id", Integer, ForeignKey("RankObjects.rankobj_id")),
    Column("id_type", Text),
    Column("bio_idstr", Text, nullable=False),
    Column("rank", Integer),
    Column("value", Float),
    sqlite_autoincrement=True,
)


collections = Table(
    "Collections",
    metadata,
    Column("collection_id", Integer, primary_key=True),
    Column("name", Text),
    sqlite_autoincrement=True,
)


pathways = Table(
    "Pathways",
    metadata,
    Column("pathway_id", Integer, primary_key=True),
    Column("name", Text, nullable=False),
    Column("ids", Text),
    Column("id_type", Text, server_default=text("'entrez'")),
    Column("species", Text, server_default=text("'Homo sapiens'")),
    Column("collection_id", Integer, ForeignKey("Collections.collection_id")),
    sqlite_autoincrement=True,
)


collection_results = Table(
    "CollectionResults",
    metadata,
    Column("result_id", Integer, primary_key=True),
    Column("rankobj_id", Integer, ForeignKey("RankObjects.rankobj_id")),
    Column("collection_id", Integer, ForeignKey("Collections.collection_id")),
    Column("pathway_id", Integer, ForeignKey("Pathways.pathway_id")),
    Column("pathway", Text),
    Column("pval", Float),
    Column("padj", Float),
    Column("log2err", Float),
    Column("es", Float),
    Column("nes", Float),
    Column("size", Integer),
    Column("leadingEdge", Text),
    Column("leadingEdge_genesymbol", Text),
    Column("mainpathway", Integer),
    Column("peak_rank_pct", Float),
    Column("leading_edge_fraction", Float),
    Column("leading_edge_span_pct", Float),
    Column("front_loaded_score", Float),
    Column("stored_at", Text, nullable=False, server_default=text("CURRENT_TIMESTAMP")),
    sqlite_autoincrement=True,
)


leading_edge_members = Table(
    "LeadingEdgeMembers",
    metadata,
    Column("leading_edge_id", Integer, primary_key=True),
    Column("result_id", Integer, ForeignKey("CollectionResults.result_id"), nullable=False),
    Column("member_ordinal", Integer, nullable=False),
    Column("bio_idstr", Text, nullable=False),
    Column("gene_symbol", Text),
    sqlite_autoincrement=True,
)


curves = Table(
    "Curves",
    metadata,
    Column("curve_id", Integer, primary_key=True),
    Column("rankobj_id", Integer, ForeignKey("RankObjects.rankobj_id")),
    Column("pathway_id", Integer, ForeignKey("Pathways.pathway_id")),
    Column("rank", Float),
    Column("ES", Float),
    sqlite_autoincrement=True,
)


Index("idx_ranks_rankobj_id", ranks.c.rankobj_id)
Index("idx_ranks_rankobj_rank", ranks.c.rankobj_id, ranks.c.rank)
Index("idx_collections_name", collections.c.name, unique=True)
Index("idx_pathways_collection_name", pathways.c.collection_id, pathways.c.name)
Index(
    "idx_collection_results_lookup",
    collection_results.c.collection_id,
    collection_results.c.rankobj_id,
    collection_results.c.pathway_id,
    unique=True,
)
Index(
    "idx_collection_results_collection_rank",
    collection_results.c.collection_id,
    collection_results.c.rankobj_id,
)
Index(
    "idx_leading_edge_members_result_ordinal",
    leading_edge_members.c.result_id,
    leading_edge_members.c.member_ordinal,
    unique=True,
)
Index("idx_curves_lookup", curves.c.rankobj_id, curves.c.pathway_id, curves.c.rank)


@dataclass(frozen=True)
class MigrationSpec:
    version: int
    name: str
    statements: tuple[str, ...]


MIGRATIONS: tuple[MigrationSpec, ...] = (
    MigrationSpec(
        version=2,
        name="add_run_metadata_and_leading_edges",
        statements=(
            """
CREATE TABLE IF NOT EXISTS SchemaInfo (
    singleton_id INTEGER PRIMARY KEY,
    schema_name TEXT NOT NULL DEFAULT 'bcm_gsea',
    schema_version INTEGER NOT NULL,
    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ck_schema_info_singleton CHECK (singleton_id = 1)
)
""".strip(),
            """
CREATE TABLE IF NOT EXISTS SchemaMigrations (
    version INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    applied_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
)
""".strip(),
            """
CREATE TABLE IF NOT EXISTS RunMetadata (
    run_id INTEGER PRIMARY KEY,
    savedir_name TEXT,
    savedir_path TEXT,
    config_path TEXT,
    config_sha256 TEXT,
    species TEXT DEFAULT 'Homo sapiens',
    ranks_from TEXT,
    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ck_run_metadata_singleton CHECK (run_id = 1)
)
""".strip(),
            """
CREATE TABLE IF NOT EXISTS LeadingEdgeMembers (
    leading_edge_id INTEGER PRIMARY KEY AUTOINCREMENT,
    result_id INTEGER NOT NULL,
    member_ordinal INTEGER NOT NULL,
    bio_idstr TEXT NOT NULL,
    gene_symbol TEXT,
    FOREIGN KEY (result_id) REFERENCES CollectionResults(result_id)
)
""".strip(),
            "ALTER TABLE CollectionResults ADD COLUMN leadingEdge_genesymbol TEXT",
            "ALTER TABLE CollectionResults ADD COLUMN peak_rank_pct REAL",
            "ALTER TABLE CollectionResults ADD COLUMN leading_edge_fraction REAL",
            "ALTER TABLE CollectionResults ADD COLUMN leading_edge_span_pct REAL",
            "ALTER TABLE CollectionResults ADD COLUMN front_loaded_score REAL",
            "CREATE UNIQUE INDEX IF NOT EXISTS idx_collections_name ON Collections(name)",
            "CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_id ON Ranks(rankobj_id)",
            "CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_rank ON Ranks(rankobj_id, rank)",
            "CREATE INDEX IF NOT EXISTS idx_pathways_collection_name ON Pathways(collection_id, name)",
            """
CREATE UNIQUE INDEX IF NOT EXISTS idx_collection_results_lookup
ON CollectionResults(collection_id, rankobj_id, pathway_id)
""".strip(),
            """
CREATE INDEX IF NOT EXISTS idx_collection_results_collection_rank
ON CollectionResults(collection_id, rankobj_id)
""".strip(),
            """
CREATE UNIQUE INDEX IF NOT EXISTS idx_leading_edge_members_result_ordinal
ON LeadingEdgeMembers(result_id, member_ordinal)
""".strip(),
            """
CREATE INDEX IF NOT EXISTS idx_curves_lookup
ON Curves(rankobj_id, pathway_id, rank)
""".strip(),
            """
INSERT OR REPLACE INTO SchemaInfo (singleton_id, schema_name, schema_version, updated_at)
VALUES (1, 'bcm_gsea', 2, CURRENT_TIMESTAMP)
""".strip(),
            """
INSERT OR REPLACE INTO SchemaMigrations (version, name, applied_at)
VALUES (2, 'add_run_metadata_and_leading_edges', CURRENT_TIMESTAMP)
""".strip(),
        ),
    ),
)


def _ddl_for_tables(tables: Iterable[Table]) -> list[str]:
    dialect = sqlite.dialect()
    statements: list[str] = []
    for table in tables:
        statement = str(CreateTable(table).compile(dialect=dialect)).strip()
        statement = statement.replace("CREATE TABLE ", "CREATE TABLE IF NOT EXISTS ", 1)
        statement = "\n".join(line.rstrip() for line in statement.splitlines())
        statements.append(statement)
    for table in tables:
        for index in sorted(table.indexes, key=lambda candidate: candidate.name or ""):
            statement = str(CreateIndex(index).compile(dialect=dialect)).strip()
            if statement.startswith("CREATE UNIQUE INDEX "):
                statement = statement.replace("CREATE UNIQUE INDEX ", "CREATE UNIQUE INDEX IF NOT EXISTS ", 1)
            elif statement.startswith("CREATE INDEX "):
                statement = statement.replace("CREATE INDEX ", "CREATE INDEX IF NOT EXISTS ", 1)
            statement = "\n".join(line.rstrip() for line in statement.splitlines())
            statements.append(statement)
    return statements


def render_sqlite_snapshot() -> str:
    """Render the latest SQLite bootstrap schema."""

    statements = [
        "-- Generated from python.sqlite_schema; do not edit by hand.",
        *[f"{statement};" for statement in _ddl_for_tables(metadata.sorted_tables)],
        (
            "INSERT OR REPLACE INTO SchemaInfo "
            "(singleton_id, schema_name, schema_version, updated_at) "
            f"VALUES (1, '{SCHEMA_NAME}', {LATEST_SCHEMA_VERSION}, CURRENT_TIMESTAMP);"
        ),
    ]

    for migration in MIGRATIONS:
        statements.append(
            "INSERT OR REPLACE INTO SchemaMigrations "
            "(version, name, applied_at) "
            f"VALUES ({migration.version}, '{migration.name}', CURRENT_TIMESTAMP);"
        )

    return "\n\n".join(statements) + "\n"


def render_migration_sql(migration: MigrationSpec) -> str:
    """Render a versioned migration file."""

    header = [
        f"-- Migration {migration.version}: {migration.name}",
        "-- Generated from python.sqlite_schema; do not edit by hand.",
    ]
    body = [f"{statement.strip()};" for statement in migration.statements]
    return "\n\n".join([*header, *body]) + "\n"


def migration_filename(migration: MigrationSpec) -> str:
    return f"{migration.version:03d}_{migration.name}.sql"


def migration_output_path(root: Path, migration: MigrationSpec) -> Path:
    return Path(root) / "sql" / "migrations" / migration_filename(migration)


def snapshot_output_path(root: Path) -> Path:
    return Path(root) / "sql" / "init_db.sql"


def write_schema_artifacts(root: Path) -> tuple[Path, tuple[Path, ...]]:
    """Write the current snapshot and migration SQL under ``sql/``."""

    root = Path(root)
    snapshot_path = snapshot_output_path(root)
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    snapshot_path.write_text(render_sqlite_snapshot(), encoding="utf-8")

    migration_dir = root / "sql" / "migrations"
    migration_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for migration in MIGRATIONS:
        path = migration_output_path(root, migration)
        path.write_text(render_migration_sql(migration), encoding="utf-8")
        written.append(path)
    return snapshot_path, tuple(written)


def committed_paths(root: Path) -> tuple[Path, tuple[Path, ...]]:
    root = Path(root)
    snapshot_path = snapshot_output_path(root)
    migration_paths = tuple(migration_output_path(root, migration) for migration in MIGRATIONS)
    return snapshot_path, migration_paths


def iter_migrations() -> Sequence[MigrationSpec]:
    return MIGRATIONS


__all__ = [
    "LATEST_SCHEMA_VERSION",
    "MIGRATIONS",
    "MigrationSpec",
    "SCHEMA_NAME",
    "committed_paths",
    "iter_migrations",
    "migration_filename",
    "migration_output_path",
    "render_migration_sql",
    "render_sqlite_snapshot",
    "snapshot_output_path",
    "write_schema_artifacts",
]
