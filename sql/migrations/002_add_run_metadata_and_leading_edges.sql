-- Migration 2: add_run_metadata_and_leading_edges

-- Generated from python.sqlite_schema; do not edit by hand.

CREATE TABLE IF NOT EXISTS SchemaInfo (
    singleton_id INTEGER PRIMARY KEY,
    schema_name TEXT NOT NULL DEFAULT 'bcm_gsea',
    schema_version INTEGER NOT NULL,
    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT ck_schema_info_singleton CHECK (singleton_id = 1)
);

CREATE TABLE IF NOT EXISTS SchemaMigrations (
    version INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    applied_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

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
);

CREATE TABLE IF NOT EXISTS LeadingEdgeMembers (
    leading_edge_id INTEGER PRIMARY KEY AUTOINCREMENT,
    result_id INTEGER NOT NULL,
    member_ordinal INTEGER NOT NULL,
    bio_idstr TEXT NOT NULL,
    gene_symbol TEXT,
    FOREIGN KEY (result_id) REFERENCES CollectionResults(result_id)
);

ALTER TABLE CollectionResults ADD COLUMN leadingEdge_genesymbol TEXT;

ALTER TABLE CollectionResults ADD COLUMN peak_rank_pct REAL;

ALTER TABLE CollectionResults ADD COLUMN leading_edge_fraction REAL;

ALTER TABLE CollectionResults ADD COLUMN leading_edge_span_pct REAL;

ALTER TABLE CollectionResults ADD COLUMN front_loaded_score REAL;

CREATE UNIQUE INDEX IF NOT EXISTS idx_collections_name ON Collections(name);

CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_id ON Ranks(rankobj_id);

CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_rank ON Ranks(rankobj_id, rank);

CREATE INDEX IF NOT EXISTS idx_pathways_collection_name ON Pathways(collection_id, name);

CREATE UNIQUE INDEX IF NOT EXISTS idx_collection_results_lookup
ON CollectionResults(collection_id, rankobj_id, pathway_id);

CREATE INDEX IF NOT EXISTS idx_collection_results_collection_rank
ON CollectionResults(collection_id, rankobj_id);

CREATE UNIQUE INDEX IF NOT EXISTS idx_leading_edge_members_result_ordinal
ON LeadingEdgeMembers(result_id, member_ordinal);

CREATE INDEX IF NOT EXISTS idx_curves_lookup
ON Curves(rankobj_id, pathway_id, rank);

INSERT OR REPLACE INTO SchemaInfo (singleton_id, schema_name, schema_version, updated_at)
VALUES (1, 'bcm_gsea', 2, CURRENT_TIMESTAMP);

INSERT OR REPLACE INTO SchemaMigrations (version, name, applied_at)
VALUES (2, 'add_run_metadata_and_leading_edges', CURRENT_TIMESTAMP);
