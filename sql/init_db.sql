-- Generated from python.sqlite_schema; do not edit by hand.

CREATE TABLE IF NOT EXISTS "Collections" (
	collection_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	name TEXT
);

CREATE TABLE IF NOT EXISTS "RankObjects" (
	rankobj_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	name TEXT NOT NULL,
	species TEXT DEFAULT 'Homo sapiens'
);

CREATE TABLE IF NOT EXISTS "RunMetadata" (
	run_id INTEGER NOT NULL,
	savedir_name TEXT,
	savedir_path TEXT,
	config_path TEXT,
	config_sha256 TEXT,
	species TEXT DEFAULT 'Homo sapiens',
	ranks_from TEXT,
	created_at TEXT DEFAULT CURRENT_TIMESTAMP NOT NULL,
	updated_at TEXT DEFAULT CURRENT_TIMESTAMP NOT NULL,
	PRIMARY KEY (run_id),
	CONSTRAINT ck_run_metadata_singleton CHECK (run_id = 1)
);

CREATE TABLE IF NOT EXISTS "SchemaInfo" (
	singleton_id INTEGER NOT NULL,
	schema_name TEXT DEFAULT 'bcm_gsea' NOT NULL,
	schema_version INTEGER NOT NULL,
	updated_at TEXT DEFAULT CURRENT_TIMESTAMP NOT NULL,
	PRIMARY KEY (singleton_id),
	CONSTRAINT ck_schema_info_singleton CHECK (singleton_id = 1)
);

CREATE TABLE IF NOT EXISTS "SchemaMigrations" (
	version INTEGER NOT NULL,
	name TEXT NOT NULL,
	applied_at TEXT DEFAULT CURRENT_TIMESTAMP NOT NULL,
	PRIMARY KEY (version)
);

CREATE TABLE IF NOT EXISTS "Pathways" (
	pathway_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	name TEXT NOT NULL,
	ids TEXT,
	id_type TEXT DEFAULT 'entrez',
	species TEXT DEFAULT 'Homo sapiens',
	collection_id INTEGER,
	FOREIGN KEY(collection_id) REFERENCES "Collections" (collection_id)
);

CREATE TABLE IF NOT EXISTS "Ranks" (
	rank_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	rankobj_id INTEGER,
	id_type TEXT,
	bio_idstr TEXT NOT NULL,
	rank INTEGER,
	value FLOAT,
	FOREIGN KEY(rankobj_id) REFERENCES "RankObjects" (rankobj_id)
);

CREATE TABLE IF NOT EXISTS "CollectionResults" (
	result_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	rankobj_id INTEGER,
	collection_id INTEGER,
	pathway_id INTEGER,
	pathway TEXT,
	pval FLOAT,
	padj FLOAT,
	log2err FLOAT,
	es FLOAT,
	nes FLOAT,
	size INTEGER,
	"leadingEdge" TEXT,
	"leadingEdge_genesymbol" TEXT,
	mainpathway INTEGER,
	peak_rank_pct FLOAT,
	leading_edge_fraction FLOAT,
	leading_edge_span_pct FLOAT,
	front_loaded_score FLOAT,
	stored_at TEXT DEFAULT CURRENT_TIMESTAMP NOT NULL,
	FOREIGN KEY(rankobj_id) REFERENCES "RankObjects" (rankobj_id),
	FOREIGN KEY(collection_id) REFERENCES "Collections" (collection_id),
	FOREIGN KEY(pathway_id) REFERENCES "Pathways" (pathway_id)
);

CREATE TABLE IF NOT EXISTS "Curves" (
	curve_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	rankobj_id INTEGER,
	pathway_id INTEGER,
	rank FLOAT,
	"ES" FLOAT,
	FOREIGN KEY(rankobj_id) REFERENCES "RankObjects" (rankobj_id),
	FOREIGN KEY(pathway_id) REFERENCES "Pathways" (pathway_id)
);

CREATE TABLE IF NOT EXISTS "LeadingEdgeMembers" (
	leading_edge_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	result_id INTEGER NOT NULL,
	member_ordinal INTEGER NOT NULL,
	bio_idstr TEXT NOT NULL,
	gene_symbol TEXT,
	FOREIGN KEY(result_id) REFERENCES "CollectionResults" (result_id)
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_collections_name ON "Collections" (name);

CREATE INDEX IF NOT EXISTS idx_pathways_collection_name ON "Pathways" (collection_id, name);

CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_id ON "Ranks" (rankobj_id);

CREATE INDEX IF NOT EXISTS idx_ranks_rankobj_rank ON "Ranks" (rankobj_id, rank);

CREATE INDEX IF NOT EXISTS idx_collection_results_collection_rank ON "CollectionResults" (collection_id, rankobj_id);

CREATE UNIQUE INDEX IF NOT EXISTS idx_collection_results_lookup ON "CollectionResults" (collection_id, rankobj_id, pathway_id);

CREATE INDEX IF NOT EXISTS idx_curves_lookup ON "Curves" (rankobj_id, pathway_id, rank);

CREATE UNIQUE INDEX IF NOT EXISTS idx_leading_edge_members_result_ordinal ON "LeadingEdgeMembers" (result_id, member_ordinal);

INSERT OR REPLACE INTO SchemaInfo (singleton_id, schema_name, schema_version, updated_at) VALUES (1, 'bcm_gsea', 2, CURRENT_TIMESTAMP);

INSERT OR REPLACE INTO SchemaMigrations (version, name, applied_at) VALUES (2, 'add_run_metadata_and_leading_edges', CURRENT_TIMESTAMP);
