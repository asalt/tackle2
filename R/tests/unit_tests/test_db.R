library(here)
library(testthat)
library(DBI)
library(RSQLite)
library(stringr)

TESTDB <- tempfile("tackle2-test-", fileext = ".db")

legacy_schema_sql <- paste(
  c(
    "CREATE TABLE IF NOT EXISTS RankObjects (",
    "    rankobj_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    name TEXT NOT NULL,",
    "    species TEXT DEFAULT 'Homo sapiens'",
    ");",
    "CREATE TABLE IF NOT EXISTS Ranks (",
    "    rank_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    rankobj_id INTEGER,",
    "    id_type TEXT,",
    "    bio_idstr TEXT NOT NULL,",
    "    rank INTEGER,",
    "    value REAL,",
    "    FOREIGN KEY (rankobj_id) REFERENCES RankObjects(rankobj_id),",
    "    UNIQUE (rankobj_id, bio_idstr)",
    ");",
    "CREATE TABLE IF NOT EXISTS Collections(",
    "    collection_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    name TEXT",
    ");",
    "CREATE TABLE IF NOT EXISTS CollectionResults (",
    "    result_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    rankobj_id INTEGER,",
    "    collection_id INTEGER,",
    "    pathway_id INTEGER,",
    "    pathway TEXT,",
    "    pval REAL,",
    "    padj REAL,",
    "    log2err REAL,",
    "    es REAL,",
    "    nes REAL,",
    "    size INTEGER,",
    "    leadingEdge TEXT,",
    "    mainpathway INTEGER,",
    "    stored_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,",
    "    FOREIGN KEY (rankobj_id) REFERENCES RankObjects(rankobj_id),",
    "    FOREIGN KEY (collection_id) REFERENCES Collections(collection_id),",
    "    FOREIGN KEY (pathway_id) REFERENCES Pathways(pathway_id),",
    "    UNIQUE (rankobj_id, collection_id, pathway_id)",
    ");",
    "CREATE TABLE IF NOT EXISTS Pathways (",
    "    pathway_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    name TEXT NOT NULL,",
    "    ids TEXT,",
    "    id_type TEXT DEFAULT 'entrez',",
    "    species TEXT DEFAULT 'Homo sapiens',",
    "    collection_id INTEGER,",
    "    FOREIGN KEY (collection_id) REFERENCES Collections(collection_id)",
    ");",
    "CREATE TABLE IF NOT EXISTS Curves (",
    "    curve_id INTEGER PRIMARY KEY AUTOINCREMENT,",
    "    rankobj_id INTEGER,",
    "    pathway_id INTEGER,",
    "    rank REAL,",
    "    ES REAL,",
    "    FOREIGN KEY (rankobj_id) REFERENCES RankObjects(rankobj_id),",
    "    FOREIGN KEY (pathway_id) REFERENCES Pathways(pathway_id)",
    ");"
  ),
  collapse = "\n"
)


make_temp_db <- function() {
  tempfile("tackle2-test-", fileext = ".db")
}

sim_tools <- new.env()
source(file.path(here("R"), "simulate.R"), local = sim_tools)

io_tools <- new.env()
source(file.path(here("R"), "./io.R"), local = io_tools)

db_tools <- new.env()
source(file.path(here("R"), "db.R"), local = db_tools)

fgsea_tools <- new.env()
source(file.path(here("R"), "fgsea.R"), local = fgsea_tools)

geneset_tools <- new.env()
source(file.path(here("R"), "geneset_utils.R"), local = geneset_tools)


testthat::test_that("test db table setup", {
  db_path <- make_temp_db()
  on.exit(if (file.exists(db_path)) unlink(db_path), add = TRUE)
  db_tools$initialize_db(db_path = db_path)
  con <- db_tools$get_con(db_path)
  on.exit(db_tools$close_con(con), add = TRUE)
  res <- dbGetQuery(con, "SELECT * FROM sqlite_master WHERE type='table';")
  for (table in c("Ranks", "RankObjects", "Collections", "Pathways", "CollectionResults", "Curves", "SchemaInfo", "SchemaMigrations", "RunMetadata", "LeadingEdgeMembers")){
    expect_true(table %in% res$name, info = table)
  }
  schema_info <- dbGetQuery(con, "SELECT schema_version FROM SchemaInfo WHERE singleton_id = 1")
  expect_equal(schema_info$schema_version[[1]], 2)
}
)


setup <- function() {

  print(TESTDB)
  db_tools$initialize_db(db_path=TESTDB)
  con <- db_tools$get_con(TESTDB)


  data1 <- sim_tools$simulate_preranked_data(seed = 1234, sample_frac = .4)
  data2 <- sim_tools$simulate_preranked_data(seed = 4321, sample_frac = .4)
  data <- list(first = data1, second = data2)

  test_data <- sim_tools$generate_test_data(pathways = c("H", "GO:BP"), preranked_data = data)

  genesets_info <- list(
                    "H_" = geneset_tools$get_collection(category="H", subcategory = ""),
                    "C5_GO:BP" = geneset_tools$get_collection(category="C5", subcategory = "GO:BP")
                    )
  genesets_list_of_lists <- purrr::map(genesets_info, geneset_tools$genesets_df_to_list)

  collections <- names(test_data)
  rank_names <- names(test_data[[1]])

  # INSERT rank names
  #for (rank_name in rank_names){
  for (ix in seq_along(rank_names)){
      rankobj_id <- db_tools$insert_rankobj(con = con, rank_name = rank_names[[ix]])
      db_tools$insert_ranks(con = con, rankobj_id = rankobj_id, ranks_data = data[[ix]])
  }

  # INSERT collections
  for (collection in collections){ # this is slow and could be optimized
      collection_id <- db_tools$insert_collection(con, collection)
      for (pathway_name in names(genesets_list_of_lists[[collection]])) {
        db_tools$insert_pathway(
            con,
            collection_id=collection_id,
            pathway_name = pathway_name,
            members = str_c(genesets_list_of_lists[[collection]][[pathway_name]])
                )
       }
  }


  # INSERT gsea results
  for (rank_name in rank_names){
    for (collection in collections){
      db_tools$insert_results(con, rank_name = rank_name, collection_name = collection,
      results = test_data[[collection]][[rank_name]]
      )
    }
  }


  rankobjs <- io_tools$ranks_dfs_to_lists(data)

  db_tools$insert_ranks(con, rankobj_name = names(rankobjs)[[1]], ranks_data = rankobjs[[1]])
  db_tools$insert_ranks(con, rankobj_name = names(rankobjs)[[2]], ranks_data = rankobjs[[2]])

  rankorder1 <- fgsea_tools$get_rankorder(
    genesets_list_of_lists[[1]][[1]],
    rankobjs[[1]]
  )

  rankorder2 <- fgsea_tools$get_rankorder(
    genesets_list_of_lists[[1]][[1]],
    rankobjs[[2]]
  )



  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[1]],
    pathway_name = names(genesets_list_of_lists[[1]])[[1]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[1]],
    pathway_name = names(genesets_list_of_lists[[1]])[[2]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[2]],
    pathway_name = names(genesets_list_of_lists[[1]])[[1]],
    curve_data = rankorder1$curve
  )

  # INSERT
  db_tools$insert_curve(
    con,
    rankobj_name = names(rankobjs)[[2]],
    pathway_name = names(genesets_list_of_lists[[1]])[[2]],
    curve_data = rankorder1$curve
  )





    # pathway_name = names(genesets_list_of_lists[[1]])[[1]]

}


testthat::test_that("initialize_db migrates a legacy database to the latest schema", {
  db_path <- make_temp_db()
  on.exit(if (file.exists(db_path)) unlink(db_path), add = TRUE)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  legacy_statements <- strsplit(legacy_schema_sql, ";", fixed = TRUE)[[1]]
  legacy_statements <- trimws(legacy_statements)
  legacy_statements <- legacy_statements[nzchar(legacy_statements)]
  for (statement in legacy_statements) {
    DBI::dbExecute(con, statement)
  }
  DBI::dbDisconnect(con)

  db_tools$initialize_db(db_path = db_path)
  con <- db_tools$get_con(db_path)
  on.exit(db_tools$close_con(con), add = TRUE)

  schema_version <- db_tools$get_schema_version(con)
  expect_equal(schema_version, 2)

  columns <- DBI::dbListFields(con, "CollectionResults")
  expect_true("leadingEdge_genesymbol" %in% columns)
  expect_true("peak_rank_pct" %in% columns)
  expect_true("leading_edge_fraction" %in% columns)
  expect_true("leading_edge_span_pct" %in% columns)
  expect_true("front_loaded_score" %in% columns)
  expect_true(DBI::dbExistsTable(con, "RunMetadata"))
  expect_true(DBI::dbExistsTable(con, "LeadingEdgeMembers"))
})


testthat::test_that("insert_results stores shape metrics and normalized leading edges", {
  db_path <- make_temp_db()
  on.exit(if (file.exists(db_path)) unlink(db_path), add = TRUE)
  db_tools$initialize_db(db_path = db_path)
  con <- db_tools$get_con(db_path)
  on.exit(db_tools$close_con(con), add = TRUE)

  collection_id <- db_tools$insert_collection(con, "H")
  rankobj_id <- db_tools$insert_rankobj(con, "group_A")
  pathway_id <- db_tools$insert_pathway(
    con,
    collection_id = collection_id,
    pathway_name = "HALLMARK_APOPTOSIS",
    members = c("7157", "581")
  )
  expect_false(is.null(pathway_id))

  results_df <- data.frame(
    pathway = "HALLMARK_APOPTOSIS",
    pval = 0.001,
    padj = 0.01,
    log2err = 0.2,
    es = 0.8,
    NES = 2.1,
    size = 12L,
    leadingEdge = "7157/581",
    leadingEdge_genesymbol = "TP53/BAX",
    mainpathway = TRUE,
    peak_rank_pct = 0.125,
    leading_edge_fraction = 0.33,
    leading_edge_span_pct = 0.04,
    front_loaded_score = 16,
    stringsAsFactors = FALSE
  )

  db_tools$insert_results(
    con,
    rankobj_id = rankobj_id,
    collection_id = collection_id,
    results = results_df
  )

  stored <- DBI::dbGetQuery(
    con,
    "SELECT * FROM CollectionResults WHERE rankobj_id = ? AND collection_id = ?",
    params = list(rankobj_id, collection_id)
  )
  expect_equal(nrow(stored), 1)
  expect_equal(stored$peak_rank_pct[[1]], 0.125)
  expect_equal(stored$leading_edge_fraction[[1]], 0.33)
  expect_equal(stored$leading_edge_span_pct[[1]], 0.04)
  expect_equal(stored$front_loaded_score[[1]], 16)
  expect_equal(stored$leadingEdge_genesymbol[[1]], "TP53/BAX")

  members <- DBI::dbGetQuery(
    con,
    "SELECT member_ordinal, bio_idstr, gene_symbol FROM LeadingEdgeMembers WHERE result_id = ? ORDER BY member_ordinal",
    params = list(stored$result_id[[1]])
  )
  expect_equal(nrow(members), 2)
  expect_equal(members$bio_idstr, c("7157", "581"))
  expect_equal(members$gene_symbol, c("TP53", "BAX"))

  updated_df <- results_df
  updated_df$padj <- 0.02
  updated_df$leadingEdge <- "7157"
  updated_df$leadingEdge_genesymbol <- "TP53"
  db_tools$insert_results(
    con,
    rankobj_id = rankobj_id,
    collection_id = collection_id,
    results = updated_df
  )

  stored_updated <- DBI::dbGetQuery(
    con,
    "SELECT * FROM CollectionResults WHERE rankobj_id = ? AND collection_id = ?",
    params = list(rankobj_id, collection_id)
  )
  expect_equal(nrow(stored_updated), 1)
  expect_equal(stored_updated$padj[[1]], 0.02)

  members_updated <- DBI::dbGetQuery(
    con,
    "SELECT member_ordinal, bio_idstr, gene_symbol FROM LeadingEdgeMembers WHERE result_id = ? ORDER BY member_ordinal",
    params = list(stored_updated$result_id[[1]])
  )
  expect_equal(nrow(members_updated), 1)
  expect_equal(members_updated$bio_idstr[[1]], "7157")
  expect_equal(members_updated$gene_symbol[[1]], "TP53")
})


testthat::test_that("upsert_run_metadata stores a singleton run record", {
  db_path <- make_temp_db()
  on.exit(if (file.exists(db_path)) unlink(db_path), add = TRUE)
  db_tools$initialize_db(db_path = db_path)
  con <- db_tools$get_con(db_path)
  on.exit(db_tools$close_con(con), add = TRUE)

  config_path <- tempfile("tackle2-config-", fileext = ".toml")
  on.exit(if (file.exists(config_path)) unlink(config_path), add = TRUE)
  writeLines(c("[params]", "savedir = 'demo'"), config_path)

  db_tools$upsert_run_metadata(
    con,
    savedir_path = "/tmp/demo_savedir",
    config_path = config_path,
    species = "Homo sapiens",
    ranks_from = "gct"
  )

  metadata <- DBI::dbGetQuery(con, "SELECT * FROM RunMetadata WHERE run_id = 1")
  expect_equal(nrow(metadata), 1)
  expect_equal(metadata$savedir_name[[1]], "demo_savedir")
  expect_equal(metadata$ranks_from[[1]], "gct")
  expect_true(nzchar(metadata$config_sha256[[1]]))
})
