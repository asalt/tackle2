suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(stringr))

util_tools <- new.env()
source(file.path(here("R"), "utils.R"), local = util_tools)
log_msg <- util_tools$make_partial(util_tools$log_msg)

LATEST_SCHEMA_VERSION <- 2L
BASELINE_SCHEMA_VERSION <- 1L


get_con <- function(db_path = file.path(here("sql"), "rankorder_data.db")){
  con <- dbConnect(RSQLite::SQLite(), db_path)
  dbExecute(con, "PRAGMA foreign_keys = ON")
  return(con)
}

close_con <- function(con){
  dbDisconnect(con)
}

.read_sql_statements <- function(sql_file) {
  sql_lines <- readLines(sql_file, warn = FALSE)
  sql_lines <- gsub("--.*$", "", sql_lines)
  sql_script <- paste(sql_lines, collapse = "\n")
  sql_statements <- unlist(strsplit(sql_script, ";", fixed = TRUE), use.names = FALSE)
  sql_statements <- sql_statements %>%
    trimws() %>%
    (\(x) x[nzchar(x)])()
  return(sql_statements)
}


.execute_sql_file <- function(con, sql_file) {
  for (statement in .read_sql_statements(sql_file)) {
    log_msg(debug = statement)
    dbExecute(con, statement)
  }
}


get_schema_version <- function(con) {
  if (!dbExistsTable(con, "SchemaInfo")) {
    return(BASELINE_SCHEMA_VERSION)
  }
  res <- dbGetQuery(
    con,
    "SELECT schema_version FROM SchemaInfo WHERE singleton_id = 1 LIMIT 1"
  )
  if (nrow(res) == 0 || is.na(res$schema_version[[1]])) {
    return(BASELINE_SCHEMA_VERSION)
  }
  as.integer(res$schema_version[[1]])
}


.list_migration_files <- function(migrations_dir) {
  if (!dir.exists(migrations_dir)) {
    return(character())
  }
  files <- list.files(migrations_dir, pattern = "\\.sql$", full.names = TRUE)
  versions <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", basename(files))))
  files <- files[!is.na(versions)]
  versions <- versions[!is.na(versions)]
  files[order(versions)]
}


apply_pending_migrations <- function(con, migrations_dir = file.path(here("sql"), "migrations")) {
  current_version <- get_schema_version(con)
  migration_files <- .list_migration_files(migrations_dir)
  if (length(migration_files) == 0) {
    return(invisible(current_version))
  }

  for (migration_file in migration_files) {
    version <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", basename(migration_file))))
    if (is.na(version) || version <= current_version) {
      next
    }
    log_msg(info = paste0("Applying migration ", basename(migration_file)))
    .execute_sql_file(con, migration_file)
    current_version <- version
  }

  invisible(current_version)
}


initialize_db <- function(
  db_path = file.path(here("sql"), "rankorder_data.db"),
  sql_file = file.path(here("sql"), "init_db.sql"),
  migrations_dir = file.path(here("sql"), "migrations")
) {
  log_msg(info = "Initializing database...")
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  con <- get_con(db_path)
  on.exit(dbDisconnect(con), add = TRUE)

  existing_tables <- dbListTables(con)
  is_fresh_db <- length(existing_tables) == 0

  if (is_fresh_db) {
    .execute_sql_file(con, sql_file)
  } else {
    apply_pending_migrations(con, migrations_dir = migrations_dir)
  }

  log_msg(info = paste0("Database initialized successfully at schema v", get_schema_version(con), "."))
}


# returns collection_id
insert_collection <- function(con, collection_name){
    res <- dbGetQuery(con, 'select collection_id from Collections where name = ?', params = collection_name)
    if (nrow(res) > 0) {
      warning(paste0(collection_name, " already exists, skipping"))
      return(res$collection_id[1])
    }

    dbExecute(con, "INSERT INTO Collections (name) VALUES (?)", params = collection_name)
    res <- dbGetQuery(con, 'select collection_id from Collections where name = ?', params = collection_name)
    return(res$collection_id[1])
}

get_pathway_id <- function(con, pathway_name, collection_id = NULL, collection_name = NULL){
    if (is.null(pathway_name)) return(NULL)
    if (is.null(collection_id) && !is.null(collection_name)) {
      res <- dbGetQuery(con, 'select collection_id from Collections where name = ?', params = collection_name)
      if (nrow(res) > 0) {
        collection_id <- res$collection_id[1]
      }
    }
    if (!is.null(collection_id)) {
      res <- dbGetQuery(
        con,
        'select pathway_id from Pathways where name = ? and collection_id = ? LIMIT 1',
        params = list(pathway_name, collection_id)
      )
    } else {
      res <- dbGetQuery(con, 'select pathway_id from Pathways where name = ? LIMIT 1', params = pathway_name)
    }
    if (nrow(res) == 0) return(NULL)
    return(res$pathway_id[1])
}

insert_pathway <- function(con, collection_id = NULL, collection_name = NULL, pathway_name = NULL, members = NULL, id_type = "entrez"){
    # dbExecute(con, "INSERT INTO Collection (name) VALUES (?)", params = collection_name)

    if (is.null(pathway_name)) pathway_name <- "empty"
    if (is.null(members)) members <- ""

    if (is.null(collection_id) && !is.null(collection_name)){
      res <- dbGetQuery(con, 'select collection_id from Collections where name = ?', params = collection_name)
      if (nrow(res) == 0){ # then create
          collection_id <- insert_collection(con, collection_name)
      } else{
          collection_id <- res$collection_id[1]
      }
    }

    maybe_pathway_id <- get_pathway_id(con, pathway_name, collection_id = collection_id, collection_name = collection_name)

    if (!is.null(maybe_pathway_id)) {
        warning(paste0("pathway ", pathway_name, " already present in db"))
        return()
    }

    dbExecute(con, "INSERT INTO Pathways (name, ids, id_type, collection_id) VALUES (?, ?, ?, ?)",
     params = c(pathway_name, str_c(members, collapse='/'), id_type, collection_id)
    )

    pathway_id <- get_pathway_id(con, pathway_name, collection_id = collection_id, collection_name = collection_name)
    return(pathway_id)

}


.col_or_default <- function(df, cols, default_val = NA) {
  for (col in cols) {
    if (col %in% colnames(df)) return(df[[col]])
  }
  rep(default_val, nrow(df))
}


.normalize_numeric <- function(x) {
  if (is.list(x)) {
    x <- vapply(x, function(val) as.numeric(val)[1], numeric(1), USE.NAMES = FALSE)
  }
  suppressWarnings(as.numeric(x))
}


.normalize_integer <- function(x) {
  if (is.list(x)) {
    x <- vapply(x, function(val) as.integer(val)[1], integer(1), USE.NAMES = FALSE)
  }
  suppressWarnings(as.integer(x))
}


.normalize_logical_integer <- function(x) {
  if (is.list(x)) {
    x <- vapply(x, function(val) as.logical(val)[1], logical(1), USE.NAMES = FALSE)
  }
  suppressWarnings(as.integer(as.logical(x)))
}


.normalize_text <- function(x) {
  if (is.list(x)) {
    x <- vapply(
      x,
      function(val) paste(as.character(val), collapse = "/"),
      character(1),
      USE.NAMES = FALSE
    )
  }
  as.character(x)
}


.split_joined_tokens <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(character())
  }
  value %>%
    strsplit("[/,;]") %>%
    unlist(use.names = FALSE) %>%
    trimws() %>%
    (\(x) x[nzchar(x)])()
}


.leading_edge_pairs <- function(raw_ids, gene_symbols) {
  ids <- .split_joined_tokens(raw_ids)
  symbols <- .split_joined_tokens(gene_symbols)
  count <- max(length(ids), length(symbols))
  if (count == 0) {
    return(data.frame())
  }

  if (length(ids) < count) {
    ids <- c(ids, rep(NA_character_, count - length(ids)))
  }
  if (length(symbols) < count) {
    symbols <- c(symbols, rep(NA_character_, count - length(symbols)))
  }

  data.frame(
    member_ordinal = seq_len(count),
    bio_idstr = as.character(ids),
    gene_symbol = as.character(symbols),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(bio_idstr) & nzchar(bio_idstr))
}


upsert_run_metadata <- function(
  con,
  savedir_path = NULL,
  config_path = NULL,
  species = NULL,
  ranks_from = NULL
) {
  savedir_name <- if (!is.null(savedir_path) && nzchar(savedir_path)) basename(savedir_path) else NULL
  config_sha256 <- NULL
  if (!is.null(config_path) && nzchar(config_path) && file.exists(config_path)) {
    config_sha256 <- digest::digest(config_path, algo = "sha256", file = TRUE)
  }

  metadata_params <- list(
    run_id = 1L,
    savedir_name = savedir_name,
    savedir_path = savedir_path,
    config_path = config_path,
    config_sha256 = config_sha256,
    species = species,
    ranks_from = ranks_from
  )

  dbExecute(
    conn = con,
    statement = paste(
      "INSERT INTO RunMetadata",
      "(run_id, savedir_name, savedir_path, config_path, config_sha256, species, ranks_from, updated_at)",
      "VALUES (?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)",
      "ON CONFLICT(run_id) DO UPDATE SET",
      "savedir_name=excluded.savedir_name,",
      "savedir_path=excluded.savedir_path,",
      "config_path=excluded.config_path,",
      "config_sha256=excluded.config_sha256,",
      "species=excluded.species,",
      "ranks_from=excluded.ranks_from,",
      "updated_at=CURRENT_TIMESTAMP"
    ),
    params = unname(metadata_params)
  )

  invisible(1L)
}


replace_leading_edge_members <- function(con, result_id, raw_ids = NULL, gene_symbols = NULL) {
  if (is.null(result_id) || is.na(result_id)) {
    return(invisible(NULL))
  }

  dbExecute(con, "DELETE FROM LeadingEdgeMembers WHERE result_id = ?", params = list(result_id))
  pairs <- .leading_edge_pairs(raw_ids, gene_symbols)
  if (nrow(pairs) == 0) {
    return(invisible(NULL))
  }

  stmt <- dbSendStatement(
    con,
    paste(
      "INSERT INTO LeadingEdgeMembers",
      "(result_id, member_ordinal, bio_idstr, gene_symbol)",
      "VALUES (?, ?, ?, ?)"
    )
  )
  on.exit(dbClearResult(stmt), add = TRUE)
  for (i in seq_len(nrow(pairs))) {
    dbBind(
      stmt,
      list(
        result_id,
        pairs$member_ordinal[[i]],
        pairs$bio_idstr[[i]],
        pairs$gene_symbol[[i]]
      )
    )
  }

  invisible(NULL)
}


.resolve_rankobj_id <- function(con, rankobj_id = NULL, rank_name = NULL) {
  if (!is.null(rankobj_id)) {
    return(rankobj_id)
  }
  if (is.null(rank_name)) {
    stop("both rankobj_id and rank_name cannot be null")
  }

  sql <- "SELECT rankobj_id from RankObjects where name = ? LIMIT 1"
  res <- dbGetQuery(con, sql, params = rank_name)
  if (nrow(res) == 0) {
    warning("no rank name found, creating")
    return(insert_rankobj(con = con, rank_name = rank_name))
  }

  res$rankobj_id[1]
}


.resolve_collection_id <- function(con, collection_id = NULL, collection_name = NULL) {
  if (!is.null(collection_id)) {
    return(collection_id)
  }
  if (is.null(collection_name)) {
    stop("both collection_id and collection_name cannot be null")
  }

  sql <- "SELECT collection_id from Collections where name = ? LIMIT 1"
  res <- dbGetQuery(con, sql, params = collection_name)
  if (nrow(res) == 0) {
    warning("no collection name found, creating")
    return(insert_collection(con = con, collection_name = collection_name))
  }

  res$collection_id[1]
}


.build_collection_result_vectors <- function(results) {
  list(
    pathway = results %>%
      .col_or_default(cols = c("pathway"), default_val = "") %>%
      as.character(),
    pval = results %>%
      .col_or_default(cols = c("pval", "PVAL"), default_val = NA_real_) %>%
      .normalize_numeric(),
    padj = results %>%
      .col_or_default(cols = c("padj", "PADJ"), default_val = NA_real_) %>%
      .normalize_numeric(),
    log2err = results %>%
      .col_or_default(cols = c("log2err", "LOG2ERR"), default_val = NA_real_) %>%
      .normalize_numeric(),
    es = results %>%
      .col_or_default(cols = c("es", "ES"), default_val = NA_real_) %>%
      .normalize_numeric(),
    nes = results %>%
      .col_or_default(cols = c("nes", "NES"), default_val = NA_real_) %>%
      .normalize_numeric(),
    size = results %>%
      .col_or_default(cols = c("size", "SIZE"), default_val = NA_integer_) %>%
      .normalize_integer(),
    mainpathway = results %>%
      .col_or_default(cols = c("mainpathway", "main_pathway", "mainPathway"), default_val = NA) %>%
      .normalize_logical_integer(),
    leading_edge = results %>%
      .col_or_default(cols = c("leadingEdge", "leading_edge", "leadingedge"), default_val = "") %>%
      .normalize_text(),
    leading_edge_genesymbol = results %>%
      .col_or_default(cols = c("leadingEdge_genesymbol", "leading_edge_genesymbol"), default_val = "") %>%
      .normalize_text(),
    peak_rank_pct = results %>%
      .col_or_default(cols = c("peak_rank_pct"), default_val = NA_real_) %>%
      .normalize_numeric(),
    leading_edge_fraction = results %>%
      .col_or_default(cols = c("leading_edge_fraction"), default_val = NA_real_) %>%
      .normalize_numeric(),
    leading_edge_span_pct = results %>%
      .col_or_default(cols = c("leading_edge_span_pct"), default_val = NA_real_) %>%
      .normalize_numeric(),
    front_loaded_score = results %>%
      .col_or_default(cols = c("front_loaded_score"), default_val = NA_real_) %>%
      .normalize_numeric()
  )
}


.collection_result_upsert_sql <- function() {
  paste(
    "INSERT INTO CollectionResults",
    "(rankobj_id, collection_id, pathway_id, pathway, pval, padj, log2err, es, nes, size, leadingEdge, leadingEdge_genesymbol, mainpathway, peak_rank_pct, leading_edge_fraction, leading_edge_span_pct, front_loaded_score)",
    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    "ON CONFLICT(collection_id, rankobj_id, pathway_id) DO UPDATE SET",
    "pathway=excluded.pathway,",
    "pval=excluded.pval,",
    "padj=excluded.padj,",
    "log2err=excluded.log2err,",
    "es=excluded.es,",
    "nes=excluded.nes,",
    "size=excluded.size,",
    "leadingEdge=excluded.leadingEdge,",
    "leadingEdge_genesymbol=excluded.leadingEdge_genesymbol,",
    "mainpathway=excluded.mainpathway,",
    "peak_rank_pct=excluded.peak_rank_pct,",
    "leading_edge_fraction=excluded.leading_edge_fraction,",
    "leading_edge_span_pct=excluded.leading_edge_span_pct,",
    "front_loaded_score=excluded.front_loaded_score,",
    "stored_at=CURRENT_TIMESTAMP"
  )
}


.collection_result_row_params <- function(
  rankobj_id,
  collection_id,
  pathway_id,
  row_index,
  result_vectors
) {
  list(
    rankobj_id = rankobj_id,
    collection_id = collection_id,
    pathway_id = pathway_id,
    pathway = result_vectors$pathway[[row_index]],
    pval = result_vectors$pval[[row_index]],
    padj = result_vectors$padj[[row_index]],
    log2err = result_vectors$log2err[[row_index]],
    es = result_vectors$es[[row_index]],
    nes = result_vectors$nes[[row_index]],
    size = result_vectors$size[[row_index]],
    leading_edge = result_vectors$leading_edge[[row_index]],
    leading_edge_genesymbol = result_vectors$leading_edge_genesymbol[[row_index]],
    mainpathway = result_vectors$mainpathway[[row_index]],
    peak_rank_pct = result_vectors$peak_rank_pct[[row_index]],
    leading_edge_fraction = result_vectors$leading_edge_fraction[[row_index]],
    leading_edge_span_pct = result_vectors$leading_edge_span_pct[[row_index]],
    front_loaded_score = result_vectors$front_loaded_score[[row_index]]
  )
}


.get_collection_result_id <- function(con, rankobj_id, collection_id, pathway_id) {
  dbGetQuery(
    con,
    paste(
      "SELECT result_id FROM CollectionResults",
      "WHERE collection_id = ? AND rankobj_id = ? AND pathway_id = ?",
      "LIMIT 1"
    ),
    params = list(collection_id, rankobj_id, pathway_id)
  )$result_id[[1]]
}


insert_results <- function(con, rankobj_id = NULL, rank_name = NULL, collection_id = NULL, collection_name = NULL, results = NULL){

  if (is.null(results) || !is.data.frame(results)) stop("results must be a data.frame")

  rankobj_id <- .resolve_rankobj_id(
    con = con,
    rankobj_id = rankobj_id,
    rank_name = rank_name
  )
  collection_id <- .resolve_collection_id(
    con = con,
    collection_id = collection_id,
    collection_name = collection_name
  )

  on.exit({
    dbRollback(con)
    message("Transaction rolled back due to an error.")
  }, add = TRUE)

  result_vectors <- .build_collection_result_vectors(results = results)

  dbBegin(con)
  sql_upsert <- .collection_result_upsert_sql()
  for (i in seq_len(nrow(results))) {
    pathway_name <- result_vectors$pathway[[i]]
    pathway_id <- get_pathway_id(
      con = con,
      pathway_name = pathway_name,
      collection_id = collection_id
    )
    if (is.null(pathway_id)) {
      pathway_id <- insert_pathway(
        con = con,
        collection_id = collection_id,
        pathway_name = pathway_name
      )
    }

    row_params <- .collection_result_row_params(
      rankobj_id = rankobj_id,
      collection_id = collection_id,
      pathway_id = pathway_id,
      row_index = i,
      result_vectors = result_vectors
    )
    dbExecute(
      conn = con,
      statement = sql_upsert,
      params = unname(row_params)
    )
    result_id <- .get_collection_result_id(
      con = con,
      rankobj_id = rankobj_id,
      collection_id = collection_id,
      pathway_id = pathway_id
    )
    replace_leading_edge_members(
      con,
      result_id = result_id,
      raw_ids = row_params$leading_edge,
      gene_symbols = row_params$leading_edge_genesymbol
    )
    if (i %% 1000 == 0) cat("Inserted row", i, "\n")
  }
  dbExecute(con, "COMMIT")
  message("Bulk insert with prepared statements completed.")
  on.exit(NULL, add = FALSE)

}

insert_rankobj <- function(con, rank_name){
    res <- dbGetQuery(con, 'select rankobj_id from RankObjects where name = ?', params = rank_name)
    if (nrow(res) > 0) {
        warning(paste0("rankname ", rank_name, " already exists, cannot insert"))
        return(res$rankobj_id[1])
    }
    dbExecute(con, "INSERT INTO RankObjects (name) VALUES (?)", params = rank_name)
    res <- dbGetQuery(con, 'select rankobj_id from RankObjects where name = ?', params = rank_name)
    return(res$rankobj_id[1])
}

insert_ranks <- function(con, rankobj_id = NULL, rankobj_name = NULL, ranks_data = NULL) {

  if (!"data.frame" %in% class(ranks_data)) {
    ranks_data <- data.frame(id = names(ranks_data), value=ranks_data, rank=rank(-ranks_data))
  }
  if (!"rank" %in% colnames(ranks_data)) {
    ranks_data[['rank']] <- rank(-ranks_data$value)
  }

  # Assuming ranks_data is a data frame with columns: rank, gene
  if (is.null(rankobj_id)) {
      sql <- "SELECT rankobj_id from RankObjects where name = ? LIMIT 1"
      res <- dbGetQuery(con, sql, params = rankobj_name)
      if (nrow(res) == 0){
        warning("no rank name found, creating")
        rankobj_id <- insert_rankobj(con, rankobj_name)
      } else{
        rankobj_id <- res$rankobj_id[1]
      }
  }

  on.exit({
    dbRollback(con)
    message("Transaction rolled back due to an error.")
  }, add = TRUE)

  dbBegin(con)
  # Prepare the statement once
  stmt <- dbSendStatement(con, "INSERT INTO Ranks (rankobj_id, id_type, bio_idstr, rank, value) VALUES (?, ?, ?, ?, ?)")
  # Bind and execute for each row
  for (i in 1:nrow(ranks_data)) {
    dbBind(stmt, list(rankobj_id, "entrez", ranks_data$id[i], ranks_data$rank[i], ranks_data$value[i]))
    # Optional: Show progress for large datasets
    if (i %% 1000 == 0) cat("Inserted row", i, "\n")
  }

  # Clear the statement and commit the transaction
  dbClearResult(stmt)
  dbExecute(con, "COMMIT")
  message("Bulk insert with prepared statements completed.")
  on.exit(NULL, add = FALSE)

}


insert_curve <- function(
  con,
  rankobj_id = NULL,
  rankobj_name = NULL,
  pathway_id = NULL,
  pathway_name = NULL,
  collection_id = NULL,
  collection_name = NULL,
  curve_data = NULL
){


  if (!"data.frame" %in% class(curve_data)) {
    stop("curve_data must be a data frame with columns: rank, stat")
  }
  if (!all(c("rank", "ES") %in% colnames(curve_data))) {
    stop("curve_data must have columns: rank, ES")
  }

  # get rankobj id if not provided but name is
  if (is.null(rankobj_id)) {
      sql <- "SELECT rankobj_id from RankObjects where name = ? LIMIT 1"
      res <- dbGetQuery(con, sql, params = rankobj_name)
      if (nrow(res) == 0){
        warning("no rank name found, creating")
        log_msg(warning="no rank name found, creating")
        rankobj_id <- insert_rankobj(con, rankobj_name)
      } else{
        rankobj_id <- res$rankobj_id[1]
      }
  }

  # get collection id if not provided but name is
  if (is.null(pathway_id)) {
      pathway_id <- get_pathway_id(con, pathway_name, collection_id = collection_id, collection_name = collection_name)
      if (is.null(pathway_id)){
        warning("no pathway id found, ")
        log_msg(warning="no pathway id found, ")
        return()
      }
  }


  # ==

  on.exit({
    dbRollback(con)
    message("Transaction rolled back.")
  }, add = TRUE)

  dbBegin(con)
  stmt <- dbSendStatement(con, "INSERT INTO Curves (rankobj_id, pathway_id, rank, ES) VALUES (?, ?, ?, ?)")
  # Bind and execute for each row
  for (i in seq_len(nrow(curve_data))) {
    dbBind(stmt, list(rankobj_id, pathway_id, curve_data$rank[i], curve_data$ES[i]))
    # Optional: Show progress for large datasets
    if (i %% 1000 == 0) cat("Inserted row", i, "\n")
  }

  # Clear the statement and commit the transaction
  dbClearResult(stmt)
  dbExecute(con, "COMMIT")
  message("Bulk insert with prepared statements completed.")
  on.exit(NULL, add = FALSE)

}

get_plot_enrichmentdata_by_pathway <- function(
  con,
  rankobj_id = NULL,
  rankobj_name = NULL,
  pathway_id = NULL,
  pathway_name = NULL
 ) {

  # get rankobj id if not provided but name is
  if (is.null(rankobj_id)) {
      sql <- "SELECT rankobj_id from RankObjects where name = ? LIMIT 1"
      res <- dbGetQuery(con, sql, params = rankobj_name)
      if (nrow(res) == 0){
        warning("no rank name found, creating")
        rankobj_id <- insert_rankobj(con, rankobj_name)
      } else{
        rankobj_id <- res[1,]
      }
  }

  # get collection id if not provided but name is
  if (is.null(pathway_id)) {
      sql <- "SELECT pathway_id from Pathways where name = ? LIMIT 1"
      res <- dbGetQuery(con, sql, params = pathway_name)
      if (nrow(res) == 0){
        warning("no pathway id found, ")
        return()
        # pathway_id <- insert_pathway(con, rankobj_name)
      } else{
        pathway_id <- res[1,]
      }
  }


  rankdata <- dbGetQuery(con, "SELECT * FROM Ranks WHERE rankobj_id = ?", params = list(rankobj_id))
  rankdata %<>% rename(stat = value)
  rankdata %<>% arrange(rank)

  curvedata <- dbGetQuery(con, "SELECT * FROM Curves WHERE rankobj_id = ? AND pathway_id = ?", params = list(rankobj_id, pathway_id))

  pathway <- dbGetQuery(con, "SELECT * FROM Pathways WHERE pathway_id = ?", params = list(pathway_id))
  members <- pathway$ids %>% strsplit("/") %>% unlist()

  ticks <- rankdata %>% filter(bio_idstr %in% members) %>% arrange(-stat)

  dbGetQuery(con, "SELECT * FROM Curves WHERE pathway_id = ?", params = list(pathway_id))


  maxAbsStat <- rankdata$stat %>% abs() %>% max()
  posES = max(curvedata$ES, na.rm = T)
  negES = min(curvedata$ES, na.rm = T)
  spreadES = posES - negES


  list(
    stats = rankdata %>% select(rank, stat),
    curve = curvedata %>% select(rank, ES),
    ticks = ticks %>% select(rank, stat),
    maxAbsStat = maxAbsStat,
    posES = posES,
    negES = negES,
    spreadES = spreadES
  )


}
