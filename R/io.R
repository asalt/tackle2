suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(here))
# suppressPackageStartupMessages(library(janitor))
# src_dir <- file.path(here("R"))
# source(file.path(src_dir, "utils.R"))

# util_tools <- new.env()
# source(file.path(here("R"), "./utils.R"), local = util_tools)

source(file.path(here("R", "lazyloader.R")))
util_tools <- get_tool_env("utils")

# map_tools <- new.env()
# source(file.path(here("R"), "./map.R"), local = map_tools)
map_tools <- get_tool_env("map")
model_tools <- get_tool_env("modeling")

log_msg <- util_tools$make_partial(util_tools$log_msg)

# ==

normalize_rank_names <- function(x) {
  out <- fs::path_file(as.character(x))
  is_missing <- is.na(out)

  # Some cached rank files are created from volcano TSVs as <name>.tsv.rnk.
  # Strip only known rank/table suffixes so dotted sample labels remain intact.
  for (i in seq_len(3)) {
    out <- sub("\\.(rnk|tsv|csv|txt)(\\.gz)?$", "", out, ignore.case = TRUE, perl = TRUE)
  }
  out[is_missing] <- NA_character_
  out
}

dedupe_rank_names <- function(lst, context = "rank inputs") {
  rank_names <- names(lst)
  if (is.null(rank_names) || length(rank_names) == 0) {
    return(lst)
  }

  duplicated_labels <- rank_names[duplicated(rank_names)]
  if (length(duplicated_labels) > 0) {
    log_msg(warning = paste0(
      "duplicate rank names after extension normalization in ",
      context,
      "; keeping first occurrence for: ",
      paste(unique(duplicated_labels), collapse = ", ")
    ))
    lst <- lst[!duplicated(rank_names)]
  }
  lst
}

repair_duplicate_rank_labels <- function(labels, sep = ".") {
  labels <- as.character(labels)
  for (label in labels) {
    if (is.na(label) || !nzchar(label)) {
      stop("labels must be non-missing, non-empty strings")
    }
  }

  duplicated_labels <- labels[duplicated(labels)]
  if (length(duplicated_labels) == 0) {
    return(labels)
  }

  log_msg(warning = paste0(
    "duplicate rank labels detected: ",
    paste(unique(duplicated_labels), collapse = ", "),
    ". Repairing display labels with numeric suffixes."
  ))

  repaired <- character(length(labels))
  counts <- integer(0)
  used <- character(0)
  for (ix in seq_along(labels)) {
    label <- labels[[ix]]
    # Check for label membership in counts before indexing the named integer vector.
    current_count <- if (label %in% names(counts)) counts[[label]] else 0L
    counts[[label]] <- current_count + 1L
    candidate <- if (counts[[label]] == 1L) label else paste0(label, sep, counts[[label]]) # guaranteed to be in counts
    while (candidate %in% used) {
      counts[[label]] <- counts[[label]] + 1L
      candidate <- paste0(label, sep, counts[[label]])
    }
    repaired[[ix]] <- candidate
    used <- c(used, candidate)
  }
  repaired
}

read_rank_label_mapping <- function(rankfiledir) {
  if (is.null(rankfiledir)) {
    return(NULL)
  }

  name_mapping_file <- file.path(rankfiledir, "names.txt")
  if (!fs::file_exists(name_mapping_file)) {
    return(NULL)
  }

  log_msg(msg = paste0("found rank label mapping file: ", name_mapping_file))
  mapping <- read_delim(
    name_mapping_file,
    col_names = c("rank_label", "rank_source"),
    delim = "=",
    comment = "#",
    trim_ws = TRUE,
    show_col_types = FALSE
  )
  if (!all(c("rank_label", "rank_source") %in% colnames(mapping))) {
    stop("names.txt must contain lines of the form Rank Label=rank_file.rnk")
  }

  mapping <- mapping %>%
    dplyr::mutate(
      rank_label = trimws(as.character(.data$rank_label)),
      rankname = normalize_rank_names(.data$rank_source)
    ) %>%
    dplyr::filter(!is.na(.data$rankname), nzchar(.data$rankname))

  if (nrow(mapping) == 0) {
    return(NULL)
  }

  empty_labels <- !nzchar(mapping$rank_label) | is.na(mapping$rank_label)
  if (any(empty_labels)) {
    log_msg(warning = "names.txt contains empty labels; using canonical rank names for those entries")
    mapping$rank_label[empty_labels] <- mapping$rankname[empty_labels]
  }

  duplicated_ranknames <- mapping$rankname[duplicated(mapping$rankname)]
  if (length(duplicated_ranknames) > 0) {
    log_msg(warning = paste0(
      "names.txt references the same rank more than once: ",
      paste(unique(duplicated_ranknames), collapse = ", "),
      ". Keeping the first mapping for each rank."
    ))
    mapping <- mapping[!duplicated(mapping$rankname), , drop = FALSE]
  }

  mapping
}

make_rank_metadata <- function(
    rank_names,
    rankfiledir = NULL,
    rankname_order = NULL) {
  rank_names <- normalize_rank_names(rank_names)
  rank_names <- rank_names[!is.na(rank_names) & nzchar(rank_names)]
  rank_names <- rank_names[!duplicated(rank_names)]

  if (length(rank_names) == 0) {
    return(data.frame(
      rankname = character(),
      rank_label = character(),
      rank_order = integer(),
      rank_label_source = character(),
      rankname_factor = ordered(character()),
      rank_label_factor = ordered(character()),
      stringsAsFactors = FALSE
    ))
  }

  mapping <- read_rank_label_mapping(rankfiledir)
  mapped_order <- character(0)
  label_by_rank <- stats::setNames(rank_names, rank_names)
  label_source_by_rank <- stats::setNames(rep("default", length(rank_names)), rank_names)

  if (!is.null(mapping) && nrow(mapping) > 0) {
    missing_rankfiles <- setdiff(mapping$rankname, rank_names)
    if (length(missing_rankfiles) > 0) {
      log_msg(warning = paste0(
        "names.txt references rank files that were not loaded: ",
        paste(missing_rankfiles, collapse = ", "),
        ". These entries will be skipped."
      ))
    }
    mapping <- mapping[mapping$rankname %in% rank_names, , drop = FALSE]
    if (nrow(mapping) > 0) {
      mapped_order <- mapping$rankname
      label_by_rank[mapping$rankname] <- mapping$rank_label
      label_source_by_rank[mapping$rankname] <- "names.txt"
    }
  }

  requested_order <- rankname_order %||% NULL
  if (!is.null(requested_order)) {
    requested_order <- normalize_rank_names(as.character(unlist(requested_order, use.names = FALSE)))
    requested_order <- requested_order[!is.na(requested_order) & nzchar(requested_order)]
    if (length(requested_order) == 0) {
      requested_order <- NULL
    }
  }

  if (!is.null(requested_order)) {
    duplicated_ranknames <- requested_order[duplicated(requested_order)]
    if (length(duplicated_ranknames) > 0) {
      log_msg(warning = paste0(
        "rankname_order contains duplicate labels: ",
        paste(unique(duplicated_ranknames), collapse = ", "),
        ". Keeping the first occurrence of each."
      ))
    }
    requested_order <- unique(requested_order)
    missing_ranknames <- setdiff(requested_order, rank_names)
    if (length(missing_ranknames) > 0) {
      log_msg(warning = paste0(
        "rankname_order entries not present in ranks: ",
        paste(missing_ranknames, collapse = ", "),
        ". These names will be ignored."
      ))
    }
    ordered_ranknames <- c(intersect(requested_order, rank_names), setdiff(rank_names, requested_order))
  } else if (length(mapped_order) > 0) {
    ordered_ranknames <- c(mapped_order, setdiff(rank_names, mapped_order))
  } else {
    ordered_ranknames <- rank_names
  }

  metadata <- data.frame(
    rankname = ordered_ranknames,
    rank_label = unname(label_by_rank[ordered_ranknames]),
    rank_order = seq_along(ordered_ranknames),
    rank_label_source = unname(label_source_by_rank[ordered_ranknames]),
    stringsAsFactors = FALSE
  )
  metadata$rank_label <- repair_duplicate_rank_labels(metadata$rank_label)
  metadata$rankname_factor <- factor(metadata$rankname, levels = metadata$rankname, ordered = TRUE)
  metadata$rank_label_factor <- factor(metadata$rank_label, levels = metadata$rank_label, ordered = TRUE)
  rownames(metadata) <- metadata$rankname
  metadata
}

rank_label_map <- function(rank_metadata, ranknames = NULL) {
  if (is.null(rank_metadata) || !is.data.frame(rank_metadata) ||
      !all(c("rankname", "rank_label") %in% colnames(rank_metadata))) {
    if (is.null(ranknames)) {
      return(character(0))
    }
    return(stats::setNames(as.character(ranknames), as.character(ranknames)))
  }
  mapping <- stats::setNames(as.character(rank_metadata$rank_label), as.character(rank_metadata$rankname))
  if (!is.null(ranknames)) {
    missing <- setdiff(as.character(ranknames), names(mapping))
    if (length(missing) > 0) {
      mapping <- c(mapping, stats::setNames(missing, missing))
    }
    mapping <- mapping[as.character(ranknames)]
  }
  mapping
}

has_explicit_rank_labels <- function(rank_metadata, ranknames = NULL) {
  if (is.null(rank_metadata) || !is.data.frame(rank_metadata) ||
      !"rank_label_source" %in% colnames(rank_metadata)) {
    return(FALSE)
  }
  md <- rank_metadata
  if (!is.null(ranknames) && "rankname" %in% colnames(md)) {
    md <- md[md$rankname %in% as.character(ranknames), , drop = FALSE]
  }
  any(md$rank_label_source != "default")
}

finalize_rank_inputs <- function(
    rnkdfs,
    rankfiledir = NULL,
    rankname_order = NULL,
    context = "rank inputs") {
  rank_names <- names(rnkdfs)
  if (!is.null(rank_names) && length(rank_names) == length(rnkdfs)) {
    names(rnkdfs) <- normalize_rank_names(rank_names)
  }
  rnkdfs <- dedupe_rank_names(rnkdfs, context = context)
  rank_metadata <- make_rank_metadata(
    names(rnkdfs),
    rankfiledir = rankfiledir,
    rankname_order = rankname_order
  )
  if (nrow(rank_metadata) > 0) {
    rnkdfs <- rnkdfs[rank_metadata$rankname]
  }
  list(
    ranks = rnkdfs %>% ranks_dfs_to_lists(),
    rank_metadata = rank_metadata
  )
}

rank_context_label <- function(rank_name) {
  if (is.null(rank_name) || length(rank_name) == 0) {
    return("unnamed rank")
  }
  rank_name <- as.character(rank_name[[1]])
  if (is.na(rank_name) || !nzchar(rank_name)) {
    return("unnamed rank")
  }
  rank_name
}

prepare_rank_df <- function(rank_df, rank_name = NULL, warn = TRUE) {
  label <- rank_context_label(rank_name)
  if (!is.data.frame(rank_df)) {
    stop("rank input must be a data.frame")
  }

  if ("GeneID" %in% colnames(rank_df) && !"id" %in% colnames(rank_df)) {
    rank_df$id <- rank_df$GeneID
  }
  if (!all(c("id", "value") %in% colnames(rank_df))) {
    stop("rank input must contain id and value columns")
  }

  out <- data.frame(
    id = trimws(as.character(rank_df$id)),
    value = suppressWarnings(as.numeric(rank_df$value)),
    stringsAsFactors = FALSE
  )

  valid <- !is.na(out$id) & nzchar(out$id) & is.finite(out$value)
  dropped_invalid <- sum(!valid)
  if (dropped_invalid > 0 && warn) {
    log_msg(warning = paste0(
      "rank ", label, ": dropping ", dropped_invalid,
      " row(s) with missing ids or non-finite values"
    ))
  }
  out <- out[valid, , drop = FALSE]
  if (nrow(out) == 0) {
    if (warn) log_msg(warning = paste0("rank ", label, ": no usable rows remain"))
    return(out)
  }

  duplicate_rows <- duplicated(out$id)
  if (any(duplicate_rows) && warn) {
    log_msg(warning = paste0(
      "rank ", label, ": dropping ", sum(duplicate_rows),
      " duplicate gene id row(s); keeping first occurrence"
    ))
  }
  out[!duplicate_rows, , drop = FALSE]
}

make_random_gct <- function(nrow = 10, ncol = 4) {
  set.seed(369)
  nrow <- max(nrow, 1)
  ncol <- max(ncol, 1)
  .mat <- matrix(runif(nrow * ncol), nrow = nrow, ncol = ncol)
  .rids <- seq(1, dim(.mat)[1]) %>% as.character()
  rownames(.mat) <- .rids

  .cids <- seq(1, dim(.mat)[2]) %>% as.character()
  .cids <- paste0("X", .cids)
  .cdesc <- data.frame(
    id = .cids,
    metavar1 = sample(letters[1:5], ncol, replace = T),
    metavar2 = sample(letters[1:5], ncol, replace = T)
  )
  .rdesc <- data.frame(
    id = rownames(.mat),
    rdesc = rownames(.mat)
  )
  gct <- cmapR::GCT(
    mat = .mat,
    rid = .rids,
    cid = .cids,
    cdesc = .cdesc,
    rdesc = .rdesc
  )
  #
  return(gct)
  #
}

create_rnkfiles_from_emat <- function(
    emat,
    apply_z_score = FALSE,
    zscore_groupby = FALSE,
    sample_exclude = NULL,
    exclude_samples_from_data = FALSE,
    ...) {
  gct <- cmapR::parse_gctx(emat)


  #if ((!is.null(sample_exclude)) && (sample_exclude != FALSE)){
  # if (!is.null(sample_exclude) && !isFALSE(sample_exclude) && length(sample_exclude) > 0L) {
  sample_exclude_normalized <- util_tools$normalize_sample_exclude(sample_exclude, gct@cdesc)

  if (exclude_samples_from_data && length(sample_exclude_normalized) > 0L) {
    to_keep <- setdiff(gct@cid, sample_exclude_normalized)
    removed <- setdiff(gct@cid, to_keep)
    if (length(removed) > 0) {
      log_msg(info = paste0(
        "sample_exclude: dropping ",
        paste(removed, collapse = ", "),
        " from expression matrix prior to rank creation"
      ))
      gct <- cmapR::subset_gct(gct, cid = to_keep)
    }
  }


  if (apply_z_score) {
    gct <- util_tools$scale_gct(gct, group_by = zscore_groupby)
    # .new <- gct %>% cmapR::melt_gct()
    # .new <- gct@mat %>%
    #   apply(MARGIN = 1, FUN = .GlobalEnv$myzscore) %>%
    #   t() %>%
    #   as.matrix()
    # colnames(.new) <- colnames(mat(gct))
    # gct@mat <- .new
  }
  # gct@

  # Initialize a list to hold each new matrix
  list_of_matrices <- list()

  # Loop through each column of the matrix
  for (i in seq_len(ncol(gct@mat))) {
    # Create a new matrix for each column with row names and the column of interest
    # new_mat <- cbind(id = rownames(gct@mat), value = gct@mat[, i])
    new_mat <- data.frame(id = rownames(gct@mat), value = gct@mat[, i])

    # Convert the matrix to data frame for more intuitive row and column handling (optional)
    new_df <- as.data.frame(new_mat)

    # Store the matrix in the list
    list_of_matrices[[colnames(gct@mat)[i]]] <- new_df
  }

  # Output or return the list of matrices
  return(list_of_matrices)
}



create_rnkfiles_from_volcano <- function(
    volcanodir = "./",
    id_col = "GeneID",
    value_col = "value") {
  if (is.null(volcanodir)) {
    stop("volcanodir not defined")
  }

  if (!fs::dir_exists(volcanodir)) {
    stop("volcanodir does not exist")
  }

  (volcanofiles <- fs::dir_ls(path = volcanodir, regexp = ".*tsv", recurse = TRUE))
  log_msg(msg = paste0("Found ", length(volcanofiles), " tsv files"))
  log_msg(msg = paste(volcanofiles, collapse = "\n"))

  lst <- volcanofiles %>%
    purrr::set_names(nm = ~ normalize_rank_names(.)) %>%
    purrr::map(~ {
      .table <- read_tsv(.x, show_col_types = F)
      if (value_col %in% colnames(.table)) {
        .table <- dplyr::rename(.table, value = !!rlang::sym(value_col))
      }
      if (tolower(value_col) %in% tolower(colnames(.table))) {
        .match <- colnames(.table)[ stringr::str_detect(tolower(colnames(.table)), tolower(value_col)) ]
        if (length(.match) >= 1) {
          .table <- dplyr::rename(.table, value = !!rlang::sym(.match[[1]]))
        }
      }
      if (id_col %in% colnames(.table)) {
        .table <- dplyr::rename(.table, id = !!rlang::sym(id_col))
      }
      return(.table)
    })

  dedupe_rank_names(lst, context = "volcano files")
}


write_rnkfiles <- function(
    lst,
    dir = "rnkfiles") {
  if (is.null(dir)) {
    dir <- "rnkfiles"
  }
  if (!fs::dir_exists(dir)) {
    log_msg(msg = paste0("creating ", dir))
    fs::dir_create(dir, recurse = TRUE)
  }
  rank_names <- names(lst)
  if (!is.null(rank_names) && length(rank_names) == length(lst)) {
    names(lst) <- normalize_rank_names(rank_names)
  }
  lst <- dedupe_rank_names(lst, context = "rank files to write")
  lst %>% purrr::iwalk( # .x is the value, .y is the name
    ~ {
      .outname <- fs::path_join(
        c(dir, paste0(.y, ".rnk"))
      )
      if (!fs::file_exists(.outname)) {
        rank_df <- prepare_rank_df(.x, rank_name = .y)
        if (nrow(rank_df) == 0) {
          log_msg(warning = paste0("rank ", .y, ": skipping empty rank file write"))
          return(invisible(NULL))
        }
        rank_df %>%
          write_tsv(.outname, col_names = FALSE)
        log_msg(msg = paste0("Wrote ", .outname))
      }
    }
  )
}

load_rnkfiles <- function(rnkfiles) {
  data <- map(rnkfiles, ~ readr::read_tsv(.x,
    col_names = c("id", "value"),
    show_col_types = F
  ) %>%
    mutate(
      id = as.character(id),
      value = as.numeric(value)
    ) %>%
    # arrange(value) %>% # do not change order of files here
    drop_na())
  data
}


ranks_dfs_to_lists <- function(rnkdfs) {

  if (!"list" %in% class(rnkdfs)) rnkdfs <- list(rnkdfs)

  ranks_list <- rnkdfs %>% purrr::imap(
    ~ {
      rank_df <- prepare_rank_df(.x, rank_name = .y)
      stats::setNames(rank_df$value, rank_df$id)
    }
  )
  return(ranks_list)
}

load_genesets_from_json <- function(json_str) {
  genesets_of_interest <- jsonlite::fromJSON(json_str)
  genesets_of_interest <- genesets_of_interest %>% dplyr::mutate(
    collection_name = stringr::str_c(category, subcategory, sep = "_")
  )
  return(genesets_of_interest)
}


# save_gsea_results <- function(
#     results_list,
#     savedir = NULL) {
#   if (is.null(savedir)) savedir <- "gsea_tables"
#   if (!file.exists(savedir)) fs::dir_create(savedir)
#   names(results_list) %>%
#     purrr::map(
#       ~ {
#         collection_name <- .x
#         names(results_list[[collection_name]]) %>%
#           purrr::map(
#             ~ {
#               comparison_name <- .x
#               result <- results_list[[collection_name]][[comparison_name]]
#               # print(collection_name)
#               # print(comparison_name)
#
#               outf <- paste0(
#                 make.names(collection_name),
#                 "_",
#                 make.names(comparison_name),
#                 ".tsv"
#               )
#               outf <- file.path(savedir, outf)
#               # one last check here
#               result <- result %>% mutate(leadingEdge = purrr::map_chr(leadingEdge, paste, collapse = "/"))
#               log_msg(msg = paste0("Writing: ", outf, "..."))
#               if (is.data.frame(result)) {
#                 result %>% readr::write_tsv(outf)
#                 log_msg(msg = "done")
#               } else {
#                 log_msg(msg = "Invalid result, cannot write to file.")
#               }
#               # if (!fs::file_exists(outf)) result %>% readr::write_tsv(outf)
#             }
#           )
#       }
#     )
# }


collapse_list_col_chr <- function(x, collapse = "/") {
  if (!is.list(x)) return(x)
  purrr::map_chr(x, function(el) {
    if (is.null(el) || length(el) == 0) return(NA_character_)
    flat <- unlist(el, use.names = FALSE)
    flat <- as.character(flat)
    flat <- flat[!is.na(flat)]
    if (length(flat) == 0) return(NA_character_)
    paste(flat, collapse = collapse)
  })
}

prepare_results_for_export <- function(result, collapse = "/") {
  if (!is.data.frame(result)) {
    log_msg(msg = "Invalid result, cannot write to file.")
    log_msg(msg = as.character(result))
    stop("Invalid result, cannot write to file.")
  }

  list_cols <- vapply(result, is.list, logical(1))
  if (!any(list_cols)) {
    return(result)
  }

  out <- result
  for (col in names(out)[list_cols]) {
    out[[col]] <- collapse_list_col_chr(out[[col]], collapse = collapse)
  }
  out
}

write_results <- function(result, outf, replace = FALSE) {
  if (is.null(replace)) replace <- FALSE

  if (fs::file_exists(outf) && !replace) {
    log_msg(msg = paste0("File ", outf, " already exists, skipping"))
    return(invisible(NULL))
  }

  result_export <- prepare_results_for_export(result)
  result_export %>% write_tsv(outf)
  log_msg(msg = paste0("Successfully written to ", outf))
  return(invisible(outf))
}

# Main function to save GSEA results
save_individual_gsea_results <- function(
  results_list,
  savedir = "gsea_tables",
  replace = FALSE,
  species = "Homo sapiens",
  rank_metadata = NULL) {

if (is.null(replace)) replace <- FALSE

log_msg(msg = "writing results")
log_msg(msg = paste0("names results list :", names(results_list)))
log_msg(msg = paste0("length results list :", length(results_list)))

  fs::dir_create(savedir) # Ensures directory exists, no error if it already does
  results_list_towrite <- results_list %>% purrr::imap(~{
    result_list  <- .x #%>% map_tools$add_leadingedges_to_results_list()
    collection_name <- .y
    # Construct a comparison name map to shorten labels for files within this collection
    comparison_names <- names(result_list)
    name_map <- if (has_explicit_rank_labels(rank_metadata, comparison_names)) {
      rank_label_map(rank_metadata, comparison_names)
    } else {
      util_tools$make_name_map(comparison_names)
    }
    result_list %>% purrr::imap(~{
      result <- .x
      comparison_name <- .y
      comparison_label <- name_map[[comparison_name]] %||% comparison_name
      filename <- paste0(
        util_tools$safe_filename(collection_name, comparison_label, fallback = "gsea_result"),
        ".tsv"
      )
      outf <- file.path(savedir, filename)
      if (!"data.frame" %in% class(result)) {
        log_msg(paste0("Invalid result, cannot write to file."))
        return()
      }
      if (fs::file_exists(outf) && !replace) {
        log_msg(msg = paste0("File ", outf, " already exists, skipping"))
        # return(result)
      } else {
        write_results(result, outf, replace = replace)
      }
    })
  })


  # results_list_towrite <- results_list %>% purrr::imap(~{
  #   result_list  <- .x %>% map_tools$add_leadingedges_to_results_list()
  #   collection_name <- .y
  #   result_list %>% purrr::imap(~{
  #     result <- .x
  #     comparison_name <- .y
  #     outf <- file.path(savedir, make.names(paste0(collection_name, "_", comparison_name, ".tsv")))
  #     log_msg(paste0("Writing: ", outf, "..."))
  #     result %>% write_tsv(outf)
  #   })

}

save_pivoted_gsea_results <- function(results_list, savedir = "gsea_tables", replace = FALSE, species = "Homo sapiens") {

  if (is.null(replace)) replace <- FALSE
  fs::dir_create(savedir)
  # here results list is concatenated list. one level per collection
  # names are the collection names
  # values are the fgsea concatenated tables/comparisons for given collection

  results_list %<>% map_tools$add_leadingedges_to_results_list(species = species)

  results_list_pivoted <- results_list %>% purrr::imap(
    ~{
      result_collection  <- .x
      collection_name <- .y

      wanted <- c(
      "pval", "padj", "log2err", "ES", "NES",
      "n_main", "size",
      "leadingEdge",
      "leadingEdge_genesymbol",
      "leadingEdge_entrezid",
      "mainpathway"
    )
      present <- intersect(colnames(result_collection), wanted)

      result_collection_export <- prepare_results_for_export(result_collection)

      res_pivoted <- result_collection_export %>%
        pivot_wider(
          id_cols = c("pathway"),
          names_from = rankname,
          # values_from = c(pval, padj, log2err, ES, NES, n_main, size, leadingEdge_genesymbol, leadingEdge_entrezid, mainpathway),
          values_from = all_of(present),
          names_sep = "_"
          # Alternatively, use names_glue for more complex naming
          # names_glue = "{rankname}_{.value}"
        )
      filename <- paste0(
        util_tools$safe_filename(collection_name, "all", fallback = "gsea_result"),
        ".tsv"
      )
      outf <- file.path(savedir, filename)
      if (fs::file_exists(outf) && !replace) {
        log_msg(msg = paste0("File ", outf, " already exists, skipping"))
        return(res_pivoted)
      }
      log_msg(msg=paste0("Writing: ", outf, "..."))
      res_pivoted %>%
        prepare_results_for_export() %>%
        write_tsv(outf)
      return(res_pivoted)
    }
  )
  return(results_list_pivoted)
}


load_from_cache <- function(filename, cache_dir = NULL) {
  if (is.null(cache_dir)) {
    cache_dir <- here("cache")
  }
  if (!fs::dir_exists(cache_dir)) fs::dir_create(cache_dir)
  target_file <- paste0(file.path(cache_dir, filename), ".rds")
  if (!fs::file_exists(target_file)) {
    log_msg(msg = paste0("File ", target_file, " not found in cache"))
    return(NULL)
  } else {
    log_msg(msg = paste0("File ", target_file, " found in cache"))
    return(readRDS(target_file))
  }
}

write_to_cache <- function(object, filename, cache_dir = NULL) {
  if (is.null(cache_dir)) {
    cache_dir <- here("cache")
  }
  if (!fs::dir_exists(cache_dir)) fs::dir_create(cache_dir)
  target_file <- paste0(file.path(cache_dir, filename), ".rds")
  log_msg(msg = paste0("saving ", target_file, " to cache"))
  saveRDS(object, file = target_file)
}


load_and_process_rank_inputs <- function(params) {
  rankfiledir <- params$rankfiledir
  volcanodir <- params$volcanodir
  gct_path <- params$gct_path
  ranks_from <- params$ranks_from
  sample_exclude <- params$sample_exclude %||% NULL
  zscore_emat <- params$zscore_emat %||% TRUE
  zscore_emat_groupby <- ifelse(
    (!is.null(params$zscore_emat_groupby) && !is.na(params$zscore_emat_groupby ) && is.character(params$zscore_emat_groupby)),
    params$zscore_emat_groupby,
    FALSE
  )
  exclude_samples_from_data <- params$advanced$exclude_samples_from_data %||% FALSE
  rankname_order <- params$extra$rankname_order %||% NULL


  log_msg(msg = paste0("ranks from : ", ranks_from))
  log_msg(msg = paste0("rankfiledir : ", rankfiledir))

  if (!is.null(rankfiledir) && file.exists(rankfiledir)) { #
    rnkfiles <- dir_ls(path = rankfiledir, regexp = ".*\\.rnk$", fail = FALSE)
    log_msg(msg = paste0("looking for rank files in ", rankfiledir))


    if (length(rnkfiles) > 0) {
      log_msg(msg = paste0("found ", length(rnkfiles), " rankfiles"))
      rnkdfs <- rnkfiles %>% load_rnkfiles()
      names(rnkdfs) <- normalize_rank_names(rnkfiles)
      return(finalize_rank_inputs(
        rnkdfs,
        rankfiledir = rankfiledir,
        rankname_order = rankname_order,
        context = "cached rank files"
      ))
    } # exit and we're done
    log_msg(msg = "couldn't find any previously saved rnkfiles")
  }
  # ==
  if (ranks_from == "model") {
    if (is.null(gct_path) || !nzchar(gct_path)) {
      stop("gct_path must be provided when ranks_from = 'model'")
    }
    model_specs <- params$models %||% list()
    if (length(model_specs) == 0) {
      model_specs <- list(params$model %||% list())
    }
    replace_outputs <- params$advanced$replace %||% FALSE
    rnkdfs <- list()

    for (idx in seq_along(model_specs)) {
      spec <- model_specs[[idx]]
      model_name <- spec$name %||% paste0("model", idx)
      model_type <- tolower(spec$type %||% "limma")
      model_dir <- file.path(
        params$savedir,
        "model",
        model_type,
        util_tools$safe_path_component(model_name, max_chars = 60)
      )
      log_msg(info = paste0("Generating rank files from model ", model_name, " (", model_type, ")"))
      spec_ranks <- model_tools$create_rnkfiles_from_model(
        gct_path = gct_path,
        model_spec = spec,
        sample_exclude = sample_exclude,
        exclude_samples_from_data = exclude_samples_from_data,
        output_dir = model_dir,
        replace = replace_outputs,
        model_index = idx,
        cache = params$advanced$cache %||% TRUE,
        cache_dir = params$advanced$cachedir %||% NULL
      )
      if (length(spec_ranks) == 0) {
        next
      }
      rnkdfs <- c(rnkdfs, spec_ranks)
    }

    if (length(rnkdfs) == 0) {
      stop("Model specifications yielded no contrasts to export.")
    }

    rnkdfs %>% write_rnkfiles(dir = rankfiledir)
    log_msg(msg = paste0("length of retrieved rankfiles: ", length(rnkdfs)))
    return(finalize_rank_inputs(
      rnkdfs,
      rankfiledir = rankfiledir,
      rankname_order = rankname_order,
      context = "model rank files"
    ))
  }
  if (ranks_from == "volcano") {
    if (is.null(volcanodir) || !file.exists(volcanodir)) {
      stop(paste0("improper volcanodir specification: ", volcanodir))
    }
    log_msg(msg = "saving rankfiles from volcano output. using signedlogp as value")
    rnkdfs <- create_rnkfiles_from_volcano(volcanodir, value_col = "signedlogP")
    rnkdfs %>% write_rnkfiles(dir = rankfiledir) # and save
    names(rnkdfs) <- names(rnkdfs) %>%
      normalize_rank_names()
    log_msg(paste0("length of retrieved rankfiles: ", length(rnkdfs)))
    return(finalize_rank_inputs(
      rnkdfs,
      rankfiledir = rankfiledir,
      rankname_order = rankname_order,
      context = "volcano rank files"
    ))
  }
  if (ranks_from == "gct" && !is.null(gct_path)) {
    apply_z_score <- zscore_emat
    rnkdfs <- create_rnkfiles_from_emat(
      gct_path,
      apply_z_score = apply_z_score,
      zscore_groupby = zscore_emat_groupby,
      sample_exclude = sample_exclude,
      exclude_samples_from_data = exclude_samples_from_data
    )

    names(rnkdfs) <- names(rnkdfs) %>%
      fs::path_file() %>%
      sub("\\.rnk$", "", .)
      #fs::path_ext_remove()
    rnkdfs %>% write_rnkfiles(dir = rankfiledir)
    log_msg(msg = paste0("length of retrieved rankfiles: ", length(rnkdfs)))
    rank_inputs <- finalize_rank_inputs(
      rnkdfs,
      rankfiledir = rankfiledir,
      rankname_order = rankname_order,
      context = "gct rank files"
    )
  }
  # not sure if this level of flow is relevant, refactor later
  if (!exists("rank_inputs")) {
    stop("No rankfiles found, problem loading")
  }
  return(rank_inputs)
}

load_and_process_ranks <- function(params) {
  load_and_process_rank_inputs(params)$ranks
}
