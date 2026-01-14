suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(digest))
# suppressPackageStartupMessages(library(dplyr))

APP_NAME <- "tackle2"

pkg_option_name <- function(suffix) {
  paste0(APP_NAME, "_", suffix)
}

source(file.path(here("R"), "lazyloader.R"))
# listen_tools <- get_tool_env("listen.R")
# util_tools <- get_tool_env("utils")



#
# Helper function to prepend root_dir only to valid paths
prepend_root_dir <- function(params, root_dir) {
  # Function to prepend root_dir if the constructed path exists
  modify_path <- function(x) {
    if (is.character(x)) { # Check if it's a character string
      potential_path <- file.path(root_dir, x)
      if (file.exists(potential_path)) { # Check if the constructed path exists
        return(potential_path) # Use the constructed path
      }
    }
    return(x) # Return the value as is if it's not a valid path
  }

  # Apply the modify_path function to all elements in params
  params <- lapply(params, function(param) {
    if (is.list(param)) {
      return(lapply(param, modify_path)) # Recursively apply for nested lists
    } else {
      return(modify_path(param)) # Apply to non-list items
    }
  })

  return(params)
}


is_absolute_path <- function(path) {
  grepl(paste0("^", normalizePath(root_dir)), normalizePath(path))
}

clean_args <- function(params, root_dir = "/") {
  # would be best for root_dir to be explicitly specified, which it is elsewhere

  # this doesn't work-
  # params <- prepend_root_dir(params, root_dir)

  if (is.null(params$savedir)) {
    params$savedir <- file.path("./plots")
  }
  params$savedir <- file.path(root_dir, (params$savedir %||% file.path("./plots")))

  # all top level params
  params$advanced <- params$advanced %||% list()
  params$barplot <- params$barplot %||% list()
  params$bubbleplot <- params$bubbleplot %||% list()
  params$heatmap_gsea <- params$heatmap_gsea %||% list()
  params$heatmap_gene <- params$heatmap_gene %||% list()
  params$enplot <- params$enplot %||% list()
  params$db <- params$db %||% list()

  params$barplot$do_individual <- params$barplot$do_individual %||% TRUE
  params$barplot$do_combined <- params$barplot$do_combined %||% TRUE
  params$bubbleplot$do_individual <- params$bubbleplot$do_individual %||% TRUE
  params$bubbleplot$do_combined <- params$bubbleplot$do_combined %||% TRUE
  params$barplot$advanced <- params$barplot$advanced %||% list()
  params$bubbleplot$advanced <- params$bubbleplot$advanced %||% list()
  params$db$enable <- params$db$enable %||% FALSE
  params$db$write_results <- params$db$write_results %||% TRUE
  params$db$write_ranks <- params$db$write_ranks %||% FALSE
  params$db$write_pathways <- params$db$write_pathways %||% TRUE

  default_limit_values <- c(10, 20, 30, 50)
  params$barplot$limit <- normalize_limit_vector(params$barplot$limit, default_limit_values)
  params$bubbleplot$limit <- normalize_limit_vector(params$bubbleplot$limit, params$barplot$limit)
  params$bubbleplot$glyph <- params$bubbleplot$glyph %||% "⁕"

  params$heatmap_gsea$do <- params$heatmap_gsea$do %||% TRUE
  params$heatmap_gene$do <- params$heatmap_gene$do %||% TRUE
  params$pca$do <- params$pca$do %||% TRUE
  params$pca_gene <- params$pca_gene %||% list()
  params$pca_gene$do <- params$pca_gene$do %||% FALSE
  params$pca_gene$components <- params$pca_gene$components %||% 3
  pca_gene_colors <- params$pca_gene$metadata_color %||% list()
  if (is.character(pca_gene_colors)) {
    pca_gene_colors <- pca_gene_colors[nzchar(pca_gene_colors)]
  }
  if (is.list(pca_gene_colors)) {
    pca_gene_colors <- unlist(pca_gene_colors, use.names = FALSE)
  }
  params$pca_gene$metadata_color <- unique(as.character(pca_gene_colors))
  params$pca_gene$metadata_shape <- params$pca_gene$metadata_shape %||% ""
  params$pca_gene$top_loadings <- params$pca_gene$top_loadings %||% 25
  params$pca_gene$heatmap <- params$pca_gene$heatmap %||% TRUE
  params$pca_gene$labSize <- params$pca_gene$labSize %||% 1.8
  params$pca_gene$pointSize <- params$pca_gene$pointSize %||% 4.0
  params$pca_gene$sizeLoadingsNames <- params$pca_gene$sizeLoadingsNames %||% 1.4
  cluster_rows <- params$pca_gene$cluster_rows
  if (is.null(cluster_rows)) {
    cluster_rows <- TRUE
  }
  cluster_rows <- as.logical(cluster_rows)
  cluster_rows <- cluster_rows[!is.na(cluster_rows)]
  if (length(cluster_rows) == 0) {
    cluster_rows <- TRUE
  }
  params$pca_gene$cluster_rows <- unique(cluster_rows)

  cluster_columns <- params$pca_gene$cluster_columns
  if (is.null(cluster_columns)) {
    cluster_columns <- c(FALSE, TRUE)
  }
  cluster_columns <- as.logical(cluster_columns)
  cluster_columns <- cluster_columns[!is.na(cluster_columns)]
  if (length(cluster_columns) == 0) {
    cluster_columns <- c(FALSE, TRUE)
  }
  params$pca_gene$cluster_columns <- unique(cluster_columns)

  cut_by <- params$pca_gene$cut_by %||% NULL
  if (is.character(cut_by)) {
    cut_by <- cut_by[nzchar(trimws(cut_by))]
    if (length(cut_by) > 0) {
      cut_by <- cut_by[[1]]
    } else {
      cut_by <- NULL
    }
  } else if (is.logical(cut_by) && !isTRUE(cut_by)) {
    cut_by <- NULL
  }
  params$pca_gene$cut_by <- cut_by

  params$umap_gene <- params$umap_gene %||% list()
  params$umap_gene$do <- params$umap_gene$do %||% FALSE
  params$umap_gene$width <- params$umap_gene$width %||% 7.2
  params$umap_gene$height <- params$umap_gene$height %||% 6.4
  params$umap_gene$n_neighbors <- params$umap_gene$n_neighbors %||% 15
  params$umap_gene$min_dist <- params$umap_gene$min_dist %||% 0.1
  params$umap_gene$metric <- params$umap_gene$metric %||% "euclidean"
  params$umap_gene$seed <- params$umap_gene$seed %||% 42
  params$umap_gene$scale <- params$umap_gene$scale %||% TRUE
  umap_colors <- params$umap_gene$metadata_color %||% list()
  if (is.character(umap_colors)) {
    umap_colors <- umap_colors[nzchar(umap_colors)]
  }
  if (is.list(umap_colors)) {
    umap_colors <- unlist(umap_colors, use.names = FALSE)
  }
  params$umap_gene$metadata_color <- unique(as.character(umap_colors))
  params$umap_gene$metadata_shape <- params$umap_gene$metadata_shape %||% ""
  point_type <- params$umap_gene$point_type %||% "gene"
  if (length(point_type) == 0 || is.na(point_type)) {
    point_type <- "gene"
  }
  point_type <- tolower(as.character(point_type[[1]]))
  if (!point_type %in% c("gene", "sample")) {
    log_msg(warning = paste0("params.umap_gene.point_type '", point_type, "' not recognised; defaulting to 'gene'."))
    point_type <- "gene"
  }
  params$umap_gene$point_type <- point_type
  rank_name <- params$umap_gene$rank_name %||% ""
  if (is.list(rank_name)) {
    rank_name <- unlist(rank_name, use.names = FALSE)
  }
  if (length(rank_name) > 0) {
    rank_name <- trimws(as.character(rank_name[[1]]))
  } else {
    rank_name <- ""
  }
  params$umap_gene$rank_name <- rank_name
  variants_raw <- params$umap_gene$variants %||% list()
  if (!is.list(variants_raw)) {
    variants_raw <- list(variants_raw)
  }
  normalized_variants <- list()
  if (length(variants_raw) > 0) {
    for (idx in seq_along(variants_raw)) {
      candidate <- variants_raw[[idx]]
      if (!is.list(candidate)) {
        next
      }
      variant <- candidate
      colors <- variant$metadata_color %||% list()
      if (is.character(colors)) {
        colors <- colors[nzchar(colors)]
      }
      if (is.list(colors)) {
        colors <- unlist(colors, use.names = FALSE)
      }
      variant$metadata_color <- unique(as.character(colors))
      variant$metadata_shape <- variant$metadata_shape %||% ""
      variant_point_type <- variant$point_type %||% point_type
      if (length(variant_point_type) == 0 || is.na(variant_point_type)) {
        variant_point_type <- point_type
      }
      variant_point_type <- tolower(as.character(variant_point_type[[1]]))
      if (!variant_point_type %in% c("gene", "sample")) {
        variant_point_type <- point_type
      }
      variant$point_type <- variant_point_type
      variant_rank <- variant$rank_name %||% rank_name
      if (is.list(variant_rank)) {
        variant_rank <- unlist(variant_rank, use.names = FALSE)
      }
      if (length(variant_rank) > 0) {
        variant_rank <- trimws(as.character(variant_rank[[1]]))
      } else {
        variant_rank <- ""
      }
      variant$rank_name <- variant_rank
      normalized_variants[[length(normalized_variants) + 1]] <- variant
    }
  }
  params$umap_gene$variants <- normalized_variants


  if (!is.null(params$volcanodir)) {
    params$volcanodir <- file.path(root_dir, params$volcanodir)
    if (!file.exists(params$volcanodir)) stop(paste0("volcanodir does not exist: ", params$volcanodir))
  }

  if (!is.null(params$gct_path) && params$gct_path != "") {
    params$gct_path <- file.path(root_dir, params$gct_path)
  } else {
    params$gct_path <- NULL
  }

  if (!is.null(params$gct_path) && !file.exists(params$gct_path)) {
    stop(paste0(params$gct_path, " does not exist, exiting.."))
  }

  params$model_file <- params$model_file %||% ""
  global_model_file <- params$model_file %||% ""

  load_model_spec <- function(path) {
    if (!file.exists(path)) {
      stop(paste0("Model file does not exist: ", path))
    }
    model_toml <- RcppTOML::parseTOML(path)
    model_toml$model %||% model_toml$Model %||% model_toml
  }

  # Normalize several ways the config can provide model definitions (single model,
  # list of models, or nothing) into a single list we can iterate.
  raw_models <- list()
  if (!is.null(params$model) && length(params$model) > 0) {
    raw_models <- c(raw_models, list(params$model))
  }
  if (!is.null(params$models) && length(params$models) > 0) {
    if (!is.list(params$models)) {
      stop("params$models must be a list of model definitions")
    }
    raw_models <- c(raw_models, params$models)
  }
  if (length(raw_models) == 0) {
    raw_models <- list(list())
  }

  normalized_models <- vector("list", length(raw_models))
  for (idx in seq_along(raw_models)) {
    spec <- raw_models[[idx]]
    if (is.null(spec) || (is.list(spec) && length(spec) == 0)) {
      spec <- list()
    } else if (!is.list(spec)) {
      if (is.character(spec) && length(spec) == 1) {
        if (nzchar(spec)) {
          spec <- list(model_file = spec)
        } else {
          spec <- list()
        }
      } else {
        stop("Model definitions must be provided as lists or single file paths")
      }
    }

    spec_model_file <- spec$model_file %||% global_model_file
    if (!is.null(spec_model_file) && nzchar(spec_model_file)) {
      model_path <- file.path(root_dir, spec_model_file)
      file_model <- load_model_spec(model_path)
      spec <- modifyList(file_model %||% list(), spec)
      spec$model_file <- model_path
    } else {
      spec$model_file <- NULL
    }

    spec$type <- spec$type %||% "limma"
    spec$design <- spec$design %||% ""
    spec$contrasts <- spec$contrasts %||% list()
    if (is.character(spec$contrasts)) {
      spec$contrasts <- as.list(spec$contrasts)
    }
    spec$name <- spec$name %||% spec$label %||% paste0("model", idx)
    normalized_models[[idx]] <- spec
  }
  # Each entry in raw_models may be a full spec list, a single filename, or empty.
  # Convert them into normalized spec lists with defaults applied.

  params$models <- normalized_models
  params$model <- normalized_models[[1]]
  if (!is.null(global_model_file) && nzchar(global_model_file)) {
    params$model_file <- file.path(root_dir, global_model_file)
  } else {
    params$model_file <- NULL
  }

  params$advanced$cache <- params$advanced$cache %||% TRUE

  # print(params$enplot$combine_by)
  params$enplot$combine_by <- params$enplot$combine_by %||% NULL
  # print(params$enplot$combine_by)

  cachedir <- params$advanced$cachedir %||% file.path(params$savedir, "cache")
  if (!is.null(cachedir)) {
    if (cachedir == "savedir") {
      cachedir <- file.path(params$savedir, "cache")
    }
  }
  params$advanced$cachedir <- cachedir

  db_path <- params$db$path %||% "savedir"
  if (!is.null(db_path)) {
    if (db_path == "savedir" || db_path == "") {
      db_path <- file.path(params$savedir, "gsea_results.sqlite")
    } else {
      db_path <- file.path(root_dir, db_path)
    }
  }
  params$db$path <- db_path


  # this block could be cleaned up
  if (!is.null(params$rankfiledir)) {
    if (params$rankfiledir == "savedir" || params$rankfiledir == "") {
      params$rankfiledir <- file.path(params$savedir, "ranks")
    } else {
      params$rankfiledir <- file.path(root_dir, params$rankfiledir)
    }
  } else {
      params$rankfiledir <- file.path(params$savedir, "ranks")
  }

  params$advanced$pivot_gsea_results <- params$advanced$pivot_gsea_results %||% FALSE
  #
  if (!is.null(params$extra$rankname_order)) {
    params$extra$rankname_order <- as.character(params$extra$rankname_order)
    if (length(params$extra$rankname_order) == 1 && params$extra$rankname_order == "sample_order") {
      log_msg(info = "params.extra.rankname_order references 'sample_order'; using params.extra.sample_order instead")
      params$extra$rankname_order <- params$extra$sample_order
    }
  } else {
    params$extra$rankname_order <- params$extra$sample_order
  }
  # think have to do params$pca <- list() and then can do nested assignment
  # but don't need to do this really anyway
  # if ((params$pca$do %||% FALSE) ==  TRUE){
  #   params$pca$extra <- params$extra
  # }

  if (!is.null(params$extra$sample_order)) {
    params$extra$sample_order <- as.character(params$extra$sample_order)
    if (length(params$extra$sample_order) == 1 && params$extra$sample_order == "rankname_order") {
      log_msg(info = "params.extra.sample_order references 'rankname_order'; using params.extra.rankname_order instead")
      params$extra$sample_order <- params$extra$rankname_order
    }
  } else {
    params$extra$sample_order <- params$extra$rankname_order
  }

  if (!is.null(params$extra$rankname_order) && !is.null(params$extra$sample_order)) {
    if (!identical(params$extra$rankname_order, params$extra$sample_order)) {
      log_msg(warning = "params.extra.rankname_order and params.extra.sample_order differ; continuing with rankname_order as canonical list")
    }
  }

  params$species <- params$species %||% "Homo sapiens"

  params$genesets <- params$genesets %||% list(list(category = "H", subcategory = "", collapse = FALSE))

  params$advanced$quiet <- params$advanced$quiet %||% FALSE

  params$advanced$parallel <- params$advanced$parallel %||% FALSE
  params$advanced$exclude_samples_from_data <- params$advanced$exclude_samples_from_data %||% FALSE

  logfile <- params$advanced$logfile %||% file.path(params$savedir, "run.log")
  # browser()
  loglevel <- params$advanced$loglevel
  options(structure(list(logfile), names = pkg_option_name("log_msg_filename")))
  options(structure(list(loglevel), names = pkg_option_name("loglevel")))
  params$advanced$logfile <- logfile

  print(str(params))

  # Optional user-provided color map file (CSV/TSV) for categorical annotations
  # Recognized keys: params$extra$colormap_file, params$extra$cmap_file,
  #                  params$advanced$colormap_file, params$advanced$cmap_file
  user_cmap_path <- params$extra$colormap_file %||% params$extra$cmap_file %||%
    params$advanced$colormap_file %||% params$advanced$cmap_file
  if (!is.null(user_cmap_path) && nzchar(user_cmap_path)) {
    full_cmap_path <- file.path(root_dir, user_cmap_path)
    if (file.exists(full_cmap_path)) {
      cmap <- load_user_colormap(full_cmap_path)
      if (!is.null(cmap)) {
        options(structure(list(cmap), names = pkg_option_name("user_colormap")))
        log_msg(info = paste0("loaded user colormap from ", full_cmap_path))
      }
    } else {
      log_msg(warning = paste0("colormap file not found: ", full_cmap_path))
    }
  }

  if (identical(params$ranks_from, "model")) {
    if (is.null(params$gct_path)) {
      stop("params$gct_path must be provided when ranks_from = 'model'")
    }
    if (is.null(params$model$design) || !nzchar(params$model$design)) {
      stop("params$model$design must be provided when ranks_from = 'model'")
    }
  }

  return(params)
}


# Helpers for safe filesystem naming ------------------------------------

.sanitize_component <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  candidate <- as.character(x)[1]
  if (is.na(candidate)) {
    return(NA_character_)
  }
  candidate <- trimws(candidate)
  if (!nzchar(candidate)) {
    return(NA_character_)
  }

  # transliterate to ASCII where possible, fall back to underscore substitution
  candidate_ascii <- suppressWarnings(iconv(candidate,
    from = "",
    to = "ASCII//TRANSLIT",
    sub = "_"
  ))
  if (!is.na(candidate_ascii) && nzchar(candidate_ascii)) {
    candidate <- candidate_ascii
  }

  candidate <- gsub("[^A-Za-z0-9._-]+", "_", candidate)
  candidate <- gsub("_+", "_", candidate)
  candidate <- gsub("^_+|_+$", "", candidate)

  if (!nzchar(candidate)) {
    return(NA_character_)
  }

  return(candidate)
}

safe_path_component <- function(x, fallback = "item", max_chars = 60, allow_empty = FALSE) {
  candidate <- .sanitize_component(x)

  if (is.na(candidate) || (!nzchar(candidate) && !allow_empty)) {
    candidate <- fallback
  } else if (allow_empty && (is.na(candidate) || !nzchar(candidate))) {
    candidate <- ""
  }

  if (nzchar(candidate) && nchar(candidate) > max_chars) {
    hash <- substr(digest::digest(candidate, algo = "crc32"), 1, 8)
    keep <- max(1, max_chars - nchar(hash) - 1)
    candidate <- paste0(substr(candidate, 1, keep), "_", hash)
  }

  if (!allow_empty && !nzchar(candidate)) {
    candidate <- fallback
  }

  return(candidate)
}

safe_filename <- function(..., fallback = "file", max_chars = 80, delim = "_") {
  parts <- list(...)
  if (length(parts) == 0) {
    return(fallback)
  }

  sanitized <- vapply(parts, function(part) {
    safe_path_component(part, fallback = "", max_chars = max_chars, allow_empty = TRUE)
  }, character(1), USE.NAMES = FALSE)

  sanitized <- sanitized[nzchar(sanitized)]
  if (length(sanitized) == 0) {
    candidate <- fallback
  } else {
    candidate <- paste(sanitized, collapse = delim)
  }

  if (!nzchar(candidate)) {
    candidate <- fallback
  }

  if (nchar(candidate) > max_chars) {
    hash <- substr(digest::digest(candidate, algo = "crc32"), 1, 8)
    keep <- max(1, max_chars - nchar(hash) - nchar(delim))
    candidate <- paste0(substr(candidate, 1, keep), delim, hash)
  }

  if (!nzchar(candidate)) {
    candidate <- fallback
  }

  return(candidate)
}

# --- Smart label shortening helpers ------------------------------------

# Compute the longest common prefix among a character vector
.longest_common_prefix <- function(strings) {
  strings <- as.character(strings)
  strings <- strings[!is.na(strings)]
  if (length(strings) < 2) return("")
  pref <- strings[[1]]
  for (s in strings[-1]) {
    maxlen <- min(nchar(pref), nchar(s))
    while (maxlen > 0 && substr(pref, 1, maxlen) != substr(s, 1, maxlen)) {
      maxlen <- maxlen - 1
    }
    pref <- substr(pref, 1, maxlen)
    if (maxlen == 0) break
  }
  pref
}

# Reverse a string (ASCII-safe)
.rev_str <- function(s) paste(rev(strsplit(s, NULL, fixed = FALSE)[[1]]), collapse = "")

# Compute the longest common suffix among a character vector
.longest_common_suffix <- function(strings) {
  strings <- as.character(strings)
  strings <- strings[!is.na(strings)]
  if (length(strings) < 2) return("")
  rev_strings <- vapply(strings, .rev_str, character(1))
  rev_pref <- .longest_common_prefix(rev_strings)
  if (!nzchar(rev_pref)) return("")
  .rev_str(rev_pref)
}

# Trim an affix to a sensible token boundary (prefix keeps up to last delimiter; suffix drops from first delimiter)
.trim_prefix_to_boundary <- function(prefix, delim = "[._\\-\\s]") {
  if (!nzchar(prefix)) return("")
  locs <- gregexpr(delim, prefix, perl = TRUE)[[1]]
  if (length(locs) == 1 && locs[1] == -1) return("")
  last_pos <- max(locs[locs > 0])
  if (is.finite(last_pos) && last_pos > 0) substr(prefix, 1, last_pos) else ""
}

.trim_suffix_to_boundary <- function(suffix, delim = "[._\\-\\s]") {
  if (!nzchar(suffix)) return("")
  locs <- gregexpr(delim, suffix, perl = TRUE)[[1]]
  if (length(locs) == 1 && locs[1] == -1) return("")
  first_pos <- min(locs[locs > 0])
  if (is.finite(first_pos) && first_pos > 0) substr(suffix, first_pos, nchar(suffix)) else ""
}

# Strip common leading/trailing tokens from a vector of labels
strip_common_affixes <- function(strings,
                                 min_affix_chars = 4,
                                 min_remaining = 6,
                                 delim = "[._\\-\\s]") {
  strings <- as.character(strings)
  if (length(strings) < 2) {
    return(list(values = strings, prefix = "", suffix = ""))
  }

  lcp_raw <- .longest_common_prefix(strings)
  lcs_raw <- .longest_common_suffix(strings)

  pref <- .trim_prefix_to_boundary(lcp_raw, delim = delim)
  suf  <- .trim_suffix_to_boundary(lcs_raw, delim = delim)

  # discard tiny affixes
  if (nchar(pref) < min_affix_chars) pref <- ""
  if (nchar(suf)  < min_affix_chars) suf  <- ""

  # Helper to apply selected affixes
  apply_affix <- function(x, use_pref = TRUE, use_suf = TRUE) {
    y <- x
    if (use_pref && nzchar(pref)) {
      y <- ifelse(startsWith(y, pref), substr(y, nchar(pref) + 1L, nchar(y)), y)
    }
    if (use_suf && nzchar(suf)) {
      y <- ifelse(endsWith(y, suf), substr(y, 1L, nchar(y) - nchar(suf)), y)
    }
    y
  }

  # Try both, then fall back to only one side if needed
  cand_both <- apply_affix(strings, TRUE, TRUE)
  if (all(nchar(cand_both) >= min_remaining & nzchar(cand_both))) {
    cleaned <- cand_both
    used_pref <- pref
    used_suf  <- suf
  } else {
    cand_pref <- apply_affix(strings, TRUE, FALSE)
    cand_suf  <- apply_affix(strings, FALSE, TRUE)
    if (nzchar(pref) && all(nchar(cand_pref) >= min_remaining & nzchar(cand_pref))) {
      cleaned <- cand_pref
      used_pref <- pref
      used_suf <- ""
    } else if (nzchar(suf) && all(nchar(cand_suf) >= min_remaining & nzchar(cand_suf))) {
      cleaned <- cand_suf
      used_pref <- ""
      used_suf <- suf
    } else {
      cleaned <- strings
      used_pref <- ""
      used_suf <- ""
    }
  }

  list(values = cleaned, prefix = used_pref, suffix = used_suf)
}

# Build a name mapping: original -> cleaned (with attributes noting removed parts)
make_name_map <- function(strings,
                          min_affix_chars = 4,
                          min_remaining = 6,
                          delim = "[._\\-\\s]",
                          allow_stem_stripping = TRUE,
                          min_stem_len = 4) {
  # Allow runtime toggle via option tackle2_name_map_strip_stems
  allow_stem_stripping <- getOption(pkg_option_name("name_map_strip_stems"), allow_stem_stripping)
  strings <- as.character(strings)
  res <- strip_common_affixes(strings,
                              min_affix_chars = min_affix_chars,
                              min_remaining = min_remaining,
                              delim = delim)
  cleaned <- res$values
  names(cleaned) <- strings
  # Log once if anything changed
  if (!identical(strings, cleaned)) {
    msg <- paste0(
      "shortened labels by stripping",
      if (nzchar(res$prefix)) paste0(" prefix '", res$prefix, "'") else "",
      if (nzchar(res$suffix)) paste0(" suffix '", res$suffix, "'") else "",
      "."
    )
    log_msg(info = msg)
  }
  attr(cleaned, "removed_prefix") <- res$prefix
  attr(cleaned, "removed_suffix") <- res$suffix

  # Optional: remove common lowercase stems at the start of tokens like 'cell', 'treat'
  allow_stem_stripping <- FALSE # doesn't work yet
  if (allow_stem_stripping) {
    cleaned2 <- strip_common_token_stems(cleaned, min_stem_len = min_stem_len, delim = delim)
    if (!identical(cleaned2, cleaned)) {
      cleaned <- cleaned2
    }
  }
  cleaned
}

# Remove common stems that appear as lowercase prefixes at the start of tokens across all strings.
# Example tokens affected: 'cellOCI3' -> 'OCI3', 'treatControl' -> 'Control'.
# this is not safe and shouldn't be used at the moment
strip_common_token_stems <- function(strings, min_stem_len = 4, delim = "[._\\-\\s]") {
  strings <- as.character(strings)
  if (length(strings) < 2) return(strings)

  # Tokenize
  split_tokens <- function(s) unlist(strsplit(s, delim, perl = TRUE), use.names = FALSE)
  toks_list <- lapply(strings, split_tokens)

  # Extract candidate stems present in each string: lowercase prefix before Uppercase/Digit
  extract_stems <- function(tokens) {
    stems <- vapply(tokens, function(tok) {
      m <- regexpr(paste0("^([a-z]{", min_stem_len, ",})(?=[A-Z0-9])"), tok, perl = TRUE)
      if (m[1] > 0) substr(tok, m[1], m[1] + attr(m, "match.length") - 1) else ""
    }, character(1))
    unique(stems[nzchar(stems)])
  }
  per_string_stems <- lapply(toks_list, extract_stems)
  if (any(lengths(per_string_stems) == 0)) return(strings)

  common_stems <- Reduce(intersect, per_string_stems)
  # Exclude trivial stems or ones we shouldn't touch explicitly, "vs", "minus", "dv"
  common_stems <- setdiff(common_stems, c("vs", "minus", "dv"))
  if (length(common_stems) == 0) return(strings)

  # Apply removal conservatively
  cleaned <- vapply(seq_along(strings), function(i) {
    tokens <- toks_list[[i]]
    out_tokens <- vapply(tokens, function(tok) {
      replaced <- tok
      for (stem in common_stems) {
        if (startsWith(replaced, stem)) {
          rest <- substr(replaced, nchar(stem) + 1L, nchar(replaced))
          # Only drop the stem if remainder remains meaningful
          if (nzchar(rest)) {
            replaced <- rest
          }
        }
      }
      replaced
    }, character(1))
    paste(out_tokens[out_tokens != ""], collapse = "_")
  }, character(1))

  # Log once
  if (!identical(strings, cleaned)) {
    log_msg(info = paste0("shortened labels by removing common token stems: ", paste(common_stems, collapse = ", "))) }

  names(cleaned) <- names(strings)
  cleaned
}

# this is also considered unsafe and untested 
strip_common_token_stems2 <- function(
  strings,
  min_stem_len      = 4,
  delim             = "[._\\-\\s]",
  min_frac          = 1.0,            # fraction of strings (0..1) that must share a stem/token
  min_remainder_len = 2,              # require >= this many chars after a trim
  protect           = c("vs"),        # tokens/stems to never drop
  drop_identical_tokens = FALSE       # optional; see notes
) {
  strings <- as.character(strings)
  if (length(strings) < 2L) return(strings)

  # --- tokenization ---
  split_tokens <- function(s) unlist(strsplit(s, delim, perl = TRUE), use.names = FALSE)
  toks_list <- lapply(strings, split_tokens)

  # --- helper: longest common head/tail runs of tokens ---
  common_run <- function(lst, from_end = FALSE) {
    if (from_end) lst <- lapply(lst, rev)
    maxk <- min(lengths(lst))
    out  <- character(0)
    for (k in seq_len(maxk)) {
      col <- vapply(lst, `[`, character(1), k)
      if (length(unique(col)) == 1L) out <- c(out, col[1]) else break
    }
    if (from_end) rev(out) else out
  }

  # drop common head/tail runs outright (safe shortening)
  head_run <- common_run(toks_list, FALSE)
  tail_run <- common_run(toks_list, TRUE)
  trim_runs <- function(x) {
    a <- length(head_run); b <- length(tail_run); n <- length(x)
    i <- a + 1L; j <- n - b
    if (i > j) character(0) else x[i:j]
  }
  core <- lapply(toks_list, trim_runs)

  # --- generalized stem extraction ---
  extract_stem <- function(tok) {
    # alpha -> digit (e.g., rep1, KRAS12)
    m1 <- regexpr(paste0("^([[:alpha:]]{", min_stem_len, ",})(?=[[:digit:]])"), tok, perl = TRUE)
    if (m1[1] > 0) return(substr(tok, m1[1], m1[1] + attr(m1, "match.length") - 1))
    # lower -> UPPER hump (e.g., sampleA, mapkErk)
    m2 <- regexpr(paste0("^([[:lower:]]{", min_stem_len, ",})(?=[[:upper:]])"), tok, perl = TRUE)
    if (m2[1] > 0) return(substr(tok, m2[1], m2[1] + attr(m2, "match.length") - 1))
    # UPPER -> digit (e.g., KRAS12, HIF1A)
    m3 <- regexpr(paste0("^([[:upper:]]{", min_stem_len, ",})(?=[[:digit:]])"), tok, perl = TRUE)
    if (m3[1] > 0) return(substr(tok, m3[1], m3[1] + attr(m3, "match.length") - 1))
    ""
  }

  per_string_stems <- lapply(core, function(tokens) {
    unique(Filter(nzchar, vapply(tokens, extract_stem, character(1))))
  })

  # find common stems (robust to 1-row/1-col cases)
  stem_universe <- unique(unlist(per_string_stems, use.names = FALSE))
  common_stems  <- character(0)
  if (length(stem_universe)) {
    present <- vapply(per_string_stems,
                      function(x) stem_universe %in% x,
                      FUN.VALUE = rep_len(FALSE, length(stem_universe)))
    if (is.null(dim(present))) present <- matrix(present, nrow = length(stem_universe))
    freq <- rowMeans(present)
    common_stems <- setdiff(stem_universe[freq >= min_frac], protect)
    # avoid prefix-of-prefix interference
    common_stems <- common_stems[order(-nchar(common_stems))]
  }

  # --- apply stem removals conservatively ---
  cleaned_core <- lapply(core, function(tokens) {
    if (!length(tokens)) return(character(0))
    vapply(tokens, function(tok) {
      replaced <- tok
      if (length(common_stems)) {
        for (stem in common_stems) {
          if (startsWith(replaced, stem)) {
            rest <- substr(replaced, nchar(stem) + 1L, nchar(replaced))
            # only drop if the remainder looks meaningful
            if (nchar(rest) >= min_remainder_len && grepl("[[:alnum:]]", rest, perl = TRUE)) {
              replaced <- rest
              break  # stop after first match to avoid over-trimming
            }
          }
        }
      }
      replaced
    }, character(1))
  })

  # --- optionally drop identical tokens AFTER stem trimming ---
  if (isTRUE(drop_identical_tokens)) {
    universe <- unique(unlist(lapply(cleaned_core, unique), use.names = FALSE))
    universe <- setdiff(universe, protect)
    if (length(universe)) {
      present <- vapply(cleaned_core,
                        function(x) universe %in% unique(x),
                        FUN.VALUE = rep_len(FALSE, length(universe)))
      if (is.null(dim(present))) present <- matrix(present, nrow = length(universe))
      freq <- rowMeans(present)
      common_tok <- universe[freq >= min_frac]
      if (length(common_tok)) {
        cleaned_core <- lapply(cleaned_core, function(x) {
          for (t in common_tok) {
            pos <- match(t, x)
            if (!is.na(pos)) x <- x[-pos]  # drop one occurrence
          }
          x
        })
      }
    }
  }

  # --- reassemble ---
  cleaned <- vapply(cleaned_core, function(tok) paste(tok[tok != ""], collapse = "_"), character(1))

  names(cleaned) <- names(strings)
  cleaned
}


cluster_flag_token <- function(flag, prefix) {
  value <- "NA"
  if (length(flag) > 0) {
    candidate <- flag[[1]]
    if (isTRUE(candidate)) {
      value <- "T"
    } else if (identical(candidate, FALSE)) {
      value <- "F"
    } else {
      logical_candidate <- suppressWarnings(as.logical(candidate))
      if (isTRUE(logical_candidate)) {
        value <- "T"
      } else if (identical(logical_candidate, FALSE)) {
        value <- "F"
      }
    }
  }
  paste0(prefix, value)
}

safe_subdir <- function(base, ..., max_chars = 60) {
  components <- list(...)
  if (length(components) == 0) {
    return(base)
  }

  sanitized <- vapply(components, function(part) {
    safe_path_component(part, max_chars = max_chars)
  }, character(1), USE.NAMES = FALSE)

  do.call(file.path, c(list(base), as.list(sanitized)))
}



normalize_limit_vector <- function(limit, fallback = c(10, 20, 30, 50)) {
  if (is.null(limit)) {
    return(fallback)
  }

  if (is.list(limit)) {
    limit <- unlist(limit, use.names = FALSE)
  }

  if (is.character(limit)) {
    if (length(limit) == 1L) {
      limit_candidate <- tolower(limit)
      if (limit_candidate %in% c("auto", "default", "all", "")) {
        return(fallback)
      }
    }
    limit <- suppressWarnings(as.numeric(limit))
  }

  limit <- suppressWarnings(as.numeric(limit))
  limit <- limit[is.finite(limit) & !is.na(limit) & limit > 0]

  if (length(limit) == 0) {
    return(fallback)
  }

  sort(unique(limit))
}


normalize_sample_exclude <- function(sample_exclude, metadata = NULL, warn = TRUE) {
  if (is.null(sample_exclude) ||
    (is.logical(sample_exclude) && length(sample_exclude) == 1L && isFALSE(sample_exclude))) {
    return(character(0))
  }

  spec <- sample_exclude
  if (is.list(spec)) {
    spec <- unlist(spec, use.names = FALSE)
  }

  if (is.null(spec)) {
    return(character(0))
  }

  resolved <- character(0)

  if (!is.null(metadata)) {
    subset_rows <- tryCatch(
      metadata[spec, , drop = FALSE],
      error = function(e) NULL
    )

    if (!is.null(subset_rows) && nrow(subset_rows) > 0) {
      resolved <- rownames(subset_rows)
      resolved <- resolved[!is.na(resolved)]
      resolved <- resolved[resolved != "NA"]
    }
  }

  char_candidates <- unique(na.omit(trimws(as.character(spec))))
  char_candidates <- char_candidates[nzchar(char_candidates)]
  char_candidates <- setdiff(char_candidates, c("TRUE", "FALSE"))

  if (!is.null(metadata)) {
    metadata_rows <- rownames(metadata)
    missing <- setdiff(char_candidates, metadata_rows)
    if (length(missing) > 0 && warn) {
      log_msg(warning = paste0(
        "sample_exclude entries not present in metadata: ",
        paste(missing, collapse = ", ")
      ))
    }
    resolved <- unique(c(resolved, char_candidates))
  } else {
    resolved <- unique(char_candidates)
  }

  resolved <- resolved[!is.na(resolved)]
  resolved <- resolved[resolved != "NA"]

  resolved
}

filter_results_by_rankname <- function(results_list, sample_exclude = NULL) {
  if (is.null(sample_exclude) || length(sample_exclude) == 0L) {
    return(results_list)
  }

  if (!is.list(results_list)) {
    return(results_list)
  }

  for (collection_name in names(results_list)) {
    collection_results <- results_list[[collection_name]]
    if (!is.list(collection_results) || is.null(names(collection_results))) {
      next
    }

    keep <- setdiff(names(collection_results), sample_exclude)
    removed <- setdiff(names(collection_results), keep)
    if (length(removed) > 0) {
      log_msg(info = paste0(
        "sample_exclude: omitting ",
        paste(removed, collapse = ", "),
        " from collection ", collection_name
      ))
    }
    results_list[[collection_name]] <- collection_results[keep]
  }

  results_list
}

filter_named_list <- function(named_list, sample_exclude = NULL) {
  if (is.null(sample_exclude) || length(sample_exclude) == 0L) {
    return(named_list)
  }

  if (!is.list(named_list) || is.null(names(named_list))) {
    return(named_list)
  }

  keep <- setdiff(names(named_list), sample_exclude)
  removed <- setdiff(names(named_list), keep)
  if (length(removed) > 0) {
    log_msg(info = paste0(
      "sample_exclude: omitting ",
      paste(removed, collapse = ", ")
    ))
  }

  named_list[keep]
}

filter_gsea_results <- function(results_frames, sample_exclude = NULL) {
  if (is.null(sample_exclude) || length(sample_exclude) == 0L) {
    return(results_frames)
  }

  if (!is.list(results_frames)) {
    return(results_frames)
  }

  for (collection_name in names(results_frames)) {
    df <- results_frames[[collection_name]]
    if (!is.data.frame(df) || !"rankname" %in% colnames(df)) {
      next
    }

    before_n <- nrow(df)
    df <- df[!(df$rankname %in% sample_exclude), , drop = FALSE]
    after_n <- nrow(df)
    if (before_n != after_n) {
      log_msg(info = paste0(
        "sample_exclude: removed ", before_n - after_n,
        " rows from collection ", collection_name
      ))
    }
    results_frames[[collection_name]] <- df
  }

  results_frames
}

filter_metadata_samples <- function(metadata, sample_exclude = NULL) {
  if (is.null(metadata) || !is.data.frame(metadata)) {
    return(metadata)
  }

  if (is.null(sample_exclude) || length(sample_exclude) == 0L) {
    return(metadata)
  }

  keep <- setdiff(rownames(metadata), sample_exclude)
  metadata[keep, , drop = FALSE]
}



infer_ordered_factor <- function(column) { # this probably doesn't do what you want it to do
  # Function to infer ordering for a given column vector

  # Extract numeric values from strings
  numeric_values <- as.numeric(gsub("[^0-9.-]+", "", column))

  # Find the minimum numeric value
  min_num <- min(numeric_values, na.rm = TRUE)
  if (is.infinite(min_num)) {
    min_num <- 0 # Default to 0 if no numeric values are found
  }

  # Initialize order values with numeric values
  order_values <- numeric_values

  # Indices of non-numeric values
  non_numeric_indices <- which(is.na(order_values))

  if (length(non_numeric_indices) > 0) {
    non_numeric_values <- tolower(column[non_numeric_indices])

    # Assign special order value for "veh" or "vehicle"
    is_vehicle <- grepl("^veh$|^vehicle|^ctrl|^dmso|^untreat$", non_numeric_values, ignore.case = TRUE)
    order_values[non_numeric_indices[is_vehicle]] <- min_num - 2 # Highest priority

    # Assign next priority to other non-numeric values
    order_values[non_numeric_indices[!is_vehicle]] <- min_num - 1
  }

  # Create a data frame for sorting
  df <- data.frame(
    original_value = column,
    order_value = order_values,
    stringsAsFactors = FALSE
  )

  # Remove duplicates while preserving order
  df_unique <- df[!duplicated(df$original_value), ]

  # Sort the data frame by order_value
  df_ordered <- df_unique[order(df_unique$order_value), ]

  # Create an ordered factor with levels in the sorted order
  ordered_factor <- factor(
    column,
    levels = df_ordered$original_value,
    ordered = TRUE
  )

  return(ordered_factor)
}


myzscore <- function(value, minval = NA, remask = TRUE) {
  mask <- is.na(value)
  if (is.na(minval)) minval <- min(value, na.rm = TRUE)

  if (minval == Inf) {
    minval <- 0
  }

  value[is.na(value)] <- minval
  # todo make smaller than min val
  out <- scale(value)

  # if all NA:
  if (sum(!is.finite(out)) == length(out)) {
    out[, 1] <- 0
  }

  if (remask == TRUE) {
    out[mask] <- NA
  }
  # coerce to vector
  out <- as.vector(out)

  return(out)
}

dist_no_na <- function(mat) {
  mat[is.na(mat)] <- min(mat, na.rm = TRUE)
  edist <- dist(mat)
  return(edist)
}

scale_gct <- function(gct, group_by = NULL) {
  log_msg(msg = "zscoring gct file by row")
  log_msg(msg = paste0("group_by is set to: ", group_by))
  if (!is.null(group_by) && length(group_by) == 1 && group_by == FALSE) group_by <- NULL
  group_cols <- group_by # this is a hack to get around the fact that group_by is a dplyr function
  res <- gct %>%
    melt_gct() # %>%
  # Conditionally add group_by
  if (!is.null(group_cols)) {
    if (length(group_cols) == 1) {
      res <- group_by(res, id.x, !!sym(group_cols))
    }
    if (length(group_cols) > 1) {
      res <- group_by(res, id.x, across(all_of(group_cols)))
    }
  } else {
    res <- group_by(res, id.x)
  }

  res <- res %>%
    dplyr::mutate(zscore = myzscore(value)) %>%
    dplyr::ungroup()

  # make a new gct and return
  res <- res %>%
    dplyr::select(id.x, id.y, zscore) %>%
    tidyr::pivot_wider(names_from = id.y, values_from = zscore) %>%
    as.data.frame()
  rownames(res) <- res$id.x
  res$id.x <- NULL
  res <- as.matrix(res)
  rdesc <- gct@rdesc
  cdesc <- gct@cdesc

  if (length(colnames(gct@rdesc)) == 1) {
    rdesc["dummy"] <- "X"
  }

  if (length(colnames(res)) == 1) {
    rdesc["dummy"] <- "X"
  }

  newgct <- new("GCT",
    mat = res,
    rid = rownames(res),
    cid = colnames(res),
    rdesc = rdesc[rownames(res), ],
    cdesc = cdesc[colnames(res), ]
  )

  return(newgct)
}


# # this is exploratory rewrite of plot_utils::make_partial
# get_arg <- function(f, arg, default = "") {
#   args <- if (!is.null(attr(f, "preset_args"))) attr(f, "preset_args") else list()
#   if (arg %in% names(args)) {
#     return(args[[arg]])
#   }
#   return(default)
# }

# another version, maybe easier to read
# rewrite of plot_utils::make_partial
get_arg <- function(f, arg, default = "") {
  args <- attr(f, "preset_args") # will return NULL if no attr
  if (!is.null(args) && !is.null(args[[arg]])) {
    return(args[[arg]])
  }
  return(default)
}

# another version, maybe easier to read
# rewrite of plot_utils::make_partial
get_args <- function(f, ...) {
  args <- attr(f, "preset_args") # will return NULL if no attr
  if (!is.null(args)) {
    return(args)
  }
  return(list())
}

# Revised make_partial using an environment for cleaner argument handling
make_partial <- function(.f, ...) {
  # Ensure the function is correctly resolved
  if (is.character(.f)) {
    .f <- get(.f, envir = parent.frame())
  }

  # print(.f)
  # print(is.function(.f)) # Should be TRUE if .f is correctly resolved
  if (!is.function(.f)) {
    stop("The first argument must be a function")
  }

  # Environment to store arguments

  env <- new.env()
  env$preset_args <- if (!is.null(attr(.f, "preset_args"))) attr(.f, "preset_args") else list()

  # New fixed arguments
  args_fixed <- list(...)

  # Combine old and new arguments
  if (!is.null(names(args_fixed))) {
    env$preset_args[names(args_fixed)] <- args_fixed
  }

  # Inner function using environment
  inner <- function(...) {
    new_args <- list(...)
    # combined_args <- c(env$preset_args, new_args)
    for (arg in names(new_args)) {
      env$preset_args[[arg]] <- new_args[[arg]]
    }
    combined_args <- env$preset_args
    # combined_args <- modifyList(env$preset_args, list(...))
    do.call(.f, combined_args)
  }

  # Attach environment as an attribute (optional but can be helpful for debugging)
  attr(inner, "preset_args") <- env$preset_args

  return(inner)
}

# logging.getLevelNamesMapping() # py
# {'CRITICAL': 50,
#  'FATAL': 50,
#  'ERROR': 40,
#  'WARN': 30,
#  'WARNING': 30,
#  'INFO': 20,
#  'DEBUG': 10,
#  'NOTSET': 0}

.log_levels <- list(
  "CRITICAL" = 50,
  "ERROR" = 40,
  "WARN" = 30,
  "WARNING" = 30,
  "INFO" = 20,
  "DEBUG" = 10,
  "NOTSET" = 0
)



globalloglevel <- options(pkg_option_name("loglevel"))[[1]] %||% "INFO"

log_msg <- function(msg = NULL, info = NULL, debug = NULL, warning = NULL, warn = NULL, error = NULL, filename = NULL, end = "\n", shell = T, loglevel = loglevel, send_over_socket = TRUE, socket_port = 8765, ...) {
  level <- dplyr::case_when(
    !is.null(msg) ~ "INFO",
    !is.null(info) ~ "INFO",
    !is.null(warning) ~ "WARNING",
    !is.null(warn) ~ "WARNING",
    !is.null(debug) ~ "DEBUG",
    !is.null(error) ~ "ERROR",
    TRUE ~ "INFO"
  )

  is_lvl_too_low <- .log_levels[[level]] < .log_levels[[globalloglevel]]
  if (is_lvl_too_low == TRUE) {
    return()
  }

  msg <- Filter(Negate(is.null), list(msg, info, warning, warn, debug, error))
  if (length(msg) > 0){
      msg <- msg[[1]]
  } else {
      msg <- "??"
  }

  prefix <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), level, ": ")

  maybe_filename <- options(pkg_option_name("log_msg_filename"))[[1]]
  if (!is.null(maybe_filename)) {
    filename <- maybe_filename[[1]]
  }
  if (is.null(filename)) {
    filename <- paste0(APP_NAME, ".log")
  }

  dir_path <- fs::path_dir(filename)
  if (!dir_exists(dir_path)) {
    dir_create(dir_path, recurse = TRUE)
  }
  cat(paste0(prefix, msg, end), file = filename, append = TRUE)
  if (shell) {
    cat(paste0(prefix, msg, end))
  }

  # if (send_to_websocket){
  #   listen_tools$send_to_websocket(msg, port=socket_port)
  # }
}


process_cut_by <- function(cut_by, cdesc) {
  #print("***")
  #print(cut_by)
  # Return NULL immediately if cut_by is NULL
  if (is.null(cut_by)) {
    return(NULL)
  }

  # If cut_by is a single string containing ':', split it into a vector
  if (is.character(cut_by) && length(cut_by) == 1 && grepl(":", cut_by)) {
    cut_by <- strsplit(cut_by, ":")[[1]]
  }

  # Ensure cut_by is now a character vector
  if (!is.character(cut_by)) {
    # warning("cut_by should be a character string or vector.")
    # return(NULL)
    # this is fine
    cut_by <- as.character(cut_by)
  }

  # Check if all elements in cut_by are valid column names
  invalid_cols <- setdiff(cut_by, colnames(cdesc))
  if (length(invalid_cols) > 0) {
    warning(
      "The following cut_by elements are not column names in cdesc: ",
      paste(invalid_cols, collapse = ", ")
    )
    return(NULL)
  }

  # Subset the relevant columns and create the interaction factor
  cut_by_factor <- interaction(cdesc[, cut_by, drop = FALSE], drop = TRUE)

  print("***")
  print(cut_by_factor)

  return(cut_by_factor)
}


# ===================== User colormap support =====================

.is_valid_hex <- function(x) {
  if (is.na(x) || !nzchar(x)) return(FALSE)
  grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{8})$", x)
}

.normalize_color <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  x <- trimws(as.character(x))
  if (.is_valid_hex(x)) return(toupper(x))
  # try to convert named colors to hex
  to_hex <- tryCatch({
    rgb <- grDevices::col2rgb(x)
    sprintf("#%02X%02X%02X", rgb[1, 1], rgb[2, 1], rgb[3, 1])
  }, error = function(e) NA_character_)
  if (.is_valid_hex(to_hex)) return(toupper(to_hex))
  NA_character_
}

# load_user_colormap
# Accepts a CSV/TSV with either:
#  - 2 columns: name,color
#  - 3 columns: column,name,color (applies only to the given metadata column)
# Header names are flexible: column/col/field | name/value/level | color/hex
# Lines starting with '#' are ignored.
load_user_colormap_tabular <- function(path) {
  # Helper to try different separators/headers
  parse_with <- function(sep, header) {
    utils::read.table(
      file = path,
      sep = sep,
      header = header,
      stringsAsFactors = FALSE,
      comment.char = "#",
      quote = "\"",
      strip.white = TRUE,
      fill = TRUE
    )
  }

  df <- NULL
  # Try TSV/CSV with and without headers
  specs <- list(
    list("\t", TRUE), list(",", TRUE),
    list("\t", FALSE), list(",", FALSE)
  )
  for (sp in specs) {
    df <- tryCatch(parse_with(sp[[1]], sp[[2]]), error = function(e) NULL)
    if (!is.null(df) && ncol(df) %in% c(2L, 3L)) break
  }

  if (is.null(df) || !(ncol(df) %in% c(2L, 3L))) {
    log_msg(warning = paste0("failed to parse colormap file (expected 2 or 3 columns): ", path))
    return(NULL)
  }

  # Normalize to lower-case names if present
  cn <- tolower(colnames(df))
  if (ncol(df) == 3L) {
    # Map to column,name,color
    idx_col   <- which(cn %in% c("column", "col", "field"))
    idx_name  <- which(cn %in% c("name", "value", "level"))
    idx_color <- which(cn %in% c("color", "colour", "hex"))
    if (length(idx_col) == 0L || length(idx_name) == 0L || length(idx_color) == 0L) {
      idx_col <- 1L; idx_name <- 2L; idx_color <- 3L
    }
    sel <- c(idx_col[1], idx_name[1], idx_color[1])
    df <- df[, sel, drop = FALSE]
    colnames(df) <- c("column", "name", "color")
  } else {
    # 2-column: name,color
    if (length(cn) >= 2L && all(c("name", "color") %in% cn)) {
      name_idx  <- which(cn == "name")[1]
      color_idx <- which(cn %in% c("color", "hex"))[1]
      df <- df[, c(name_idx, color_idx), drop = FALSE]
    } else if (length(cn) >= 2L && all(c("value", "hex") %in% cn)) {
      df <- df[, c(which(cn == "value")[1], which(cn == "hex")[1]), drop = FALSE]
    } else {
      df <- df[, 1:2, drop = FALSE]
    }
    colnames(df) <- c("name", "color")
  }

  # Clean/validate
  df$name  <- trimws(as.character(df$name))
  df$color <- vapply(df$color, .normalize_color, character(1))
  df <- df[nzchar(df$name) & !is.na(df$color), , drop = FALSE]

  # Build result mapping
  mapping <- list(global = character(0), by_column = list())
  if (ncol(df) == 2L) {
    keep <- !duplicated(df$name)
    mapping$global <- stats::setNames(df$color[keep], df$name[keep])
  } else {
    df$column <- trimws(as.character(df$column))
    df <- df[nzchar(df$column), , drop = FALSE]
    split_list <- split(df, df$column)
    mapping$by_column <- lapply(split_list, function(d) {
      keep <- !duplicated(d$name)
      stats::setNames(d$color[keep], d$name[keep])
    })
  }

  mapping
}

# JSON loader: supports either
# 1) {"global": {name: color, ...}, "by_column": { column: {name: color, ...}, ... }}
# 2) {name: color, ...}  (treated as global)
# 3) [ {column?: string, name: string, color: string}, ... ]
load_user_colormap_json <- function(path) {
  obj <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(obj)) {
    log_msg(warning = paste0("failed to parse colormap JSON: ", path))
    return(NULL)
  }

  mapping <- list(global = character(0), by_column = list())

  # Case 1: object with global/by_column
  if (is.list(obj) && !is.null(names(obj))) {
    nm <- tolower(names(obj))
    if ("global" %in% nm || "by_column" %in% nm) {
      # global
      if ("global" %in% nm) {
        g <- obj[[which(nm == "global")[1]]]
        if (is.list(g) && !is.null(names(g))) {
          names_vec <- names(g)
          cols_vec <- vapply(g, .normalize_color, character(1))
          keep <- nzchar(names_vec) & !is.na(cols_vec)
          mapping$global <- setNames(cols_vec[keep], names_vec[keep])
        }
      }
      # by_column
      if ("by_column" %in% nm) {
        bc <- obj[[which(nm == "by_column")[1]]]
        if (is.list(bc) && !is.null(names(bc))) {
          for (colname in names(bc)) {
            if (!nzchar(colname)) next
            m <- bc[[colname]]
            if (is.list(m) && !is.null(names(m))) {
              nms <- names(m)
              cols <- vapply(m, .normalize_color, character(1))
              keep <- nzchar(nms) & !is.na(cols)
              mapping$by_column[[colname]] <- setNames(cols[keep], nms[keep])
            }
          }
        }
      }
      return(mapping)
    }
    # Case 2: object with name:color pairs (global)
    if (length(obj) > 0 && all(vapply(obj, is.character, logical(1)))) {
      names_vec <- names(obj)
      cols_vec <- vapply(obj, .normalize_color, character(1))
      keep <- nzchar(names_vec) & !is.na(cols_vec)
      mapping$global <- setNames(cols_vec[keep], names_vec[keep])
      return(mapping)
    }
  }

  # Case 3: array of records
  if (is.list(obj) && is.null(names(obj))) {
    entries <- obj
    # each entry expected to be list with name/color and optional column
    for (e in entries) {
      if (!is.list(e)) next
      n <- e[["name"]] %||% e[["value"]] %||% e[["level"]]
      col <- e[["color"]] %||% e[["hex"]]
      colname <- e[["column"]] %||% e[["col"]] %||% e[["field"]]
      n <- if (is.null(n)) NA_character_ else as.character(n)
      col <- .normalize_color(col)
      if (is.null(colname)) {
        if (!is.na(n) && !is.na(col)) {
          if (is.null(mapping$global[[n]])) mapping$global[[n]] <- col
        }
      } else {
        colname <- as.character(colname)
        if (!nzchar(colname)) next
        if (!is.na(n) && !is.na(col)) {
          if (is.null(mapping$by_column[[colname]])) mapping$by_column[[colname]] <- character(0)
          if (is.null(mapping$by_column[[colname]][[n]])) mapping$by_column[[colname]][[n]] <- col
        }
      }
    }
    return(mapping)
  }

  log_msg(warning = paste0("unrecognized JSON structure in colormap: ", path))
  NULL
}

# Wrapper that chooses JSON vs tabular parser
load_user_colormap <- function(path) {
  ext <- tolower(fs::path_ext(path))
  # Prefer JSON if extension suggests it, otherwise sniff first non-space char
  if (ext == "json") {
    return(load_user_colormap_json(path))
  }
  # sniff
  first <- tryCatch({
    con <- file(path, "r"); on.exit(close(con))
    repeat {
      line <- readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) break
      trimmed <- trimws(line)
      if (nzchar(trimmed) && substr(trimmed, 1, 1) != "#") {
        return(if (grepl("^[\\[{]", trimmed)) load_user_colormap_json(path) else load_user_colormap_tabular(path))
      }
    }
    load_user_colormap_tabular(path)
  }, error = function(e) load_user_colormap_tabular(path))
  first
}
