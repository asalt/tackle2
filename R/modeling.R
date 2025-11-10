suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(ggrepel))

source(file.path(here("R"), "lazyloader.R"))

util_tools <- get_tool_env("utils")
log_msg <- util_tools$make_partial(util_tools$log_msg)

# Custom MNAR imputation function: normal draw shifted down by 2 SD
# Usage from mice: set method = "mnar_shift_norm" for a variable
mice.impute.mnar_shift_norm <- function(y, ry, x, wy = NULL, ...) {
  imp <- mice::mice.impute.norm(y, ry, x, wy, ...)
  shift <- 2 * stats::sd(y[ry], na.rm = TRUE)
  imp_shifted <- imp - shift
  imp_shifted
}

# Ensure the custom imputation function is discoverable by mice by placing it in the
# global environment as `mice.impute.mnar_shift_norm` when not already defined there.
tryCatch({
  if (!exists("mice.impute.mnar_shift_norm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("mice.impute.mnar_shift_norm", mice.impute.mnar_shift_norm, envir = .GlobalEnv)
  }
}, error = function(e) {
  if (exists("log_msg")) log_msg(warning = paste0("Failed to register custom imputer globally: ", conditionMessage(e)))
})

# Plot overlay histogram/density of observed vs imputed values
plot_imputation_diagnostics <- function(imp, original_metadata, vars, outdir, model_label = NULL, logger = NULL) {
  if (is.null(outdir) || !nzchar(outdir) || is.null(imp)) return(invisible(NULL))
  diag_dir <- fs::path(outdir, "diagnostics", "imputation")
  fs::dir_create(diag_dir, recurse = TRUE)
  m <- tryCatch(imp$m, error = function(e) 1L)
  for (var in vars) {
    if (!var %in% colnames(original_metadata)) next
    orig <- original_metadata[[var]]
    miss_idx <- which(is.na(orig))
    if (length(miss_idx) == 0) next
    observed <- orig[!is.na(orig)]

    # Collect imputed values across all m imputations
    imputed_vals <- c()
    for (k in seq_len(max(1L, m))) {
      completed <- mice::complete(imp, k)
      if (var %in% colnames(completed)) {
        imputed_vals <- c(imputed_vals, completed[miss_idx, var, drop = TRUE])
      }
    }

    # Skip plotting if we have no observed or no imputed values to compare
    if (length(observed) == 0 || length(imputed_vals) == 0) {
      if (!is.null(logger)) logger(warning = paste0("Skipping imputation diagnostic for ", var, ": insufficient data (observed=", length(observed), ", imputed=", length(imputed_vals), ")"))
      next
    }

    fname <- util_tools$safe_filename(model_label %||% "model", var, fallback = var)
    pdf_path <- fs::path(diag_dir, paste0(fname, ".pdf"))

    # Choose numeric vs categorical plotting
    if (is.numeric(observed) || is.integer(observed)) {
      # Drop non-finite values to avoid ggplot errors
      observed_num <- as.numeric(observed)
      imputed_num <- suppressWarnings(as.numeric(imputed_vals))
      observed_num <- observed_num[is.finite(observed_num)]
      imputed_num <- imputed_num[is.finite(imputed_num)]
      if (length(observed_num) == 0 || length(imputed_num) == 0) {
        if (!is.null(logger)) logger(warning = paste0("Skipping numeric imputation diagnostic for ", var, ": non-finite values only"))
        next
      }
      df_obs <- data.frame(value = observed_num, status = "observed")
      df_imp <- data.frame(value = imputed_num, status = "imputed")
      df <- rbind(df_obs, df_imp)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = status, color = status)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ..density..), bins = 30, alpha = 0.35, position = "identity") +
        ggplot2::geom_density(alpha = 0.7) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::labs(
          title = paste0("Imputation diagnostics: ", var, if (!is.null(model_label)) paste0(" (", model_label, ")") else ""),
          x = var,
          y = "Density"
        )
      ggplot2::ggsave(pdf_path, plot = p, width = 7.2, height = 5.0)
      if (!is.null(logger)) logger(info = paste0("Wrote imputation diagnostic plot: ", pdf_path))
    } else {
      df_obs <- data.frame(value = as.character(observed), status = "observed")
      df_imp <- data.frame(value = as.character(imputed_vals), status = "imputed")
      df <- rbind(df_obs, df_imp)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = status)) +
        ggplot2::geom_bar(position = "dodge") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::labs(
          title = paste0("Imputation diagnostics: ", var, if (!is.null(model_label)) paste0(" (", model_label, ")") else ""),
          x = var,
          y = "Count"
        )
      ggplot2::ggsave(pdf_path, plot = p, width = 7.2, height = 5.0)
      if (!is.null(logger)) logger(info = paste0("Wrote imputation diagnostic plot: ", pdf_path))
    }
  }
  invisible(NULL)
}

# Univariate expression imputation via mice's imputation functions, per gene row.
# Applies the given method (e.g., "pmm", "norm", "mnar_shift_norm").
# - expression_mat: numeric matrix [genes x samples]
# - meta_df: data.frame of sample metadata aligned to columns in expression_mat
# Returns an imputed expression matrix (same dimensions) with NAs filled where possible.
impute_expression_with_mice <- function(expression_mat, meta_df, method = "pmm", logger = NULL) {
  if (!is.matrix(expression_mat)) return(expression_mat)
  if (!requireNamespace("mice", quietly = TRUE)) {
    if (!is.null(logger)) logger(warning = "mice not installed; skipping expression imputation")
    return(expression_mat)
  }
  # Identify genes (rows) with any NA
  na_rows <- which(rowSums(is.na(expression_mat)) > 0)
  if (length(na_rows) == 0) return(expression_mat)

  # Build numeric predictor matrix from metadata (one-hot encode factors/characters)
  meta_df2 <- meta_df
  if (".sample_id" %in% colnames(meta_df2)) meta_df2$.sample_id <- NULL
  # Drop columns with all NA (unusable as predictors)
  drop_all_na <- vapply(meta_df2, function(v) all(is.na(v)), logical(1))
  meta_df2 <- meta_df2[, !drop_all_na, drop = FALSE]
  # If empty, fall back to intercept-only predictor
  if (ncol(meta_df2) == 0) {
    X <- matrix(1, nrow = nrow(meta_df), ncol = 1)
  } else {
    # Characters to factor
    for (cn in colnames(meta_df2)) {
      if (is.character(meta_df2[[cn]])) meta_df2[[cn]] <- factor(meta_df2[[cn]])
    }
    # model.matrix will drop rows with NA in predictors; replace remaining NA in predictors with explicit level for factors
    # For numeric predictors with NA, set to column mean (na.rm) so imputation can still proceed
    for (cn in colnames(meta_df2)) {
      v <- meta_df2[[cn]]
      if (is.factor(v)) {
        levels(v) <- c(levels(v), "NA_level")
        v[is.na(v)] <- "NA_level"
        meta_df2[[cn]] <- v
      } else if (is.numeric(v)) {
        mu <- suppressWarnings(mean(v, na.rm = TRUE))
        if (!is.finite(mu)) mu <- 0
        v[is.na(v)] <- mu
        meta_df2[[cn]] <- v
      }
    }
    X <- stats::model.matrix(~ . , data = meta_df2)
  }

  method_fun <- switch(tolower(method),
    pmm = mice::mice.impute.pmm,
    norm = mice::mice.impute.norm,
    mnar_shift_norm = mice.impute.mnar_shift_norm,
    # default fallback
    mice::mice.impute.pmm
  )

  # Iterate over NA rows and impute missing entries using the selected method.
  # We keep predictors constant across genes for speed.
  for (ri in na_rows) {
    y <- expression_mat[ri, ]
    ry <- !is.na(y)
    if (sum(ry) < 2) {
      # Not enough observed values; fill with median or zero
      fill <- suppressWarnings(median(y, na.rm = TRUE))
      if (!is.finite(fill)) fill <- 0
      y[!ry] <- fill
      expression_mat[ri, ] <- y
      next
    }
    # Call the univariate imputer; wy defaults to missing positions
    imp_vals <- tryCatch({
      method_fun(y = as.numeric(y), ry = ry, x = X, wy = NULL)
    }, error = function(e) {
      if (!is.null(logger)) logger(warning = paste0("Expression impute failed for row ", ri, ": ", conditionMessage(e), "; using median"))
      NA_real_
    })
    if (length(imp_vals) != sum(!ry) || any(!is.finite(imp_vals))) {
      # Fallback
      fill <- suppressWarnings(median(y, na.rm = TRUE))
      if (!is.finite(fill)) fill <- 0
      y[!ry] <- fill
    } else {
      y[!ry] <- imp_vals
    }
    expression_mat[ri, ] <- y
  }

  if (!is.null(logger)) logger(info = paste0("Imputed expression matrix rows with NA: ", length(na_rows)))
  expression_mat
}

compute_model_cache_key <- function(
    gct_path,
    design,
    contrasts,
    sample_exclude,
    exclude_samples_from_data,
    model_name,
    model_index,
    model_file = NULL,
    volcano_cutoff = 0.05,
    volcano_top_n = 35,
    imputation = NULL) {
  info <- tryCatch(fs::file_info(gct_path), error = function(e) NULL)
  sample_exclude <- sort(unique(as.character(sample_exclude %||% character(0))))
  digest::digest(
    list(
      version = 1L,
      gct_path = normalizePath(gct_path, winslash = "/", mustWork = FALSE),
      gct_size = if (!is.null(info)) as.numeric(info$size) else NA_real_,
      gct_mtime = if (!is.null(info)) as.numeric(info$modification_time) else NA_real_,
      design = design,
      contrasts = unname(unlist(contrasts)),
      sample_exclude = sample_exclude,
      exclude_samples_from_data = isTRUE(exclude_samples_from_data),
      model_name = model_name,
      model_index = model_index,
      model_file = model_file %||% "",
      volcano_cutoff = volcano_cutoff,
      volcano_top_n = volcano_top_n,
      imputation = imputation %||% list()
    ),
    algo = "xxhash64"
  )
}

load_cached_model_results <- function(cache_dir, cache_key, logger = NULL) {
  if (is.null(cache_dir) || !nzchar(cache_dir) || is.null(cache_key)) {
    return(NULL)
  }
  target <- fs::path(cache_dir, paste0("model_", cache_key, ".rds"))
  if (!fs::file_exists(target)) {
    if (!is.null(logger)) {
      logger(debug = paste0("Model cache miss: ", target))
    }
    return(NULL)
  }
  if (!is.null(logger)) {
    logger(info = paste0("Model cache hit: ", target))
  }
  readRDS(target)
}

write_cached_model_results <- function(cache_dir, cache_key, results, logger = NULL) {
  if (is.null(cache_dir) || !nzchar(cache_dir) || is.null(cache_key)) {
    return(invisible(NULL))
  }
  if (!fs::dir_exists(cache_dir)) {
    fs::dir_create(cache_dir, recurse = TRUE)
  }
  target <- fs::path(cache_dir, paste0("model_", cache_key, ".rds"))
  saveRDS(results, target)
  if (!is.null(logger)) {
    logger(info = paste0("Saved model cache: ", target))
  }
}

prettify_term_label <- function(original_term, sanitized_term) {
  sanitized_term <- sanitized_term %||% ""
  candidate <- original_term %||% ""
  candidate <- gsub("`", "", candidate)
  candidate <- trimws(candidate)
  if (!nzchar(candidate)) {
    candidate <- sanitized_term
  }

  hi_pattern <- "^factor\\(I\\(([A-Za-z0-9._-]+)\\s*>\\s*median\\(([^)]+)\\)\\)\\)(TRUE|1)$"
  hi_match <- regmatches(candidate, regexec(hi_pattern, candidate))[[1]]
  if (length(hi_match) >= 2) {
    gene <- hi_match[2]
    return(paste0(gene, "_hi_indicator"))
  }

  lo_pattern <- "^factor\\(I\\(([A-Za-z0-9._-]+)\\s*>\\s*median\\(([^)]+)\\)\\)\\)(FALSE|0)$"
  lo_match <- regmatches(candidate, regexec(lo_pattern, candidate))[[1]]
  if (length(lo_match) >= 2) {
    gene <- lo_match[2]
    return(paste0(gene, "_lo_indicator"))
  }

  factor_pattern <- "^factor\\(([^)]+)\\)([A-Za-z0-9._-]+)$"
  factor_match <- regmatches(candidate, regexec(factor_pattern, candidate))[[1]]
  if (length(factor_match) >= 3) {
    base <- gsub("\\s+", "_", factor_match[2])
    level <- gsub("\\s+", "_", factor_match[3])
    return(paste(base, level, sep = "_"))
  }

  label <- candidate
  replacements <- c(
    "factor\\(" = "",
    "I\\(" = "",
    "\\)" = "",
    ":" = "_by_",
    "~" = "_tilde_",
    "," = "_",
    "\\+" = "_plus_",
    "-" = "_minus_",
    "\\*" = "_times_",
    ">" = "_gt_",
    "<" = "_lt_",
    "=" = "_eq_"
  )
  for (pattern in names(replacements)) {
    label <- gsub(pattern, replacements[[pattern]], label)
  }
  label <- gsub("\\s+", "_", label)
  label <- gsub("__+", "_", label)
  label <- gsub("^_", "", label)
  label <- gsub("_$", "", label)
  if (!nzchar(label)) {
    label <- sanitized_term
  }
  label
}

append_metadata_to_gct <- function(gct, metadata_values, output_dir, model_label, replace = FALSE, logger = NULL) {
  if (!length(metadata_values) || is.null(output_dir) || !nzchar(output_dir)) {
    return(invisible(NULL))
  }

  sample_ids <- gct@cid
  augmented <- gct
  for (name in names(metadata_values)) {
    values <- metadata_values[[name]]
    aligned <- rep(NA_real_, length(sample_ids))
    names(aligned) <- sample_ids
    if (!is.null(names(values))) {
      overlap <- intersect(names(values), sample_ids)
      aligned[overlap] <- as.numeric(values[overlap])
    } else if (length(values) == length(sample_ids)) {
      aligned <- as.numeric(values)
    } else {
      aligned[] <- as.numeric(values)
    }
    augmented@cdesc[[name]] <- aligned
  }

  metadata_dir <- fs::path(output_dir, "metadata")
  fs::dir_create(metadata_dir, recurse = TRUE)
  filename <- paste0(util_tools$safe_filename(model_label, fallback = "model"), "_annotated.gct")
  out_path <- fs::path(metadata_dir, filename)
  if (!replace && fs::file_exists(out_path)) {
    if (!is.null(logger)) {
      logger(info = paste0("Metadata GCT exists, skipping write: ", out_path))
    }
    return(invisible(out_path))
  }
  cmapR::write_gct(augmented, out_path, appenddim = FALSE)
  if (!is.null(logger)) {
    logger(info = paste0("Wrote metadata-augmented GCT: ", out_path))
  }
  invisible(out_path)
}

persist_model_outputs <- function(
    results,
    output_dir,
    replace,
    model_label,
    logger = NULL,
    sig_cutoff = 0.05,
    label_top_n = 35) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    return(invisible(NULL))
  }

  tables <- attr(results, "tables") %||% list()
  if (!length(tables)) {
    return(invisible(NULL))
  }

  volcano <- attr(results, "volcano") %||% vector("list", length(tables))
  file_stubs <- attr(results, "file_stubs") %||% setNames(names(tables), names(tables))
  aliases <- attr(results, "aliases") %||% setNames(names(tables), names(tables))
  sig_cutoff <- attr(results, "volcano_sig_cutoff") %||% sig_cutoff
  label_top_n <- attr(results, "volcano_top_n") %||% label_top_n

  tables_dir <- fs::path(output_dir, "tables")
  volcano_plot_dir <- fs::path(output_dir, "volcano_plots")

  fs::dir_create(tables_dir, recurse = TRUE)
  fs::dir_create(volcano_plot_dir, recurse = TRUE)

  for (name in names(tables)) {
    stub <- file_stubs[[name]] %||% util_tools$safe_filename(name, fallback = name)
    alias <- aliases[[name]] %||% name
    table <- tables[[name]]
    volcano_df <- volcano[[name]]

    if (!is.data.frame(table)) {
      next
    }
    if (!is.data.frame(volcano_df)) {
      next
    }

    table_path <- fs::path(tables_dir, paste0(stub, ".tsv"))
    if (replace || !fs::file_exists(table_path)) {
      readr::write_tsv(table, table_path)
      if (!is.null(logger)) {
        logger(info = paste0("Wrote limma table: ", table_path))
      }
    }

    volcano_plot <- volcano_df
    sig_vec <- volcano_plot$`adj.P.Val`
    volcano_plot$significant <- !is.na(sig_vec) & sig_vec < sig_cutoff
    volcano_plot$direction <- ifelse(
      volcano_plot$significant,
      ifelse(volcano_plot$logFC >= 0, "up", "down"),
      "ns"
    )
    plot_df <- volcano_plot
    plot_df$.label <- if ("GeneSymbol" %in% names(plot_df)) {
      labs <- as.character(plot_df$GeneSymbol)
      fallback <- as.character(plot_df$GeneID)
      labs[is.na(labs) | labs == ""] <- fallback[is.na(labs) | labs == ""]
      labs
    } else {
      as.character(plot_df$GeneID)
    }

    volcano_pdf <- fs::path(volcano_plot_dir, paste0(stub, ".pdf"))
    if (replace || !fs::file_exists(volcano_pdf)) {
      label_top_n_val <- max(label_top_n, 0)
      label_df <- NULL
      if (label_top_n_val > 0) {
        ord <- order(plot_df$P.Value, na.last = NA)
        ord <- ord[seq_len(min(label_top_n_val, length(ord)))]
        if (length(ord) > 0) {
          label_df <- plot_df[ord, , drop = FALSE]
        }
      }

      plot_title <- paste(model_label, alias, sep = " - ")
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = logFC, y = neg_log10_p, color = direction)
      ) +
        ggplot2::geom_point(
          alpha = 0.7,
          size = 1.3,
          na.rm = TRUE
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed",
          colour = "#bbbbbb"
        ) +
        ggplot2::geom_hline(
          yintercept = -log10(0.05),
          linetype = "dashed",
          colour = "#bbbbbb"
        ) +
        ggplot2::labs(
          title = plot_title,
          x = "log2 Fold Change",
          y = "-log10(P value)"
        ) +
        ggplot2::theme_minimal(base_size = 11)

      p <- p +
        ggplot2::scale_color_manual(
          values = c("up" = "#d73027", "down" = "#1f77b4", "ns" = "#bdbdbd"),
          breaks = c("up", "down", "ns"),
          guide = "none"
        )

      if (!is.null(label_df) && nrow(label_df) > 0) {
        p <- p +
          ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = .label),
            size = 3.0,
            max.overlaps = Inf,
            min.segment.length = 0,
            box.padding = 0.25,
            point.padding = 0.2,
            na.rm = TRUE
          )
      }

      sig_count <- sum(plot_df$significant, na.rm = TRUE)
      sig_total <- nrow(plot_df)
      sig_label <- sprintf(
        "%d / %d genes padj < %.3g",
        sig_count,
        sig_total,
        sig_cutoff
      )
      x_pos <- max(plot_df$logFC, na.rm = TRUE)
      y_pos <- min(plot_df$neg_log10_p, na.rm = TRUE)
      if (is.finite(x_pos) && is.finite(y_pos)) {
        p <- p +
          ggplot2::annotate(
            "text",
            x = x_pos,
            y = y_pos,
            label = sig_label,
            hjust = 1,
            vjust = -0.3,
            size = 3.2
          )
      }

      ggplot2::ggsave(filename = volcano_pdf, plot = p, width = 7.2, height = 5.0)
    }
  }
  invisible(NULL)
}

sanitize_model_terms <- function(terms) {
  if (length(terms) == 0) {
    return(character(0))
  }

  sanitized <- gsub("[^A-Za-z0-9]+", "", terms)
  empty_idx <- which(!nzchar(sanitized))
  if (length(empty_idx) > 0) {
    sanitized[empty_idx] <- paste0("coef", empty_idx)
  }
  make.unique(sanitized, sep = "_")
}

sanitize_predictor_column_name <- function(token, existing_names = character(0)) {
  base <- util_tools$safe_filename("expr", token, fallback = "expr")
  if (!(base %in% existing_names)) {
    return(base)
  }
  candidate <- base
  counter <- 2
  while (candidate %in% existing_names) {
    candidate <- paste0(base, "_", counter)
    counter <- counter + 1
  }
  candidate
}

parse_contrast_specs <- function(contrast_strings) {
  if (length(contrast_strings) == 0) {
    return(list())
  }

  specs <- vector("list", length(contrast_strings))
  for (idx in seq_along(contrast_strings)) {
    spec <- contrast_strings[[idx]]
    parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      label <- trimws(parts[1])
      expression <- trimws(paste(parts[-1], collapse = "="))
    } else {
      label <- ""
      expression <- trimws(spec)
    }

    if (!nzchar(expression)) {
      stop("Contrast definition is empty after parsing: '", spec, "'")
    }

    display <- if (nzchar(label)) label else expression
    alias_candidate <- if (nzchar(label)) label else expression
    file_candidate <- if (nzchar(display)) display else expression

    specs[[idx]] <- list(
      label = display,
      expression = expression,
      alias_candidate = alias_candidate,
      fallback = paste0("contrast_", idx),
      file_candidate = file_candidate
    )
  }

  alias_raw <- vapply(
    specs,
    function(item) {
      candidate <- item$alias_candidate
      if (!nzchar(candidate)) {
        item$fallback
      } else {
        candidate
      }
    },
    character(1)
  )
  alias_unique <- make.unique(make.names(alias_raw), sep = "_")

  stub_candidates <- vapply(
    specs,
    function(item) {
      util_tools$safe_filename(
        if (nzchar(item$file_candidate)) item$file_candidate else item$fallback,
        fallback = item$fallback
      )
    },
    character(1)
  )
  stub_unique <- make.unique(stub_candidates, sep = "_")

  for (idx in seq_along(specs)) {
    specs[[idx]]$alias <- alias_unique[[idx]]
    specs[[idx]]$file_stub <- stub_unique[[idx]]
  }

  specs
}

resolve_expression_predictor <- function(token, gct) {
  rid_matches <- which(gct@rid == token)
  if (length(rid_matches) == 1) {
    values <- gct@mat[rid_matches, , drop = TRUE]
    return(setNames(as.numeric(values), gct@cid))
  }

  for (column_name in colnames(gct@rdesc)) {
    column <- gct@rdesc[[column_name]]
    if (is.null(column)) next
    column <- as.character(column)
    matches <- which(!is.na(column) & column == token)
    if (length(matches) == 1) {
      values <- gct@mat[matches, , drop = TRUE]
      return(setNames(as.numeric(values), gct@cid))
    }
    if (length(matches) > 1) {
      stop("Multiple rows matched predictor '", token, "' in rdesc column '", column_name, "'. Please disambiguate.")
    }
  }

  NULL
}

create_rnkfiles_from_model <- function(
    gct_path,
    model_spec,
    sample_exclude = NULL,
    exclude_samples_from_data = FALSE,
    output_dir = NULL,
    replace = FALSE,
    model_index = 1,
    cache = TRUE,
    cache_dir = NULL) {
  if (is.null(gct_path) || !nzchar(gct_path)) {
    stop("gct_path must be provided when ranks_from='model'")
  }
  if (!file.exists(gct_path)) {
    stop("gct_path '", gct_path, "' does not exist")
  }

  spec <- model_spec %||% list()
  model_name <- spec$name %||% paste0("model", model_index)
  model_label <- model_name
  model_type <- tolower(spec$type %||% "limma")
  if (!identical(model_type, "limma")) {
    stop("Unsupported model type '", spec$type, "'. Currently only 'limma' is supported.")
  }

  design_str <- spec$design %||% ""
  if (!nzchar(design_str)) {
    stop("Model design formula must be provided when ranks_from='model'")
  }

  contrasts <- spec$contrasts %||% list()
  if (is.character(contrasts)) {
    contrasts <- as.list(contrasts)
  }

  # Optional automatic factor-style contrast generation
  factor_spec <- spec$factor_contrasts %||% spec$glm_factor %||% NULL
  factor_var <- NULL
  factor_modes <- NULL
  factor_baseline <- NULL
  if (!is.null(factor_spec)) {
    factor_var <- as.character(factor_spec$variable %||% factor_spec$var %||% factor_spec$factor %||% "")
    if (nzchar(factor_var)) {
      # modes: baseline, rest (one-vs-rest), pairwise
      factor_modes <- tolower(as.character(unlist(factor_spec$modes %||% list("rest")) ))
      factor_baseline <- as.character(factor_spec$baseline %||% factor_spec$reference %||% "")
    } else {
      factor_var <- NULL
    }
  }

  sig_cutoff <- spec$volcano_padj_cutoff %||% 0.05
  label_top_n <- spec$volcano_top_n %||% 35

  # Optional missing value imputation configuration (e.g., via mice)
  impute_spec <- spec$imputation %||% spec$missing_imputation %||% list()
  do_impute <- impute_spec$do %||% FALSE
  impute_m <- as.integer(impute_spec$m %||% 1L)
  impute_engine <- tolower(as.character(impute_spec$engine %||% impute_spec$method %||% "mice"))
  impute_maxit <- as.integer(impute_spec$maxit %||% 5L)
  impute_seed <- impute_spec$seed %||% NULL
  impute_vars <- impute_spec$variables %||% NULL
  impute_defaultMethod <- impute_spec$defaultMethod %||% NULL
  impute_methods <- impute_spec$methods %||% NULL  # mice methods vector or single string
  impute_blocks <- impute_spec$blocks %||% NULL    # optional blocks definition
  impute_apply_to <- tolower(as.character(impute_spec$apply_to %||% "metadata"))
  impute_expr_method <- impute_spec$expression_method %||% impute_spec$expr_method %||% (if (is.character(impute_methods) && length(impute_methods) == 1 && is.null(names(impute_methods))) impute_methods else "pmm")
  # Prepare a light-weight summary of imputation settings for caching
  imputation_cache_info <- list(
    engine = impute_engine,
    m = impute_m,
    maxit = impute_maxit,
    seed = if (!is.null(impute_seed)) as.character(impute_seed) else NULL,
    variables = if (!is.null(impute_vars)) as.character(impute_vars) else NULL,
    defaultMethod = impute_defaultMethod,
    methods = if (is.character(impute_methods)) impute_methods else NULL,
    blocks = if (!is.null(impute_blocks)) "custom" else NULL,
    apply_to = impute_apply_to,
    expression_method = impute_expr_method
  )

  cache_enabled <- isTRUE(cache) && !is.null(cache_dir) && nzchar(cache_dir)
  cache_key <- NULL

  design_formula <- tryCatch(
    stats::as.formula(design_str),
    error = function(e) {
      stop("Failed to parse model formula '", design_str, "': ", e$message)
    }
  )

  gct <- cmapR::parse_gctx(gct_path)
  metadata <- as.data.frame(gct@cdesc, stringsAsFactors = FALSE)
  metadata$.sample_id <- gct@cid
  rownames(metadata) <- gct@cid

  symbol_map <- NULL
  symbol_candidates <- c(
    "pr_gene_symbol",
    "GeneSymbol",
    "gene_symbol",
    "Symbol",
    "symbol"
  )
  for (candidate in symbol_candidates) {
    if (candidate %in% colnames(gct@rdesc)) {
      symbol_map <- gct@rdesc[[candidate]]
      break
    }
  }
  if (!is.null(symbol_map)) {
    symbol_map <- as.character(symbol_map)
    names(symbol_map) <- as.character(gct@rid)
  }

  normalized_exclude <- util_tools$normalize_sample_exclude(sample_exclude, metadata)
  if (length(normalized_exclude) > 0) {
    keep <- setdiff(rownames(metadata), normalized_exclude)
    if (length(keep) == 0) {
      stop("All samples were removed by sample_exclude; cannot fit model.")
    }
    removed <- setdiff(rownames(metadata), keep)
    if (length(removed) > 0) {
      log_msg(info = paste0("Excluding ", length(removed), " samples from model fitting: ", paste(removed, collapse = ", ")))
    }
    metadata <- metadata[keep, , drop = FALSE]
  }

  if (cache_enabled) {
    cache_key <- compute_model_cache_key(
      gct_path = gct_path,
      design = design_str,
      contrasts = contrasts,
      sample_exclude = normalized_exclude,
      exclude_samples_from_data = exclude_samples_from_data,
      model_name = model_name,
      model_index = model_index,
      model_file = spec$model_file %||% NULL,
      volcano_cutoff = sig_cutoff,
      volcano_top_n = label_top_n,
      imputation = imputation_cache_info
    )
    cached <- load_cached_model_results(cache_dir, cache_key, logger = log_msg)
    if (!is.null(cached) && !isTRUE(replace)) {
      cached_metadata <- attr(cached, "metadata_values") %||% list()
      if (length(cached_metadata) > 0) {
        annotated_path <- append_metadata_to_gct(
          gct = gct,
          metadata_values = cached_metadata,
          output_dir = output_dir,
          model_label = model_label,
          replace = FALSE,
          logger = log_msg
        )
        if (!is.null(annotated_path)) {
          attr(cached, "annotated_gct_path") <- annotated_path
        }
      }
      persist_model_outputs(
        cached,
        output_dir,
        replace = FALSE,
        model_label = attr(cached, "model_name") %||% model_label,
        logger = log_msg
      )
      return(cached)
    }
  }

  expression_mat <- gct@mat[, rownames(metadata), drop = FALSE]

  predictor_columns <- list()

  vars <- all.vars(design_formula)
  missing_vars <- setdiff(vars, colnames(metadata))
  if (length(missing_vars) > 0) {
    for (token in missing_vars) {
      predictor <- resolve_expression_predictor(token, gct)
      if (is.null(predictor)) {
        stop("Formula variable '", token, "' not found in metadata or expression matrix.")
      }
      metadata[[token]] <- as.numeric(predictor[rownames(metadata)])
      safe_name <- sanitize_predictor_column_name(
        token,
        existing_names = names(predictor_columns)
      )
      predictor_columns[[safe_name]] <- list(values = predictor, token = token)
    }
  }

  # Optionally perform missing value imputation on covariates used by the model
  # before constructing the model frame. Defaults to the previous behaviour
  # (omit rows with missing values) when no imputation configuration is provided.
  #do_impute <- identical(impute_engine, "mice")

  # Helper: run fitting from a given metadata (possibly imputed), with optional label suffix
  fit_one <- function(current_metadata, label_suffix = NULL) {
    model_frame <- stats::model.frame(
      design_formula,
      current_metadata,
      na.action = stats::na.omit,
      drop.unused.levels = TRUE
    )

    sample_ids <- rownames(model_frame)
    current_metadata <- current_metadata[sample_ids, , drop = FALSE]
    if (length(sample_ids) == 0) {
      stop("No samples remain after constructing the model frame (check missing values).")
    }

    expr_mat <- expression_mat[, sample_ids, drop = FALSE]
    if (!is.numeric(expr_mat)) {
      storage.mode(expr_mat) <- "double"
    }

    design_matrix <- stats::model.matrix(design_formula, model_frame)
    original_terms <- colnames(design_matrix)
    sanitized_terms <- sanitize_model_terms(original_terms)
    colnames(design_matrix) <- sanitized_terms
    term_lookup <- stats::setNames(original_terms, sanitized_terms)

    if (ncol(design_matrix) == 0) {
      stop("Design matrix contains no columns after sanitisation.")
    }

    keep_cols <- vapply(
      seq_len(ncol(design_matrix)),
      function(i) any(abs(design_matrix[, i]) > .Machine$double.eps),
      logical(1)
    )
    if (!all(keep_cols)) {
      dropped <- sanitized_terms[!keep_cols]
      log_msg(warning = paste0("Dropping zero-variance design columns: ", paste(dropped, collapse = ", ")))
      design_matrix <- design_matrix[, keep_cols, drop = FALSE]
      sanitized_terms <- colnames(design_matrix)
      original_terms <- term_lookup[sanitized_terms]
      term_lookup <- stats::setNames(original_terms, sanitized_terms)
    }

    if (isTRUE(spec$print_model_matrix)) {
      preview <- utils::capture.output(print(head(as.data.frame(design_matrix), n = 10)))
      log_msg(debug = paste0("Limma model matrix preview (first 10 rows):\n", paste(preview, collapse = "\n")))
    }

    if (ncol(design_matrix) == 0) {
      stop("All design columns were dropped; check model specification.")
    }

  apply_contrasts <- FALSE
  pretty_alias_map <- stats::setNames(
    vapply(sanitized_terms, function(term) prettify_term_label(term_lookup[[term]], term), character(1)),
    sanitized_terms
  )
  alias_to_term <- stats::setNames(names(pretty_alias_map), pretty_alias_map)

    normalize_contrast_expression <- function(expr) {
      updated <- expr
      for (alias_label in names(alias_to_term)) {
        target <- alias_to_term[[alias_label]]
        pattern <- sprintf("(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])", alias_label)
        updated <- gsub(pattern, target, updated, perl = TRUE)
      }
      updated
    }

  # Build user-specified contrast defs (if any)
  contrast_defs <- list()
  if (length(contrasts) > 0) {
    user_defs <- parse_contrast_specs(contrasts)
    for (idx in seq_along(user_defs)) {
      raw_expr <- user_defs[[idx]]$expression
      converted <- normalize_contrast_expression(raw_expr)
      user_defs[[idx]]$expression <- converted
      if (converted %in% sanitized_terms) {
        user_defs[[idx]]$coef_name <- converted
      } else if (raw_expr %in% alias_to_term) {
        user_defs[[idx]]$coef_name <- alias_to_term[[raw_expr]]
      }
      if (is.null(user_defs[[idx]]$coef_name)) {
        user_defs[[idx]]$coef_name <- converted
      }
    }
    contrast_defs <- c(contrast_defs, user_defs)
  }

  # Auto-generate factor contrasts if requested
  if (!is.null(factor_var)) {
    # Identify columns for this factor in original terms
    factor_cols_idx <- which(
      grepl(paste0("^", factor_var), original_terms) |
        grepl(paste0("^factor\\(", factor_var, "\\)"), original_terms)
    )
    if (length(factor_cols_idx) >= 2) {
      factor_terms_orig <- original_terms[factor_cols_idx]
      factor_terms_sanitized <- sanitized_terms[factor_cols_idx]
      # Extract levels from original names by stripping prefix
      levels_from_orig <- vapply(seq_along(factor_terms_orig), function(i) {
        o <- factor_terms_orig[[i]]
        o <- sub(paste0("^factor\\(", factor_var, "\\)"), "", o)
        o <- sub(paste0("^", factor_var), "", o)
        o
      }, character(1))

      # Helper: safe label for file names
      mk_stub <- function(...) util_tools$safe_filename(..., fallback = paste0(factor_var, "_contrast"))

      # Modes
      modes <- unique(factor_modes %||% character(0))
      if (length(modes) == 0) modes <- c("rest")

      # baseline contrasts
      if ("baseline" %in% modes) {
        base <- factor_baseline
        if (!nzchar(base)) base <- levels_from_orig[[1]]
        base_idx <- match(base, levels_from_orig)
        if (!is.na(base_idx)) {
          base_term <- factor_terms_sanitized[[base_idx]]
          for (i in seq_along(levels_from_orig)) {
            if (i == base_idx) next
            li <- factor_terms_sanitized[[i]]
            lname <- levels_from_orig[[i]]
            expr <- sprintf("%s - %s", li, base_term)
            alias <- paste0(factor_var, "_", lname, "_vs_", base)
            contrast_defs[[length(contrast_defs) + 1]] <- list(
              alias = alias,
              expression = expr,
              file_stub = mk_stub(alias),
              original_term = alias,
              coef_name = NULL
            )
          }
        }
      }

      # one-vs-rest contrasts
      if ("rest" %in% modes || "one-vs-rest" %in% modes) {
        nlev <- length(factor_terms_sanitized)
        for (i in seq_along(factor_terms_sanitized)) {
          li <- factor_terms_sanitized[[i]]
          lname <- levels_from_orig[[i]]
          others <- setdiff(factor_terms_sanitized, li)
          if (length(others) > 0) {
            mean_other <- paste("(", paste(others, collapse = "+"), ")/", length(others), sep = "")
            expr <- sprintf("%s - %s", li, mean_other)
            alias <- paste0(factor_var, "_", lname, "_vs_rest")
            contrast_defs[[length(contrast_defs) + 1]] <- list(
              alias = alias,
              expression = expr,
              file_stub = mk_stub(alias),
              original_term = alias,
              coef_name = NULL
            )
          }
        }
      }

      # pairwise contrasts
      if ("pairwise" %in% modes) {
        k <- length(factor_terms_sanitized)
        if (k >= 2) {
          for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
              li <- factor_terms_sanitized[[i]]; lj <- factor_terms_sanitized[[j]]
              lname_i <- levels_from_orig[[i]]; lname_j <- levels_from_orig[[j]]
              expr <- sprintf("%s - %s", li, lj)
              alias <- paste0(factor_var, "_", lname_i, "_vs_", lname_j)
              contrast_defs[[length(contrast_defs) + 1]] <- list(
                alias = alias,
                expression = expr,
                file_stub = mk_stub(alias),
                original_term = alias,
                coef_name = NULL
              )
            }
          }
        }
      }
    } else {
      log_msg(warning = paste0("factor_contrasts: variable '", factor_var, "' has <2 levels in design; skipping auto contrasts"))
    }
  }

  if (length(contrast_defs) > 0) {
    apply_contrasts <- TRUE
    contrast_args <- setNames(
      lapply(contrast_defs, function(def) def$expression),
      vapply(contrast_defs, function(def) def$alias, character(1))
    )
    contrast_matrix <- do.call(
      limma::makeContrasts,
      c(contrast_args, list(levels = sanitized_terms))
    )
    if (is.null(dim(contrast_matrix))) {
      contrast_matrix <- matrix(
        contrast_matrix,
        ncol = 1,
        dimnames = list(sanitized_terms, names(contrast_args))
      )
    }
  } else {
    is_intercept <- sanitized_terms %in% c("(Intercept)", "Intercept", "XIntercept", "X.Intercept.")
    base_terms <- sanitized_terms[!is_intercept]
    if (length(base_terms) == 0) {
      base_terms <- sanitized_terms
    }
      if (length(base_terms) == 0) {
        stop("Design matrix contains no estimable terms after removing intercept.")
      }
      contrast_defs <- lapply(seq_along(base_terms), function(idx) {
        term <- base_terms[[idx]]
        original_term <- term_lookup[[term]] %||% term
        pretty_alias <- prettify_term_label(original_term, term)
        list(
          alias = pretty_alias,
          expression = term,
          file_stub = util_tools$safe_filename(
            pretty_alias,
            fallback = paste0("coef_", idx)
          ),
          original_term = original_term,
          coef_name = term
        )
      })
      alias_unique <- make.unique(vapply(contrast_defs, function(def) def$alias, character(1)), sep = "_")
      stub_unique <- make.unique(vapply(contrast_defs, function(def) def$file_stub, character(1)), sep = "_")
      for (idx in seq_along(contrast_defs)) {
        contrast_defs[[idx]]$alias <- alias_unique[[idx]]
        contrast_defs[[idx]]$file_stub <- stub_unique[[idx]]
      }
    }

    # Fit limma
    print("Design matrix is :")
    print(design_matrix)
    fit <- limma::lmFit(expr_mat, design_matrix)
    if (apply_contrasts) {
      fit <- limma::contrasts.fit(fit, contrast_matrix)
    }
    fit <- limma::eBayes(fit, trend=TRUE, robust=TRUE)

    t_stats <- as.matrix(fit$t)
    if (is.null(colnames(t_stats))) {
      colnames(t_stats) <- vapply(contrast_defs, function(def) def$alias, character(1))
    }

    # Allow distinct naming/labeling when running multiple imputations
    run_model_name <- model_name
    run_model_label <- model_label
    if (!is.null(label_suffix) && nzchar(label_suffix)) {
      run_model_name <- paste0(model_name, "_", label_suffix)
      run_model_label <- paste0(model_label, " (", label_suffix, ")")
    }

    prefix <- util_tools$safe_filename(run_model_name, fallback = paste0("model", model_index))

    output <- vector("list", length(contrast_defs))
    result_tables <- vector("list", length(contrast_defs))
    volcano_tables <- vector("list", length(contrast_defs))

    names(output) <- vapply(
      contrast_defs,
      function(def) util_tools$safe_filename(prefix, def$file_stub),
      character(1)
    )

    contrast_aliases <- vapply(contrast_defs, function(def) def$alias, character(1))
    contrast_file_stubs <- vapply(contrast_defs, function(def) def$file_stub, character(1))
    names(contrast_aliases) <- names(output)
    names(contrast_file_stubs) <- names(output)

    lookup_term <- function(key) {
      if (!is.character(key) || length(key) == 0) {
        return(NULL)
      }
      candidate <- key[[1]]
      if (!nzchar(candidate) || is.null(names(term_lookup))) {
        return(NULL)
      }
      if (!(candidate %in% names(term_lookup))) {
        return(NULL)
      }
      term_lookup[[candidate]]
    }

    for (idx in seq_along(contrast_defs)) {
      alias <- contrast_defs[[idx]]$alias
      candidates <- c(
        contrast_defs[[idx]]$coef_name,
        contrast_defs[[idx]]$expression,
        alias,
        lookup_term(contrast_defs[[idx]]$coef_name),
        lookup_term(contrast_defs[[idx]]$expression)
      )
      candidates <- unique(Filter(function(x) is.character(x) && nzchar(x), candidates))
      available_cols <- colnames(t_stats)
      hit <- candidates[candidates %in% available_cols]
      if (length(hit) == 0) {
        stop("Contrast '", alias, "' was not found in the fitted model output.")
      }
      coef_name <- hit[[1]]
      output_name <- names(output)[[idx]]
      values <- t_stats[, coef_name, drop = TRUE]
      df <- data.frame(
        id = rownames(t_stats),
        value = as.numeric(values),
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$value), , drop = FALSE]
      output[[output_name]] <- df

      top_table <- limma::topTable(
        fit,
        coef = coef_name,
        number = Inf,
        sort.by = "none"
      )
      top_table$GeneID <- rownames(top_table)
      if (!is.null(symbol_map)) {
        gene_symbols <- symbol_map[as.character(top_table$GeneID)]
        top_table$GeneSymbol <- ifelse(
          is.na(gene_symbols) | gene_symbols == "",
          NA_character_,
          gene_symbols
        )
        top_table <- top_table[, c("GeneID", "GeneSymbol", setdiff(colnames(top_table), c("GeneID", "GeneSymbol")))]
      } else {
        top_table <- top_table[, c("GeneID", setdiff(colnames(top_table), "GeneID"))]
      }
      result_tables[[idx]] <- top_table

      volcano_plot <- top_table
      volcano_plot$signedlogP <- sign(volcano_plot$logFC) * -log10(pmax(volcano_plot$P.Value, .Machine$double.eps))
      volcano_plot$neg_log10_p <- -log10(pmax(volcano_plot$P.Value, .Machine$double.eps))
      volcano_plot$significant <- !is.na(volcano_plot$`adj.P.Val`) & volcano_plot$`adj.P.Val` < sig_cutoff
      volcano_plot$direction <- ifelse(
        volcano_plot$significant,
        ifelse(volcano_plot$logFC >= 0, "up", "down"),
        "ns"
      )
      volcano_tables[[idx]] <- volcano_plot
    }

    updated_metadata <- list()
    annotated_gct_path <- NULL
    if (length(predictor_columns) > 0) {
      for (safe_name in names(predictor_columns)) {
        predictor <- predictor_columns[[safe_name]]$values
        aligned <- rep(NA_real_, length(current_metadata$.sample_id))
        names(aligned) <- current_metadata$.sample_id
        matches <- intersect(names(predictor), names(aligned))
        aligned[matches] <- as.numeric(predictor[matches])
        current_metadata[[safe_name]] <- aligned
        updated_metadata[[safe_name]] <- aligned
      }
      attr(output, "metadata_columns") <- names(updated_metadata)
      attr(output, "metadata_values") <- updated_metadata
      log_msg(info = paste0(
        "Captured predictor metadata columns for model ",
        run_model_name, ": ",
        paste(names(updated_metadata), collapse = ", ")
      ))
      annotated_gct_path <- append_metadata_to_gct(
        gct = gct,
        metadata_values = updated_metadata,
        output_dir = output_dir,
        model_label = run_model_label,
        replace = replace,
        logger = log_msg
      )
    }

    attr(output, "tables") <- setNames(result_tables, names(output))
    attr(output, "volcano") <- setNames(volcano_tables, names(output))
    attr(output, "model_name") <- run_model_name
    attr(output, "model_type") <- model_type
    attr(output, "aliases") <- contrast_aliases
    attr(output, "file_stubs") <- contrast_file_stubs
    attr(output, "volcano_sig_cutoff") <- sig_cutoff
    attr(output, "volcano_top_n") <- label_top_n
    if (!is.null(annotated_gct_path)) {
      attr(output, "annotated_gct_path") <- annotated_gct_path
    }

    persist_model_outputs(
      output,
      output_dir,
      replace = replace,
      model_label = run_model_label,
      logger = log_msg,
      sig_cutoff = sig_cutoff,
      label_top_n = label_top_n
    )

    output
  }

  if (do_impute) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("The 'mice' package is required for imputation but is not installed.")
    }
    # If expression has any NA, widen scope to include expression even if user
    # only requested metadata (to satisfy "any NA in whole dataset" semantics)
    expr_na_total <- sum(is.na(expression_mat))
    impute_apply_to_eff <- impute_apply_to
    if (identical(impute_apply_to, "metadata") && expr_na_total > 0) {
      impute_apply_to_eff <- "both"
      log_msg(info = paste0("Detected NAs in expression matrix (", expr_na_total, "); expanding imputation scope to 'both'."))
    }
    # Determine variables eligible for imputation across the entire metadata
    all_cols <- setdiff(colnames(metadata), c(".sample_id"))
    if (!is.null(impute_vars)) {
      impute_cols <- intersect(as.character(impute_vars), all_cols)
    } else {
      impute_cols <- all_cols
    }
    # keep only atomic (non-list) columns
    if (length(impute_cols) > 0) {
      keep_mask <- vapply(metadata[, impute_cols, drop = FALSE], function(v) is.atomic(v) && !is.list(v), logical(1))
      impute_cols <- impute_cols[keep_mask]
    }

    if (length(impute_cols) == 0) {
      log_msg(info = "mice imputation requested but no eligible metadata columns found; proceeding without imputation.")
      do_impute <- FALSE
    } else {
      # Count missing before
      na_counts <- vapply(impute_cols, function(v) sum(is.na(metadata[[v]])), integer(1))
      total_na <- sum(na_counts)
      impute_cols_with_na <- impute_cols[na_counts > 0]
      if (total_na == 0 && !(impute_apply_to_eff %in% c("both", "expression") && expr_na_total > 0)) {
        log_msg(info = "No missing values detected in metadata; skipping imputation.")
        do_impute <- FALSE
      } else {
        log_msg(info = paste0(
          "Running mice imputation on ", length(impute_cols),
          " metadata columns (", total_na, " total missing values)"
        ))

        # Prepare data for mice: selected metadata; convert characters to factors for mice
        mice_data <- metadata[, impute_cols, drop = FALSE]
        orig_types <- vapply(mice_data, function(x) class(x)[1], character(1))
        char_cols <- names(orig_types)[orig_types == "character"]
        if (length(char_cols) > 0) {
          for (cn in char_cols) mice_data[[cn]] <- factor(mice_data[[cn]])
        }

        # Respect user-provided defaultMethod if supplied
        args <- list(
          data = mice_data,
          m = max(1L, impute_m),
          maxit = max(1L, impute_maxit),
          printFlag = FALSE
        )
        # Optional: methods vector mapping variables to specific imputation functions
        if (!is.null(impute_methods)) {
          cols <- colnames(mice_data)
          if (is.character(impute_methods) && length(impute_methods) == 1 && is.null(names(impute_methods))) {
            # Single method string applied to all columns with missingness
            miss_cols <- cols[vapply(cols, function(cn) any(is.na(mice_data[[cn]])), logical(1))]
            method_vec <- rep("", length(cols))
            names(method_vec) <- cols
            method_vec[miss_cols] <- impute_methods
            args$method <- method_vec
          } else {
            # Named vector or list mapping variable -> method
            method_vec <- rep("", length(cols))
            names(method_vec) <- cols
            if (is.list(impute_methods)) impute_methods <- unlist(impute_methods, use.names = TRUE)
            for (nm in intersect(names(impute_methods), cols)) {
              if (any(is.na(mice_data[[nm]]))) method_vec[[nm]] <- as.character(impute_methods[[nm]])
            }
            args$method <- method_vec
          }
        }
        # Optional: blocks grouping
        if (!is.null(impute_blocks)) {
          args$blocks <- impute_blocks
        }
        if (!is.null(impute_defaultMethod)) {
          args$defaultMethod <- impute_defaultMethod
        }
        if (!is.null(impute_seed)) {
          # do not set global seed; pass to mice if provided
          args$seed <- impute_seed
        }

        imp <- do.call(mice::mice, args)

        # Create overlay density/hist diagnostics for imputed variables
        tryCatch({
          plot_imputation_diagnostics(
            imp = imp,
            original_metadata = metadata,
            vars = impute_cols_with_na,
            outdir = output_dir,
            model_label = model_label,
            logger = log_msg
          )
        }, error = function(e) {
          if (!is.null(log_msg)) log_msg(warning = paste0("Imputation diagnostics failed: ", conditionMessage(e)))
        }, warning = function(w) {
          if (!is.null(log_msg)) log_msg(warning = paste0("Imputation diagnostics warning: ", conditionMessage(w)))
        })

        # Single imputation: fit once on completed data
        if (imp$m <= 1L) {
          completed <- mice::complete(imp, 1)
          # Convert back character columns that were coerced to factor
          for (cn in intersect(char_cols, colnames(completed))) {
            completed[[cn]] <- as.character(completed[[cn]])
          }
          metadata[, colnames(completed)] <- completed
          # Optional expression matrix imputation
          if (impute_apply_to_eff %in% c("both", "expression")) {
            expression_mat <<- impute_expression_with_mice(
              expression_mat,
              meta_df = metadata,
              method = impute_expr_method,
              logger = log_msg
            )
          }
          output <- fit_one(metadata, label_suffix = NULL)
          # attach imputation meta
          attr(output, "imputation") <- imputation_cache_info

          if (cache_enabled && !is.null(cache_key)) {
            write_cached_model_results(cache_dir, cache_key, output, logger = log_msg)
          }
          return(output)
        }

        # Multiple imputations: fit per completed dataset and combine outputs
        all_outputs <- list()
        for (k in seq_len(imp$m)) {
          completed <- mice::complete(imp, k)
          md_k <- metadata
          # Convert back character columns that were coerced to factor
          for (cn in intersect(char_cols, colnames(completed))) {
            completed[[cn]] <- as.character(completed[[cn]])
          }
          md_k[, colnames(completed)] <- completed
          # Optional expression matrix imputation per imputation k
          if (impute_apply_to_eff %in% c("both", "expression")) {
            expression_mat_k <- impute_expression_with_mice(
              expression_mat,
              meta_df = md_k,
              method = impute_expr_method,
              logger = log_msg
            )
            # replace for this run
            expression_mat_backup <- expression_mat
            expression_mat <<- expression_mat_k
          } else {
            expression_mat_backup <- NULL
          }
          suffix <- paste0("imp", k)
          out_k <- fit_one(md_k, label_suffix = suffix)
          # restore expression matrix if modified
          if (!is.null(expression_mat_backup)) expression_mat <<- expression_mat_backup
          attr(out_k, "imputation") <- modifyList(imputation_cache_info, list(imputation_index = k))
          all_outputs <- c(all_outputs, out_k)
        }

        # For caching, store the combined list
        if (cache_enabled && !is.null(cache_key)) {
          write_cached_model_results(cache_dir, cache_key, all_outputs, logger = log_msg)
        }
        return(all_outputs)
      }
    }
  }

  # No imputation requested or needed: fit once on the original data
  output <- fit_one(metadata, label_suffix = NULL)
  attr(output, "imputation") <- list(engine = "none")
  if (cache_enabled && !is.null(cache_key)) {
    write_cached_model_results(cache_dir, cache_key, output, logger = log_msg)
  }
  output
}
