suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(colorspace)

suppressPackageStartupMessages(library(here))

basedir <- file.path(here())
src_dir <- file.path(here("R"))

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)
make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)

PATHWAY_COUNTS_FILENAME <- "pathway_counts.json"
PATHWAY_COUNTS_SCHEMA_VERSION <- 1L

pathway_summary_ranknames <- function(results_by_rank, rank_metadata = NULL) {
  result_names <- names(results_by_rank)
  result_names <- as.character(result_names[!is.na(result_names) & nzchar(result_names)])

  metadata_names <- character(0)
  if (is.data.frame(rank_metadata) && "rankname" %in% colnames(rank_metadata)) {
    metadata_names <- as.character(rank_metadata$rankname)
    metadata_names <- metadata_names[!is.na(metadata_names) & nzchar(metadata_names)]
  }

  unique(c(metadata_names, result_names))
}

pathway_count_statistics <- function(df) {
  empty_counts <- list(
    n_pathways = 0L,
    n_main_pathways = 0L,
    n_main_padj_lt_0_25 = 0L,
    n_main_padj_lt_0_05 = 0L
  )
  if (is.null(df) || !is.data.frame(df) || !"pathway" %in% colnames(df) || nrow(df) == 0) {
    return(empty_counts)
  }

  pathway <- as.character(df$pathway)
  valid_pathway <- !is.na(pathway) & nzchar(pathway)

  mainpathway <- if ("mainpathway" %in% colnames(df)) {
    suppressWarnings(as.logical(df$mainpathway))
  } else {
    rep(TRUE, nrow(df))
  }
  mainpathway[is.na(mainpathway)] <- FALSE

  padj <- if ("padj" %in% colnames(df)) {
    suppressWarnings(as.numeric(df$padj))
  } else {
    rep(NA_real_, nrow(df))
  }

  list(
    n_pathways = as.integer(length(unique(pathway[valid_pathway]))),
    n_main_pathways = as.integer(length(unique(pathway[valid_pathway & mainpathway]))),
    n_main_padj_lt_0_25 = as.integer(length(unique(
      pathway[valid_pathway & mainpathway & !is.na(padj) & padj < 0.25]
    ))),
    n_main_padj_lt_0_05 = as.integer(length(unique(
      pathway[valid_pathway & mainpathway & !is.na(padj) & padj < 0.05]
    )))
  )
}

.pathway_summary_frame <- function(df, rankname) {
  if (is.null(df) || !is.data.frame(df)) {
    return(data.frame(
      pathway = character(),
      padj = numeric(),
      mainpathway = logical(),
      rankname = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- df
  if (!"pathway" %in% colnames(out)) {
    out$pathway <- rep(NA_character_, nrow(out))
  }
  if (!"padj" %in% colnames(out)) {
    out$padj <- rep(NA_real_, nrow(out))
  }
  if (!"mainpathway" %in% colnames(out)) {
    out$mainpathway <- rep(TRUE, nrow(out))
  }
  out$rankname <- rep(rankname, nrow(out))
  out
}

make_pathway_count_summary <- function(
    results_by_rank,
    rank_labels = NULL,
    expected_ranknames = NULL,
    collapse = FALSE,
    combined = FALSE,
    main_pathway_ratio = 0.1,
    generated_at = Sys.time()) {
  if (is.null(results_by_rank)) {
    results_by_rank <- list()
  }
  if (is.null(expected_ranknames)) {
    expected_ranknames <- names(results_by_rank)
  }
  expected_ranknames <- unique(as.character(expected_ranknames))
  expected_ranknames <- expected_ranknames[
    !is.na(expected_ranknames) & nzchar(expected_ranknames)
  ]

  rank_entries <- stats::setNames(vector("list", length(expected_ranknames)), expected_ranknames)
  rank_frames <- stats::setNames(vector("list", length(expected_ranknames)), expected_ranknames)

  for (rankname in expected_ranknames) {
    rank_result <- if (rankname %in% names(results_by_rank)) {
      results_by_rank[[rankname]]
    } else {
      NULL
    }
    label <- unname(rank_labels[rankname])
    if (length(label) == 0 || is.na(label[[1]]) || !nzchar(as.character(label[[1]]))) {
      label <- rankname
    } else {
      label <- as.character(label[[1]])
    }

    rank_entries[[rankname]] <- c(
      list(label = label),
      pathway_count_statistics(rank_result)
    )
    rank_frames[[rankname]] <- .pathway_summary_frame(rank_result, rankname)
  }

  retained_df <- if (length(rank_frames) == 0) {
    .pathway_summary_frame(NULL, "")
  } else if (isTRUE(combined)) {
    combined_df <- dplyr::bind_rows(rank_frames)
    fgsea_tools$filter_plot_candidates(
      combined_df,
      collapse = collapse,
      combined = TRUE,
      main_pathway_ratio = main_pathway_ratio
    )
  } else {
    rank_frames %>%
      purrr::map(~ fgsea_tools$filter_plot_candidates(
        .x,
        collapse = collapse,
        combined = FALSE
      )) %>%
      dplyr::bind_rows()
  }

  retained_pathways <- if ("pathway" %in% colnames(retained_df)) {
    pathway <- as.character(retained_df$pathway)
    unique(pathway[!is.na(pathway) & nzchar(pathway)])
  } else {
    character(0)
  }

  generated_at <- if (inherits(generated_at, "POSIXt")) {
    format(generated_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  } else {
    as.character(generated_at[[1]])
  }

  list(
    statistics = list(
      n_retained_pathways = as.integer(length(retained_pathways)),
      ranks = rank_entries
    ),
    schema_version = PATHWAY_COUNTS_SCHEMA_VERSION,
    generated_at = generated_at
  )
}

write_pathway_count_summary <- function(
    summary,
    path,
    replace = TRUE,
    filename = PATHWAY_COUNTS_FILENAME) {
  target <- file.path(path, filename)
  if (fs::file_exists(target) && !isTRUE(replace)) {
    log_msg(msg = paste0("pathway counts: skipping existing file ", target))
    return(invisible(target))
  }

  fs::dir_create(path, recurse = TRUE)
  jsonlite::write_json(
    summary,
    target,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    na = "null"
  )
  log_msg(msg = paste0("pathway counts: wrote ", target))
  invisible(target)
}

write_pathway_count_sidecars <- function(
    results_by_rank,
    output_root,
    rank_labels = NULL,
    expected_ranknames = NULL,
    collapse = FALSE,
    combined = FALSE,
    main_pathway_ratio = 0.1,
    replace = TRUE) {
  if (is.null(output_root) || length(output_root) == 0 || !nzchar(output_root[[1]])) {
    return(invisible(character(0)))
  }
  if (is.null(expected_ranknames)) {
    expected_ranknames <- names(results_by_rank)
  }
  expected_ranknames <- unique(as.character(expected_ranknames))
  expected_ranknames <- expected_ranknames[
    !is.na(expected_ranknames) & nzchar(expected_ranknames)
  ]

  if (isTRUE(combined)) {
    summary <- make_pathway_count_summary(
      results_by_rank = results_by_rank,
      rank_labels = rank_labels,
      expected_ranknames = expected_ranknames,
      collapse = collapse,
      combined = TRUE,
      main_pathway_ratio = main_pathway_ratio
    )
    return(write_pathway_count_summary(summary, output_root, replace = replace))
  }

  targets <- stats::setNames(character(length(expected_ranknames)), expected_ranknames)
  for (rankname in expected_ranknames) {
    label <- unname(rank_labels[rankname])
    if (length(label) == 0 || is.na(label[[1]]) || !nzchar(as.character(label[[1]]))) {
      label <- rankname
    }
    rank_result <- if (rankname %in% names(results_by_rank)) {
      results_by_rank[[rankname]]
    } else {
      NULL
    }
    one_result <- stats::setNames(list(rank_result), rankname)
    summary <- make_pathway_count_summary(
      results_by_rank = one_result,
      rank_labels = rank_labels,
      expected_ranknames = rankname,
      collapse = collapse,
      combined = FALSE,
      main_pathway_ratio = main_pathway_ratio
    )
    rank_dir <- util_tools$safe_subdir(
      output_root,
      util_tools$safe_path_component(label)
    )
    targets[[rankname]] <- write_pathway_count_summary(
      summary,
      rank_dir,
      replace = replace
    )
  }
  invisible(targets)
}

# moved to utils.R
# # Helper function to get current preset arguments or an empty list if none are set
# get_args <- function(f) {
#   if (!is.null(attr(f, "preset_args"))) {
#     return(attr(f, "preset_args"))
#   } else {
#     return(list()) # Return an empty list for easier handling
#   }
# }
# get_arg <- function(f, arg) {
#   args <- get_args(f)
#   val <- args[[arg]]
#   if (is.null(val)) {
#     return("")
#   }
#   return(val)
# }

# # Custom partial function with dynamic argument handling
# make_partial <- function(.f, ...) {
#   # Retrieve current preset arguments, if any
#   current_args <- get_args(.f)

#   # New fixed arguments
#   args_fixed <- list(...)

#   # Ensure that named arguments are handled properly
#   if (!is.null(names(args_fixed))) {
#     # Overwrite or add new arguments
#     current_args[names(args_fixed)] <- args_fixed
#   }

#   # Inner function to call .f with the correct arguments
#   inner <- function(...) {
#     # Combine fixed arguments with any new ones provided at call time
#     args <- modifyList(current_args, list(...))
#     do.call(.f, args)
#   }

#   # Attach updated preset arguments to the inner function for later inspection
#   attr(inner, "preset_args") <- current_args

#   return(inner)
# }


plot_and_save_unsafe <- function(
    plot_code,
    filename,
    path = file.path(
      basedir,
      "plots/"
    ),
    type = "pdf",
    width = 8,
    height = 6,
    replace = T,
    ...) {
  # Setup: Open the appropriate graphics device

  #log_msg(msg = "plot_and_save")
  if (is.null(filename)){
      log_msg(warning = "filename is null")
  }

  filename <- util_tools$safe_path_component(filename, fallback = "plot")

  full_path <- file.path(path, paste0(filename, ".", type))
  log_msg(msg = paste0("plot_and_save: target ", full_path))

  max_total <- 240
  if (nchar(full_path) > max_total) {
    base_len <- nchar(file.path(path, ""))
    extension_len <- nchar(paste0(".", type))
    allowed <- max_total - base_len - extension_len
    allowed <- max(16, allowed)
    filename <- util_tools$safe_path_component(filename, fallback = "plot", max_chars = allowed)
    full_path <- file.path(path, paste0(filename, ".", type))
  }
  # log_msg(msg = paste0("full_path: ", full_path))

  if (!fs::dir_exists(path)) fs::dir_create(path)

  if (file.exists(full_path) && replace == FALSE) {
    log_msg(msg = paste0("plot_and_save: skipping existing file ", full_path))
    graphics.off() # turn off anything that opened
    return()
  }
  # ??
  on.exit(dev.off(), add = TRUE)

  if (type == "pdf") {
    # pdf(full_path, width = width, height = height)
    cairo_pdf(full_path, width = width, height = height)
  } else if (type == "png") {
    png(full_path, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file type")
  }


  # Execute the plot code (passed as a function)
  h <- plot_code()
  # if (inherits(h, "ggplot")) {  # this may be included already in the plot_code closure and unnecessary and not worth pytting here
  #   print(h)  # Required for ggplot rendering inside functions
  # }

  # Teardown: Close the graphics device
  # dev.off()
  # Ensure device closes even if an error occurs by putting the command here
  # is this the proper way to do this to ensure avoiding error of having too many open devices?
  # check all open devices with dev.list()
  # on.exit(dev.off(), add = TRUE)  moved to above

  # log_msg(msg = paste0("done"))

  log_msg(msg = paste0("plot_and_save: wrote ", full_path))
  return(h)
}

safe_plot_and_save <- function(...) {
  tryCatch(
    plot_and_save_unsafe(...),
    error = function(e) {
      message("\n❗️ Error caught: ", conditionMessage(e))
      message("\nLast traceback:")
      print(rlang::last_trace(drop = FALSE))

      message("\nOpen graphics devices:")
      print(dev.list())

      # Force close them
      if (length(dev.list()) > 0) {
        message("\nClosing all devices...")
        graphics.off()
      }

      # Optional: Write to log file
      log_msg(warning = paste("Plotting failed:", conditionMessage(e)))

      # Optionally re-throw if you want upstream failure
      # stop(e)
    }
  )
}
plot_and_save <- safe_plot_and_save #


# Matplotlib default colors
#matplotlib_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
#                       "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
#                       "#bcbd22", "#17becf")


is_numericish <- function(x) {
  if (is.numeric(x)) return(TRUE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]                       # drop empty strings
  if (length(x) == 0) return(FALSE)
  suppressWarnings(all(!is.na(as.numeric(x))))
}

create_named_color_list <- function(df, columns, c=80) {

  # Initialize an empty list to store the result
  color_list <- list()

  # Iterate over each column
  user_map <- getOption(util_tools$pkg_option_name("user_colormap"), NULL)
  for (col_name in columns) {
    col_data <- df[[col_name]]

    # If numeric-like, build a continuous color mapping;
    # otherwise create a discrete palette with explicit NA color.
    if (is_numericish(col_data)) {
      suppressWarnings({ vals <- as.numeric(as.character(col_data)) })
      vals <- vals[is.finite(vals)]
      n_unique <- length(unique(vals))
      tol <- .Machine$double.eps^0.5
      is_integerish <- length(vals) > 0 && all(abs(vals - round(vals)) < tol)

      # Treat small-cardinality integer-like as categorical; otherwise continuous
      if (n_unique <= 5 && is_integerish) {
        # Discrete palette for small integer sets (e.g., 0/1 groups)
        unique_vals <- sort(unique(df[[col_name]]))
        n_vals <- length(unique_vals)
        colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
        assigned <- setNames(colors_assigned, as.character(unique_vals))
        # Apply user overrides (column-specific first, then global)
        if (!is.null(user_map)) {
          col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
          if (!is.null(col_override) && length(col_override) > 0) {
            matches <- intersect(names(assigned), names(col_override))
            assigned[matches] <- col_override[matches]
          }
          global_override <- tryCatch(user_map$global, error = function(e) NULL)
          if (!is.null(global_override) && length(global_override) > 0) {
            # apply only for values still not overridden
            remaining <- setdiff(names(assigned), names(col_override))
            matches2 <- intersect(remaining, names(global_override))
            assigned[matches2] <- global_override[matches2]
          }
        }
        if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
        if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
        color_list[[col_name]] <- assigned
      } else if (length(vals) >= 2) {
        minv <- min(vals, na.rm = TRUE)
        maxv <- max(vals, na.rm = TRUE)

        if (is.finite(minv) && is.finite(maxv) && minv != maxv) {
          # Diverging if range spans zero; otherwise sequential
          if (minv < 0 && maxv > 0) {
            brks <- c(minv, 0, maxv)
            cols <- c("#2166ac", "#f7f7f7", "#b2182b")
          } else {
            brks <- c(minv, maxv)
            cols <- c("#eff3ff", "#08519c")
          }
          color_list[[col_name]] <- circlize::colorRamp2(brks, cols)
          next
        } else {
          # Fallback for constant range
          unique_vals <- sort(unique(df[[col_name]]))
          n_vals <- length(unique_vals)
          colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
          assigned <- setNames(colors_assigned, as.character(unique_vals))
          if (!is.null(user_map)) {
            col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
            if (!is.null(col_override) && length(col_override) > 0) {
              matches <- intersect(names(assigned), names(col_override))
              assigned[matches] <- col_override[matches]
            }
            global_override <- tryCatch(user_map$global, error = function(e) NULL)
            if (!is.null(global_override) && length(global_override) > 0) {
              remaining <- setdiff(names(assigned), names(col_override))
              matches2 <- intersect(remaining, names(global_override))
              assigned[matches2] <- global_override[matches2]
            }
          }
          if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
          if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
          color_list[[col_name]] <- assigned
        }
      } else {
        # Not enough values to infer; treat as discrete
        unique_vals <- sort(unique(df[[col_name]]))
        n_vals <- length(unique_vals)
        colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
        assigned <- setNames(colors_assigned, as.character(unique_vals))
        if (!is.null(user_map)) {
          col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
          if (!is.null(col_override) && length(col_override) > 0) {
            matches <- intersect(names(assigned), names(col_override))
            assigned[matches] <- col_override[matches]
          }
          global_override <- tryCatch(user_map$global, error = function(e) NULL)
          if (!is.null(global_override) && length(global_override) > 0) {
            remaining <- setdiff(names(assigned), names(col_override))
            matches2 <- intersect(remaining, names(global_override))
            assigned[matches2] <- global_override[matches2]
          }
        }
        if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
        if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
        color_list[[col_name]] <- assigned
      }
    } else {
      # Discrete categorical mapping
      col_data[is.na(col_data)] <- "NA"
      unique_vals <- sort(unique(col_data))
      n_vals <- length(unique_vals)
      colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
      assigned <- setNames(colors_assigned, as.character(unique_vals))
      # Apply user overrides (column-specific first, then global)
      if (!is.null(user_map)) {
        col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
        if (!is.null(col_override) && length(col_override) > 0) {
          matches <- intersect(names(assigned), names(col_override))
          assigned[matches] <- col_override[matches]
        }
        global_override <- tryCatch(user_map$global, error = function(e) NULL)
        if (!is.null(global_override) && length(global_override) > 0) {
          remaining <- setdiff(names(assigned), names(col_override))
          matches2 <- intersect(remaining, names(global_override))
          assigned[matches2] <- global_override[matches2]
        }
      }
      if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
      if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
      color_list[[col_name]] <- assigned
    }
  }

  return(color_list)
}



process_cut_by <- function(cut_by, cdesc) {
  # Return NULL immediately if cut_by is NULL
  if (is.null(cut_by) || length(cut_by) == 0) {
    return(NULL)
  }

  # Treat NA and FALSE as "no cut"
  if (is.logical(cut_by) && length(cut_by) == 1 && isFALSE(cut_by)) {
    return(NULL)
  }
  if (all(is.na(cut_by))) {
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

  # Drop empty/NA values after coercion
  cut_by <- cut_by[!is.na(cut_by)]
  cut_by <- cut_by[nzchar(trimws(cut_by))]
  if (length(cut_by) == 0) {
    return(NULL)
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

  return(cut_by_factor)
}
