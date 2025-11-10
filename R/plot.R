suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(assertthat))

suppressPackageStartupMessages(library(here))
source(file.path(here("R", "lazyloader.R")))


basedir <- file.path(here())
src_dir <- file.path(here("R"))

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)


plot_utils <- new.env()
source(file.path(src_dir, "./plot_utils.R"), local = plot_utils)

util_tools <- get_tool_env("utils")
make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg

log_msg <- util_tools$make_partial(util_tools$log_msg)

select <- dplyr::select # always

handle_save_func <- function(save_func, path, filename, width = NULL, height = NULL) {
  if (!is.null(save_func)) {
    if (!is.null(path)) {
      save_func <- make_partial(save_func, path = path)
    }
    if (!is.null(filename)) {
      save_func <- make_partial(save_func, filename = filename)
    }
    if (!is.null(width) && !is.null(height)) {
      save_func <- make_partial(save_func, width = width, height = height)
    }
  }
  return(save_func)
}


make_heatmap_fromgct <- function(
    gct,
    row_title = "",
    column_title = "",
    save_func = NULL,
    cluster_rows = T,
    cluster_columns = F,
    color_mapper = NULL,
    sample_order = NULL,
    cut_by = NULL,
    autoscale = TRUE,
    meta_to_include = NULL, # if NULL uses all
    meta_to_exclude = NULL, # if NULL filters nothing out
    colorbar_title = "zscore",
    cluster_column_slices = TRUE,
    sample_exclude = NULL,
    row_annotation = NULL,
    # scale = T
    ...) {
  # gct <- subgct
  # gct@cdesc$treat <-
  #   factor(.gct@cdesc$treat , levels = c("untreated", "carboplatin", "IMT", "carboplatin_IMT"), ordered = T)

  # print(paste0("nrow of gct in inner func: ", nrow(gct@mat)))
  # print(paste0("sample order is ", sample_order))
  # print(paste0("sample order is null? ", is.null(sample_order)))

  if (!is.null(sample_order)) {
    missing_samples <- setdiff(sample_order, gct@cid)
    if (length(missing_samples) > 0) {
      log_msg(warning = paste0(
        "sample_order entries not present in expression matrix: ",
        paste(missing_samples, collapse = ", "),
        ". They will be ignored."
      ))
    }
    sample_order <- intersect(sample_order, gct@cid)
    if (length(sample_order) == 0) {
      log_msg(info = "sample_order did not match any samples; falling back to matrix order")
    } else {
      gct <- cmapR::subset_gct(gct, cid = sample_order)
    }
  }

  if (!is.null(sample_exclude)) {
    keep_cid <- setdiff(gct@cid, sample_exclude)
    removed <- setdiff(gct@cid, keep_cid)
    if (length(removed) > 0) {
      log_msg(debug = paste0(
        "sample_exclude: removing ", paste(removed, collapse = ", "), " from heatmap"
      ))
    }
    if (length(keep_cid) == 0) {
      log_msg(warning = "sample_exclude removed all samples; skipping heatmap")
      return(NULL)
    }
    gct <- cmapR::subset_gct(gct, cid = keep_cid)
  }

  if (ncol(gct@mat) == 0L) {
    log_msg(warning = "heatmap skipped because no columns remain after filtering")
    return(NULL)
  }

  default_meta_to_exclude <- c("recno", "runno", "searchno", "label", "expr_col", "expr_file", "assay", "tube", "Tube", "TubeLabel", "id")
  # If user specifies additional meta_to_exclude, merge with defaults
  if (!is.null(meta_to_exclude)) {
    meta_to_exclude <- union(default_meta_to_exclude, meta_to_exclude)
  } else {
    meta_to_exclude <- default_meta_to_exclude
  }

  # cat("make_heatmap\n")
  ca <- NULL

  cdesc_for_annot <- gct@cdesc %>% select(-any_of(meta_to_exclude))

  # would be better to do this earlier
  # Prepare metadata columns for annotation: keep numeric as numeric (continuous),
  # and coerce non-numeric to ordered factors.
  for (col in colnames(cdesc_for_annot)) {
    .col <- cdesc_for_annot[[col]]
    # Detect numeric-ish columns (including decimals stored as character)
    if (plot_utils$is_numericish(.col)) {
      # Preserve as numeric for continuous annotations; keep NA as NA
      suppressWarnings({ cdesc_for_annot[[col]] <- as.numeric(as.character(.col)) })
    } else {
      # For categorical, replace NA with explicit label for consistent coloring
      .col <- replace(.col, is.na(.col), "NA")
      if (!is.factor(.col)) {
        .col <- util_tools$infer_ordered_factor(.col)
      }
      cdesc_for_annot[[col]] <- .col
    }
  }
  cdesc_for_annot <- cdesc_for_annot %>%
    arrange(across(everything()))
  .mat <- gct@mat[, rownames(cdesc_for_annot)]

  .colors <- plot_utils$create_named_color_list(cdesc_for_annot, colnames(cdesc_for_annot), c = 124)
  # if ("group" %in% colnames(gct@cdesc)) {
  ca <- ComplexHeatmap::columnAnnotation(
    df = cdesc_for_annot,
    col = .colors
  )
  # }



  #.title <- ifelse(linear == TRUE, "iBAQ", "log(iBAQ)")
  #if (is.null(z_score)) z_score <- FALSE
  #if (z_score == TRUE | z_score == "0") .title <- paste0(.title, " zscore")
  #if (!is.null(z_score_by)) .title <- paste0(.title, " by ", z_score_by)
  # {
  #   heatmap_legend_param$title <- paste0("zscore ", heatmap_legend_param$title)
  # }
  # if (!is.null(standard_scale) && standard_scale == TRUE) .title <- paste0(.title, " (standardized)")

  .legend_width <- unit(6, "cm")
  #.cbar_title <- "zscore"
  heatmap_legend_param <- list(
    title = colorbar_title,
    direction = "horizontal",
    just = "bottom",
    legend_width = .legend_width
  )

  # .note <- paste0(description, '\nNES ', sprintf("%.2f", NES), '  pvalue: ', pval)

  if ("rdesc" %in% colnames(gct@rdesc)) {  # clean this
    row_labels <- gct@rdesc$rdesc
  } else if ("genesymbol" %in% colnames(gct@rdesc)) {
    row_labels <- gct@rdesc$genesymbol
  } else if ("gene_symbol" %in% colnames(gct@rdesc)) {
    row_labels <- gct@rdesc$gene_symbol
  } else if ("GeneSymbol" %in% colnames(gct@rdesc)) {
    row_labels <- gct@rdesc$GeneSymbol
  } else {
    row_labels <- gct@rid
  }


  # gct@mat %>% apply( 1, function(x) scale(x, center = T, scale = scale)) %>% t() %>% as.matrix()),

  if (autoscale == TRUE) {
    heatmap_matrix_width <- unit(ncol(.mat) * .2, "in")
    heatmap_matrix_height <- unit(nrow(.mat) * .2, "in")
  } else {
    heatmap_matrix_width <- NULL
    heatmap_matrix_height <- unit(nrow(.mat) * .2, "in")
  }

  cut_by <- plot_utils$process_cut_by(cut_by, cdesc_for_annot)

  # print(paste0("cutby : ", cut_by))

  # if (!is.null(cut_by) && cut_by %in% colnames(gct@cdesc)) {
  #   cut_by <- gct@cdesc[[cut_by]]
  #   cut_by <- factor(cut_by, levels = unique(cut_by))
  # } else {
  #   cut_by <- NULL
  # }

  # print(paste0("nrow of mat: ", nrow(.mat)))

  if (!is.null(row_annotation) && !"HeatmapAnnotation" %in% class(row_annotation)) {
    if (is.data.frame(row_annotation)) {
      row_annotation <- do.call(ComplexHeatmap::rowAnnotation, c(list(df = row_annotation), list(show_annotation_name = FALSE)))
    } else if (is.vector(row_annotation) || is.factor(row_annotation)) {
      row_annotation <- ComplexHeatmap::rowAnnotation(annotation = row_annotation, show_annotation_name = FALSE)
    } else {
      log_msg(warning = "row_annotation must be a ComplexHeatmap::HeatmapAnnotation, data.frame, or vector; ignoring.")
      row_annotation <- NULL
    }
  }

  ht <- ComplexHeatmap::Heatmap(
    .mat,
    width = heatmap_matrix_width,
    height = heatmap_matrix_height,
    # TODO use z_score_withna or other custom func for handling nas when scaling
    # This appears to be done

    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    #

    row_title = row_title,
    row_labels = row_labels,
    column_labels = colnames(.mat), # id should always be here

    row_title_gp = grid::gpar(fontsize = 10.4, fontface="bold", hjust=1, vjust=1), #just="left"),  # remove just='left' to restore default behavior, defualt just = "center"
    column_title_gp = grid::gpar(fontsize = 7.4, fontface="bold", hjust=1.4), #just="left"),  # remove just='left' to restore default behavior, defualt just = "center"
    column_title_rot = 30,

    column_split = cut_by,
    top_annotation = ca,
    heatmap_legend_param = heatmap_legend_param,
    row_names_gp = grid::gpar(fontsize = 7, fontface="bold"),
    column_names_gp = grid::gpar(fontsize = 7, fontface="bold"),
    cluster_column_slices = cluster_column_slices,
    column_names_side = "top",
    left_annotation = row_annotation
  )



  # log_msg(debug = paste0("defining draw func"))
  do_draw <- function() {
    draw(ht,
      heatmap_legend_side = "bottom",
      column_title = column_title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold", just = "left"),

      padding = unit(c(2, 24, 2, 24), "mm"), # top, left, bottom, right
    )
  }

  # log_msg(debug = paste0("save func: ", class(save_func) %>% as.character()))
  # log_msg(debug = paste0("is null save func: ", is.null(save_func)))

  height <- 5 + (nrow(.mat) * .20) + (nrow(cdesc_for_annot) * 0.2)  # add space for row labels
  width <- 8 + (ncol(.mat) * .26)

  if (!is.null(save_func)) {
    # log_msg(debug = paste0("save func attrs before: "))
    # log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))

    save_func <- make_partial(save_func, height = height, width = width)

    # log_msg(debug = paste0("save func attrs after: "))
    # log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))

    save_func(plot_code = do_draw)
  }


  return(ht)
}

plot_heatmap_of_edges <- function(
    gct,
    results_list,
    scale = T,
    scale_by = NULL,
    save_func = NULL,
    cluster_rows = TRUE,
    cluster_columns = c(FALSE, TRUE),
    limit = 20,
    sample_order = NULL,
    to_include = NULL,
    cut_by = NULL,
    pstat_cutoff = 1,
    replace = TRUE,
    combine_by = NULL,
    meta_to_include = NULL,
    meta_to_exclude = NULL,
    parallel = FALSE,
    pathways_of_interest = NULL, # TODO: implement this. this has started somewhere
    sample_exclude = NULL,
    ...) {
  # transform for plotting
  .gct <- if (scale) util_tools$scale_gct(gct, group_by = scale_by) else gct

  colorbar_title <- if (scale_by %in% colnames(.gct@cdesc)) {
    paste0("zscore by ", scale_by)
  } else {
    "zscore"
  }


  if (!is.null(combine_by)) { # check to  ensure is a dataframe with columns facet and rankname
    if (!all(c("facet", "rankname") %in% colnames(combine_by))) {
      log_msg(error = "combine_by must have columns facet and rankname")
      return(NULL)
    }
    # one additional check
    assertthat::assert_that(all(rownames(combine_by) == combine_by$rankname))
  }

  # Process each collection
  process_collection <- function(collection_name, results_list) {
    collection_results <- results_list[[collection_name]]

    # Build a friendly mapping for comparison names to strip shared affixes
    comparison_names <- names(collection_results)
    name_map <- util_tools$make_name_map(comparison_names)

    if (!is.null(combine_by)) {
      ranknames <- names(results_list[[collection_name]])
      facets <- combine_by[ranknames, "facet"] %>% unique()
      aggregated_results <- facets %>%
        purrr::map(~ {
          facet <- .x
          .ranknames <- combine_by[combine_by$facet == facet, "rankname"] # this needs more error checking
          collection_results[.ranknames] %>%
            purrr::imap(~ {
              .x %>% dplyr::mutate(comparison_name = .y)
            }) %>%
            bind_rows() %>%
            group_by(pathway) %>%
            arrange(abs(pval)) %>%
            slice_head(n = 1) %>%
            ungroup()
        }) %>%
        setNames(facets)
      collection_results <- aggregated_results
    }

    purrr::map(names(collection_results), function(comparison_name) {
      comparison_label <- name_map[[comparison_name]] %||% comparison_name
      process_comparison(collection_name, comparison_name, comparison_label, collection_results[[comparison_name]], .gct)
    })
    #
  }





  # Process each row
  process_row <- function(row, collection_name, comparison_name, comparison_label, gct) {
    # expected row names: pathway, pval, padj, NES, leadingEdge
    expected_row_names <- c("pathway", "pval", "padj", "NES", "leadingEdge")
    if (!all(expected_row_names %in% colnames(row))) {
      log_msg(error = paste0("row does not contain all expected columns: ", row))
      return(NULL)
    }

    .id <- row$pathway
    .row_title <- paste(
      paste0("pval: ", round(row$pval, 4)),
      paste0("padj: ", round(row$padj, 4)),
      paste0("NES: ", round(row$NES, 4)),
      sep = "\n"
    )

      subgct <- gct %>% cmapR::subset_gct(row$leadingEdge[[1]])



    param_grid <- expand.grid(
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      cut_by = cut_by,
      stringsAsFactors = FALSE
    )

    #param_grid %>% purrr::pmap(~ {
    get_arg <- util_tools$get_arg
    param_grid %>% purrr::pmap( ~{
      params <- list(...)
      path <- NULL
      newpath <- NULL
      local_save_func <- save_func
      if (!is.null(local_save_func)) {
        path <- get_arg(local_save_func, "path")
        collection_dir <- util_tools$safe_path_component(strip_pathway_leader(collection_name))
        pathway_dir <- util_tools$safe_path_component(strip_pathway_leader(row$pathway), max_chars = 40)
        newpath <- util_tools$safe_subdir(path, collection_dir, "heatmaps_gene", pathway_dir)
      }
      cluster_rows <- params$cluster_rows %||% FALSE
      cluster_columns <- params$cluster_columns %||% FALSE
      cut_by_val <- params$cut_by

      .cut_by_val <- plot_utils$process_cut_by(cut_by_val, subgct@cdesc)
      cut_by_label <- if (!is.null(.cut_by_val)) {
        paste0("cut_", util_tools$safe_path_component(cut_by_val, max_chars = 32))
      } else {
        ""
      }

      filename <- util_tools$safe_filename(
        get_arg(local_save_func, "filename"),
        util_tools$cluster_flag_token(cluster_rows, "rc"),
        util_tools$cluster_flag_token(cluster_columns, "cc"),
        comparison_label,
        cut_by_label,
        paste0(nrow(subgct@mat), "x", ncol(subgct@mat)),
        fallback = "gene_heatmap"
      )
      if (!is.null(local_save_func)) {
        local_save_func <- make_partial(local_save_func, filename = filename, path = newpath)
      }
      filepath <- if (!is.null(newpath)) file.path(newpath, paste0(filename, ".pdf")) else NULL

      if (!replace && !is.null(filepath) && file.exists(filepath)) {
        log_msg(debug = paste0("skipping ", filename, " as it already exists"))
        return(NULL)
      }
      else {
      log_msg(msg = paste0("plotting gene heatmap for ", .id, " ", comparison_name))

      tryCatch(
        make_heatmap_fromgct(
          subgct,
          save_func = local_save_func,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          sample_order = sample_order,
          row_title = .row_title,
          column_title = paste0(.id, "\n", str_wrap(str_replace_all(comparison_label, "_", " "), width=64)),
          cut_by = cut_by_val,
          meta_to_include = meta_to_include,
          meta_to_exclude = meta_to_exclude,
          sample_exclude = sample_exclude,
          colorbar_title = colorbar_title
          # ...
        ),
        error = function(e) log_msg(error = paste("failed to produce heatmap:", e$message))
      )
      }
    })
  }

  # Process each comparison
  process_comparison <- function(collection_name, comparison_name, comparison_label, result, gct) {

    xtra <- NULL
    if (!is.null(pathways_of_interest)){ # an input argument list, or NULL
      xtra <- result %>% dplyr::filter(pathway %in% pathways_of_interest)
    }

    forplot <- result %>% fgsea_tools$select_topn(
      limit = limit,
      to_include = to_include,
      pstat_cutoff = pstat_cutoff
    )

    if (!is.null(xtra)){
      forplot <- bind_rows(forplot, xtra)
    }

    log_msg(info = paste0("plotting heatmaps for ", collection_name, " ", comparison_name))
    log_msg(info = paste0("n forplot: ", forplot %>% nrow()))

    globals_list <- list(
      gct = gct,
      collection_name = collection_name,
      comparison_name = comparison_name,
      forplot = forplot
    )

    map_func <- purrr::map
    if (parallel == TRUE) { # this is unstable
      workers <- future::availableCores() - 2
      options(future.globals.maxSize = 200 * 1024^2)
      log_msg(msg = paste0("using ", workers, " workers"))
      map_func <- function(...) furrr::future_map(..., .options = furrr::furrr_options(globals = globals_list))
    } else {
      map_func <- purrr::map
    }

    map_func(seq_len(nrow(forplot)), function(row_index) {
      tryCatch({
      row = forplot[row_index, ]
      process_row(row, collection_name, comparison_name, comparison_label, gct)
      }, error = function(e) {
      message("Error at row_index=", row_index, ": ", e$message)
        NULL
      })
    })

    # purrr::map(seq_len(nrow(forplot)), function(row_index) {
    #   row = forplot[row_index, ]
    #   process_row(row, collection_name, comparison_name, gct)
    # })

  }

  # Main logic: Iterate through collections
  res <- purrr::map(names(results_list), function(collection_name) {
    process_collection(collection_name, results_list)
  })
  return(res)
}

  # hts <- names(results_list) %>%
  #   purrr::map(
  #     ~ {
  #       collection_name <- .x
  #       hts <- names(results_list[[collection_name]]) %>%
  #         purrr::map(
  #           ~ {
  #             comparison_name <- .x
  #             result <- results_list[[collection_name]][[comparison_name]]
  #             # print(collection_name)
  #             # print(comparison_name)
  #             forplot <- result %>% fgsea_tools$select_topn(
  #               limit = limit,
  #               to_include = to_include, # extra pathways to expilcitly include
  #               pstat_cutoff = pstat_cutoff
  #             )
  #             log_msg(info = paste0("plotting heatmaps for ", collection_name, " ", comparison_name))
  #             log_msg(info = paste0("n forplot: ", forplot %>% nrow()))

  #             # result %>%
  #             # forplot <- result %>%
  #             #   arrange(pval) %>%
  #             #   head(10) %>%
  #             #   mutate(leadingedgelist = stringr::str_split(leadingEdge, ","))
  #             # for (i in seq_len(nrow(forplot))) {
  #             hts <- seq_len(nrow(forplot)) %>% purrr::map(~{
  #               row <- forplot[.x, ]
  #               .id <- row$pathway
  #               .row_title <- row$pval
  #               # .leading_edge <- unlist(lapply(row$leadingedgelist, function(x) gsub('[c\\(\\)" ]', "", x)))
  #               print(paste0("length of leading edge: ", length(row$leadingEdge[[1]])))
  #               print(paste0("head of leading edge: ", head(row$leadingEdge[[1]])))
  #               subgct <- .gct %>% cmapR::subset_gct(row$leadingEdge[[1]])
  #               print(paste0("nrow of subgct: ", nrow(subgct@mat)))

  #               # print(.id)
  #               .row_title <- paste(
  #                 paste0("pval: ", row$pval %>% round(4) %>% as.character()),
  #                 paste0("padj: ", row$padj %>% round(4) %>% as.character()),
  #                 paste0("NES: ", row$NES %>% round(4) %>% as.character()),
  #                 sep = "\n")
  #               log_msg(debug = paste0("plotting gene heatmap for ", .id, " ", comparison_name))

  #               path <- get_arg(save_func, "path")
  #               newpath <- file.path(path, paste0(make.names(collection_name)), make.names("heatmaps_gene"))
  #               filename <- paste0(get_arg(save_func, "filename"), make.names(row$pathway), make.names(comparison_name), nrow(subgct@mat))
  #               save_func <- make_partial(save_func, filename = filename, path = newpath)


  #               ht <- NULL
  #               ht <- tryCatch(
  #                 {
  #                   #
  #                   # save_func <- make_partial(save_func, filename = filename)
  #                   ht <- make_heatmap_fromgct(subgct,
  #                     save_func = save_func,
  #                     cluster_rows = F,
  #                     cluster_columns = F,
  #                     sample_order = sample_order,
  #                     row_title = .row_title,
  #                     column_title = paste0(.id, "\n", comparison_name),
  #                     cut_by = cut_by,
  #                     ...
  #                   )
  #                   # print(ht)
  #                 },
  #                 error = function(e) print(sprintf("failed to produce heatmap, %s", e$message))
  #               )
  #               return(ht)
  #             }) # end of one heatmap
  #           }) # end of all heatmaps for one collection
  #       return(hts)
  #     }
  #   )
  # return(hts)



# do_plot_edge_heatmaps_onecollection <- function(
#     df,
#     cut_by = NULL,
#     limit = 120,
#     pstat_cutoff = 1,
#     pstat_usetype = "padj",
#     main_pathway_ratio = 0.1,
#     cluster_rows = F,
#     cluster_columns = F,
#     ...) {
#   kwargs <- list(...)
#   sampleresults_names <- names(sampleresults_list)
#   if (!"rankname" %in% colnames(df)) {
#     warning("rankname not in colnames df, returning")
#     return(NULL)
#   }
#
#   df <- fgsea_tools$filter_on_mainpathway(df, main_pathway_ratio = main_pathway_ratio)
#   # Limit the number of pathways if necessary
#   top_pathways <- df %>%
#     arrange(-abs(NES)) %>%
#     filter(!!as.symbol(pstat_usetype) < pstat_cutoff) %>%
#     slice_head(n = limit) %>%
#     pull(pathway)
#   df <- df %>% filter(pathway %in% top_pathways)
# }
#
# do_plot_edge_heatmaps_allcollections <- function(results_list, ...) {
#   # input is named list
#   # names are collection names
#   # values are gsea dataframe results, long form
#   # with comparison/sample names inside
#   log_msg(msg = "starting all edge plots")
#   kwargs <- list(...)
#   collections <- names(results_list)
#   # if is null then exit
#   purrr::map(~ {
#     collection <- .x
#     do_plot_edge_heatmaps_onecollection(results_list[[collection]], ...)
#   })
# }






strip_pathway_leader <- function(x) {
  # Accepts either a character vector of pathway names or a data.frame with a 'pathway' column.
  remove_tokens <- function(v) {
    v %>%
      stringr::str_remove("^HALLMARK_") %>%
      stringr::str_remove("^KEGG_") %>%
      stringr::str_remove("^REACTOME_") %>%
      stringr::str_remove("^MEDICUS_") %>%
      stringr::str_remove("^GOMF_") %>%
      stringr::str_remove("^GOBP_") %>%
      stringr::str_remove("^GOCC_")
  }

  if (is.character(x)) {
    return(remove_tokens(x))
  }
  if (is.data.frame(x)) {
    if (!"pathway" %in% colnames(x)) return(x)
    x <- dplyr::mutate(x, pathway = remove_tokens(.data$pathway))
    return(x)
  }
  x
}

# Function: prepare_data_for_barplot
# Description: This function prepares the data for a barplot by taking in a dataframe as input.
# Parameters:
#   - df: The dataframe containing the data for the barplot.
# Returns: None
prepare_data_for_barplot <- function(df) {
  df_renamed <- df %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "HALLMARK_")) %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "KEGG_")) %>%
    # Strip vendor-specific prefixes if present
    dplyr::mutate(pathway = stringr::str_remove_all(pathway, stringr::regex("(?i)MEDICUS_?"))) %>%
    dplyr::mutate(pathway = stringr::str_remove_all(pathway, stringr::regex("(?i)REFERENCE_?"))) %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "GOMF_")) %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "REACTOME_")) %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "GOBP_")) %>%
    dplyr::mutate(pathway = stringr::str_remove(pathway, "GOCC_"))
  df <- df_renamed

  # is across necessary?
  #   df_renamed <- df %>%
  #   mutate(across(starts_with("pathway"), ~str_remove(., "HALLMARK_"))) %>%
  #   mutate(across(starts_with("pathway"), ~str_remove(., "KEGG_"))) %>%
  #   mutate(across(starts_with("pathway"), ~str_remove(., "GOMF_"))) %>%
  #   mutate(across(starts_with("pathway"), ~str_remove(., "REACTOME_")))
  # df <- df_renamed


  sel <- df %>%
    arrange(-abs(NES)) %>%
    arrange(-NES) %>%
    mutate(pathway = str_replace_all(pathway, "_", " ") %>% str_wrap(width = 52)) %>%
    mutate(pathway = factor(pathway, levels = unique(pathway), ordered = T)) %>%
    arrange(pathway) # %>%

  sel <- sel %>%
    rowwise() %>%
    mutate(leadingEdgeNum = length(leadingEdge)) %>%
    mutate(leadingEdgeFrac = paste0(leadingEdgeNum, "/", size)) %>%
    ungroup() %>%
    mutate(outline_val = dplyr::if_else(padj < .05, "black", NA))

  return(sel)
}


#' Barplot with Numbers
#'
#' This function creates a barplot with numbers using the provided dataframe.
#'
#' @param df dataframe with columns pathway, NES, padj, leadingEdge, size
#' @param title title of the plot
#' @param save_func function to save the plot
#' @param ... additional arguments to be passed to the function
#'
#' @return a ggplot object representing the barplot with numbers
#' automatic wrapping of pathway names
#' facet wrap by rankname if present
#' @examples
#' df <- data.frame(
#'   pathway = c("Pathway 1", "Pathway 2", "Pathway 3"),
#'   NES = c(1.5, 2.0, 1.8),
#'   padj = c(1.05, 01, 0.001),
#'   leadingEdge = c("A", "B", "C"),
#'   size = c(10, 15, 20)
#' )
#' barplot_with_numbers(df, title = "Pathway Analysis")
#'
#' @export
barplot_with_numbers <- function(
    df,
    title = "",
    subtitle = NULL,
    save_func = NULL,
    facet_order = NULL,
    nes_range = NULL, # or a list of len 2
    ...) {

  if (!is.null(nes_range)){
    if (length(nes_range) != 2){
       stop("nes_range should be a list of length 2")
    }
  }


  sel <- prepare_data_for_barplot(df)


  custom_labeller <- function(value) {
    wrapped_labels <- sapply(value, function(label) {
      label %>%
        str_replace_all("_", " ") %>%
        str_wrap(width = 54)
    })
    return(wrapped_labels)
  }
  labeller_func <- custom_labeller

  # Format title/subtitle with less aggressive wrapping
  formatted_title <- title %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_wrap(width = 54)
  formatted_subtitle <- if (is.null(subtitle)) NULL else subtitle %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_wrap(width = 72)


  if (!is.null(facet_order)) {
    sel <- sel %>%
      mutate(rankname = factor(rankname, levels = facet_order, ordered = T)) %>%
      arrange(rankname)
  }

  if ("rankname" %in% names(sel)) {
    # right here we make a mapping and rename for display
    # for facetted plots
    name_map <- util_tools$make_name_map(unique(sel$rankname))
    sel <- sel %>% dplyr::mutate(
                          rankname=dplyr::recode(
                                                 rankname,
                                                 !!!name_map  # named vector
                                                 )
                          )
  }

  # Dynamically choose a single y-axis font size based on longest label
  pathway_chars <- sel$pathway %>% as.character() %>%
    stringr::str_replace_all("\n", " ") %>%
    stringr::str_squish()
  max_label_chars <- if (length(pathway_chars) > 0) max(nchar(pathway_chars)) else 0
  axis_text_y_size <- dplyr::case_when(
    max_label_chars < 36 ~ 7.6,
    max_label_chars < 64 ~ 6.6,
    max_label_chars < 84 ~ 6.2,
    TRUE ~ 5.2
  )

  # Compute bubble-like fill value: sign(NES) * (1 - pval)
  if (!"pval" %in% colnames(sel)) {
    sel <- sel %>% mutate(pval = padj)
  }
  sel <- sel %>% mutate(
    pval = {
      tmp <- pval
      if (is.list(tmp)) {
        tmp <- vapply(tmp, function(x) as.numeric(x)[1], numeric(1), USE.NAMES = FALSE)
      }
      tmp <- suppressWarnings(as.numeric(tmp))
      ifelse(is.na(tmp), 1, tmp)
    },
    fill_value = sign(NES) * (1 - pmin(pval, 1))
  )
  # rewrite below, check for equivalence 
  # get_size <- function(x) {
  #   # font size for geneset names
  #   x <- as.character(x)
  #   n <- nchar(x)

  #   cut(n,
  #       breaks = c(-Inf, 27, 64, 84, Inf),
  #       labels = c(7.6, 6.6, 6.2, 5.2),
  #       right = FALSE
  #   ) |> as.numeric()
  # }


  p <- sel %>%
    ggplot2::ggplot(
      aes(
        y = pathway,
        x = NES,
        # size = leadingEdgeNum,
        fill = fill_value,
        # fill = rankname
        # color = outline_val,
        # col = id
      )
    ) +
    # Reference line at x=0 behind bars
    geom_vline(xintercept = 0, colour = scales::alpha("#555555", 0.6), size = 0.4, show.legend = FALSE) +
    # scale_color_gradient(low = "#0000ffee", high = "#ff0000ee") +  # Adjust colors to represent p-values
    # geom_point() +
    scale_fill_gradient2(
      low = "#084594",
      mid = "#ffffff",
      high = "#b30000",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "#f7f7f7",
      guide = guide_colourbar(title = "1 - pval")
    ) +
    geom_col(linewidth = 0.8, aes(color = outline_val)) +
    scale_color_identity() +
    scale_y_discrete(position = "right") +
    labs(title = formatted_title, subtitle = formatted_subtitle) +
    # scale_color_manual(values=c("black", 'blue'))+
    # scale_size(range = c(4, 12)) +  # Adjust point sizes
    geom_text(
      aes(
        label = leadingEdgeFrac,
        x = sign(NES) * .46,
      ),
      color = "white",
      fontface = "bold",
      size = 2.2,
      # position = position_dodge(width = 0.8),
      vjust = 0.5, hjust = 0.5
    ) +
    theme_bw() + theme(
      axis.text.y = element_text(size = axis_text_y_size, face = "bold"),
      axis.text.x = element_text(size = 7.0),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      strip.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0)
    )

  if (!is.null(nes_range)){
    nes_max <- max(sel$NES, na.rm=T) * 1.05
    nes_min <- min(sel$NES, na.rm=T) * 1.05
    # Expand user range if needed
    nes_range[1] <- min(nes_range[1], nes_min)
    nes_range[2] <- max(nes_range[2], nes_max)
    # use symmetric range only if both sides have signal
    if (nes_min < 0 && nes_max > 0) {
      max_abs <- max(abs(nes_min), abs(nes_max))
      nes_range <- c(-max_abs, max_abs)
    } else {
      nes_range <- c(min(nes_min, 0, na.rm=T), max(nes_max, 0, na.rm=T))
    }
    p <- p + xlim(nes_range[1], nes_range[2])

  }

  if ("rankname" %in% colnames(df) && (length(unique(df$rankname)) > 1)) {
    log_msg(info = "facet wrapping by rankname")
    p <- p + facet_wrap(~rankname, labeller = as_labeller(labeller_func))
    num_panels <- length(unique(df$rankname))
    ncol <- ceiling(sqrt(num_panels)) # ggplot2 default behavior if ncol is not specified
    nrow <- ceiling(num_panels / ncol) # Calculate rows
  } else {
    ncol <- 1
    nrow <- 1
  }
  panel_width_in <- 4.0
  # Height scales linearly with number of rows and respects facet rows
  n_pathways <- dplyr::n_distinct(sel$pathway)
  text_scale <- axis_text_y_size / 6.6
  per_row_in <- dplyr::case_when(
    n_pathways <= 20 ~ 0.28,
    n_pathways <= 60 ~ 0.25,
    TRUE ~ 0.26
  ) * text_scale
  per_row_in <- max(min(per_row_in, 0.34), 0.18)
  min_panel_height_in <- 3.2
  facet_strip_pad_in <- 0.35
  panel_height_in <- max(min_panel_height_in, per_row_in * n_pathways + facet_strip_pad_in)

  # Calculate total figure size
  total_width_in <- 2.2 + (panel_width_in * ncol)
  total_height_in <- panel_height_in * nrow
  log_msg(msg = paste0("total width: ", total_width_in, " total height: ", total_height_in))


  if (!is.null(save_func)) {
    save_func(
      plot_code = function() {
        print(p)
      },
      width = total_width_in,
      height = total_height_in
    )
  }
  return(p)
}

# Format a readable source/collection label for subtitles
format_source_label <- function(name) {
  name %>%
    stringr::str_replace_all("[:._]", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_squish()
}

all_barplots_with_numbers <- function(
    results_list,
    save_func = NULL,
    facet_order = NULL,
    limit = 20,
    ...) {
  if (!is.null(save_func)) {
    filename <- util_tools$safe_filename(
      get_arg(save_func, "filename"),
      "bar",
      fallback = "bar"
    )
    save_func <- make_partial(save_func, filename = filename)
  }

  plts <- results_list %>% purrr::imap(
    ~ {
      collection_name <- .y
      collection_label <- format_source_label(collection_name)
      list_of_comparisons <- .x
      # Shorten comparison names by stripping shared affixes for path components
      comparison_names <- names(list_of_comparisons)
      name_map <- util_tools$make_name_map(comparison_names)
      plts <- list_of_comparisons %>% imap(
        ~ {
          dataframe <- .x
          comparison_name <- .y
          comparison_label <- name_map[[comparison_name]] %||% comparison_name
          .plts <- limit %>% purrr::map(
            ~ {

              nes_max <- max(abs(dataframe$NES), na.rm=T)
              nes_range <- c(-nes_max, nes_max)
              .limit <- .x
              sel <- fgsea_tools$select_topn(dataframe, limit = .limit, pstat_cutoff=1)
              if (nrow(sel) == 0) {
                log_msg(warning = paste0(
                  "no pathways available for barplot (requested top ", .limit, ") in ",
                  collection_name, " / ", comparison_name
                ))
                return(NULL)
              }
              pathway_count <- dplyr::n_distinct(sel$pathway)
              effective_limit <- min(.limit, pathway_count)
              .title <- comparison_label %||% comparison_name # %>% fs::path_file() %>% fs::path_ext_remove() #%>% gsub(pattern="_", replacement=" ", x=.)

              rank_label <- sel$rankname %>% na.omit() %>% unique()
              rank_label <- if (length(rank_label) == 0) NULL else rank_label %>%
                stringr::str_replace_all("_", " ") %>%
                paste(collapse = ", ")
              base_subtitle <- if (!is.null(rank_label)) {
                paste0("rank: ", rank_label, " • top ", effective_limit)
              } else {
                paste0("top ", effective_limit, " pathways")
              }
              subtitle_text <- paste0(base_subtitle, " • source: ", collection_label)

              local_save_func <- save_func
              save_path <- NULL
              if (!is.null(local_save_func)) {
                collection_dir <- util_tools$safe_path_component(collection_name)
                comparison_dir <- util_tools$safe_path_component(comparison_label)
                save_path <- util_tools$safe_subdir(
                  get_arg(local_save_func, "path"),
                  collection_dir,
                  "bar",
                  comparison_dir
                )
                filename <- util_tools$safe_filename(
                  "bar",
                  paste0("top", effective_limit),
                  paste0("n", pathway_count),
                  fallback = "bar"
                )
                local_save_func <- make_partial(local_save_func, filename = filename, path = save_path)
                log_msg(msg = paste0(
                  "plotting bar target ", file.path(save_path, paste0(filename, ".pdf")), " ", collection_name
                ))
              }

              p <- barplot_with_numbers(sel,
                title = .title,
                subtitle = subtitle_text,
                save_func = local_save_func,
                facet_order = facet_order,
                nes_range = nes_range,
                ...
              )

              return(p)
            }
          )
          return(.plts)
        }
      )
    }
  )
  return(plts)
}


do_combined_barplots <- function(
    results_list,
    save_func = NULL,
    facet_order = NULL,
    limit = 20,
    ...) {
  genesets <- names(results_list)

  # args <- list(...)
  # if ("save_func" %in% names(args)) {
  #   save_func <- args$save_func
  # } else {
  #   save_func <- NULL
  # }

  purrr::map(genesets, ~ {
    geneset_name <- .x
    fgsea_res_list <- results_list[[geneset_name]]
    collection_label <- format_source_label(geneset_name)
    # geneset <- genesets_list[[geneset_name]]

    plts <- limit %>% purrr::map(~ {
      .limit <- .x
      res <- fgsea_res_list %>% bind_rows(.id = "rankname") # all comparisons 1 gene set
      # res <- res %>% mutate(rankname = rankname %>% fs::path_file() %>% fs::path_ext_remove())
      res <- fgsea_tools$select_topn(res, limit = .limit, pstat_cutoff=1)
      n_sel <- res %>%
        distinct(pathway) %>%
        nrow()
      effective_limit <- min(.limit, n_sel)
      if (n_sel == 0) {
        log_msg(warning = paste0(
          "no pathways available for combined barplot (requested top ", .limit, ") for ",
          geneset_name
        ))
        return(NULL)
      }
      log_msg(msg = paste0(
        "bar combined: geneset=", geneset_name,
        " limit=", .limit,
        " selected_rows=", nrow(res),
        " distinct_pathways=", n_sel
      ))
      # Symmetric NES range across all facets for consistency
      nes_max <- suppressWarnings(max(abs(res$NES), na.rm = TRUE))
      nes_range <- if (is.finite(nes_max)) c(-nes_max, nes_max) else NULL
      # pathway_df <- get_pathway_info(geneset_name)
      # .merge <- left_join(res , pathway_df, by= )
      local_save_func <- save_func
      if (!is.null(local_save_func)) {
        geneset_dir <- util_tools$safe_path_component(geneset_name)
        filename <- util_tools$safe_filename(
          get_arg(local_save_func, "filename"),
          "bar",
          paste0("top", effective_limit),
          paste0("n", n_sel),
          "all",
          fallback = "bar_all"
        )
        path <- util_tools$safe_subdir(get_arg(local_save_func, "path"), geneset_dir, "bar")
        log_msg(msg = paste0(
          "plotting bar target ", file.path(path, paste0(filename, ".pdf")), " ", geneset_name, " "
        ))
        local_save_func <- make_partial(local_save_func, filename = filename, path = path)
      }
      p <- res %>% barplot_with_numbers(
        title = geneset_name,
        subtitle = paste0("top ", effective_limit, " pathways • source: ", collection_label),
        nes_range = nes_range,
        save_func = local_save_func,
        facet_order = facet_order,
        ...
      )
      return(p)
    })
    return(plts)
  })
}



plot_results_all_collections <- function(
    list_of_lists,
    metadata = NULL,
    cut_by = NULL,
    limit = 120,
    pstat_cutoff = 1,
    pstat_usetype = "padj",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    main_pathway_ratio = 0.1,
    save_func = NULL,
    rankname_order = NULL,
    meta_to_include = NULL,
    meta_to_exclude = NULL,
    # full_path = NULL,
    # sample_order = NULL,
    ...) {
  xtra_args <- list(...) # dunno what to do with these

  # if (length(limit) == 1) {
  #   limit <- rep(limit, length(list_of_lists))
  # }
  res <- list_of_lists %>% purrr::imap(
    ~ {
      .data <- .x
      .title <- .y
      param_grid <- expand.grid(
        limit = limit,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        cut_by = cut_by %||% NA,
        stringsAsFactors = F
      )
      plts <- param_grid %>% purrr::pmap(
        ~ {
          params <- list(...)
          .cut_by <- params$cut_by
          limit <- params$limit
          cluster_rows <- params$cluster_rows
          cluster_columns <- params$cluster_columns
          .cut_by_val <- plot_utils$process_cut_by(.cut_by, metadata)
          cut_by_label <- if (!is.null(.cut_by_val)) {
            paste0("cut_", util_tools$safe_path_component(.cut_by, max_chars = 32))
          } else {
            ""
          }

          local_save_func <- save_func
          if (!is.null(local_save_func)) {
            title_dir <- util_tools$safe_path_component(.title)
            filename <- util_tools$safe_filename(
              get_arg(local_save_func, "filename"),
              "gsea_heatmap",
              util_tools$cluster_flag_token(cluster_rows, "rc"),
              util_tools$cluster_flag_token(cluster_columns, "cc"),
              cut_by_label,
              title_dir,
              fallback = "gsea_heatmap"
            )
            path <- util_tools$safe_subdir(get_arg(local_save_func, "path"), title_dir, "heatmaps_gsea")
            local_save_func <- make_partial(local_save_func, filename = filename, path = path)
          }

          log_msg(info = paste0("cutting by ", .cut_by))
          plot_results_one_collection(.data,
            title = .title,
            metadata = metadata,
            cut_by = .cut_by,
            limit = limit,
            pstat_cutoff = pstat_cutoff,
            pstat_usetype = pstat_usetype,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns,
            save_func = local_save_func,
            rankname_order = rankname_order,
            meta_to_include = meta_to_include,
            meta_to_exclude = meta_to_exclude,
          )
        }
      )
      return(plts)
    }
  )
  # res <- purrr::map(list_of_lists, function(item) {
  #   do.call("plot_results_one_collection", c(list(df = item), args))
  # })
  return(res)
}

plot_results_one_collection <- function(
    df,
    metadata = NULL,
    pathway_metadata = NULL,
    title = "",
    cut_by = NULL,
    limit = 120,
    pstat_cutoff = 1,
    pstat_usetype = "padj",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    sample_order = NULL,
    rankname_order = NULL,
    meta_to_exclude = c("recno", "runno", "searchno", "label"),
    meta_to_include = NULL,
    row_annotation = NULL, # can be passed if made somewhere else
    save_func = NULL,
    ...) {
  if (nrow(df) == 0){
   log_msg(msg = paste0("no results!"))
   return()
  }
  log_msg(msg = paste0("calling plot results one collection"))


  default_meta_to_exclude <- c(
    "recno", "runno", "searchno", "label", "expr_col", "expr_file",
    "assay", "rec_run_search",
    "SampleName",
    "SampleName_old",
    "RawSampleID",
    "SampleID",
    "SampleLabel",
    "SampleTube",
    "RNA_GARP_SampleID",
    "upper_quartile",
    "median_upper_quartile",
    "Unnamed",
    "QC_recno",
    "Tube",
    "TubeLabel"
  )
  # If user specifies additional meta_to_exclude, merge with defaults
  if (!is.null(meta_to_exclude)) {
    meta_to_exclude <- union(default_meta_to_exclude, meta_to_exclude)
  } else {
    meta_to_exclude <- default_meta_to_exclude
  }

  # log_msg(msg = paste0("rankname order is ", rankname_order))

  # Ensure necessary columns are present
  required_cols <- c("pathway", "NES")
  if (!all(required_cols %in% colnames(df))) {
    stop("Required columns not found in dataframe")
  }

  # if (is.null(sample_order)) {
  #   sample_order <- unique(df$rankname)
  # } else {
  #   sample_order <- union(sample_order, unique(df$rankname))
  # }

  requested_rank_order <- rankname_order
  available_ranknames <- unique(df$rankname)

  if (!is.null(requested_rank_order)) {
    missing_ranknames <- setdiff(requested_rank_order, available_ranknames)
    if (length(missing_ranknames) > 0) {
      log_msg(warning = paste0(
        "rankname_order entries not present in results: ",
        paste(missing_ranknames, collapse = ", "),
        ". These names will be ignored."
      ))
    }

    duplicated_ranknames <- requested_rank_order[duplicated(requested_rank_order)]
    if (length(duplicated_ranknames) > 0) {
      log_msg(warning = paste0(
        "rankname_order contains duplicate labels: ",
        paste(unique(duplicated_ranknames), collapse = ", "),
        ". Keeping the first occurrence of each."
      ))
    }
  }

  if (is.null(rankname_order)) {
    rankname_order <- available_ranknames
    log_msg(debug = "rankname_order not supplied; using detected order from results")
  } else {
    rankname_order <- intersect(rankname_order, available_ranknames)
    if (length(rankname_order) == 0) {
      warning("rankname_order is empty, using default")
      rankname_order <- available_ranknames # default backup
    } else if (!is.null(requested_rank_order)) {
      if (length(rankname_order) < length(unique(requested_rank_order))) {
        log_msg(info = paste0(
          "Applied rankname_order for ",
          length(rankname_order),
          " comparisons (", length(unique(requested_rank_order)), " requested)."
        ))
      }
    }
  }


  # Align metadata with dataframe
  if (!is.null(metadata) && length(metadata) > 0) {
    if (!all(rownames(metadata) %in% df$rankname)) {
      cat("Metadata not aligned with df")
      log_msg(msg = "Metadata not aligned with df")
      metadata <- NULL
    }
  }

  if (is.null(metadata)) {
    metadata <- data.frame(id = unique(df$rankname), dummy = "X") # necessary for some reason
    rownames(metadata) <- metadata$id
  }
  metadata <- metadata[rankname_order, ] # order metadata by rankname_order

  metadata %<>% dplyr::select(-any_of(meta_to_exclude)) # remove columns

  # Handling cut_by parameter
  # print(paste0("cutting by :", cut_by))
  if (!is.null(cut_by)) {
    if (length(cut_by) > 1) {
      cut_by <- cut_by[1]
      warning("cut_by is a vector, using first element")
    }
    if (!is.null(cut_by) && cut_by %in% colnames(metadata)) {
      cut_by <- metadata[, cut_by]
      cut_by <- factor(cut_by, levels = unique(cut_by))
    } else {
      warning("cut_by not found in metadata, setting to NULL")
      cut_by <- NULL
    }
  }

  # calculate quantile before selection
  # Set up color scale
  q01 <- quantile(abs(df$NES), 0.99, na.rm = TRUE)
  q01 <- max(q01, 2.8)
  # print(paste0("q01 is ", q01))
  num_colors <- 9
  my_colors <- colorRampPalette(c(
    "#0000eebb",
    # "#0000bbbb",
    # "#8888ffbb",
    # "#ddddff77",
    "#ffffff",
    # "#ffdddd77",
    # "#ff8888bb"
    # "#bb0000bb",
    "#ee0000bb"
  ))(num_colors)
  break_points <- seq(-q01, q01, length.out = num_colors)
  col <- colorRamp2(breaks = break_points, colors = my_colors)

  df <- fgsea_tools$select_topn(
    df,
    pstat_cutoff = pstat_cutoff,
    pstat_usetype = pstat_usetype,
    limit = limit
  )

  if (nrow(df) == 0) {
    log_msg(msg = "No pathways to plot")
    return(NULL)
  }

  if (!is.null(row_annotation)) {
    if (!"HeatmapAnnotation" %in% class(row_annotation)) {
      warning("row_annotation is not a HeatmapAnnotation object")
      row_annotation <- NULL
    }
  }

  # Prepare data for heatmap
  dfp <- df %>%
    pivot_wider(id_cols = pathway, values_from = NES, names_from = rankname) %>%
    dplyr::mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "KEGG_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOMF_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOBP_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOCC_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "REACTOME_")) %>%
    as.data.frame()
  # dfp %<>% as.data.frame() # because we will be assigning row names


  rownames(dfp) <- dfp$pathway
  dfp <- dfp[, metadata$id, drop = FALSE] # select ranknames as ordered inside metadata dataframe

  # dfp["pathway"] <- NULL # now remove this column to exclude from heatmap

  dfp_padj <- df %>%
    pivot_wider(id_cols = pathway, values_from = padj, names_from = rankname) %>%
    as.data.frame()
  rownames(dfp_padj) <- dfp_padj$pathway
  dfp_padj["pathway"] <- NULL
  dfp_padj <- dfp_padj[, metadata$id, drop = FALSE]
  # %>% select(-pathway, all_of(metadata$id))
  logical_matrix <- dfp_padj < 0.25
  star_matrix <- ifelse(logical_matrix, "*", "")
  star_matrix <- star_matrix[, metadata$id, drop = FALSE] # shouldn't be necessary as comes from dfp_padj
  star_matrix[is.na(star_matrix)] <- ""


  # Prepare column annotation
  # log_msg(msg = paste0("metadata: ", head(metadata)))

  # it isn't strictly necessary to exclude these colum s here,
  # as they will be excluded upon creation of the column_annotation object
  metadata_for_colors <- metadata %>% dplyr::select(-any_of(c("id", "dummy")))
  # Coerce numeric-like metadata (incl. decimals) to numeric for continuous annotations
  if (ncol(metadata_for_colors) > 0) {
    for (.colname in colnames(metadata_for_colors)) {
      .col <- metadata_for_colors[[.colname]]
      if (plot_utils$is_numericish(.col)) {
        suppressWarnings({ metadata_for_colors[[.colname]] <- as.numeric(as.character(.col)) })
      }
    }
  }
  if (ncol(metadata_for_colors) > 0) {
    # colors_list <- metadata_for_colors %>% colnames() %>%
    #   map(~{
    #     # Extract unique values for the current column
    #     .unique_vals <- unique(metadata[[.x]])

    #     # Generate a qualitative color palette
    #     .colors <- colorspace::qualitative_hcl(n = length(.unique_vals))

    #     # Assign unique values as names to the colors
    #     names(.colors) <- .unique_vals

    #     # Return the named color vector
    #     .colors
    #   }) %>% set_names(colnames(metadata_for_colors))
    colors_list <- plot_utils$create_named_color_list(metadata_for_colors,
      colnames(metadata_for_colors),
      c = 124
    )

    # log_msg(msg = paste0("color list ", as.character(colors_list)))

    column_annotation <- ComplexHeatmap::columnAnnotation(
      df = metadata_for_colors,
      col = colors_list
    )
  } else {
    column_annotation <- NULL
  }

  # Construct heatmap

  heatmap_legend_param <- list(
    title = "NES",
    direction = "horizontal",
    just = "bottom",
    legend_width = unit(6, "cm"),
    at = break_points %>% round(1)
  )

  # cell_fun <- function(j, i, x, y, width, height, fill) {
  #   # Retrieve the value that indicates whether to draw an asterisk
  #   value <- star_matrix[i, j]
  #   if (value == "*") {
  #     # Draw asterisk
  #     grid::grid.text(value, x, y, gp = grid::gpar(fontsize = 12, col = "black"))
  #     # Draw border around the cell
  #     grid::grid.rect(x, y,
  #       width = width, height = height,
  #       gp = grid::gpar(col = "black", fill = NA, lwd = 1)
  #     )
  #   }
  # }

  cell_fun <- function(j, i, x, y, width, height, fill) {
    # Ensure value is not NA before comparison
    # value <- star_matrix[i, j]
    value <- dfp_padj[i, j]
    # print(paste("Processing cell:", i, j, "Value:", value))
    # cat(sprintf("NES Cell [%d, %d] with value '%s'\n", i, j, dfp[i, j]))
    # cat(sprintf("star_matrix Cell [%d, %d] with value '%s'\n", i, j, star_matrix    [i, j]))

    if (!is.na(value) && value < 0.05) {
      # Draw asterisk
      grid::grid.text(
        "\u2B51",
        # "*",
        x, y,
        gp = grid::gpar(fontsize = 12, col = "black")
      )
    }
    if (!is.na(value) && value < 0.25) {
      # Draw border around the cell
      grid::grid.rect(x, y,
        width = width, height = height,
        gp = grid::gpar(col = "black", fill = NA, lwd = 1)
      )
    }
  }

  # height <- 6 + (nrows(dfp) * .16)
  .ncol <- ncol(dfp)
  heatmap_matrix_width <- unit(ncol(dfp) * .2, "in")
  heatmap_matrix_height <- unit(nrow(dfp) * .2, "in")

  # .row_fontsizes <- ifelse(nchar(rownames(dfp) < 20), 9.8, 5.8)
  # make the below better, later
  .column_fontsizes <- lapply(colnames(dfp), function(x) ifelse(nchar(x) < 22, 9.6, ifelse(nchar(x) < 28, 7.6, ifelse(nchar(x) < 54, 6.2, 5.2))))
  .row_fontsizes <- lapply(rownames(dfp), function(x) ifelse(nchar(x) < 36, 7.6, ifelse(nchar(x) < 64, 6.6, ifelse(nchar(x) < 84, 6.2, 5.2))))
  # rownames_gp <- ifelse(nchar(rownames(dfp) < 20), c(grid::gpar(fontsize=9.8, lineheight=.8)), c(grid::gpar(fontsize=5.8, lineheight=.8)))

  ht <- ComplexHeatmap::Heatmap(
    dfp %>% as.matrix(),
    name = "mat",
    col = col,
    na_col = "#444444",
    top_annotation = column_annotation,
    left_annotation = row_annotation,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    heatmap_legend_param = heatmap_legend_param,
    column_gap = unit(1, "mm"),
    width = heatmap_matrix_width,
    height = heatmap_matrix_height,
    # width = ncol(dfp) * 9,
    # height = nrow(dfp) * 5,
    border = T,
    column_split = cut_by,
    row_labels = rownames(dfp) %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 42),
    # row_names_gp = grid::gpar(fontsize = 6.8, lineheight=.8),
    row_names_gp = grid::gpar(fontsize = .row_fontsizes, lineheight = .8),
    column_names_gp = grid::gpar(fontsize = .column_fontsizes, lineheight = .8),
    column_labels = colnames(dfp) %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 36),
    row_names_side = "right",
    # row_names_rot=(pi/24)*180,
    column_title_gp = grid::gpar(fontsize = 12, hjust = 2),
    clustering_distance_rows = util_tools$dist_no_na,
    clustering_distance_columns = util_tools$dist_no_na,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    column_names_side = "top",
    column_title = title,
    cell_fun = cell_fun
  )


  # log_msg(msg = paste0("defining draw func"))
  do_draw <- function() {
    ht <- draw(ht,
      heatmap_legend_side = "bottom",
      padding = unit(c(2, 24, 2, 24), "mm"), # top, left, bottom, right
    )
    xunit <- ifelse(cluster_rows == TRUE, 1, 2.4)
    decorate_heatmap_body("mat", {
      grid.text(
        paste0(
          "\u25A0  padj < .25",
          "\n",
          "\u2605  padj < .05"
        ),
        unit(xunit, "cm"), unit(-4, "cm"),
        gp = gpar(fontsize = 9)
      )
    })
    return(ht)
  }

  ht <- do_draw()



  log_msg(debug = paste0("save func: ", class(save_func) %>% as.character()))
  log_msg(debug = paste0("is null save func: ", is.null(save_func)))

  height <- 4 + (nrow(dfp) * .20) + (ncol(metadata_for_colors) / 3)
  width <- 8 + (ncol(dfp) * .26)


  if (!is.null(save_func)) {
    log_msg(debug = paste0("save func attrs before: "))
    log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))
    filename <- paste0(get_arg(save_func, "filename"), nrow(dfp), "x", ncol(dfp))

    save_func <- make_partial(save_func, height = height, width = width, filename = filename)

    log_msg(debug = paste0("save func attrs after: "))
    log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))

    save_func(plot_code = do_draw)
  }

  return(ht)
}

make_selection <- function(x) {
  .first <- which(x$ES == es_min)
  .second <- which(x$ES == es_max)
  .top <- x %>%
    arrange(rank) %>%
    head(.first)
  .bot <- x %>%
    arrange(rank) %>%
    tail(dim(x)[1] - .second)
  return(bind_rows(.top, .bot))
}


edgeplot1 <- function(rankorder_object, ...) {
  posES <- rankorder_object$posES
  negES <- rankorder_object$negES
  rankorder_edge <- rankorder_object$edge
  # rankorder_edge %>%
  #   filter(!is.na(stat)) %>%
  #   dim()
  # ggplot(aes(x = rank, y = ES)) +
  #   geom_point()

  p <- rankorder_edge %>%
    # filter(!is.na(stat_tick)) %>%
    ggplot(aes(x = stat_tick, y = ES, col = rank)) +
    geom_point() +
    scale_color_continuous(type = "viridis", option = "H") +
    geom_hline(yintercept = posES, colour = "red", linetype = "dashed") +
    geom_hline(yintercept = negES, colour = "red", linetype = "dashed")

  p
}


# across all pathways
plot_top_ES_across <- function(
    gsea_results,
    ranks_list,
    geneset_collections,
    limit = 30,
    do_individual = TRUE,
    do_combined = TRUE,
    combine_by = NULL,
    combine_by_name = NULL,
    save_func = NULL,
    combined_show_ticks = FALSE,
    width = 3.4,
    height = 4,
    pathways_of_interest = NULL,
    ...) {
  if (!"list" %in% class(gsea_results)) {
    stop(cat("gsea_results should be a list of data frames"))
  }

  if (!"list" %in% class(geneset_collections)) {
    stop(cat("geneset_collections should be a list of data frames"))
  }

  missing <- setdiff(
    names(gsea_results),
    names(geneset_collections)
  )
  if (length(missing) > 0) {
    cat(paste0("missing some genesets..."))
  }

  list_of_plts <- gsea_results %>%
    purrr::imap(~ {
      df <- .x
      collection_name <- .y
      geneset_collection <- geneset_collections[[collection_name]]
      if (!is.null(save_func)) {
        path <- get_arg(save_func, "path")
        collection_dir <- util_tools$safe_path_component(collection_name)
        newpath <- util_tools$safe_subdir(path, collection_dir, "enrichplots")
        # Avoid duplicating collection info in filename; directory already encodes it
        filename <- util_tools$safe_filename(
          get_arg(save_func, "filename"),
          fallback = "enrichplots"
        )
        save_func <- make_partial(save_func, path = newpath, filename = filename)
      }
      plts <- plot_top_ES(df, ranks_list, geneset_collection,
        limit = limit,
        do_individual = do_individual,
        do_combined = do_combined,
        combine_by = combine_by,
        combine_by_name = combine_by_name,
        save_func = save_func,
        panel_width = width,
        panel_height = height,
        combined_show_ticks = combined_show_ticks,
        pathways_of_interest = pathways_of_interest
      )
      return(plts)
    })
  return(list_of_plts)
}

# across all comparisons in one pathway
plot_top_ES <- function(
    df,
    ranks_list,
    geneset_collection,
    limit = 30,
    do_individual = T,
    do_combined = T,
    combine_by = NULL,
    combine_by_name = NULL,
    save_func = NULL,
    panel_width = 3.6,
    panel_height = 3.4,
    combined_show_ticks = FALSE,
    combined_label_size = 1.75,
    filter_on_mainpathway = TRUE,
    pathways_of_interest = NULL,
    ...) {
  #

  xtra <- NULL
  if (!is.null(pathways_of_interest)){
    xtra <- df %>% dplyr::filter(pathway %in% pathways_of_interest)
  }

  if (filter_on_mainpathway == TRUE) {
    df <- fgsea_tools$filter_on_mainpathway(df)
  }

  df %<>% fgsea_tools$select_topn(limit = limit)

  if (!is.null(xtra)){
      df %<>% dplyr::bind_rows(xtra)
  }

  if (nrow(df) == 0) {
    return(NULL)
  }
  pathways <- df$pathway %>% unique()
  geneset_lists <- geneset_collection[pathways]

  # print(paste0("limit: ", limit))
  # print(nrow(df))
  # print(paste0(" do individual: ", do_individual))
  # print(paste0(" do combined: ", do_combined))
  # print(paste0("pathways ", pathways))
  # print(paste0("geneset lists ", geneset_lists))

  rankorder_by_pw <- fgsea_tools$get_rankorder_across(
    df,
    ranks_list,
    geneset_lists,
    limit = max(120, limit), # why do we even need to set this?
    pathways_of_interest = pathways_of_interest
  ) # this yields a named list
  # names are pathways
  # values are named lists of rankorder for each sample


  plts <- list()

  if (do_combined) {
    if (is.null(combine_by)) {
      combine_by <- list() ##
    }
    # everything needed for plotting gets assigned here
    rankorder_samplepathway <- fgsea_tools$combine_rankorders_on_sample(
      rankorder_by_pw,
      metadata = combine_by
    )

    #combine_by <- plot_utils$process_cut_by(cut_by, cdesc_for_annot)
    combine_by_name <- combine_by_name %||% "all"
    combine_by_name <- ifelse(combine_by_name == FALSE, "all", combine_by_name)



    .plts <- rankorder_samplepathway %>% purrr::imap(~ {
      rankorder <- .x # rankorder is a list of dataframes, each df has the facet info
      pathway_name <- .y

      log_msg(msg = paste0("plotting combined ", pathway_name))
      log_msg(msg = paste0("show ticks: ", combined_show_ticks))

      wrap_width <- getOption(util_tools$pkg_option_name("enrichplot_title_wrap"), 40)
      .plt <- plotES_combined(rankorder,
        title = pathway_name %>% str_replace_all("_", " ") %>% str_wrap(width = wrap_width),
        show_ticks = combined_show_ticks,
        label_size = combined_label_size,
      ) # will be faceted if "facet" in .x


      if (!is.null(save_func)) {
        pathway_label <- strip_pathway_leader(pathway_name)
        newdir <- util_tools$safe_filename(
          get_arg(save_func, "filename"),
          pathway_label,
          fallback = "combined_dir",
          max_chars = 60
        )
        newpath <- util_tools$safe_subdir(get_arg(save_func, "path"), newdir)
        newfilename <- util_tools$safe_filename(
          "combined",
          combine_by_name,
          fallback = "combined"
        )

        if (!"facet" %in% names(rankorder$edge)) { # this should be true if do_combined is true
          rankorder_samplepathway$facet <- "x"
        }
        num_panels <- max(length(unique(rankorder$edge$facet)), 1)
        # more scaling challenges
        if (num_panels > 1) {
          panel_width <- panel_width * .68
          panel_height <- panel_height * .75
        }
        .ncol <- ceiling(sqrt(num_panels)) # ggplot2 default behavior if ncol is not specified
        .nrow <- ceiling(num_panels / .ncol) # Calculate rows
        .width <- panel_width * .ncol
        .height <- panel_height * .8 * .nrow

        save_func <- make_partial(save_func, path = newpath, filename = newfilename, width = .width, height = .height)
        # and now call it
        save_func(plot_code = function() print(.plt))
      }
      return(.plt)
    })
    plts <- c(plts, .plts)
  }

  # invidual plots
  if (!do_individual) {
    return(plts)
  }


  # individual plots
  .plts <- rankorder_by_pw %>%
    purrr::imap(~ {
      rankorders <- .x
      pathway_name <- .y
      # Build name mapping for comparison labels within this pathway's rankorders
      name_map <- util_tools$make_name_map(names(rankorders))

      rankorders %>% purrr::imap(~ {
        rankorder <- .x
        comparison <- .y
        comparison_label <- name_map[[comparison]] %||% comparison
        .stats <- df %>% dplyr::filter(pathway == pathway_name, rankname == comparison)
        if (nrow(.stats) == 0 ) {
            # problem
            return(NULL)
        }
        wrap_width <- getOption(util_tools$pkg_option_name("enrichplot_title_wrap"), 54)
        comp_title <- comparison_label %>% stringr::str_replace_all("_", " ") %>% stringr::str_wrap(width = wrap_width)
        path_title <- pathway_name %>% stringr::str_replace_all("_", " ") %>% stringr::str_wrap(width = wrap_width)
        .title <- paste0(comp_title, "\n", path_title)
        .subtitle <- ""
        if (nrow(.stats) == 1) {
          .subtitle <- paste0(
            "ES: ", .stats[["ES"]] %>% round(2) %>% as.character(),
            "\tNES: ", .stats[["NES"]] %>% round(2) %>% as.character(),
            "\tpval: ", .stats[["pval"]] %>% round(2) %>% as.character(),
            "\tpadj: ", .stats[["padj"]] %>% round(2) %>% as.character()
          )
        }
        plt <- plotES(rankorder, title = .title, subtitle = .subtitle)


        if (!is.null(save_func)) {
          newdir <- util_tools$safe_filename(
            get_arg(save_func, "filename"),
            pathway_name,
            fallback = "enrich",
            max_chars = 60
          )
          newpath <- util_tools$safe_subdir(get_arg(save_func, "path"), newdir)
          newfilename <- util_tools$safe_filename(comparison_label, fallback = "comparison")
          save_func <- make_partial(save_func, path = newpath, filename = newfilename)
          # and now call it
          save_func(
            plot_code = function() print(plt),
            width = panel_width,
            height = panel_height
          )
        }
        return(plt)
      })
    })
  plts <- c(plts, .plts)
  return(plts)
}

plotES_combined <- function(enplot_data, label_size = 1.85, title = "", show_ticks = F, ...) {
  # plot ES for multiple samples/comparisons.
  spreadES <- max(enplot_data$curve$ES, na.rm = T) - min(enplot_data$curve$ES, na.rm = T)
  # print(spreadES)
  # look at add_color_mappings and similar called from combine_rankorders_on_sample for color assignments
  curve <- enplot_data$curve
  ticks <- enplot_data$ticks
  p <- ggplot(data = curve) +
    geom_line(aes(x = rank, y = ES, color = rankname), alpha = .6, show.legend = F) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_label_repel(
      data = curve %>% group_by(rankname) %>% arrange(-abs(ES)) %>% slice_head(n = 1) %>% ungroup(),
      # what is this doing? above?
      aes(label = rankname %>% stringr::str_replace_all("_", " ") %>% stringr::str_wrap(width = 28),
        x = rank, y = ES, color = rankname),
      # nudge_x = 0.5,
      # nudge_y = 0.5,
      size = label_size,
      show.legend = FALSE,
      max.overlaps = Inf,
      fill = "#FFFFFFcc",
      #xlim = c(-Inf, Inf), ylim = c(-Inf, Inf) # experimental

    )
  if (show_ticks == TRUE) {
    p <- p + geom_segment(
      data = ticks,
      mapping = aes(
        x = rank, y = -spreadES / 32,
        xend = rank, yend = spreadES / 32,
        color = rankname
      ),
      alpha = .4,
      show.legend = F
    )
  }
  p <- p + theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    title = element_text(size = 8),
    subtitle = element_text(size = 5)
  ) +
    labs(x = "rank", y = "enrichment score", title = title)
  if ("facet" %in% names(enplot_data$curve)) {
    p <- p + facet_wrap(~facet)
  }
  p
}

plotES <- function(enplot_data, ticksSize = 4, title = "", subtitle = "") {
  # this is directly from the example in fgsea plotEnrichmentData

  expected_names <- c("curve", "stats", "ticks", "maxAbsStat", "spreadES", "posES", "negES")
  for (name in expected_names) {
    if (!name %in% names(enplot_data)) {
      stop(paste0("missing ", name, " in enplot_data"))
    }
  }

  maxAbsStat <- enplot_data$maxAbsStat
  spreadES <- enplot_data$spreadES
  posES <- enplot_data$posES
  negES <- enplot_data$negES

  p <- with(
    enplot_data,
    ggplot(data = curve) +
      geom_line(aes(x = rank, y = ES), color = "green") +
      geom_ribbon(
        data = stats,
        mapping = aes(
          x = rank, ymin = 0,
          ymax = stat / maxAbsStat * (spreadES / 4)
        ),
        fill = "grey"
      ) +
      geom_segment(
        data = ticks,
        mapping = aes(
          x = rank, y = -spreadES / 16,
          xend = rank, yend = spreadES / 16
        ),
        size = 0.2
      ) +
      geom_hline(yintercept = posES, colour = "red", linetype = "dashed") +
      geom_hline(yintercept = negES, colour = "red", linetype = "dashed") +
      geom_hline(yintercept = 0, colour = "black") +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey92")
      ) +
      labs(x = "rank", y = "enrichment score")
  )
  p <- p + labs(title = title, subtitle = subtitle) +
    theme(
      title = element_text(size = 8),
      subtitle = element_text(size = 5)
    )
  return(p)
}



plot_table <- function(fgsea_res,
                       ranks,
                       pathways,
                       gsea_param = 1.0,
                       savefunc = NULL,
                       ...) {
  # top_pathways_u <- fgsea_res[ES > 0][head(order(pval), n=10), pathway]
  # top_pathways_d <- fgsea_res[ES < 0][head(order(pval), n=10), pathway]
  top_pathways_u <- fgsea_res %>%
    filter(ES > 0) %>%
    arrange(pval) %>%
    head(10) %>%
    pull(pathway)
  top_pathways_d <- fgsea_res %>%
    filter(ES < 0) %>%
    arrange(pval) %>%
    head(10) %>%
    pull(pathway)
  top_pathways <- c(top_pathways_u, rev(top_pathways_d)) #  rev reverse order
  tableobject <- pathways[top_pathways] %>%
    plotGseaTable(ranks,
      fgsea_res,
      gseaParam = gsea_param
      # gseaParam=1.0
    )
  return(tableobject)
}


plot_tables <- function(
    results_list,
    ranks_list,
    pathways_list) {
  # Ensure names are aligned and iterate over them
  pathway_names <- names(results_list)

  map(pathway_names, ~ {
    geneset_name <- .x
    fgsea_res_list <- results_list[[geneset_name]]
    genesets <- pathways_list[[geneset_name]]

    # Now, iterate over each rank file within this pathway category
    map(names(ranks_list), function(rank_name) {
      rankobj <- ranks_list[[rank_name]]
      # Generate a unique identifier for this plot/table, e.g., combining pathway name and rank file name
      rank_name_nice <- rank_name %>%
        fs::path_file() %>%
        fs::path_ext_remove()
      plot_id <- paste(geneset_name, rank_name_nice, sep = "_")


      .fgsea_res <- fgsea_res_list[[rank_name]]
      p <- plot_table(
        fgsea_res = .fgsea_res,
        ranks = rankobj,
        pathways = genesets,
        gsea_param = 0.5,
      )

      # outpath <- file.path(ENPLOT_BASEDIR, make.names(geneset_name), make.names(rank_name_nice))
      # if (!fs::dir_exists(outpath)) fs::dir_create(outpath)

      # print(p)
      return(p)

      # c(".png", ".pdf") %>% purrr::walk(
      # ~ggsave(filename = paste0(make.names(geneset_name), "_toptable", "_gseaparam.5", .x),
      #        path = outpath,
      #        dpi = 300,
      #        width = 19,
      #        height=12,
      #        )
      # )
    })
  })
}

plot_tables_faster <- function(
    results_list,
    ranks_list,
    pathways_list) {
  # Ensure names are aligned and iterate over them
  pathway_names <- names(results_list)

  # furrr::future_map(pathway_names, ~{
  purrr::map(pathway_names, ~ {
    geneset_name <- .x
    fgsea_res_list <- results_list[[geneset_name]]
    genesets <- pathways_list[[geneset_name]]

    # Now, iterate over each rank file within this pathway category
    # map(names(ranks_list), function(rank_name) {
    # furrr::future_map(names(ranks_list), function(rank_name) {
    furrr::future_imap(ranks_list, function(rankobj, rank_name) {
      # rankobj <- ranks_list[[rank_name]]
      # Generate a unique identifier for this plot/table, e.g., combining pathway name and rank file name
      rank_name_nice <- rank_name %>%
        fs::path_file() %>%
        fs::path_ext_remove()
      plot_id <- paste(geneset_name, rank_name_nice, sep = "_")


      .fgsea_res <- fgsea_res_list[[rank_name]]
      p <- plot_table(
        fgsea_res = .fgsea_res,
        ranks = rankobj,
        pathways = genesets,
        gsea_param = 0.5,
      )

      outpath <- util_tools$safe_subdir(
        ENPLOT_BASEDIR,
        util_tools$safe_path_component(geneset_name),
        util_tools$safe_path_component(rank_name_nice)
      )
      if (!fs::dir_exists(outpath)) fs::dir_create(outpath)

      # print(p)
      return(p)
      # c(".png", ".pdf") %>% purrr::walk(
      # ~ggsave(filename = paste0(make.names(geneset_name), "_toptable", "_gseaparam.5", .x),
      #        path = outpath,
      #        dpi = 300,
      #        width = 19,
      #        height=12,
      #        )
      # )
    })
  })
}
