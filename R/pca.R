# pca.R
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

src_dir <- file.path(here("R"))

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)


plot_tools <- new.env()
source(file.path(src_dir, "./plot.R"), local = plot_tools)


plot_utils <- new.env()
source(file.path(src_dir, "./plot_utils.R"), local = plot_utils)

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)

make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)

name_cleaner <- function(df){
  df %>%
    dplyr::mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "KEGG_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOMF_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "REACTOME_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOBP_")) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOCC_")) %>%
    dplyr::mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
    dplyr::mutate(pathway = str_wrap(pathway, 20))
}

# ==

#' handle 1 long form gsea result table
#'
#' This function runs PCA on a single GSEA result table
#' with columns "pathway", "pval", "padj", "ES", "NES", size, leadingEdge, mainpathway (logical), and rankname
#' pivot on rankname, using NES as the value
#'
#' @param gsea_object
#' @param metadata The metadata to use for the PCA
#' @return The mean of the numeric vector. If the input vector has zero length, the result is NA.
#' @examples
#' # Basic usage
#' gsea_object <- data.frame()
#'
#' # Handling NA values
#' mean_numeric(c(1, 2, NA, 4, 5), na.rm = TRUE)
#'
#' @export
do_one <- function(
    gsea_object,
    metadata = NULL,
    main_pathway_ratio = 0.1,
    ...) {
  # if ("mainpathway" %in% colnames(gsea_object)) {
  #   gsea_object <- gsea_object %>% filter(mainpathway == TRUE)
  # }

  log_msg(msg = "running pca")

  required_cols <- c("pathway", "NES", "rankname")
  for (col in required_cols) {
    if (!(col %in% colnames(gsea_object))) {
      stop(paste0(col, " column not found in the input data"))
    }
  }

  if (is.null(metadata)) {
    log_msg(msg = "metadata is null, making standin")
    metadata <- data.frame(id = unique(gsea_object$rankname))
    rownames(metadata) <- metadata$id
    metadata$dummy <- "a" ## ?
    log_msg(msg = paste0("metadata is \n", metadata))
  }

  gsea_object <- fgsea_tools$filter_on_mainpathway(gsea_object,
    main_pathway_ratio = main_pathway_ratio
  )
  if (nrow(gsea_object) == 0){return()}

  # clean names
  gsea_object <- gsea_object %>% name_cleaner()
    # dplyr::mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
    # dplyr::mutate(pathway = str_remove(pathway, "KEGG_")) %>%
    # dplyr::mutate(pathway = str_remove(pathway, "GOMF_")) %>%
    # dplyr::mutate(pathway = str_remove(pathway, "REACTOME_")) %>%
    # dplyr::mutate(pathway = str_remove(pathway, "GOBP_")) %>%
    # dplyr::mutate(pathway = str_remove(pathway, "GOCC_")) %>%
    # dplyr::mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
    # dplyr::mutate(pathway = str_wrap(pathway, 20))

  wide_df <- gsea_object %>%
    pivot_wider(id_cols = pathway, values_from = NES, names_from = rankname) %>%
    as.data.frame()
  # set pathway as rowname, remove from columns
  rownames(wide_df) <- wide_df$pathway
  wide_df$pathway <- NULL
  wide_df[is.na(wide_df)] <- 0 #min(wide_df, na.rm = TRUE) # zero because its an NES of zero, not down

  print(colnames(wide_df))
  print(rownames(metadata))
  pca_res <- wide_df %>% PCAtools::pca(metadata = metadata[colnames(wide_df), ])

  return(pca_res)
}

do_all <- function(
    gsea_objects,
    metadata = NULL) {
  pca_objects <- gsea_objects %>%
    purrr::map(~ do_one(.x, metadata = metadata))
  names(pca_objects) <- names(gsea_objects)
  return(pca_objects)
}

plot_biplot <- function(
    pca_object,
    top_pc = 3,
    showLoadings = T,
    labSize = 1.8,
    pointSize = 3,
    sizeLoadingsNames = 2,
    colby = NULL, # or a string like 'group'
    shape = NULL, # or a string like 'group'
    encircle = !is.null(colby),
    title = "",
    ...) {
  args <- list(...)
  if ("save_func" %in% names(args)) {
    save_func <- args$save_func
  } else {
    save_func <- NULL
  }
  if (!is.logical(encircle) || length(encircle) != 1 || is.na(encircle)) {
    stop("`encircle` must be a single logical value (TRUE or FALSE). Received: ", encircle)
  }

  vec <- paste0("PC", 1:top_pc)
  # vec <- c("PC1", "PC2", "PC3") # "PC4")
  pcs <- combn(vec, 2) %>%
    as.data.frame() %>%
    as.list()
  #

    log_msg(info=paste0('colby is: ', colby))
    log_msg(info=paste0('metadata is : ', pca_object$metadata))
  if (!is.null(colby) &&
    !is.null(pca_object$metadata) &&
    !colby %in% colnames(pca_object$metadata)) {
    warning(paste0(colby, " not found in metadata"))
    colby <- NULL
    encircle <- F
  }

  if (!is.null(shape) &&
    !is.null(pca_object$metadata) &&
    !shape %in% colnames(pca_object$metadata)) {
    warning(paste0(shape, " not found in metadata"))
    shape <- NULL
  }

  n_features <- nrow(pca_object$loadings)
  footnote <- paste0("n = ", as.character(n_features))


  plts <- pcs %>%
    purrr::map(
      ~ {
        # stopifnot(~COLBY%in%colnames(.metadata))
        .x1 <- .x[[1]]
        .x2 <- .x[[2]]

        if (!.x1 %in% names(pca_object$rotated)) {
          warning("not enough PCs")
          return()
        }

        if (!.x2 %in% names(pca_object$rotated)) {
          warning("not enough PCs")
          return()
        }

        .max_x <- max(
          pca_object$rotated[[.x1]] %>% abs(),
          pca_object$rotated[[.x2]] %>% abs()
        )
        .max_y <- .max_x

        plt <- PCAtools::biplot(
          pca_object,
          x = .x1,
          y = .x2,
          showLoadings = showLoadings,
          labSize = labSize,
          pointSize = pointSize,
          alphaLoadingsArrow = 0.5,
          sizeLoadingsNames = sizeLoadingsNames,
          colby = colby,
          shape = shape,
          drawConnectorsLoadings = TRUE,
          boxedLoadingsNames = TRUE,
          fillBoxedLoadings = "#ededed99",
          lengthLoadingsArrowsFactor = 1.1,
          # shape="source",
          legendPosition = "right",
          encircle = encircle,
          title = title,
          max.overlaps = Inf,
          maxoverlapsConnectors = Inf,
          ntopLoadings = 5,
        ) +
          coord_fixed(ratio = 1) +
          xlim(-.max_x, .max_x) +
          ylim(-.max_y, .max_y) +
          labs(caption = footnote) +
          geom_hline(yintercept = 0, color = "grey50", show.legend = NA) +
          geom_vline(xintercept = 0, color = "grey50", show.legend = NA) +
          colorspace::scale_color_discrete_qualitative(palette = "Dynamic") +
          colorspace::scale_fill_discrete_qualitative(palette = "Dynamic") +
          scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14)) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )

        if (!is.null(save_func)) {
          current_args <- get_args(save_func)
          filename <- current_args$filename
          if (is.null(filename)) {
            filename <- paste0("pca_biplot_")
          }
          filename <- paste0(filename, "_", .x1, "_", .x2)
          save_func(plot_code = function() print(plt), filename = filename)
        } else {
          print(plt)
        }
        return(plt)
      } # exit inner
    ) # exit outer
  return(plts)
}

plot_all_biplots <- function(
    pca_objects,
    top_pc = 3,
    showLoadings = T,
    labSize = 3,
    pointSize = 4,
    sizeLoadingsNames = 1.75,
    colby = "group",
    shape = NULL,
    fig_width = 8.4,
    fig_height = 7.6,
    save_func = NULL,
    ...) {
  log_msg(debug=paste0('colby equal to : ', colby))
  pca_objects %>%
    purrr::imap(
      ~ {
        pca_object <- .x
        title <- .y
        collection_name <- .y

        plts <- colby %>% purrr::map(
            ~{
              .colby <- .x

              # Additional validation
              if (!is.character(.colby) || length(.colby) != 1) {
                stop(paste0("Invalid .colby value: ", .colby))
              }

              log_msg(msg=paste0('.colby equal to : ', .colby))
              if (!is.null(save_func)) {
              collection_dir <- util_tools$safe_path_component(collection_name)
              save_func <- make_partial(save_func,
                filename = util_tools$safe_filename(
                  "pca_biplot",
                  title,
                  paste0("col", .colby),
                  fallback = "pca_biplot"
                ),
                path = util_tools$safe_subdir(get_arg(save_func, "path"), collection_dir, "pca"),
                width = fig_width, height = fig_height
              )
              }
              p <- tryCatch(
              {
                  plot_biplot(
                    pca_object,
                    top_pc = top_pc,
                    showLoadings = showLoadings,
                    labSize = labSize,
                    pointSize = pointSize,
                    sizeLoadingsNames = sizeLoadingsNames,
                    colby = .colby,
                    shape = shape,
                    title = title,
                    save_func = save_func
                    # ...
                )
              }, error = function(e) {
                log_msg(msg = paste0("error in plot_all_biplots: ", e))
              })
              return(p)
          })
          return(plts)
        }
    )
}

get_top_loadings <- function(pcaobj,
                             components = c("PC1", "PC2", "PC3", "PC4"),
                             rangeRetain = 0.05,
                             limit = Inf
                             ) { # directly lifted from https://github.com/kevinblighe/PCAtools/blob/50f240ba76fe07e2a546998e76dc8aa1c3570693/R/plotloadings.R#L7 and edited
  components <- intersect(components, pcaobj$components)
  x <- pcaobj$loadings[, components, drop = FALSE]
  membership <- list()
  retain <- c()
  # could this be rewritten with tidy group_by and arrange and slice top largest diff?
  for (i in seq_along(components)) {
    # for each PC, based on the loadings range, calculate the rangeRetain
    # fraction of the range
    offset <- (max(x[, i]) - min(x[, i])) * rangeRetain

    # to obtain upper and lower cut-offs, offset max and min by the offset
    uppercutoff <- max(x[, i]) - offset
    lowercutoff <- min(x[, i]) + offset

    # build a consensus list of variables that will be included
    retain_vals <- c(
      which(x[, i] >= uppercutoff),
      which(x[, i] <= lowercutoff)
    ) %>% unique()
    if (limit < Inf) {
      retain_vals <- order(abs(x[retain_vals, i]), decreasing=T) %>% head(n=limit)
    }

    retain <- unique(c(
      retain,
      retain_vals
    ))

    membership[[paste0("PC", i, "_member")]] <- retain_vals
  }
  membership_df <- stack(membership)
  membership_df[["var"]] <- rownames(x)[membership_df$values]
  membership_df[["boolvalue"]] <- !is.na(membership_df$var)
  membership_dfw <- membership_df %>% pivot_wider(id_cols = c("values", "var"), names_from = "ind", values_from = "boolvalue")
  membership_dfw[is.na(membership_dfw)] <- FALSE

  message("-- variables retained:")
  message(paste0(rownames(x)[retain], collapse = ", "))
  x_final <- x[retain, , drop = FALSE]
  # create a plot object (2-col df) of PC names and
  # explained variance of each
  x_final <- data.frame(rownames(x_final), x_final[, components, drop = FALSE])
  colnames(x_final)[1] <- "var"

  x_final <- x_final %>% left_join(membership_dfw)
  x_final$values <- NULL
  x_final %<>% as.data.frame
  rownames(x_final) <- x_final$var
  return(x_final)
}


make_enplots_from_loadings <- function(){
  stop("not implemented")
}

make_heatmap_from_loadings <- function(
    gsea_object,
    pca_object,
    components = c("PC1", "PC2", "PC3", "PC4"),
    save_func = NULL,
    cut_by = NULL,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    meta_to_include = NULL,
    meta_to_exclude = NULL,
    rankname_order = NULL,
    ...){


  components <- intersect(components, pca_object$components) # there may be less than 4 components
  top_loadings <- get_top_loadings(pca_object, components = components, rangeRetain = 0.05, limit = Inf)
  submat <- gsea_object %>% name_cleaner() %>% dplyr::filter(pathway %in% rownames(top_loadings))

  ra <- ComplexHeatmap::rowAnnotation(
    PC1 = top_loadings$PC1_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC2 = top_loadings$PC2_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC3 = top_loadings$PC3_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC4 = top_loadings$PC4_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white"))
  )


  maybe_metadata <- pca_object$metadata
  if ("dummy" %in% colnames(maybe_metadata)) {
    maybe_metadata <- NULL
  }

  param_grid <- expand.grid(
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    cut_by = unique(c(cut_by, NA)),
    stringsAsFactors = FALSE
  )

  param_grid %>% purrr::pmap(~{

    params = list(...)
    .cut_by <- params$cut_by
    .cut_by_val <- plot_utils$process_cut_by(.cut_by, pca_object$metadata)
    cut_by_label <- if (!is.null(.cut_by_val)) {
      paste0("cut_", util_tools$safe_path_component(.cut_by, max_chars = 32))
    } else {
      ""
    }
    cluster_rows <- params$cluster_rows
    cluster_columns <- params$cluster_columns

    .save_func <- NULL
    if (!is.null(save_func)) {
      filename <- util_tools$safe_filename(
        get_arg(save_func, "filename"),
        "gsea_heatmap_5pct",
        util_tools$cluster_flag_token(cluster_rows, "rc"),
        util_tools$cluster_flag_token(cluster_columns, "cc"),
        cut_by_label,
        fallback = "gsea_heatmap"
      )
        .save_func <- make_partial(save_func, filename = filename)
      }

    tryCatch({
      ht <- plot_tools$plot_results_one_collection(
        df = submat,
        metadata = pca_object$metadata, # guaranteed to match the colnames of df which was used originally to make the pca obj
        # pathway_metadata = loadings,
        row_annotation = ra,
        limit = Inf,
        pstat_cutoff = 1, # we want no filtering
        title = "Top 5% Loadings",
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        rankname_order = rankname_order,
        cut_by = .cut_by,
        save_func = .save_func,
      )}, error = function(msg) { log_msg(error=msg) }
    ) # end of tryCatch
  }) # end of pmap


  top_loadings <- get_top_loadings(pca_object, components = components, rangeRetain = 1, limit = 5)
  submat <- gsea_object %>% name_cleaner() %>% dplyr::filter(pathway %in% rownames(top_loadings))
  submat %<>% mutate(pathway = factor(pathway, levels = rownames(top_loadings), ordered=TRUE)) %>% arrange(pathway)

  ra <- ComplexHeatmap::rowAnnotation(
    PC1 = top_loadings$PC1_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC2 = top_loadings$PC2_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC3 = top_loadings$PC3_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC4 = top_loadings$PC4_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white"))
  )


  param_grid %>% purrr::pmap(~{

    params = list(...)
    .cut_by <- params$cut_by
    .cut_by_val <- plot_utils$process_cut_by(.cut_by, pca_object$metadata)
    cut_by_label <- if (!is.null(.cut_by_val)) {
      paste0("cut_", util_tools$safe_path_component(.cut_by, max_chars = 32))
    } else {
      ""
    }
    cluster_rows <- params$cluster_rows
    cluster_columns <- params$cluster_columns

    .save_func <- NULL
    if (!is.null(save_func)) {
      filename <- util_tools$safe_filename(
        get_arg(save_func, "filename"),
        "gsea_heatmap_top5",
        util_tools$cluster_flag_token(cluster_rows, "rc"),
        util_tools$cluster_flag_token(cluster_columns, "cc"),
        cut_by_label,
        fallback = "gsea_heatmap"
      )
        .save_func <- make_partial(save_func, filename = filename)
      }


    tryCatch({
      ht <- plot_tools$plot_results_one_collection(
        df = submat,
        metadata = pca_object$metadata, # guaranteed to match the colnames of df which was used originally to make the pca obj
        # pathway_metadata = loadings,
        row_annotation = ra,
        limit = Inf,
        pstat_cutoff = 1, # we want no filtering
        title = "Top 5 Loadings Per PC",
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        cut_by = .cut_by,
        save_func = .save_func,
        # rankname_order = params$extra$rankname_order
        rankname_order = rankname_order,
      )}, error = function(msg) { log_msg(error=msg) }
    ) # end of tryCatch
  }) # end of pmap



  # param_grid %>% purrr:pmap(~{
  #   params <- list(...)
  #   tryCatch({
  #     ht <- plot_tools$plot_results_one_collection(
  #       df = submat,
  #       metadata = maybe_metadata, # guaranteed to match the colnames of df which was used originally to make the pca obj
  #       # pathway_metadata = loadings,
  #       row_annotation = ra,
  #       limit = Inf,
  #       pstat_cutoff = 1, # we want no filtering
  #       save_func = save_func,
  #       title = "Top 5 Loadings Per PC",
  #       ...
  #       )}, error = function(msg){ log_msg(error=msg) }
  #   ) # end of tryCatch
  # })

}

plot_gene_loadings_heatmap <- function(
    pca_object,
    components,
    top_n = 25,
    save_func = NULL,
    gct = NULL,
    cluster_rows = TRUE,
    cluster_columns = c(FALSE, TRUE),
    cut_by = NULL) {
  available <- intersect(components, colnames(pca_object$loadings))
  if (length(available) == 0) {
    log_msg(warning = "Gene PCA heatmap skipped: requested components missing from loadings matrix.")
    return(NULL)
  }

  loadings_mat <- pca_object$loadings[, available, drop = FALSE]
  loadings_mat <- loadings_mat[apply(loadings_mat, 1, function(x) any(is.finite(x))), , drop = FALSE]
  if (nrow(loadings_mat) == 0) {
    log_msg(warning = "Gene PCA heatmap skipped: no finite loadings available.")
    return(NULL)
  }

  gene_subset <- unique(unlist(lapply(seq_along(available), function(idx) {
    pc <- available[[idx]]
    ord <- order(abs(loadings_mat[, pc]), decreasing = TRUE)
    head(rownames(loadings_mat)[ord], min(top_n, length(ord)))
  })))

  if (length(gene_subset) == 0) {
    log_msg(warning = "Gene PCA heatmap skipped: no genes selected for loadings heatmap.")
    return(NULL)
  }

  if (is.null(gct)) {
    log_msg(warning = "Gene PCA heatmap skipped: expression GCT not supplied.")
    return(NULL)
  }

  missing_genes <- setdiff(gene_subset, gct@rid)
  if (length(missing_genes) > 0) {
    log_msg(info = paste0(
      "Gene PCA heatmap: dropping ", length(missing_genes),
      " genes absent from expression matrix: ",
      paste(missing_genes, collapse = ", ")
    ))
  }

  ordered_genes <- gene_subset[gene_subset %in% gct@rid]
  if (length(ordered_genes) == 0) {
    log_msg(warning = "Gene PCA heatmap skipped: selected genes missing from expression matrix.")
    return(NULL)
  }

  gct_subset <- tryCatch(
    {
      cmapR::subset_gct(gct, rid = ordered_genes)
    },
    error = function(e) {
      log_msg(error = paste0("Gene PCA heatmap subset failed: ", conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(gct_subset)) {
    return(NULL)
  }

  loadings_subset <- loadings_mat[ordered_genes, available, drop = FALSE]
  membership_raw <- apply(loadings_subset, 1, function(row) {
    if (all(!is.finite(row))) {
      return(NA_character_)
    }
    available[[which.max(abs(row))]]
  })
  membership_levels <- available
  include_none <- any(is.na(membership_raw))
  membership_factor <- if (include_none) {
    factor(ifelse(is.na(membership_raw), "none", membership_raw), levels = c(membership_levels, "none"))
  } else {
    factor(membership_raw, levels = membership_levels)
  }

  palette <- colorspace::qualitative_hcl(length(membership_levels), palette = "Dark 3")
  if (length(palette) < length(membership_levels)) {
    palette <- grDevices::rainbow(length(membership_levels))
  }
  names(palette) <- membership_levels
  if (include_none) {
    palette <- c(palette, none = "#bfbfbf")
  }

  row_annotation <- ComplexHeatmap::rowAnnotation(
    PC = membership_factor,
    col = list(PC = palette),
    show_annotation_name = FALSE
  )

  gct_subset@rdesc$gene_pca_component <- as.character(membership_factor)

  cluster_rows_vals <- as.logical(cluster_rows)
  cluster_rows_vals <- cluster_rows_vals[!is.na(cluster_rows_vals)]
  if (length(cluster_rows_vals) == 0) {
    cluster_rows_vals <- TRUE
  }
  cluster_columns_vals <- as.logical(cluster_columns)
  cluster_columns_vals <- cluster_columns_vals[!is.na(cluster_columns_vals)]
  if (length(cluster_columns_vals) == 0) {
    cluster_columns_vals <- FALSE
  }

  cut_by_value <- cut_by
  if (is.character(cut_by_value)) {
    cut_by_value <- cut_by_value[nzchar(cut_by_value)]
  }
  if (length(cut_by_value) > 0) {
    cut_by_value <- cut_by_value[[1]]
  } else {
    cut_by_value <- NULL
  }
  if (is.logical(cut_by_value) && !isTRUE(cut_by_value)) {
    cut_by_value <- NULL
  }

  row_title <- paste0("Top Loadings (", paste(available, collapse = ", "), ")")
  cut_label <- if (!is.null(cut_by_value)) {
    paste0("cut_", util_tools$safe_path_component(cut_by_value, max_chars = 32))
  } else {
    NULL
  }

  param_grid <- expand.grid(
    cluster_rows = unique(cluster_rows_vals),
    cluster_columns = unique(cluster_columns_vals),
    stringsAsFactors = FALSE
  )

  for (idx in seq_len(nrow(param_grid))) {
    cluster_rows_flag <- isTRUE(param_grid$cluster_rows[[idx]])
    cluster_columns_flag <- isTRUE(param_grid$cluster_columns[[idx]])
    filename <- util_tools$safe_filename(
      "gene_loadings_heatmap",
      util_tools$cluster_flag_token(cluster_rows_flag, "rc"),
      util_tools$cluster_flag_token(cluster_columns_flag, "cc"),
      cut_label
    )
    save_target <- if (!is.null(save_func)) make_partial(save_func, filename = filename) else NULL

    plot_tools$make_heatmap_fromgct(
      gct = gct_subset,
      row_title = row_title,
      column_title = "Scaled expression",
      save_func = save_target,
      cluster_rows = cluster_rows_flag,
      cluster_columns = cluster_columns_flag,
      colorbar_title = "zscore",
      cut_by = cut_by_value,
      row_annotation = row_annotation
    )
  }

  invisible(gct_subset)
}

run_gene_pca_pipeline <- function(gct, params, savedir, replace = TRUE) {
  if (is.null(gct) || ncol(gct@mat) < 2) {
    log_msg(warning = "Gene PCA skipped: GCT is NULL or has fewer than 2 samples.")
    return(NULL)
  }

  expr <- gct@mat
  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }

  sample_ids <- colnames(expr)
  metadata <- as.data.frame(gct@cdesc, stringsAsFactors = FALSE)
  metadata <- metadata[sample_ids, , drop = FALSE]
  rownames(metadata) <- sample_ids

  variance <- apply(expr, 1, stats::var, na.rm = TRUE)
  keep_genes <- which(!is.na(variance) & variance > 0)
  if (length(keep_genes) < 2) {
    log_msg(warning = "Gene PCA skipped: fewer than two genes with non-zero variance.")
    return(NULL)
  }
  expr <- expr[keep_genes, , drop = FALSE]
  gct_filtered <- cmapR::subset_gct(gct, rid = rownames(expr))

  expr_scaled <- t(scale(t(expr), center = TRUE, scale = TRUE))
  expr_scaled <- as.matrix(expr_scaled)
  mask_non_finite <- !is.finite(expr_scaled)
  if (any(mask_non_finite)) {
    log_msg(info = paste0("Gene PCA: replacing ", sum(mask_non_finite), " non-finite scaled values with 0."))
    expr_scaled[mask_non_finite] <- 0
  }
  gct_scaled <- gct_filtered
  gct_scaled@mat <- expr_scaled

  metadata_colors <- params$metadata_color
  valid_colors <- metadata_colors[metadata_colors %in% colnames(metadata)]
  if (length(valid_colors) < length(metadata_colors)) {
    missing <- setdiff(metadata_colors, valid_colors)
    if (length(missing) > 0) {
      log_msg(warning = paste0("Gene PCA: metadata_color entries not found; ignoring: ", paste(missing, collapse = ", ") ))
    }
  }
  if (length(valid_colors) == 0) {
    valid_colors <- NA_character_
  }

  metadata_shape <- params$metadata_shape
  if (!nzchar(metadata_shape) || !metadata_shape %in% colnames(metadata)) {
    if (nzchar(params$metadata_shape)) {
      log_msg(warning = paste0("Gene PCA: metadata_shape '", params$metadata_shape, "' not found; ignoring."))
    }
    metadata_shape <- NULL
  }

  pca_res <- PCAtools::pca(expr_scaled,
    metadata = metadata,
    center = FALSE,
    scale = FALSE
  )

  pc_count <- max(2, min(params$components, length(pca_res$components)))
  components <- paste0("PC", seq_len(pc_count))

  base_dir <- util_tools$safe_subdir(savedir, "pca_gene")

  biplot_save <- util_tools$make_partial(
    plot_utils$plot_and_save,
    path = util_tools$safe_subdir(base_dir, "biplots"),
    replace = replace,
    width = params$width,
    height = params$height
  )

  for (color_col in valid_colors) {
    colour_label <- if (is.na(color_col)) "none" else color_col
    biplot_save_col <- make_partial(
      biplot_save,
      filename = util_tools$safe_filename("gene_pca_biplot", colour_label)
    )
    plot_biplot(
      pca_res,
      top_pc = pc_count,
      showLoadings = TRUE,
      labSize = params$labSize,
      pointSize = params$pointSize,
      sizeLoadingsNames = params$sizeLoadingsNames,
      colby = if (is.na(color_col)) NULL else color_col,
      shape = metadata_shape,
      title = paste0("Gene Expression PCA", if (is.na(color_col)) "" else paste0(" - ", color_col)),
      save_func = biplot_save_col
    )
  }

  if (isTRUE(params$heatmap)) {
    heatmap_cluster_rows <- params$cluster_rows %||% TRUE
    heatmap_cluster_columns <- params$cluster_columns %||% c(FALSE, TRUE)
    heatmap_cut_by <- util_tools$process_cut_by(params$cut_by %||% NULL)

   # requested_rank_order <- params$extra$rankname_order
   # available_ranknames <- unique(df$rankname)

   # if (!is.null(requested_rank_order)) {
   #   missing_ranknames <- setdiff(requested_rank_order, available_ranknames)
   #   if (length(missing_ranknames) > 0) {
   #     log_msg(warning = paste0(
   #       "rankname_order entries not present in results: ",
   #       paste(missing_ranknames, collapse = ", "),
   #       ". These names will be ignored."
   #     ))
   #   }

   #   duplicated_ranknames <- requested_rank_order[duplicated(requested_rank_order)]
   #   if (length(duplicated_ranknames) > 0) {
   #     log_msg(warning = paste0(
   #       "rankname_order contains duplicate labels: ",
   #       paste(unique(duplicated_ranknames), collapse = ", "),
   #       ". Keeping the first occurrence of each."
   #     ))
   #   }
   #   rankname_order <- available_ranknames
   #   log_msg(debug = "rankname_order not supplied; using detected order from results")
   # } else {
   #   rankname_order <- intersect(rankname_order, available_ranknames)
   #   if (length(rankname_order) == 0) {
   #     warning("rankname_order is empty, using default")
   #     rankname_order <- available_ranknames # default backup
   #   } else if (!is.null(requested_rank_order)) {
   #     if (length(rankname_order) < length(unique(requested_rank_order))) {
   #       log_msg(info = paste0(
   #         "Applied rankname_order for ",
   #         length(rankname_order),
   #         " comparisons (", length(unique(requested_rank_order)), " requested)."
   #       ))
   #     }
   #   }
   #   pca_res <- pca_res[, rankname_order]

    # if (is.character(heatmap_cut_by)) {
    #   heatmap_cut_by <- heatmap_cut_by[nzchar(heatmap_cut_by)]
    # }
    # if (length(heatmap_cut_by) > 0) {
    #   heatmap_cut_by <- heatmap_cut_by[[1]]
    # } else {
    #   heatmap_cut_by <- NULL
    # }
    # if (is.logical(heatmap_cut_by) && !isTRUE(heatmap_cut_by)) {
    #   heatmap_cut_by <- NULL
    # }

    heatmap_save <- util_tools$make_partial(
      plot_utils$plot_and_save,
      path = util_tools$safe_subdir(base_dir, "heatmaps"),
      replace = replace,
      width = 7,
      height = max(4, params$top_loadings * 0.22 + 2)
    )
    plot_gene_loadings_heatmap(
      pca_res,
      components = components,
      top_n = params$top_loadings,
      save_func = heatmap_save,
      gct = gct_scaled,
      cluster_rows = heatmap_cluster_rows,
      cluster_columns = heatmap_cluster_columns,
      cut_by = heatmap_cut_by
    )
  }

  tables_dir <- util_tools$safe_subdir(base_dir, "tables")
  fs::dir_create(tables_dir, recurse = TRUE)
  loadings_tbl <- as.data.frame(pca_res$loadings[, components, drop = FALSE])
  loadings_tbl <- cbind(gene = rownames(loadings_tbl), loadings_tbl)
  scores_tbl <- as.data.frame(pca_res$rotated[, components, drop = FALSE])
  scores_tbl <- cbind(sample = rownames(scores_tbl), scores_tbl)
  readr::write_tsv(loadings_tbl, fs::path(tables_dir, "gene_pca_loadings.tsv"))
  readr::write_tsv(scores_tbl, fs::path(tables_dir, "gene_pca_scores.tsv"))

  log_msg(info = paste0("Gene PCA outputs written to ", base_dir))

  invisible(pca_res)
}
