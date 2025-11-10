suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(furrr))

src_dir <- file.path(here("R"))

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)

geneset_tools <- new.env()
source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)

db_tools <- new.env()
source(file.path(src_dir, "./db.R"), local = db_tools)

map_tools <- new.env()
source(file.path(here("R"), "./map.R"), local = map_tools)

util_tools <- new.env()
source(file.path(src_dir, "utils.R"), local = util_tools)
log_msg <- util_tools$make_partial(util_tools$log_msg)



filter_on_mainpathway <- function(
    pathway_object,
    main_pathway_ratio = .1) {
  if (!"rankname" %in% colnames(pathway_object)) {
    stop("rankname column not found in the input data")
  }

  if (!"mainpathway" %in% colnames(pathway_object)) {
    pathway_object$mainpathway <- TRUE # will simply keep all
  }


  pathway_object <- pathway_object %>%
    group_by(rankname) %>%
    mutate(n_main = sum(mainpathway == T)) %>%
    mutate(ratio_main = n_main / n()) %>%
    ungroup()

  # Keep rows where the per-rankname fraction of TRUE mainpathway
  # meets the provided threshold.
  filtered_pathway_object <- pathway_object %>%
    dplyr::filter(ratio_main >= main_pathway_ratio)

  return(filtered_pathway_object)
}

#' Select top pathways based on statistical cutoff
#'
#' This function selects the top pathways from a data frame based on a statistical cutoff.
#' The statistical cutoff is determined by the p-value or adjusted p-value, depending on the value of the `pstat_usetype` parameter.
#' The function arranges the data frame in descending order of the absolute value of the NES column,
#' filters the rows where the p-value or adjusted p-value is less than the specified cutoff,
#' selects the specified number of top pathways, and returns a subset of the data frame containing only the top pathways.
#'
#' @param df A data frame containing the pathways and statistical values.
#' @param pstat_cutoff The cutoff value for the p-value or adjusted p-value.
#' @param limit The maximum number of top pathways to select.
#' @param pstat_usetype The type of statistical value to use for the cutoff (either "pval" or "padj").
#'
#' @return A subset of the input data frame containing only the top pathways.
#'
#' @examples
#' df <- data.frame(
#'   pathway = c("Pathway A", "Pathway B", "Pathway C"),
#'   NES = c(1.5, -2.3, 0.8),
#'   pval = c(0.01, 0.05, 0.001),
#'   padj = c(0.05, 0.1, 0.01)
#' )
#' select_topn(df, pstat_cutoff = 0.05, limit = 2, pstat_usetype = "pval")
#' # Returns a data frame with Pathway A and Pathway C
#' select_topn(df, pstat_cutoff = 0.05, limit = 2, pstat_usetype = "padj")
#' # Returns a data frame with Pathway A and Pathway B
#'
#' @seealso \code{\link{arrange}}, \code{\link{filter}}, \code{\link{slice_head}}, \code{\link{pull}}
#'
#' @import dplyr
#'
#' @export
select_topn <- function(df,
                        pstat_cutoff = 1,
                        limit = 120,
                        pstat_usetype = c("padj", "pval"),
                        to_include = NULL, # extra pathways to explictly include
                        ...) {
  pstat_usetype <- match.arg(pstat_usetype)

  if (!"data.frame" %in% class(df)) {
    stop(
      cat("df should be a data frame")
    )
  }

  if (!"NES" %in% colnames(df)) {
    stop(
      cat("NES should be a column in df")
    )
  }

  if (!pstat_usetype %in% colnames(df)) {
    stop(
      cat(paste0(pstat_usetype, " should be a column in df"))
    )
  }

  top_pathways <- df %>%
    arrange(-abs(NES)) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    filter(!!as.symbol(pstat_usetype) < pstat_cutoff) %>%
    slice_head(n = limit) %>%
    pull(pathway)
  if (!is.null(to_include)) {
    top_pathways <- union(top_pathways, to_include)
  }
  subdf <- df %>% filter(pathway %in% top_pathways)
  return(subdf)
}

fgsea_cache_manager <- function(
    rankobj,
    geneset,
    minSize = 15,
    maxSize = 500,
    collapse = FALSE,
    cache = TRUE,
    cache_dir = NULL,
    logger = NULL,
    final_result = NULL,
    do_load = is.null(final_result),
    save = !is.null(final_result),
    ...) {
  if (is.null(logger)) logger <- log_msg
  logger(msg=paste0("cache dir set to ", cache_dir))

  get_hash_val <- function() {
    rlang::hash(
      c(
        as.character(rankobj),
        as.character(geneset),
        minSize,
        maxSize,
        collapse
      )
    )
  }
  # print(cache_dir)
  hashval <- get_hash_val()
  if (do_load) {
    cache_load <- io_tools$load_from_cache(hashval, cache_dir = cache_dir)
    if (!is.null(cache_load)) {
      logger(msg = paste0("cache hit: ", hashval))
      return(cache_load)
    } else {
      return(NULL)
    }
  } else if (save) {
    io_tools$write_to_cache(object = final_result, filename = hashval, cache_dir = cache_dir)
  }
}

run_one <- function(
    rankobj = NULL,
    geneset = NULL,
    minSize = 15,
    maxSize = 500,
    collapse = FALSE,
    logger = NULL,
    ...) {
  if (is.null(logger)) logger <- log_msg
  logger(msg = paste0("starting run_one"))

  # to look for duplicate gene names
  # rankobj %>% names %>% table %>% as.data.frame %>% pull(Freq) %>% max
  .maxcounts <- rankobj %>%
    names() %>%
    table() %>%
    as.data.frame() %>%
    pull(Freq) %>%
    max()
  assertthat::are_equal(.maxcounts, 1)
  # this doesn't actually error out??
  # set.seed(789)

  nperm <- 5000
  #nperm_max <- 100000

  fgseaRes <- tryCatch(
    {
      # Attempt to run fgsea
      do_run <- TRUE
      #while (do_run){
      fgseaRes <- fgsea(
        pathways = geneset,
        stats = rankobj,
        minSize = minSize,
        maxSize = maxSize,
        nPermSimple=nperm,
      )
      n_fail <- fgseaRes %>% dplyr::filter(is.na(NES))
      #if ((nrow(n_fail) == 0) | (nperm > nperm_max)) do_run <- FALSE
      #nperm <- nperm * 2
      #}
      return(fgseaRes)  # Return the result if successful
    },
    error = function(e) {
      cat("Error in FGSEA: ", e$message, "\n")
      return(NULL)
    }
  )
  if (is.null(fgseaRes)) {
    return(NULL)
  }

  if (length(collapse) != 1) { # ??
    # stop("Expected a single logical value for 'collapse'")
    collapse <- collapse[[1]]
  }

  fgseaRes$mainpathway <- TRUE
  if (!is.null(collapse) && collapse) {
    cat("finding main pathways")
    logger(msg = "finding main pathways")
    collapse_results <- fgseaRes %>%
      fgsea::collapsePathways(
        pathways = geneset,
        stats = rankobj,
        pval.threshold = 0.05,
        gseaParam = 1
      )
    fgseaRes <- fgseaRes %>% dplyr::mutate(
      mainpathway = pathway %in% collapse_results$mainPathways
    )
  } 

  return(fgseaRes)
}

run_all_rankobjs <- function(
    pathway,
    rankobjs,
    parallel = F,
    minSize = 15,
    maxSize = 500,
    collapse = FALSE,
    cache = TRUE,
    cache_dir = NULL,
    logger = NULL,
    species = "Homo sapiens",
    ...) {
  if (is.null(logger)) logger <- log_msg
  logger(msg = paste0("starting run_all_rankobjs"))
  logger(msg = paste0("parallel is ", parallel))
  # rankobjs %>% furrr::future_map( # maybe later
  if (parallel == TRUE) {
    workers <- future::availableCores() - 1
    options(future.globals.maxSize = 20000 * 1024^2)
    logger(msg = paste0("using ", workers, " workers"))
    .map_func <- furrr::future_map
  } else {
    .map_func <- purrr::map
  }

  fgsea_args <- list(
    geneset = pathway,
    minSize = minSize,
    maxSize = maxSize,
    collapse = collapse,
    cache = cache,
    cache_dir = cache_dir,
    logger = logger
  )


  if (!is.null(cache) && cache == TRUE) {
    logger(msg = "caching is enabled")
    cache_results <- rankobjs %>%
      purrr::map(~ do.call(fgsea_cache_manager, c(list(rankobj = .x), fgsea_args)) ) %>%
      Filter(Negate(is.null), .)

    .to_do <- setdiff(names(rankobjs), names(cache_results))
    rankobjs <- rankobjs[.to_do]

  } else {
    cache_results <- NULL
  }

  results <- rankobjs %>% .map_func(~ do.call(run_one, c(list(rankobj = .), fgsea_args)))
  results <- Filter(Negate(is.null), results)
  # print(paste0('length rankobjs : ', length(rankobjs)))
  # print(paste0('length results : ', length(results)))

  # species <- species %||% "Homo sapiens"
  if (length(results) > 0){ # only trigger for new results
    tryCatch(
      {
        results <- results %>% map_tools$add_leadingedges_to_results_list(., species=species)
        results <- purrr::map(results, ~{
          if (!is.null(.) && !"mainpathway" %in% colnames(.)) {
            dplyr::mutate(., mainpathway = TRUE)
          } else {
            .
          }
        })
      },
      error = function(e){
        print(paste0('error mapping names to genes, ', e))
      }
    )
  }


  # results <- rankobjs %>% .map_func(
  #   ~ run_one(.,
  #     fgsea_args
  #     # geneset = pathway,
  #     # minSize = minSize,
  #     # maxSize = maxSize,
  #     # collapse = collapse,
  #     # cache = cache,
  #     # cache_dir = cache_dir,
  #     # logger = logger
  #   )
  # )

  if (!is.null(cache) && cache == TRUE) {
    # Save only successfully computed (non-NULL) results.
    # Align by name to avoid positional length mismatches when some results are filtered out.
    if (length(results) > 0) {
      .common_names <- intersect(names(rankobjs), names(results))
      if (length(.common_names) > 0) {
        purrr::walk2(
          rankobjs[.common_names],
          results[.common_names],
          ~ do.call(
            fgsea_cache_manager,
            c(list(rankobj = .x, final_result = .y), fgsea_args)
          )
        )
      }
    }
  }

  final_results <- c(
    Filter(function(x) !is.null(x), cache_results),
    Filter(function(x) !is.null(x), results)
  )

  return(final_results)
}

run_all_pathways <- function(
    geneset_lists,
    ranks,
    parallel = FALSE,
    minSize = 15,
    maxSize = 500,
    genesets_additional_info = NULL,
    collapse = FALSE,
    cache = TRUE,
    cache_dir = NULL,
    logger = NULL,
    species = "Homo sapiens",
    ...) {
  if (any(is.null(names(geneset_lists)))) {
    stop(
      cat(
        "each geneset must be named",
        "e.g. geneset_lists <- list('H_' = geneset1, 'C5_GO:BP' = geneset2)"
      ),
      call. = FALSE
    )
  }

  if (is.null(logger)) logger <- log_msg

  if (is.null(cache_dir)) cache_dir <- here("cache")

  if (any(is.null(names(ranks)))) {
    stop(
      cat(
        "each rank must be named",
        "e.g. rank <- list(comparison1 = rank1, comparison2 = rank2)"
      ),
      call. = FALSE
    )
  }

  out <- geneset_lists %>% purrr::imap(
    ~ {
      geneset_list <- .x
      geneset_name <- .y

      current_collapse <- collapse # Use local variable for clarity

      # Check for non-null and appropriate columns of genesets_additional_info
      # this is messy but works
      if (!is.null(genesets_additional_info) && !collapse) {
        # Check for necessary columns in genesets_additional_info

        if (!"collection_name" %in% colnames(genesets_additional_info)) {
          genesets_additional_info <- genesets_additional_info %>%
            dplyr::mutate(collection_name = stringr::str_c(category, subcategory, sep = "_"))
        }

        if (!"collapse" %in% colnames(genesets_additional_info)) {
          genesets_additional_info <- genesets_additional_info %>%
            dplyr::mutate(collapse = F)
        }

        # Extract additional info specific to the current pathway
        geneset_additional_info <- genesets_additional_info[genesets_additional_info$collection_name == geneset_name, ]

        # Check if specific geneset info was found and update collapse if so
        if (nrow(geneset_additional_info) > 0) {
          current_collapse <- geneset_additional_info$collapse
        } else {
          logger(msg = paste0("No matching geneset info found for pathway: ", geneset_name))
          warning("No matching geneset info found for pathway: ", geneset_name)
        }
      }

      logger(msg = paste0("starting ", geneset_name))
      # Pass the potentially updated collapse value to the next function
      # could add cache check here instead
      results <- geneset_list %>% run_all_rankobjs(.,
        rankobjs = ranks,
        parallel = parallel,
        minSize = minSize,
        maxSize = maxSize,
        collapse = current_collapse,
        cache = cache,
        cache_dir = cache_dir,
        logger = logger,
        species = species,
        ...
      )
      return(results)
    }
  )
  return(out)
}


#' Get Rank Order Data for GSEA Pathways
#'
#' This function generates rank order data for a specified gene set and its
#' corresponding rank object. The function computes the rank order of each gene
#' in the `rankobj` and joins it with enrichment plot data (e.g., curve, ticks,
#' statistics). Optionally, additional gene information from `geneset_df` can be
#' joined to the result.
#'
#' @param geneset A character vector specifying the gene set of interest. It should
#'   contain gene IDs or symbols. A warning is issued if the input is of the wrong
#'   type.
#' @param rankobj A named numeric vector representing the ranking values for the
#'   genes. If a list is passed, the function will attempt to handle it but will
#'   issue a warning if it contains more than one element.
#' @param geneset_df Optional. A data frame that contains additional information
#'   about the genes in the gene set. It must have an `entrez_gene` column, which
#'   will be matched to the `id` column (the entrez ids of that pathway) in the rank order data.
#'
#' @details
#' This function performs the following tasks:
#' - Uses the `plotEnrichmentData()` function to generate enrichment data (curve,
#'   ticks, and statistics) for the specified `geneset` and `rankobj`.
#' - Computes rank order of genes by ranking the values in `rankobj`.
#' - Merges the rank order with enrichment data (curve and ticks).
#' - Optionally joins the resulting rank order data with `geneset_df` if provided,
#'   adding additional gene-level information.
#'
#' The function returns a list containing the rank order (`edge`), the enrichment
#' curve, ticks, and other GSEA statistics.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{`edge`}{Data frame with the rank order of genes, enrichment curve, and
#'     tick data.}
#'   \item{`curve`}{Data frame representing the enrichment score curve.}
#'   \item{`ticks`}{Data frame representing the tick marks on the enrichment curve.}
#'   \item{`stats`}{Data frame with additional statistics for the enrichment analysis.}
#'   \item{`posES`, `negES`, `spreadES`, `maxAbsStat`}{Various statistics for the
#'     enrichment analysis (positive ES, negative ES, spread ES, and max absolute
#'     statistic).}
#' }
#'
#' @example
#' # Example usage:
#' rankorder <- get_rankorder(geneset = c("gene1", "gene2"), rankobj = named_vector)
#'
#' @seealso `fgsea::plotEnrichmentData`
#'
get_rankorder <- function(
    geneset_ids,
    rankobj,
    geneset_df = NULL) {
  # geneset_df has addl info
  # if (!class(geneset) == "character") {
  #   stop("geneset must be a character")
  # }
  # if (!class(rankobj) == "numeric") {
  #   stop("rankobj must be a named numeric vector")
  # }

  # if (!"character" %in% class(geneset)){
  #   warning("geneset may be of wrong type")
  #   }
  # Coerce/validate inputs
  if (is.data.frame(geneset_ids)) {
    warning("geneset_ids should be a character vector of IDs; received data.frame")
  }
  if (is.list(geneset_ids)) {
    warning("geneset_ids is a list; using the first element")
    if (length(geneset_ids) >= 1) geneset_ids <- geneset_ids[[1]]
  }
  if (!is.character(geneset_ids)) {
    geneset_ids <- as.character(geneset_ids)
  }

  # type checking for rankobj
  if (is.list(rankobj)) {
    message("rankobj is a list; using the first element")
    if (length(rankobj) >= 1) {
      rankobj <- rankobj[[1]]
    }
  }
  if (!is.numeric(rankobj) || is.null(names(rankobj))) {
    stop("rankobj must be a named numeric vector (names are gene IDs)")
  }

  # browser()
  enplot_data <- plotEnrichmentData(geneset_ids, rankobj)
  rnkorder <- -rankobj %>% rank()
  rankorder_df <- data.frame(
    id = names(rnkorder),
    rank = rnkorder,
    stat = rankobj
  )

  rankorder_edge <- rankorder_df %>% left_join(enplot_data$curve, by = "rank")
  rankorder_edge %<>% left_join(
    dplyr::rename(enplot_data$ticks, stat_tick = stat),
    by = "rank"
  )
  # rankorder_edge %<>% left_join(rename(enplot_data$stats, stat_stat = stat, by = "rank"))
  rankorder_edge$stat == rankorder_edge$stat_tick

  if (!is.null(geneset_df)) {
    rankorder_edge %<>% left_join(
      geneset_df %>%
        mutate(entrez_gene = as.character(entrez_gene)) %>%
        distinct(entrez_gene, .keep_all = TRUE),
      by = c("id" = "entrez_gene")
    )
  }

  return(
    list(
      edge = rankorder_edge,
      curve = enplot_data$curve,
      ticks = enplot_data$ticks,
      stats = enplot_data$stats,
      posES = enplot_data$posES,
      negES = enplot_data$negES,
      spreadES = enplot_data$spreadES,
      maxAbsStat = enplot_data$maxAbsStat
    )
  )
}



get_rankorder_db <- function(
  ranobj_name = NULL,
  pathway_name = NULL,
  db = NULL,
  ...
){

  if (is.null(db)) {
    log_msg(warning="db must be provided")
    return(NULL)
  }

  rankorder <- db$get_plot_enrichmentdata_by_pathway(
    rankobj_name = ranobj_name,
    pathway_name = pathway_name
  )

  return(rankorder)

}

get_rankorder_across <- function(
    df, # long form data frame with rankname col
    ranks_list,
    geneset_lists,
    collection_name = "",
    topn = 25,
    limit = 120,
    title = "",
    pstat_cutoff = 1,
    pstat_usetype = "padj",
    filter_on_mainpathway = TRUE,
    main_pathway_ratio = 0.1,
    db = NULL,
    genesets_additional_info = NULL,
    pathways_of_interest = NULL, # force fetch pathways
    ...) {

  if (!"rankname" %in% colnames(df)) {
    stop("rankname not in fgesa_longdf")
  }

  xtra <- NULL
  if (!is.null(pathways_of_interest)){
    xtra <- df %>% dplyr::filter(pathway %in% pathways_of_interest) %>%
        arrange(pval)
  }

  if (isTRUE(filter_on_mainpathway)){
    # Avoid name collision with logical parameter by fetching the function from the parent env
    filter_on_mainpathway_fn <- get("filter_on_mainpathway", envir = parent.env(environment()), mode = "function")
    df <- filter_on_mainpathway_fn(df, main_pathway_ratio = main_pathway_ratio)
  }


  df <- df %>%
    filter(!!as.symbol(pstat_usetype) < pstat_cutoff) %>%
    arrange(pval) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    head(n = limit)

  if (!is.null(xtra)){
    df <- bind_rows(df, distinct(xtra, .keep_all = TRUE))
  }


  pathways_to_plot <- df$pathway
  rank_ids <- names(ranks_list)

  # if (!is.null(db)) {
  #   con <-db_tools$get_con(db)
  #   on.exit(db_tools$close_con(con), add = TRUE)
  #   rankorders_cached <- pathways_to_plot %>% purrr::map(~{
  #     pathway_name <- .x
  #     rank_ids %>% purrr::map(~{
  #       db_tools$insert_rankorder(
  #         con = con,
  #         rankobj_name = .x,
  #         pathway_name = pathway_name,
  #         rankorder = get_rankorder_db(
  #           ranobj_name = .x,
  #           pathway_name = pathway_name,
  #           con = con
  #         )
  #       )
  #     })
  #   })
  #    #%>%
  #     #set_names(pathways_to_plot)
  #   #rankorders
  #   missing <- Filter(Negate(is.null), rankorders_cached)
  # } # check which are missing for below

  # print(paste0("missing: ", missing))
  # print(paste0("pathways to plot: ", pathways_to_plot))
  # print(paste0("rank ids ", rank_ids))



  # rankorders <- pathways_to_plot %>%
  #   purrr::map(~ {
  #     if (!.x %in% names(geneset_lists)) {
  #       cat(paste0("does not have access to geneset : ", .x))
  #       return()
  #     }
  #     geneset <- geneset_lists[[.x]]
  #     rank_ids %>%
  #       purrr::map(~ 
  #         tryCatch({
  #           rank_order <- fgsea_tools$get_rankorder(
  #             geneset,
  #             ranks_list[[.x]],
  #             geneset_df = genesets_additional_info
  #           ) #set_names(rank_ids)}, 
  #           # add check the dims are the same 
  #          names(rank_order) <- rank_ids 
  #          return(rank_order)
  #         }
  #      error = function(e) {print(e)})

  #   )}
  #     ) %>%
  #   set_names(pathways_to_plot)

  # this is all broken
  # it was an attem pted refactor that has gone wrong
  # rankorders <- pathways_to_plot %>%
  #   purrr::map(~ {
  #     pathway <- .x
  #     if (!pathway %in% names(geneset_lists)) {
  #       cat(paste0("❗️ Missing geneset: ", pathway, "\n"))
  #       return(NULL)
  #     }
  #     geneset <- geneset_lists[[pathway]]
  
    #   ranks_for_pathway <- rank_ids %>%
    #     purrr::map(~ tryCatch({
    #       rank_order <- get_rankorder(
    #         geneset,
    #         ranks_list[[.x]],
    #         geneset_df = genesets_additional_info
    #       )
    #       #names(rank_order) <- rank_ids
    #       if (length(rank_order) == length(rank_ids)) {
    #         names(rank_order) <- rank_ids
    #       } else {
    #         message("❗️ Length mismatch in rank_order and rank_ids for pathway: ", pathway)
    #         message("length(rank_order) = ", length(rank_order), ", length(rank_ids) = ", length(rank_ids))
    #         # browser()
    #         names(rank_order) <- NULL  # optional: clear names to avoid misleading. this may have unexpected downstream consequences
    #       }  # a better fallback - and figure why this can happen - would help
    #
    #
    #       return(rank_order)
    #     },
    #     error = function(e) {
    #       message("\n⚠️ Error in get_rankorder for pathway: ", pathway)
    #       message(conditionMessage(e))
    #       NULL
    #     }))
    #
    #   ranks_for_pathway
    # }) %>%
    # set_names(pathways_to_plot)
  # this is all broken


  rankorders <- pathways_to_plot %>%
    purrr::map(~ {
      if (!.x %in% names(geneset_lists)) {
        cat(paste0("does not have access to geneset : ", .x))
        return()
      }
      geneset <- geneset_lists[[.x]]
      rank_ids %>%
        purrr::map(~ {
          fgsea_tools$get_rankorder(
            geneset,
            ranks_list[[.x]],
            geneset_df = genesets_additional_info
          )
        }) %>%
        set_names(rank_ids)
    }) %>%
    set_names(pathways_to_plot)




  return(rankorders)
}


# combine_rankorders_on_sample <- function(
#     rankorders,
#     metadata = NULL,
#     ...) {
#   if (!is.null(metadata)) {
#     if (!"rankname" %in% colnames(metadata)) {
#       warn("metadata must have a 'rankname' column")
#       metadata <- NULL
#     }
#     if (!"facet" %in% colnames(metadata)) {
#       metadata$facet <- metadata$rankname
#     }
#   }
#
#   res_list <- rankorders %>%
#     purrr::imap(~ {
#       pw_name <- .y
#       list_of_rankorders <- .x # this is a list of the "rankorder" info
#       # names incldue:
#       # edge (df), curve (df), ticks, stats, posES, negES, spreadES, maxAbsStat
#
#       curves <- list_of_rankorders %>%
#         purrr::imap(~ {
#           .x$curve %>% mutate(pathway = pw_name, rankname = .y)
#         }) %>%
#         bind_rows()
#       if (!is.null(metadata)) curves %<>% left_join(metadata, by = "rankname") # by = "rankname"
#
#       edges <- list_of_rankorders %>%
#         purrr::imap(~ {
#           .x$edge %>% mutate(pathway = pw_name, rankname = .y)
#         }) %>%
#         bind_rows()
#       if (!is.null(metadata)) edges %<>% left_join(metadata, by = "rankname")
#
#       ticks <- list_of_rankorders %>%
#         purrr::imap(~ {
#           .x$ticks %>% mutate(pathway = pw_name, rankname = .y)
#         }) %>%
#         bind_rows()
#       if (!is.null(metadata)) ticks %<>% left_join(metadata, by = "rankname")
#
#       stats <- list_of_rankorders %>%
#         purrr::imap(~ {
#           .x$stats %>% mutate(pathway = pw_name, rankname = .y)
#         }) %>%
#         bind_rows()
#       if (!is.null(metadata)) stats %<>% left_join(metadata, by = "rankname")
#
#       res <- list(
#         curve = curves,
#         edge = edges,
#         ticks = ticks,
#         stats = stats
#       )
#
#       return(res)
#     })
#   return(res_list)
# }


concat_results_one_collection <- function(list_of_dfs, main_pathway_ratio=0.1) {
  if (length(list_of_dfs) == 0) {
    return(data.frame())
  }
  res <- list_of_dfs %>%
    purrr::imap(~ .x %>% dplyr::mutate(rankname = .y)) %>%
    dplyr::bind_rows()
  res <- filter_on_mainpathway(res, main_pathway_ratio = main_pathway_ratio)
  return(res)
}

concat_results_all_collections <- function(list_of_lists, main_pathway_ratio = .1, ...) {
  .dotargs <- list(...) ## this is not used nor passed to inner func

  res <- list_of_lists %>%
    purrr::map(
      ~ {
        concat_results_one_collection(.x)
      }
    )
  return(res)
}



# new test better refactor

process_rankorder_list <- function(rankorders, pw_name, metadata = NULL, key) {
  rankorders %>%
    purrr::imap(~ {
      .x[[key]] %>% dplyr::mutate(pathway = pw_name, rankname = .y)
    }) %>%
    bind_rows() %>%
    { if (!is.null(metadata)) left_join(., metadata, by = "rankname") else . }
}


combine_rankorders_on_sample <- function(rankorders, metadata = NULL, ...) {
  if (!is.null(metadata)) {
    if (!"rankname" %in% colnames(metadata)) {
      warn("metadata must have a 'rankname' column")
      metadata <- NULL
    }
    if (!"facet" %in% colnames(metadata)) {
      metadata$facet <- metadata$rankname
    }
  }

  res_list <- rankorders %>%
    purrr::imap(~ {
      pw_name <- .y
      list_of_rankorders <- .x

      curves <- process_rankorder_list(list_of_rankorders, pw_name, metadata, "curve")
      edges <- process_rankorder_list(list_of_rankorders, pw_name, metadata, "edge")
      ticks <- process_rankorder_list(list_of_rankorders, pw_name, metadata, "ticks")
      stats <- process_rankorder_list(list_of_rankorders, pw_name, metadata, "stats")

      list(
        curve = curves,
        edge = edges,
        ticks = ticks,
        stats = stats
      )
    })

  return(res_list)
}
