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
    purrr::set_names(nm = ~ basename(.) %>%
        sub("\\.rnk$", "", .)) %>%  # fs::path_ext_remove()) %>% # set names first
    purrr::map(~ {
      .table <- read_tsv(.x, show_col_types = F)
      if (value_col %in% colnames(.table)) {
        .table <- .table %>% rename(value = !!value_col)
      }
      if (tolower(value_col) %in% tolower(colnames(.table))) {
        .match <- colnames(.table)[ stringr::str_detect(tolower(colnames(.table)), tolower(value_col)) ]
        .table <- dplyr::rename(.table, value = !!rlang::sym(.match))
      }
      if (id_col %in% colnames(.table)) {
        .table <- .table %>% rename(id = !!id_col)
      }
      return(.table)
    })


  log_msg(msg = "trying to shorten names")

  shorternames <- names(lst) %>%
    stringr::str_extract(., pattern = "(?<=group_)([^.*]*)$")
  log_msg(msg = paste0("shorter names are ", paste0(shorternames, "\n")))
  if (!all(is.na(shorternames)) && length(unique(shorternames)) == length(names(lst))) {
    names(lst) <- shorternames
  } else {
    log_msg(msg = "nas in shorter names, not reassigning")
  }
  lst
}


write_rnkfiles <- function(
    lst,
    dir = "rnkfiles") {
  if (is.null(dir)) {
    dir <- "rnkfiles"
  }
  if (!fs::dir_exists(dir)) {
    log_msg(msg = paste0("creating ", dir))
    fs::dir_create(dir)
  }
  lst %>% purrr::iwalk( # .x is the value, .y is the name
    ~ {
      .outname <- fs::path_join(
        c(dir, paste0(.y, ".rnk"))
      )
      if (!fs::file_exists(.outname)) {
        if (("GeneID" %in% colnames(.x)) && (!"id" %in% colnames(.x))) .x %<>% dplyr::rename(id = GeneID)
        .x %>%
          dplyr::select(id, value) %>%
          write_tsv(.outname, col_names = FALSE)
        print(paste0("Wrote ", .outname))
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

  ranks_list <- rnkdfs %>% purrr::map(
    ~ with(.x, setNames(value, id))
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


write_results <- function(result, outf, replace = FALSE) {

  if (is.null(replace)) replace <- FALSE
  if (!is.data.frame(result)) {
    log_msg(msg = "Invalid result, cannot write to file.")
    log_msg(msg = as.character(result))
    stop("Invalid result, cannot write to file.")
  }

  if (!"leadingEdge" %in% colnames(result)) {
    log_msg(msg = "leadingEdge column not found in the input data")
    stop("leadingEdge column not found in the input data")
  }

  if (fs::file_exists(outf) && !replace) {
    log_msg(msg = paste0("File ", outf, " already exists, skipping"))
    return()
  }

  result %>%
    mutate(leadingEdge = purrr::map_chr(leadingEdge, paste, collapse = "/")) %>%
    write_tsv(outf)
  log_msg(msg = paste0("Successfully written to ", outf))
}

# Main function to save GSEA results
save_individual_gsea_results <- function(
  results_list,
  savedir = "gsea_tables",
  replace = FALSE,
  species = "Homo sapiens") {

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
    name_map <- util_tools$make_name_map(comparison_names)
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
        result %>% write_tsv(outf)
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

save_pivoted_gsea_results <- function(results_list, savedir = "gsea_tables", replace = FALSE, species = species) {

  if (is.null(replace)) replace <- FALSE
  # here results list is concatenated list. one level per collection
  # names are the collection names
  # values are the fgsea concatenated tables/comparisons for given collection

  results_list %<>% map_tools$add_leadingedges_to_results_list(species = species)

  results_list_pivoted <- results_list %>% purrr::imap(
    ~{
      result_collection  <- .x
      collection_name <- .y

      res_pivoted <- result_collection %>%
        pivot_wider(
          id_cols = c("pathway"),
          names_from = rankname,
          values_from = c(pval, padj, log2err, ES, NES, n_main, size, leadingEdge_genesymbol, leadingEdge_entrezid, mainpathway),
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
      res_pivoted %>% write_tsv(outf)
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


load_and_process_ranks <- function(params) {
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


  log_msg(msg = paste0("ranks from : ", ranks_from))
  log_msg(msg = paste0("rankfiledir : ", rankfiledir))

  if (!is.null(rankfiledir) && file.exists(rankfiledir)) { #
    rnkfiles <- dir_ls(path = rankfiledir, regexp = ".*\\.rnk$", fail = FALSE)
    log_msg(msg = paste0("looking for rank files in ", rankfiledir))


    if (length(rnkfiles) > 0) {
      log_msg(msg = paste0("found ", length(rnkfiles), " rankfiles"))
      rnkdfs <- rnkfiles %>% load_rnkfiles()
      names(rnkdfs) <- names(rnkdfs) %>%
        fs::path_file() %>%
        sub("\\.rnk$", "", .)

        #fs::path_ext_remove() # this is no good


      name_mapping_file <- file.path(rankfiledir, 'names.txt')
      if (fs::file_exists(name_mapping_file)) {
        log_msg(msg = paste0("found name mapping file: ", name_mapping_file))
        name_mapping <- read_delim(name_mapping_file,
          col_names = c("new", "old"),
          delim = '=',
           comment = '#',
           show_col_types = F
        ) %>% mutate(old = sub("\\.rnk$", "", old))

        # fs::path_ext_remove(old))

        missing_rankfiles <- setdiff(name_mapping$old, names(rnkdfs))
        if (length(missing_rankfiles) > 0) {
          log_msg(warning = paste0(
            "names.txt references rank files that were not loaded: ",
            paste(missing_rankfiles, collapse = ", "),
            ". These entries will be skipped."
          ))
        }

        duplicated_new_labels <- name_mapping$new[duplicated(name_mapping$new)]
        if (length(duplicated_new_labels) > 0) {
          log_msg(warning = paste0(
            "names.txt contains duplicate target labels: ",
            paste(unique(duplicated_new_labels), collapse = ", "),
            ". Later mappings will overwrite earlier ones."
          ))
        }

        for (ix in seq_len(length(rnkdfs))) {
          .new <- name_mapping[ix, ]$new
          .old <- name_mapping[ix, ]$old
          if (.old %in% names(rnkdfs)) {
            rnkdfs[.new] <- rnkdfs[.old]
            rnkdfs[.old] <- NULL
          } else {
            log_msg(debug = paste0("Skipping mapping for missing rank file: ", .old))
          }
        }

      }

      ranks_list <- rnkdfs %>% ranks_dfs_to_lists()

      return(ranks_list)
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
    ranks_list <- rnkdfs %>% ranks_dfs_to_lists()
    return(ranks_list)
  }
  if (ranks_from == "volcano") {
    if (is.null(volcanodir) || !file.exists(volcanodir)) {
      stop(paste0("improper volcanodir specification: ", volcanodir))
    }
    log_msg(msg = "saving rankfiles from volcano output. using signedlogp as value")
    rnkdfs <- create_rnkfiles_from_volcano(volcanodir, value_col = "signedlogP")
    rnkdfs %>% write_rnkfiles(dir = rankfiledir) # and save
    names(rnkdfs) <- names(rnkdfs) %>%
      fs::path_file() %>%
      fs::path_ext_remove()
    log_msg(paste0("length of retrieved rankfiles: ", length(rnkdfs)))
    ranks_list <- rnkdfs %>% ranks_dfs_to_lists()
    return(ranks_list)
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
    ranks_list <- rnkdfs %>% ranks_dfs_to_lists()
  }
  # not sure if this level of flow is relevant, refactor later
  if (!exists("ranks_list")) {
    stop("No rankfiles found, problem loading")
  }
  return(ranks_list)
}
