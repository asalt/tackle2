# geneset_utils.R
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(memoise))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tibble))

# util_tools <- new.env()
# source(file.path(here("R"), "./utils.R"), local = util_tools)

source(file.path(here("R", "lazyloader.R")))
util_tools <- get_tool_env("utils")
log_msg <- util_tools$make_partial(util_tools$log_msg)


c1_fix <- function(df){
    df %>% mutate( db_ncbi_gene = ensembl_gene ) %>% mutate( ncbi_gene = ensembl_gene)
}

get_collection_raw <- function(
    category,
    subcategory,
    species = "Homo sapiens",
    cache = TRUE,
    logger = NULL) {
  #
  if (is.null(logger)) logger <- log_msg
  collection_id <- paste(make.names(category), make.names(subcategory), make.names(species), sep = "_")

  post_func <- if (category == "C1") c1_fix else identity

  cache_dir <- here("cache")
  collection_id_path <- file.path(cache_dir, collection_id)

  if (cache && fs::file_exists(collection_id_path)) {
    cat(paste0("reading ", collection_id, " from ", cache_dir, "\n"))
    logger(msg = paste0("reading ", collection_id, " from ", cache_dir, "\n"))
    log_msg(info = paste0("reading ", collection_id, " from ", cache_dir, "\n"))
    df <- readr::read_tsv(collection_id_path, show_col_types = FALSE)
    return(df %>% post_func() )
  }
  cat(paste0(collection_id, " not found in ", cache_dir, "\n"))


  msig <- function(dbsp) {
  msigdbr::msigdbr(
    db_species   = dbsp,
    species      = species,
    collection   = category,
    subcollection = if (nzchar(subcategory)) subcategory else NULL
    )
  }
  db_species_primary <- if (species != "Mus musculus") "HS" else "MM"

  # Data fetching

  df <- tryCatch(
  msig(db_species_primary),
  error = function(e_primary) {
    # only fall back if the primary was "MM" (your stated requirement)
    if (db_species_primary == "MM") {
      message(sprintf(
        "msigdbr failed with db_species='%s': %s\nRetrying with db_species='HS'...",
        db_species_primary, conditionMessage(e_primary)
      ))
      tryCatch(
        msig("HS"),
        error = function(e_alt) {
          stop(sprintf(
            "msigdbr failed with db_species='MM' and fallback 'HS'.\nPrimary error: %s\nFallback error: %s",
            conditionMessage(e_primary), conditionMessage(e_alt)
          ))
        }
      )} else {
      # primary wasn't MM → per spec, do not flip; propagate the original error
      stop(e_primary)
      }
    }
  )

  # df <- msigdbr::msigdbr(
  #   db_species = if (species != "Mus musculus") "HS" else "MM",
  #   species = species,
  #   collection = category,
  #   subcollection = if (nchar(subcategory)>0) subcategory else NULL
  #   # category = category,
  #   # subcategory = subcategory
  # )

  if (cache == TRUE) {
    if (!fs::dir_exists(cache_dir)) fs::dir_create(cache_dir)
    if (!fs::file_exists(collection_id_path)) {
      readr::write_tsv(df, collection_id_path)
      cat(paste0("writing ", collection_id, " to ", cache_dir, "\n"))
    }
  }

  # special case if C1 must be indexed on ensembl

  return(df %>% post_func())
}

get_collection <- memoise::memoise(get_collection_raw) # this is convienent for interactive sessions

# get_collections <- function(list_of_collections, species = "Homo sapiens") {
#   res <- list_of_collections %>% purrr::map(
#     ~ {
#       list_name <- paste(.x$category, .x$subcategory, sep = "_")
#       collection <- get_collection(
#         category = .x$category,
#         subcategory = .x$subcategory,
#         species = species,
#       )
#       setNames(list(collection), list_name)
#     }
#   )
#   res_reduced <- purrr::reduce(res, c) # 1 level list names set appropriately
#   return(res_reduced)
# }

get_collections <- function(dataframe_obj, species = "Homo sapiens") {
  if (!"collection_name" %in% colnames(dataframe_obj)) {
    dataframe_obj <- dplyr::mutate(dataframe_obj,
      collection_name = stringr::str_c(category, subcategory, sep = "_")
    )
  }
  res <- dataframe_obj %>% purrr::pmap(
    function(category, subcategory, ...) {
      list_name <- paste(category, subcategory, sep = "_")
      if (category == "C5" & subcategory == "All") {
        .collections <- c("GO:BP", "GO:CC", "GO:MF") %>%
          purrr::map(~ {
            get_collection(
              category = category,
              subcategory = .x,
              species = species
            )
          })
        collection <- dplyr::bind_rows(.collections)
      } else {
        collection <- get_collection(
          category = category,
          subcategory = subcategory,
          species = species,
        )
      }
      setNames(list(collection), list_name)
    }
  )
  res_reduced <- purrr::reduce(res, c) # 1 level list names set appropriately
  return(res_reduced)
}



# list_of_geneset_dfs <- pathways_of_interest %>%
#   purrr::map(  # we use pmap for rowwise access to multiple values (?)
#     ~ {
#     # Here you can customize how you want to name each list element based on the inputs
#     list_name <- paste(.x$category, .x$subcategory, sep = "_")
#     list_data <- get_genesets(.x$category, .x$subcategory)
#     # Return a named list for each row
#     setNames(list(list_data), list_name)
#   }) %>%
#   purrr::reduce(c) # reduce a list to a single value by applying a binary function, in this case c.

get_pathway_info <- function(gsname) {
  pathways_dfs[[gsname]]
}

# now turn each pathway dataframe into a named list
# Transform each DataFrame and convert it into a named list of gene ids
genesets_df_to_list <- function(list_of_geneset_dfs) { # the input here is a dataframe, consider renaming variable
  genesets_list <- list_of_geneset_dfs %>%
    # mutate(gs_fullname = str_c(gs_exact_source, gs_name, sep=" ")) %>%
    # group_by(gs_fullname) %>%
    group_by(gs_name) %>%
    #summarise(entrez_gene_ids = list(as.character(db_ncbi_gene)), .groups = "drop") %>% # not db_ncbi_gene, but ncbi_gene
    summarise(entrez_gene_ids = list(as.character(ncbi_gene)), .groups = "drop") %>%
    tibble::deframe()
    # summarise(entrez_gene_ids = list(as.character(entrez_gene)), .groups = "drop") %>%
  genesets_list
}


geneset_array_to_df <- function(gs) {
  # normalise a mixed list of user-supplied entries into a rectangular data frame
  # (one row per geneset declaration). Each field is collapsed to a single scalar
  # so that `bind_rows()` never encounters vectors of differing lengths.
  first_or <- function(value, default) {
    if (is.null(value) || length(value) == 0) default else value[[1]]
  }

  entries <- gs %||% list()
  if (length(entries) == 0) {
    return(tibble::tibble(
      category = character(0),
      subcategory = character(0),
      collapse = logical(0)
    ))
  }

  rows <- purrr::imap(entries, function(entry, idx) {
    entry <- entry %||% list()

    category_raw <- first_or(entry$category, NA_character_)
    category <- as.character(category_raw)[[1]]
    if (is.na(category)) category <- NA_character_
    category <- trimws(category)
    if (is.na(category) || !nzchar(category)) {
      rlang::warn(sprintf("Dropping geneset entry %s with missing category", idx))
      return(tibble::tibble())
    }

    subcategory_raw <- first_or(entry$subcategory, "")
    subcategory <- as.character(subcategory_raw)[[1]]
    if (is.na(subcategory)) subcategory <- ""
    subcategory <- trimws(subcategory)

    collapse_raw <- first_or(entry$collapse, FALSE)
    collapse <- if (is.logical(collapse_raw)) {
      isTRUE(collapse_raw)
    } else if (is.numeric(collapse_raw)) {
      collapse_raw != 0
    } else if (is.character(collapse_raw)) {
      tolower(collapse_raw) %in% c("true", "t", "1", "yes")
    } else {
      FALSE
    }

    tibble::tibble(
      category = category,
      subcategory = subcategory,
      collapse = collapse
    )
  })

  rows <- purrr::discard(rows, ~ nrow(.x) == 0)
  if (length(rows) == 0) {
    return(tibble::tibble(
      category = character(0),
      subcategory = character(0),
      collapse = logical(0)
    ))
  }

  dplyr::bind_rows(rows)
}
