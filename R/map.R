# map.R
# Load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
# library(writexl)
# library(org.Hs.eg.db) # For human gene annotations

# required_bioc <- c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db") # Add more as needed
# installed_bioc <- rownames(installed.packages())
# for (pkg in required_bioc) {
#   if (!pkg %in% installed_bioc) {
#     BiocManager::install(pkg)
#   }
# }


# in the future we will need to support others as well
# Define the gene mapping function
normalize_ensembl_ids <- function(ids) {
  ids <- as.character(ids)
  sub("\\.[0-9]+$", "", ids)
}

infer_gene_id_keytype <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids)]
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0) return(NA_character_)

  ids_norm <- normalize_ensembl_ids(ids)
  is_entrez <- stringr::str_detect(ids_norm, "^[0-9]+$")
  # Conservative "ENSEMBL gene id" pattern, with fallback to any ENS* id.
  is_ensembl <- stringr::str_detect(ids_norm, "^ENS[A-Z]*G[0-9]+$") |
    stringr::str_detect(ids_norm, "^ENS[0-9A-Z]+$")

  if (all(is_entrez)) return("ENTREZID")
  if (all(is_ensembl)) return("ENSEMBL")

  if (mean(is_entrez) >= 0.8) return("ENTREZID")
  if (mean(is_ensembl) >= 0.8) return("ENSEMBL")

  NA_character_
}

map_entrez_to_symbol <- function(entrez_list, species = "Homo sapiens", keytype = "ENTREZID") {

  # Define species-specific annotation packages
  species_db <- list(
    "Homo sapiens" = "org.Hs.eg.db",
    "Mus musculus" = "org.Mm.eg.db",
    "Rattus norvegicus" = "org.Rn.eg.db"
    # Add more species and their corresponding annotation packages here
  )

  # Validate species input
  if (!(species %in% names(species_db))) {
    stop("Unsupported species. Please choose from: ", paste(names(species_db), collapse = ", "))
  }

  # Load the appropriate annotation package
  db_package <- species_db[[species]]
  if (!requireNamespace(db_package, quietly = TRUE)) {
    stop("Annotation package ", db_package, " not found. Please install it.")
  }
  suppressPackageStartupMessages(library(db_package, character.only = TRUE))

  entrez_list <- as.character(entrez_list)
  entrez_list <- entrez_list[!is.na(entrez_list)]
  entrez_list <- entrez_list[nzchar(entrez_list)]
  if (length(entrez_list) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (identical(keytype, "ENSEMBL")) {
    entrez_list <- normalize_ensembl_ids(entrez_list)
  }
  entrez_list <- unique(entrez_list)

  collapse_symbols <- function(values) {
    values <- as.character(values)
    values <- values[!is.na(values)]
    values <- values[nzchar(values)]
    if (length(values) == 0) return(NA_character_)
    paste(unique(values), collapse = "|")
  }

  # Map gene IDs to Gene Symbols
  quiet <- isTRUE(getOption("tackle2_quiet", FALSE))
  map_call <- function() {
    mapIds(
      get(db_package),
      keys = entrez_list,
      column = "SYMBOL",
      keytype = keytype,
      multiVals = collapse_symbols
    )
  }
  gene_symbols <- if (quiet) suppressMessages(map_call()) else map_call()

  return(gene_symbols)
}

extract_entrezids <- function(fgsea_res){

  # print(head(fgsea_res))
  # print(head(fgsea_res$leadingEdge))
  # print(class(fgsea_res$leadingEdge))

  if (!is.data.frame(fgsea_res)){
    stop("fgsea_res must be a data.frame")
  }

  if (!"leadingEdge" %in% colnames(fgsea_res)){
    stop("leadingEdge not in fgsea results")
  }

  le <- fgsea_res$leadingEdge
  ids <- NULL

  if (is.list(le)) {
    ids <- unlist(le, use.names = FALSE)
  } else if (is.character(le)) {
    ids <- stringr::str_split(le, pattern = "/") %>% unlist(use.names = FALSE)
  } else {
    ids <- le
  }

  ids <- as.character(ids)
  ids <- ids[!is.na(ids)]
  ids <- stringr::str_trim(ids)
  ids <- ids[nzchar(ids)]
  unique(ids)
}

format_entrezids <- function(fgsea_res, mapping, keytype = "ENTREZID"){

  if (is.null(mapping)) stop("must provide mapping")
  if (class(fgsea_res$leadingEdge) == "character")  {
    fgsea_res <- fgsea_res %>% mutate( leadingEdge = str_split(leadingEdge, '/'))
  }

  normalize_id <- if (identical(keytype, "ENSEMBL")) normalize_ensembl_ids else identity

  # Function to replace Entrez IDs with symbols or keep original if no mapping
  replace_ids <- function(ids) {
    ids <- Filter(Negate(is.na), ids) # Remove NAs
    ids_norm <- normalize_id(ids)
    symbols <- sapply(seq_along(ids_norm), function(ix) {
      id_norm <- ids_norm[[ix]]
      id_raw <- ids[[ix]]
      if (id_norm %in% names(mapping)) {
        symbol <- mapping[[id_norm]]
        if (!is.na(symbol)) {
          return(symbol)
        } else {
          return(id_raw) # Keep original ID if no symbol
        }
      } else {
        return(id_raw) # Keep original ID if not found
      }
    })
    paste(symbols, collapse = "/")
  }

  fgsea_res <- fgsea_res %>% mutate(
    leadingEdge_entrezid = map_chr(leadingEdge, ~ paste(.x, collapse = "/")),
    leadingEdge_genesymbol = map_chr(leadingEdge, ~replace_ids(.x))
  )

  return(fgsea_res)
}

add_leadingedges_to_results_list <- function(fgsea_res_list, species = "Homo sapiens"){
  if (!is.list(fgsea_res_list)){
    stop("fgsea_res_list should be a list")
  }

  gene_ids <- fgsea_res_list %>%
    purrr::map(~extract_entrezids(.x)) %>%
    purrr::flatten_chr() %>%
    unique()

  keytype <- infer_gene_id_keytype(gene_ids)
  if (is.na(keytype)) {
    return(fgsea_res_list)
  }

  if (!isTRUE(getOption("tackle2_quiet", FALSE))) {
    message(paste0("mapping gene ids (", keytype, ") to symbols for ", species))
  }

  mapping <- tryCatch(
    map_entrez_to_symbol(gene_ids, species = species, keytype = keytype),
    error = function(e) {
      message("[map] failed to map gene ids to symbols: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(mapping)) {
    return(fgsea_res_list)
  }

  output <- fgsea_res_list %>% purrr::map(~{
    format_entrezids(.x, mapping = mapping, keytype = keytype)
  })

  return(output)
}
