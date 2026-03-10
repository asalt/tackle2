suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))
map_tools <- new.env()
source(file.path(src_dir, "./map.R"), local = map_tools)

test_that("infer_gene_id_keytype detects ENSEMBL and ENTREZID", {
  expect_equal(map_tools$infer_gene_id_keytype(c("1", "2", "42")), "ENTREZID")
  expect_equal(map_tools$infer_gene_id_keytype(c("ENSG00000141510", "ENSG00000171862")), "ENSEMBL")
  expect_true(is.na(map_tools$infer_gene_id_keytype(c("TP53", "PTEN"))))
})

test_that("add_leadingedges_to_results_list maps ENSEMBL ids to symbols", {
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE))

  df <- data.frame(
    pathway = "demo",
    pval = 0.1,
    padj = 0.2,
    log2err = 0,
    ES = 1,
    NES = 1,
    size = 2,
    rankname = "rank1",
    leadingEdge = I(list(c("ENSG00000141510", "ENSG00000171862")))
  )

  res <- map_tools$add_leadingedges_to_results_list(list(df), species = "Homo sapiens")
  expect_true("leadingEdge_genesymbol" %in% colnames(res[[1]]))
  expect_equal(res[[1]]$leadingEdge_genesymbol[[1]], "TP53/PTEN")
})

test_that("map_entrez_to_symbol collapses multi-mappings with |", {
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  skip_if_not(requireNamespace("AnnotationDbi", quietly = TRUE))

  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))

  # Find a key with multiple SYMBOL mappings (if any exist for this DB version).
  all_keys <- head(keys(org.Hs.eg.db, keytype = "ENSEMBL"), 2000)
  sel <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_keys,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  ))
  counts <- table(sel$ENSEMBL)
  multi <- names(counts[counts > 1])
  if (length(multi) == 0) {
    skip("No multi-mapping ENSEMBL keys found in sampled set")
  }

  key <- multi[[1]]
  mapping <- map_tools$map_entrez_to_symbol(key, species = "Homo sapiens", keytype = "ENSEMBL")
  expect_true(grepl("\\|", unname(mapping[[1]])))
})
