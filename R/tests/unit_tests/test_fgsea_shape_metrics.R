suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)


test_that("add_pathway_shape_metrics computes pathway-level shape metrics", {
  fgsea_res <- data.frame(
    pathway = c("positive_pathway", "negative_pathway"),
    ES = c(2.0, -1.5),
    size = c(6, 5),
    stringsAsFactors = FALSE
  )
  fgsea_res$leadingEdge <- list(c(2L, 5L, 6L), c(18L, 15L, 12L))

  rankobj <- setNames(seq(20, 1), paste0("gene", seq_len(20)))
  out <- fgsea_tools$add_pathway_shape_metrics(fgsea_res, rankobj)

  expect_equal(out$peak_rank_pct[[1]], 6 / 20, tolerance = 1e-9)
  expect_equal(out$leading_edge_fraction[[1]], 3 / 6, tolerance = 1e-9)
  expect_equal(out$leading_edge_span_pct[[1]], (6 - 2) / 20, tolerance = 1e-9)
  expect_equal(out$front_loaded_score[[1]], 2 / (6 / 20), tolerance = 1e-9)

  expect_equal(out$peak_rank_pct[[2]], 12 / 20, tolerance = 1e-9)
  expect_equal(out$leading_edge_fraction[[2]], 3 / 5, tolerance = 1e-9)
  expect_equal(out$leading_edge_span_pct[[2]], (18 - 12) / 20, tolerance = 1e-9)
  expect_equal(out$front_loaded_score[[2]], -1.5 / (12 / 20), tolerance = 1e-9)
})


test_that("prepare_results_for_export preserves numeric precision", {
  dataframe <- data.frame(
    NES = c(2.0186087352224678, -1.2345678901),
    padj = c(3.312002011811893e-7, 0.004856511886780497),
    stringsAsFactors = FALSE
  )
  dataframe$leadingEdge <- list(c("TP53", "BAX"), c("SMAD2"))

  out <- io_tools$prepare_results_for_export(dataframe)

  expect_equal(out$NES[[1]], 2.0186087352224678)
  expect_equal(out$NES[[2]], -1.2345678901)
  expect_equal(out$padj[[1]], 3.312002011811893e-7)
  expect_equal(out$padj[[2]], 0.004856511886780497)
  expect_equal(out$leadingEdge[[1]], "TP53/BAX")
  expect_equal(out$leadingEdge[[2]], "SMAD2")
})
