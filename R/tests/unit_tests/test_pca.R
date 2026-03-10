suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(testthat))


src_dir <- file.path(here("R"))
pca_tools <- new.env()
source(file.path(src_dir, "pca.R"), local = pca_tools)


make_fake_gsea <- function(n_pathways = 5, n_ranks = 3) {
  grid <- expand.grid(
    pathway = paste0("path", seq_len(n_pathways)),
    rankname = paste0("rank", seq_len(n_ranks)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$NES <- rnorm(nrow(grid))
  grid$var <- paste0("var", seq_len(nrow(grid)))
  grid$mainpathway <- sample(c(TRUE, FALSE), nrow(grid), replace = TRUE)
  grid
}



context("Testing do_one function")

test_that("do_one stops with non-existent required columns", {
  fake_data <- data.frame(NES = rnorm(10), var = letters[1:10])
  expect_error(pca_tools$do_one(fake_data), "pathway column not found in the input data")
})

test_that("do_one computes PCA correctly with minimal inputs", {
  # Assuming a minimal correct dataset
  fake_data <- make_fake_gsea(n_pathways = 5, n_ranks = 2)
  result <- pca_tools$do_one(fake_data)
  expect_true("pca" %in% class(result)) # Check if PCA result is returned, it should be pca obj from PCAtools::pca
})

test_that("do_one handles different main_pathway_ratios correctly", {
  # Make ratios deterministic:
  # - rank1, rank2: all main pathways (ratio 1.0)
  # - rank3, rank4: partial main pathways (ratio < 0.95)
  fake_data <- make_fake_gsea(n_pathways = 5, n_ranks = 4)
  fake_data$mainpathway <- FALSE
  fake_data$mainpathway[fake_data$rankname %in% c("rank1", "rank2")] <- TRUE
  fake_data$mainpathway[fake_data$rankname == "rank3" & fake_data$pathway == "path1"] <- TRUE
  fake_data$mainpathway[fake_data$rankname == "rank4" & fake_data$pathway %in% c("path1", "path2")] <- TRUE

  low_ratio_result <- pca_tools$do_one(fake_data, main_pathway_ratio = 0.05)
  high_ratio_result <- pca_tools$do_one(fake_data, main_pathway_ratio = 0.95)
  expect_true(nrow(high_ratio_result$rotated) <= nrow(low_ratio_result$rotated))
})


context("Testing do_all function")

test_that("do_all applies do_one to multiple GSEA objects correctly", {
  fake_list <- list(
    data1 = make_fake_gsea(n_pathways = 4, n_ranks = 3),
    data2 = make_fake_gsea(n_pathways = 3, n_ranks = 2)
  )
  results <- pca_tools$do_all(fake_list)
  expect_equal(length(results), 2)
  expect_true(all(sapply(results, function(x) "pca" %in% class(x))))
})

# test_that("do_all handles empty lists gracefully", {
#   expect_error(pca_tools$do_all(list()), "Error in do_one")  # Depending on how do_one handles empty data
# })



context("Testing plot_biplot function")

test_that("plot_biplot doesnot crash if not enough pcs", {
  fake_data <- make_fake_gsea(n_pathways = 5, n_ranks = 2)
  pca_object <- pca_tools$do_one(fake_data)
  expect_warning(pca_tools$plot_biplot(pca_object, colby = "nonexistent"))
})

test_that("plot_biplot works", {
  fake_data <- make_fake_gsea(n_pathways = 6, n_ranks = 5)
  pca_object <- pca_tools$do_one(fake_data)
  plts <- pca_tools$plot_biplot(pca_object)
  for (plt in plts) {
    testthat::expect_true(
      "ggplot" %in% class(plt)
      # all(
      #   "gg" %in% class(plt),
      #   "ggplot" %in% class(plt)
      # )
    )
  }
})

test_that("plot_biplot handles non-existent color by metadata gracefully", {
  fake_data <- make_fake_gsea(n_pathways = 6, n_ranks = 5)
  pca_object <- pca_tools$do_one(fake_data)
  # pca_tools$plot_biplot(pca_object, colby = "nonexistent")
  expect_warning(pca_tools$plot_biplot(pca_object, colby = "nonexistent"), regexp = "nonexistent not found in metadata")
})

# test_that("plot_biplot creates a plot with the correct dimensions", {
#
#   fake_data <- data.frame(pathway = rep("path", 10), NES = rnorm(10), var = letters[1:10], mainpathway = sample(c(TRUE, FALSE), 10, replace = TRUE))
#   pca_object <- pca_tools$do_one(fake_data)
#   plot_result <- pca_tools$plot_biplot(pca_object, top_pc = 2)
#   expect_equal(length(plot_result), choose(2, 2))  # 2 choose 2 = 1 plot for PC1 vs PC2
# })
