# test_plot.R
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))

io_tools <- new.env()
source(file.path(file.path(src_dir, "./io.R")), local = io_tools)

geneset_tools <- new.env()
source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)

plot_tools <- new.env()
source(file.path(src_dir, "./plot.R"), local = plot_tools)

sim_tools <- new.env()
source(file.path(src_dir, "./simulate.R"), local = sim_tools)
generate_test_data <- sim_tools$generate_test_data


# ================================ test data ================================
TEST_DATA <- generate_test_data()
TEST_DATA_withCollapse <- generate_test_data(collapse = TRUE)

# ================================ tests ================================


test_that("test formatting for barplot", {
  data <- TEST_DATA[[1]][[1]]
  out <- plot_tools$prepare_data_for_barplot(data)

  testthat::expect_true(
    all(c("leadingEdgeNum", "leadingEdgeFrac") %in% colnames(out))
  )
})

test_that("enrichplot combined and individual directories share the same pathway slug", {
  combined_dir <- plot_tools$make_enrichplot_dirname(
    "enrichplots",
    "KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_SPINDLE_ASSEMBLY",
    fallback = "combined_dir"
  )
  individual_dir <- plot_tools$make_enrichplot_dirname(
    "enrichplots",
    "KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_SPINDLE_ASSEMBLY",
    fallback = "enrich"
  )

  testthat::expect_equal(combined_dir, individual_dir)
  testthat::expect_match(combined_dir, "^enrichplots_PATHOGEN_HTLV_1_TAX_TO_SPINDLE_ASSEMBLY")
})


test_that("test plot a single barplot", {
  plt <- TEST_DATA$H_[[1]] %>%
    plot_tools$barplot_with_numbers()
  testthat::expect_true(
    all(
      "gg" %in% class(plt),
      "ggplot" %in% class(plt)
    )
  )
  plt_b <- ggplot2::ggplot_build(plt)
  testthat::expect_true(
    "facet" %in% names(plt_b$layout)
  )
  testthat::expect_true(
    length(plt_b$layout$facet$params) == 0
  )
})

test_that("test plot a faceted barplot", {
  df <- TEST_DATA$H_ %>% fgsea_tools$concat_results_one_collection()
  # this should be tested in test_fgsea but we can test it here too
  testthat::expect_true(
    "data.frame" %in% class(df),
    info = "df is not a data.frame. this is an fgsea_tools problem"
  )
  testthat::expect_true(
    "rankname" %in% colnames(df),
    info = "rankname not in colnames of df. this is an fgsea_tools problem"
  )

  plt <- df %>% plot_tools$barplot_with_numbers()
  testthat::expect_true(
    all(
      "gg" %in% class(plt),
      "ggplot" %in% class(plt)
    ),
    info = "plt is not a ggplot object"
  )
  plt_b <- ggplot2::ggplot_build(plt)

  testthat::expect_true(
    "facet" %in% names(plt_b$layout)
  )
  testthat::expect_true(
    length(plt_b$layout$facet$params) > 0,
    info = "the plot is not faceted, but it should be as there are multiple ranknames (so it is a long form table)"
  )
})


test_that("test heatmap of NES", {
  gsea_res <- TEST_DATA
  res_c <- gsea_res %>% fgsea_tools$concat_results_all_collections()

  # ht <- res_c[[1]] %>% suppressWarnings({plot_tools$plot_results_one_collection()})
  suppressWarnings({
    ht <- res_c[[1]] %>% plot_tools$plot_results_one_collection()
  })

  testthat::expect_true(
    "HeatmapList" %in% class(ht)
  )
})

test_that("test run though all heatmaps of NES", {
  gsea_res <- TEST_DATA
  res_c <- gsea_res %>% fgsea_tools$concat_results_all_collections()
  suppressWarnings({
    ht_list <- plot_tools$plot_results_all_collections(res_c)
  })
  testthat::expect_true(
    "list" %in% class(ht_list)
  )


  testthat::expect_equal(
    names(ht_list), names(gsea_res)
  )


  for (collection in ht_list) {
    for (ht in collection) {
      testthat::expect_true(
        "HeatmapList" %in% class(ht)
      )
    }
  }
})
