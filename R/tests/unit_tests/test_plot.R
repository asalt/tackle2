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

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)

plot_tools <- new.env()
source(file.path(src_dir, "./plot.R"), local = plot_tools)

bubble_tools <- new.env()
source(file.path(src_dir, "./plot_bubble.R"), local = bubble_tools)

sim_tools <- new.env()
source(file.path(src_dir, "./simulate.R"), local = sim_tools)
generate_test_data <- sim_tools$generate_test_data


# ================================ test data ================================
TEST_DATA <- generate_test_data()
TEST_DATA_withCollapse <- generate_test_data(collapse = TRUE)

make_main_pathway_plot_results <- function() {
  make_comparison <- function(pathway, NES, mainpathway) {
    tibble::tibble(
      pathway = pathway,
      pval = rep(0.01, length(pathway)),
      padj = rep(0.02, length(pathway)),
      ES = NES / 10,
      NES = NES,
      size = rep(20L, length(pathway)),
      leadingEdge = rep(list(c("GENE1", "GENE2")), length(pathway)),
      mainpathway = mainpathway
    )
  }

  list(
    TEST_COLLECTION = list(
      comparison_A = make_comparison(
        c("REDUNDANT_STRONG", "SHARED_MAIN", "MAIN_A_ONLY"),
        c(9, 4, 3),
        c(FALSE, TRUE, TRUE)
      ),
      comparison_B = make_comparison(
        c("REDUNDANT_STRONG", "SHARED_MAIN", "MAIN_B_ONLY"),
        c(8, 3, 2),
        c(FALSE, FALSE, TRUE)
      )
    )
  )
}

# ================================ tests ================================


test_that("test formatting for barplot", {
  data <- TEST_DATA[[1]][[1]]
  out <- plot_tools$prepare_data_for_barplot(data)

  testthat::expect_true(
    all(c("leadingEdgeNum", "leadingEdgeFrac") %in% colnames(out))
  )
})

test_that("barplot selection variants suffix saved filenames", {
  saved_files <- character(0)
  fake_save <- function(plot_code, path = "", filename = "", width = NULL, height = NULL, ...) {
    saved_files <<- c(saved_files, file.path(path, paste0(filename, ".pdf")))
    TRUE
  }
  fake_save <- util_tools$make_partial(fake_save, path = tempdir(), filename = "base")

  plot_tools$all_barplots_with_numbers(
    list(H_ = list(sample_A = TEST_DATA$H_[[1]])),
    save_func = fake_save,
    limit = 1,
    sort_by = "pval",
    variant_name = "fdr25_pval"
  )

  testthat::expect_true(any(grepl("fdr25_pval", basename(saved_files), fixed = TRUE)))
})

test_that("individual bar selection honors the collection collapse setting", {
  results <- make_main_pathway_plot_results()

  collapsed <- plot_tools$all_barplots_with_numbers(
    results,
    limit = 1,
    collapse = TRUE
  )$TEST_COLLECTION$comparison_A[[1]]
  uncollapsed <- plot_tools$all_barplots_with_numbers(
    results,
    limit = 1,
    collapse = FALSE
  )$TEST_COLLECTION$comparison_A[[1]]

  expect_equal(unique(as.character(collapsed$data$pathway)), "SHARED MAIN")
  expect_true(all(collapsed$data$mainpathway))
  expect_equal(unique(as.character(uncollapsed$data$pathway)), "REDUNDANT STRONG")
  expect_false(all(uncollapsed$data$mainpathway))
})

test_that("individual bubble selection honors the collection collapse setting", {
  results <- make_main_pathway_plot_results()

  collapsed <- bubble_tools$all_bubble_plots(
    results,
    limit = 1,
    collapse = TRUE
  )$TEST_COLLECTION$comparison_A[[1]]
  uncollapsed <- bubble_tools$all_bubble_plots(
    results,
    limit = 1,
    collapse = FALSE
  )$TEST_COLLECTION$comparison_A[[1]]

  expect_equal(unique(as.character(collapsed$data$pathway)), "SHARED MAIN")
  expect_true(all(collapsed$data$mainpathway))
  expect_equal(unique(as.character(uncollapsed$data$pathway)), "REDUNDANT STRONG")
  expect_false(all(uncollapsed$data$mainpathway))
})

test_that("combined bar selection retains a main pathway across every facet", {
  results <- make_main_pathway_plot_results()

  collapsed <- plot_tools$do_combined_barplots(
    results,
    limit = 1,
    collapse = TRUE,
    main_pathway_ratio = 0.5
  )[[1]][[1]]
  uncollapsed <- plot_tools$do_combined_barplots(
    results,
    limit = 1,
    collapse = FALSE,
    main_pathway_ratio = 0.5
  )[[1]][[1]]

  expect_equal(unique(as.character(collapsed$data$pathway)), "SHARED MAIN")
  expect_equal(dplyr::n_distinct(collapsed$data$rankname), 2L)
  expect_true(any(!collapsed$data$mainpathway))
  expect_equal(unique(as.character(uncollapsed$data$pathway)), "REDUNDANT STRONG")
})

test_that("combined bubble selection retains a main pathway across every facet", {
  results <- make_main_pathway_plot_results()

  collapsed <- bubble_tools$do_combined_bubble_plots(
    results,
    limit = 1,
    collapse = TRUE,
    main_pathway_ratio = 0.5
  )[[1]][[1]]
  uncollapsed <- bubble_tools$do_combined_bubble_plots(
    results,
    limit = 1,
    collapse = FALSE,
    main_pathway_ratio = 0.5
  )[[1]][[1]]

  expect_equal(unique(as.character(collapsed$data$pathway)), "SHARED MAIN")
  expect_equal(dplyr::n_distinct(collapsed$data$rankname), 2L)
  expect_true(any(!collapsed$data$mainpathway))
  expect_equal(unique(as.character(uncollapsed$data$pathway)), "REDUNDANT STRONG")
})

test_that("bar and bubble plots write individual and combined pathway count sidecars", {
  output_dir <- withr::local_tempdir()
  results <- make_main_pathway_plot_results()
  rank_metadata <- data.frame(
    rankname = c("comparison_A", "comparison_B", "comparison_missing"),
    rank_label = c("Rank A", "Rank B", "Missing Rank"),
    rank_order = 1:3,
    rank_label_source = "names.txt",
    stringsAsFactors = FALSE
  )
  fake_save <- function(plot_code, path = "", filename = "", ...) {
    TRUE
  }
  save_func <- util_tools$make_partial(
    fake_save,
    path = output_dir,
    filename = "base",
    replace = TRUE
  )

  plot_tools$all_barplots_with_numbers(
    results,
    save_func = save_func,
    limit = 1,
    collapse = TRUE,
    rank_metadata = rank_metadata
  )
  plot_tools$do_combined_barplots(
    results,
    save_func = save_func,
    limit = 1,
    collapse = TRUE,
    main_pathway_ratio = 0.5,
    rank_metadata = rank_metadata
  )
  bubble_tools$all_bubble_plots(
    results,
    save_func = save_func,
    limit = 1,
    collapse = TRUE,
    rank_metadata = rank_metadata
  )
  bubble_tools$do_combined_bubble_plots(
    results,
    save_func = save_func,
    limit = 1,
    collapse = TRUE,
    main_pathway_ratio = 0.5,
    rank_metadata = rank_metadata
  )

  bar_individual <- file.path(
    output_dir, "TEST_COLLECTION", "bar", "Rank_A", "pathway_counts.json"
  )
  bar_missing <- file.path(
    output_dir, "TEST_COLLECTION", "bar", "Missing_Rank", "pathway_counts.json"
  )
  bar_combined <- file.path(
    output_dir, "TEST_COLLECTION", "bar", "pathway_counts.json"
  )
  bubble_individual <- file.path(
    output_dir, "TEST_COLLECTION", "bubble", "Rank_A", "pathway_counts.json"
  )
  bubble_missing <- file.path(
    output_dir, "TEST_COLLECTION", "bubble", "Missing_Rank", "pathway_counts.json"
  )
  bubble_combined <- file.path(
    output_dir, "TEST_COLLECTION", "bubble", "pathway_counts.json"
  )

  expect_true(all(file.exists(c(
    bar_individual,
    bar_missing,
    bar_combined,
    bubble_individual,
    bubble_missing,
    bubble_combined
  ))))

  individual <- jsonlite::read_json(bar_individual, simplifyVector = FALSE)
  missing <- jsonlite::read_json(bar_missing, simplifyVector = FALSE)
  combined_bar <- jsonlite::read_json(bar_combined, simplifyVector = FALSE)
  combined_bubble <- jsonlite::read_json(bubble_combined, simplifyVector = FALSE)

  expect_equal(names(individual$statistics$ranks), "comparison_A")
  expect_equal(individual$statistics$n_retained_pathways, 2L)
  expect_equal(
    missing$statistics$ranks$comparison_missing$n_pathways,
    0L
  )
  expect_equal(
    names(combined_bar$statistics$ranks),
    rank_metadata$rankname
  )
  expect_equal(combined_bar$statistics$n_retained_pathways, 3L)
  expect_equal(combined_bubble$statistics, combined_bar$statistics)
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

test_that("enrichplot limit zero still plots pathways of interest", {
  ranks_list <- list(
    sample_A = stats::setNames(c(3, 2, 1, -1, -2, -3), paste0("gene", 1:6)),
    sample_B = stats::setNames(c(-3, -2, -1, 1, 2, 3), paste0("gene", 1:6))
  )
  geneset_collection <- list(
    TOP_PATHWAY = c("gene1", "gene2"),
    POI_PATHWAY = c("gene5", "gene6")
  )
  df <- tidyr::expand_grid(
    pathway = names(geneset_collection),
    rankname = names(ranks_list)
  ) %>%
    dplyr::mutate(
      pval = 0.01,
      padj = 0.02,
      ES = dplyr::if_else(.data$pathway == "TOP_PATHWAY", 0.8, 0.2),
      NES = dplyr::if_else(.data$pathway == "TOP_PATHWAY", 3.0, 0.2),
      size = 2L,
      mainpathway = .data$pathway == "TOP_PATHWAY"
    )

  plts <- plot_tools$plot_top_ES(
    df = df,
    ranks_list = ranks_list,
    geneset_collection = geneset_collection,
    limit = 0,
    do_individual = TRUE,
    do_combined = FALSE,
    pathways_of_interest = "POI_PATHWAY"
  )

  testthat::expect_equal(names(plts), "POI_PATHWAY")
  testthat::expect_equal(names(plts$POI_PATHWAY), names(ranks_list))
  testthat::expect_true(all(vapply(plts$POI_PATHWAY, inherits, logical(1), "ggplot")))
})

test_that("enrichplot skips missing pathway-rank stats rows before dplyr filtering", {
  testthat::skip_if_not_installed("data.table")
  ranks_list <- list(
    sample_A = stats::setNames(c(3, 2, 1, -1, -2, -3), paste0("gene", 1:6)),
    sample_B = stats::setNames(c(-3, -2, -1, 1, 2, 3), paste0("gene", 1:6))
  )
  geneset_collection <- list(ONLY_PATHWAY = c("gene1", "gene2"))
  df <- data.table::data.table(
    pathway = "ONLY_PATHWAY",
    rankname = "sample_A",
    pval = 0.01,
    padj = 0.02,
    ES = 0.8,
    NES = 3.0,
    size = 2L,
    mainpathway = TRUE
  )

  expect_equal(
    plot_tools$match_enrichplot_stats_rows(df, "ONLY_PATHWAY", "sample_B"),
    integer(0)
  )
  expect_no_error(
    plts <- plot_tools$plot_top_ES(
      df = df,
      ranks_list = ranks_list,
      geneset_collection = geneset_collection,
      limit = 1,
      do_individual = TRUE,
      do_combined = FALSE
    )
  )

  expect_true(inherits(plts$ONLY_PATHWAY$sample_A, "ggplot"))
  expect_null(plts$ONLY_PATHWAY$sample_B)
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
