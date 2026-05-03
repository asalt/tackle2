# test_fgsea.R
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))

# source("../fgsea.R")
# source("../io.R")
# source("../geneset_utils.R")

src_dir <- file.path(here("R"))

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)

geneset_tools <- new.env()
source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)

fgsea_tools_monkeypatch <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools_monkeypatch)
fgsea_tools_monkeypatch$run_one <- function(...) {
  return("success")
}

plot_tools <- new.env()
source(file.path(src_dir, "./plot.R"), local = plot_tools)

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)

sim_tools <- new.env()
source(file.path(src_dir, "./simulate.R"), local = sim_tools)




testthat::test_that("test fgsea parse additional_info", {
  #
  geneset <- geneset_tools$get_collection("H", "")
  geneset1 <- geneset_tools$get_collection("C5", "GO:MF")

  data <- sim_tools$simulate_preranked_data()
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  geneset_list1 <- geneset_tools$genesets_df_to_list(geneset1)
  geneset_lists <- list("H_" = geneset_list, "C5_GO:MF" = geneset_list1)

  ranks <- io_tools$ranks_dfs_to_lists(list(test = data))

  # genesets_info <- tibble::tribble(
  #   ~category, ~subcategory, ~collapse,
  #   "H", "", FALSE,
  #   "C3", "GO:MF", TRUE
  # ) %>% dplyr::mutate(collection_name = stringr::str_c(category, subcategory, sep = "_")) # todo put this somewhere else


  genesets_info <- tibble::tribble(
    ~category, ~subcategory,
    "H", "",
    "C5", "GO:MF",
  ) %>% dplyr::mutate(collection_name = stringr::str_c(category, subcategory, sep = "_")) # todo put this somewhere else

  testthat::expect_no_error(
    .res <- fgsea_tools_monkeypatch$run_all_pathways(
      geneset_lists = geneset_lists,
      ranks = ranks,
      parallel = FALSE,
    )
  )
  #

  #
})


test_fgsea_runone <- function() {
  # data <- .GlobalEnv$simulate_preranked_data()
  data <- sim_tools$simulate_preranked_data()
  geneset <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  rankobjs <- io_tools$ranks_dfs_to_lists(list(data))
  rankobj <- rankobjs[[1]]
  res <- rankobj %>% fgsea_tools$run_one(geneset_list)
  return(res)
}

test_that("test fgsea runone", {
  res <- test_fgsea_runone()
  expect_true(
    "data.frame" %in% class(res)
  )
})



test_that("test run one collapse", {
  #
  geneset <- geneset_tools$get_collection("C5", "GO:BP")
  spike_terms <- c("CYCLE", "CHECKPOINT")
  data <- sim_tools$simulate_preranked_data(geneset = geneset, sample_frac = .75)

  #
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  rankobjs <- io_tools$ranks_dfs_to_lists(list(data))
  rankobj <- rankobjs[[1]]

  res_withcollapse <- rankobj %>% fgsea_tools$run_one(geneset_list, collapse = TRUE)
  res_all <- rankobj %>% fgsea_tools$run_one(geneset_list, collapse = FALSE)

  testthat::expect_true( # ES should be true, NES can vary based on, .. rng?
    all(res_withcollapse$ES == res_all$ES)
  )

  .nrow1 <- res_withcollapse %>%
    dplyr::filter(mainpathway == TRUE) %>%
    nrow()
  .nrow2 <- res_all %>%
    dplyr::filter(mainpathway == TRUE) %>%
    nrow()

  testthat::expect_true(
    .nrow1 <= .nrow2
  )

  testthat::expect_true(
    res_all %>%
      dplyr::filter(mainpathway == TRUE) %>%
      nrow() == res_all %>% nrow()
  )

  testthat::expect_true( # this one should obviously be true
    res_withcollapse %>%
      dplyr::filter(mainpathway == FALSE) %>%
      nrow() > 1
  )
})

test_that("signed Hallmark scenarios recover their intended dominant programs", {
  catalog <- sim_tools$get_hallmark_scenario_catalog()
  sim <- sim_tools$simulate_signed_hallmark_scenarios(
    seed = 5150,
    n_replicates = 1
  )

  geneset <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  rankobjs <- io_tools$ranks_dfs_to_lists(sim$rank_dfs)

  expected_negative <- list(
    oxphos_metabolic = "HALLMARK_GLYCOLYSIS",
    interferon_inflammatory = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    emt_hypoxia = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  )

  for (scenario_id in names(catalog)) {
    rank_name <- sprintf("%s_rep01", scenario_id)
    scenario_res <- suppressWarnings(
      fgsea_tools$run_one(rankobjs[[rank_name]], geneset_list, collapse = FALSE)
    ) %>%
      dplyr::filter(!is.na(NES))

    top_positive <- scenario_res %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(dplyr::desc(NES)) %>%
      utils::head(5) %>%
      dplyr::pull(pathway)

    expect_true(
      any(top_positive %in% catalog[[scenario_id]]$focus_pathways),
      info = paste("Expected a focus pathway among the top positive NES results for", scenario_id)
    )

    if (scenario_id %in% names(expected_negative)) {
      top_negative <- scenario_res %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES) %>%
        utils::head(5) %>%
        dplyr::pull(pathway)

      expect_true(
        expected_negative[[scenario_id]] %in% top_negative,
        info = paste("Expected a contrast pathway among the top negative NES results for", scenario_id)
      )
    }
  }
})

test_that("signed Hallmark expression simulator preserves Hallmark themes after limma back-calculation", {
  catalog <- sim_tools$get_hallmark_scenario_catalog()
  sim <- sim_tools$simulate_signed_hallmark_expression(
    seed = 5150,
    target_source = "signal_mean",
    n_samples_per_group = 5
  )

  geneset <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  rankobjs <- io_tools$ranks_dfs_to_lists(sim$recovered_rank_dfs)

  expected_negative <- list(
    oxphos_metabolic = "HALLMARK_GLYCOLYSIS",
    interferon_inflammatory = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    emt_hypoxia = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  )

  for (scenario_id in names(catalog)) {
    scenario_res <- suppressWarnings(
      fgsea_tools$run_one(rankobjs[[scenario_id]], geneset_list, collapse = FALSE)
    ) %>%
      dplyr::filter(!is.na(NES))

    top_positive <- scenario_res %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(dplyr::desc(NES)) %>%
      utils::head(10) %>%
      dplyr::pull(pathway)

    expect_true(
      any(top_positive %in% catalog[[scenario_id]]$focus_pathways),
      info = paste("Expected a focus pathway among the recovered positive NES results for", scenario_id)
    )

    if (scenario_id %in% names(expected_negative)) {
      top_negative <- scenario_res %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES) %>%
        utils::head(10) %>%
        dplyr::pull(pathway)

      expect_true(
        expected_negative[[scenario_id]] %in% top_negative,
        info = paste("Expected a contrast pathway among the recovered negative NES results for", scenario_id)
      )
    }
  }
})

test_that("teaching Hallmark datasets preserve their intended treatment programs after limma and GSEA", {
  testthat::skip_if_not_installed("limma")

  catalog <- sim_tools$get_hallmark_teaching_dataset_catalog()
  sim <- sim_tools$simulate_hallmark_teaching_datasets(
    seed = 6404,
    n_samples_per_group_batch = 3
  )

  geneset <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)

  for (dataset_id in names(sim$datasets)) {
    ds <- sim$datasets[[dataset_id]]
    meta <- ds$sample_metadata
    meta$group <- factor(meta$group, levels = c("A", "B", "C", "D"))
    meta$batch <- factor(meta$batch, levels = unique(meta$batch))

    design <- stats::model.matrix(~ 0 + group + batch, data = meta)
    fit <- limma::lmFit(ds$expression_matrix, design)
    contrast_matrix <- limma::makeContrasts(
      D_vs_B = groupD - groupB,
      levels = design
    )
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

    top_table <- limma::topTable(
      fit,
      coef = "D_vs_B",
      number = Inf,
      sort.by = "none"
    )
    top_table$GeneID <- rownames(top_table)
    rank_df <- top_table %>%
      dplyr::transmute(
        id = GeneID,
        value = sign(logFC) * -log10(pmax(P.Value, .Machine$double.eps))
      )
    ranks <- stats::setNames(rank_df$value, rank_df$id)

    scenario_res <- suppressWarnings(
      fgsea_tools$run_one(ranks, geneset_list, collapse = FALSE)
    ) %>%
      dplyr::filter(!is.na(NES))

    expected_up <- catalog[[dataset_id]]$expected_up_pathways
    expected_down <- catalog[[dataset_id]]$expected_down_pathways

    if (length(expected_up) > 0) {
      top_positive <- scenario_res %>%
        dplyr::filter(NES > 0) %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        utils::head(10) %>%
        dplyr::pull(pathway)

      expect_true(
        any(top_positive %in% expected_up),
        info = paste("Expected an intended up pathway among recovered positive NES results for", dataset_id)
      )
    }

    if (length(expected_down) > 0) {
      top_negative <- scenario_res %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES) %>%
        utils::head(10) %>%
        dplyr::pull(pathway)

      expect_true(
        any(top_negative %in% expected_down),
        info = paste("Expected an intended down pathway among recovered negative NES results for", dataset_id)
      )
    }
  }
})

test_that("filter_on_mainpathway retains all comparison rows for globally retained pathways", {
  df <- tibble::tibble(
    pathway = c("p1", "p1", "p2", "p2"),
    rankname = c("rank1", "rank2", "rank1", "rank2"),
    NES = rnorm(4),
    mainpathway = c(TRUE, TRUE, TRUE, FALSE)
  )

  # With ratio=1, keep only pathways that are main in 100% of comparisons.
  # p2 is excluded entirely because it is not a main pathway in rank2.
  filtered_strict <- fgsea_tools$filter_on_mainpathway(df, main_pathway_ratio = 1)
  testthat::expect_true(all(filtered_strict$pathway == "p1"))
  testthat::expect_equal(sort(unique(filtered_strict$rankname)), c("rank1", "rank2"))
  testthat::expect_true(all(filtered_strict$mainpathway))

  # With ratio=0.1, p2 is retained because it is a main pathway in 1/2 comparisons.
  # Once a pathway is retained, we intentionally keep *all* of its comparison rows,
  # including rows where mainpathway == FALSE in a specific comparison.
  filtered_loose <- fgsea_tools$filter_on_mainpathway(df, main_pathway_ratio = 0.1)
  testthat::expect_equal(sort(unique(filtered_loose$pathway)), c("p1", "p2"))
  testthat::expect_equal(
    filtered_loose %>%
      dplyr::filter(pathway == "p2") %>%
      dplyr::arrange(rankname) %>%
      dplyr::pull(rankname),
    c("rank1", "rank2")
  )
  testthat::expect_true(any(filtered_loose$pathway == "p2" & filtered_loose$mainpathway == FALSE))
  testthat::expect_equal(nrow(filtered_loose), 4)
})


# ================================ test data ================================
#
TEST_DATA <- sim_tools$generate_test_data() # this uses from above like geneset_tools, io_tools, so if the above passed this should too

test_that("test concat results one collection", {
  res <- TEST_DATA
  res1 <- res[[names(res)[1]]]

  testthat::expect_true(
    "list" %in% class(res1)
  )

  res1_c <- res1 %>% fgsea_tools$concat_results_one_collection()

  testthat::expect_true(
    "data.frame" %in% class(res1_c)
  )

  testthat::expect_true(
    "rankname" %in% colnames(res1_c)
  )

  testthat::expect_true(
    all(sort(unique(res1_c$rankname)) == c("first", "second"))
  )
})


test_that("test concat results all collections", {
  res <- TEST_DATA
  res_c <- res %>% fgsea_tools$concat_results_all_collections()
  testthat::expect_true(
    "list" %in% class(res_c)
  )
})


# ===

test_that("test run all geneset lists not named.", { # this will take a while. testing if can set collapse. var
  geneset <- geneset_tools$get_collection("C5", "GO:BP")
  data <- sim_tools$simulate_preranked_data(geneset = geneset)
  data %<>% dplyr::sample_frac(size = .25)

  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  geneset_lists <- list(geneset_list)

  rankobjs <- io_tools$ranks_dfs_to_lists(list(test = data))
  # rankobj <- rankobjs[[1]]

  testthat::expect_error(
    fgsea_tools$run_all_pathways(
      geneset_lists = geneset_lsits,
      ranks = rankobjs,
      parallel = FALSE,
    )
  )
})


test_that("test run all ranks lists not named.", { # this will take a while. testing if can set collapse. var
  geneset <- geneset_tools$get_collection("C5", "GO:BP")
  data <- sim_tools$simulate_preranked_data(geneset = geneset)
  data %<>% dplyr::sample_frac(size = .25)

  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  geneset_lists <- list("C5_GO:BP" = geneset_list)

  rankobjs <- io_tools$ranks_dfs_to_lists(list(data))
  # rankobj <- rankobjs[[1]]

  testthat::expect_error(
    fgsea_tools$run_all_pathways(
      geneset_lists = geneset_lsits,
      ranks = rankobjs,
      parallel = FALSE,
    )
  )
})


#
# tests if the "collapse" argument is working as expected
test_that("test run all collapse.", { # this will take a while. testing if can set collapse. var
  geneset <- geneset_tools$get_collection("C5", "GO:BP")
  spike_terms <- c("CYCLE", "CHECKPOINT")
  data <- sim_tools$simulate_preranked_data(geneset = geneset)
  data %<>% dplyr::sample_frac(size = .25)

  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  geneset_lists <- list("C5_GO:BP" = geneset_list)

  rankobjs <- io_tools$ranks_dfs_to_lists(list(test = data))
  # rankobj <- rankobjs[[1]]

  res <- fgsea_tools$run_all_pathways(
    geneset_lists,
    rankobjs,
    collapse = TRUE
  )

  res_all <- fgsea_tools$run_all_pathways(
    geneset_lists,
    rankobjs,
    collapse = FALSE
  )

  testthat::expect_true(
    all(res[[1]]$NES == res_all[[1]]$NES)
  )

  testthat::expect_true(
    res[[1]][[1]] %>% dplyr::filter(mainpathway == TRUE) %>% nrow() <=
      res_all[[1]][[1]] %>%
        dplyr::filter(mainpathway == TRUE) %>%
        nrow()
  )

  testthat::expect_true(
    res[[1]][[1]] %>%
      dplyr::filter(mainpathway == TRUE) %>%
      nrow() > 1
  )
  #
})

# ==

test_get_edge <- function() {
  data <- sim_tools$simulate_preranked_data()
  geneset <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset)
  rankobjs <- io_tools$ranks_dfs_to_lists(list(data))
  rankobj <- rankobjs[[1]]

  res <- rankobj %>% fgsea_tools$run_one(geneset_list) # we aren't actually using this result
  # all we need for this test is the rankobj and gene list

  geneset_name <- names(geneset_list)[1]
  geneset_collection_ids <- geneset_list[[geneset_name]]

  rankorder_edge_list <- fgsea_tools$get_rankorder(geneset_collection_ids, rankobj)

  assertthat::assert_that(
    "list" %in% class(rankorder_edge_list)
  )

  assertthat::assert_that(
    "edge" %in% names(rankorder_edge_list)
  )
  # now look at edge

  edge <- rankorder_edge_list$edge

  assertthat::assert_that(
    "data.frame" %in% class(edge)
  )

  assertthat::assert_that(
    "id" %in% names(edge)
  )

  .expected_names <- c("id", "rank", "stat", "ES", "stat_tick") # , "stat_stat")
  for (name in .expected_names) {
    assertthat::assert_that(
      name %in% names(edge)
    )
  }



  edge_specific <- edge %>% filter(id %in% geneset_collection_ids)

  assertthat::noNA(edge_specific$stat_tick)
  assertthat::assert_that(
    all(edge_specific$stat == edge_specific$stat_tick),
    TRUE
  )
  return("Success")
}

test_that("test get edge", {
  expect_equal(test_get_edge(), "Success")
})
