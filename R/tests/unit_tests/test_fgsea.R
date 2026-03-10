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

test_that("filter_on_mainpathway filters per-pathway across comparisons", {
  df <- tibble::tibble(
    pathway = c("p1", "p1", "p2", "p2"),
    rankname = c("rank1", "rank2", "rank1", "rank2"),
    NES = rnorm(4),
    mainpathway = c(TRUE, TRUE, TRUE, FALSE)
  )

  # With ratio=1, keep only pathways that are main in 100% of ranknames.
  filtered_strict <- fgsea_tools$filter_on_mainpathway(df, main_pathway_ratio = 1)
  testthat::expect_true(all(filtered_strict$pathway == "p1"))
  testthat::expect_equal(sort(unique(filtered_strict$rankname)), c("rank1", "rank2"))
  testthat::expect_true(all(filtered_strict$mainpathway))

  # With ratio=0.1, keep any pathway that is main in >=10% of ranknames,
  # but only retain rows where mainpathway==TRUE.
  filtered_loose <- fgsea_tools$filter_on_mainpathway(df, main_pathway_ratio = 0.1)
  testthat::expect_true(all(filtered_loose$mainpathway))
  testthat::expect_true(all(c("p1", "p2") %in% unique(filtered_loose$pathway)))
  testthat::expect_true(!any(filtered_loose$pathway == "p2" & filtered_loose$rankname == "rank2"))
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
