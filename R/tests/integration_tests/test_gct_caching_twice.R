suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))

source(file.path(here("R"), "lazyloader.R"))
sim_tools <- get_tool_env("simulate")
run_tools <- get_tool_env("run")

make_twice_test_gct <- function(data_dir) {
  fs::dir_create(data_dir, recurse = TRUE)

  # Small synthetic GCT with repeated group structure to exercise GCT rank loading.
  datas1 <- purrr::map(1:2, ~ sim_tools$simulate_preranked_data(seed = 4321))
  datas2 <- purrr::map(1:2, ~ sim_tools$simulate_preranked_data(seed = 1234))
  datas3 <- purrr::map(1:2, ~ sim_tools$simulate_preranked_data(seed = 9999))
  datas <- c(datas1, datas2, datas3)
  names(datas) <- c(paste0("group_A_", seq_along(datas1)),
                    paste0("group_B_", seq_along(datas2)),
                    paste0("group_C_", seq_along(datas3)))

  .mat <- Reduce(function(...) full_join(..., by = "id"), datas) %>% as.data.frame()
  colnames(.mat) <- c("id", names(datas))
  rownames(.mat) <- .mat$id
  .mat$id <- NULL
  .meta <- data.frame(id = colnames(.mat), group = c(rep("A", 2), rep("B", 2), rep("C", 2)))
  rownames(.meta) <- .meta$id
  .rdesc <- data.frame(id = rownames(.mat), dummy = "X")
  rownames(.rdesc) <- rownames(.mat)
  gct <- new("GCT", mat = as.matrix(.mat), cdesc = .meta, rdesc = .rdesc)
  gct_path <- file.path(data_dir, "twice_test.gct")
  cmapR::write_gct(gct, gct_path, appenddim = FALSE)
  gct_path
}

list_cache_rds <- function(cache_dir) {
  fs::dir_ls(cache_dir, recurse = TRUE, type = "file", regexp = "\\.rds$", fail = FALSE)
}

test_that("gct run populates fgsea cache and reruns preserve cache files", {
  root_dir <- fs::path(
    tempdir(),
    sprintf("gsea-digest-gct-cache-%s-%04d", Sys.getpid(), sample.int(9999, 1))
  )
  data_dir <- fs::path(root_dir, "data")
  output_dir <- fs::path(root_dir, "output")
  cache_dir <- fs::path(output_dir, "cache_tmp")
  on.exit({
    if (fs::dir_exists(root_dir)) {
      fs::dir_delete(root_dir)
    }
  }, add = TRUE)

  fs::dir_create(output_dir, recurse = TRUE)
  fs::dir_create(cache_dir, recurse = TRUE)
  gct_path <- make_twice_test_gct(data_dir)

  params <- list(
    rankfiledir = data_dir,
    savedir = output_dir,
    gct_path = gct_path,
    ranks_from = "gct",
    zscore_emat = TRUE,
    combine_by = "group",
    genesets = list(list(category = "H", subcategory = "", collapse = FALSE)),
    barplot = list(do_individual = FALSE, do_combined = FALSE),
    enplot = list(do_individual = FALSE, do_combined = FALSE),
    heatmap_gene = list(do = FALSE),
    heatmap_gsea = list(do = FALSE),
    pca = list(do = FALSE),
    bubbleplot = list(do_individual = FALSE, do_combined = FALSE),
    advanced = list(
      cache = TRUE,
      cachedir = cache_dir,
      replace = TRUE,
      quiet = TRUE,
      parallel = FALSE,
      logfile = fs::path(output_dir, "run.log")
    )
  )

  expect_no_error(run_tools$run(params))
  cache_files_1 <- list_cache_rds(cache_dir)
  expect_gt(length(cache_files_1), 0L)

  expect_no_error(run_tools$run(params))
  cache_files_2 <- list_cache_rds(cache_dir)
  expect_true(all(basename(cache_files_1) %in% basename(cache_files_2)))
  expect_gte(length(cache_files_2), length(cache_files_1))

  expect_no_error(run_tools$run(params))
  cache_files_3 <- list_cache_rds(cache_dir)
  expect_setequal(basename(cache_files_3), basename(cache_files_2))
})
