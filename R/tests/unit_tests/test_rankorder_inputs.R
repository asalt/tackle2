suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))

io_tools <- new.env(); source(file.path(src_dir, "./io.R"), local = io_tools)
geneset_tools <- new.env(); source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)
fgsea_tools <- new.env(); source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)
sim_tools <- new.env(); source(file.path(src_dir, "./simulate.R"), local = sim_tools)

test_that("get_rankorder handles list-wrapped inputs and returns expected structure", {
  # geneset
  gs_df <- geneset_tools$get_collection("H", "")
  gs_list <- geneset_tools$genesets_df_to_list(gs_df)
  gs_name <- names(gs_list)[1]
  gs_ids <- gs_list[[gs_name]]

  # rankobj
  d <- sim_tools$simulate_preranked_data(seed = 2025)
  ranks <- io_tools$ranks_dfs_to_lists(list(one = d))
  rankvec <- ranks[[1]]

  # Wrap both as lists to exercise coercion
  expect_no_error({
    out <- fgsea_tools$get_rankorder(list(gs_ids), list(rankvec))
    expect_true(is.list(out))
    expect_true(all(c("edge", "curve", "ticks", "stats") %in% names(out)))
    expect_true(is.data.frame(out$edge))
  })
})

test_that("get_rankorder errors for unnamed numeric rankobj", {
  # unnamed numeric vector should fail
  rankvec <- stats::rnorm(10)
  gs_ids <- as.character(seq_len(5))
  expect_error(
    fgsea_tools$get_rankorder(gs_ids, rankvec),
    regexp = "named numeric vector"
  )
})

