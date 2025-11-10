suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(here))

# Load helpers similar to other unit tests
src_dir <- file.path(here("R"))

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)

geneset_tools <- new.env()
source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)

# Create a fresh fgsea tools env and monkeypatch run_one to always return NULL
fgsea_tools_mp <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools_mp)
fgsea_tools_mp$run_one <- function(...) {
  return(NULL)
}

sim_tools <- new.env()
source(file.path(src_dir, "./simulate.R"), local = sim_tools)

test_that("run_all_rankobjs does not crash when some/all results are NULL (cache save alignment)", {
  # Prepare a small geneset list (any content is fine since run_one is mocked)
  geneset_df <- geneset_tools$get_collection("H", "")
  geneset_list <- geneset_tools$genesets_df_to_list(geneset_df)

  # Two small rank objects (named) to ensure length > 0
  d1 <- sim_tools$simulate_preranked_data(seed = 101)
  d2 <- sim_tools$simulate_preranked_data(seed = 202)
  ranks <- io_tools$ranks_dfs_to_lists(list(a = d1, b = d2))

  # Use a temp cache dir so the test is isolated
  tmp_cache <- tempdir(check = TRUE)

  expect_no_error({
    res <- fgsea_tools_mp$run_all_rankobjs(
      pathway = geneset_list,
      rankobjs = ranks,
      cache = TRUE,
      cache_dir = tmp_cache
    )
    # With NULL results across the board, we should get an empty list
    expect_true(is.list(res))
    expect_equal(length(res), 0)
  })
})

