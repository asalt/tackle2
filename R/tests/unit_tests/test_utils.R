suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))
# ==

testthat::test_that("test_create_env", {
  # Create a new environment

  src_dir <- file.path(here("R"))
  util_tools <- new.env()
  source(file.path(src_dir, "./utils.R"), local = util_tools)

  testthat::expect_true(
    !is.null(util_tools$make_partial) && ("function" %in% class(util_tools$make_partial))
  )

  testthat::expect_true(
    !is.null(util_tools$get_arg) && ("function" %in% class(util_tools$get_arg))
  )

  make_partial <- util_tools$make_partial
})

# ==

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)

make_partial <- util_tools$make_partial
get_arg <- util_tools$get_arg
get_args <- util_tools$get_args




testthat::test_that("test get_args", {
  f <- function() {}
  attr(f, "preset_args") <- list(a = 1, b = 2)

  expect_equal(get_args(f), list(a = 1, b = 2))
  expect_equal(get_args(function() {}), list()) # Should return an empty list for functions without preset_args
})



testthat::test_that("test scale gct", {
  gct <- io_tools$make_random_gct(10, 6)
  newgct <- util_tools$scale_gct(gct)
  testthat::expect_equal(
    gct@mat %>% dim(),
    newgct@mat %>% dim()
  )
})

testthat::test_that("test scale gct group_by NA does not error", {
  gct <- io_tools$make_random_gct(10, 6)
  newgct <- util_tools$scale_gct(gct, group_by = NA)
  testthat::expect_equal(
    gct@mat %>% dim(),
    newgct@mat %>% dim()
  )
})

testthat::test_that("process_cut_by treats NA as NULL", {
  cdesc <- data.frame(group = c("A", "B"), stringsAsFactors = FALSE)
  rownames(cdesc) <- c("s1", "s2")
  testthat::expect_warning(
    util_tools$process_cut_by(NA, cdesc),
    NA
  )
  testthat::expect_null(util_tools$process_cut_by(NA, cdesc))
  testthat::expect_null(util_tools$process_cut_by(FALSE, cdesc))
})

testthat::test_that("clean_args normalizes numeric do flags", {
  params <- list(
    savedir = "plots",
    genesets = list(list(category = "H", subcategory = "", collapse = FALSE)),
    barplot = list(do_individual = 2, do_combined = 0),
    bubbleplot = list(do_individual = 0, do_combined = 2),
    heatmap_gsea = list(do = 0),
    heatmap_gene = list(do = 2),
    enplot = list(do_individual = 2, do_combined = 0),
    pca = list(do = 2)
  )

  cleaned <- util_tools$clean_args(params, root_dir = tempdir())

  testthat::expect_true(isTRUE(cleaned$barplot$do_individual))
  testthat::expect_false(isTRUE(cleaned$barplot$do_combined))
  testthat::expect_false(isTRUE(cleaned$bubbleplot$do_individual))
  testthat::expect_true(isTRUE(cleaned$bubbleplot$do_combined))
  testthat::expect_false(isTRUE(cleaned$heatmap_gsea$do))
  testthat::expect_true(isTRUE(cleaned$heatmap_gene$do))
  testthat::expect_true(isTRUE(cleaned$pca$do))
  testthat::expect_true(isTRUE(cleaned$enplot$do_individual))
  testthat::expect_false(isTRUE(cleaned$enplot$do_combined))
})

testthat::test_that("clean_args defaults enplot toggles to TRUE and honors enplot.do fallback", {
  base_params <- list(
    savedir = "plots",
    genesets = list(list(category = "H", subcategory = "", collapse = FALSE))
  )

  cleaned_default <- util_tools$clean_args(base_params, root_dir = tempdir())
  testthat::expect_true(isTRUE(cleaned_default$enplot$do_individual))
  testthat::expect_true(isTRUE(cleaned_default$enplot$do_combined))

  cleaned_do_false <- util_tools$clean_args(
    c(base_params, list(enplot = list(do = FALSE))),
    root_dir = tempdir()
  )
  testthat::expect_false(isTRUE(cleaned_do_false$enplot$do_individual))
  testthat::expect_false(isTRUE(cleaned_do_false$enplot$do_combined))

  cleaned_override <- util_tools$clean_args(
    c(base_params, list(enplot = list(do = FALSE, do_individual = TRUE, do_combined = FALSE))),
    root_dir = tempdir()
  )
  testthat::expect_true(isTRUE(cleaned_override$enplot$do_individual))
  testthat::expect_false(isTRUE(cleaned_override$enplot$do_combined))
})

# logging

testthat::test_that("test log msg", {
  .f <- "test.log"
  testthat::expect_error(util_tools$log_msg(info = "test", filename = .f), NA)
  # This is the standard way to assert that a block of code should execute without any errors in testthat.
  testthat::expect_true(fs::file_exists(.f))
  testthat::expect_true(stringr::str_detect(readLines(.f), "INFO"))
  fs::file_delete(.f)
})


testthat::test_that("test log msg levels", {
  .f <- "test.log"
  if (fs::file_exists(.f)) fs::file_delete(.f)
  util_tools$log_msg(debug = "test", filename = .f, loglevel = "DEBUG")
  testthat::expect_true(
    stringr::str_detect(readLines(.f), "DEBUG")
  )
  fs::file_delete(.f)
})
