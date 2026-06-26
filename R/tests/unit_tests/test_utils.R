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

testthat::test_that("make_name_map preserves common suffixes by default", {
  labels <- c(
    "groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
    "groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
  )

  name_map <- util_tools$make_name_map(labels, delim = "[._\\s]")

  testthat::expect_equal(unname(as.character(name_map)), labels)
  testthat::expect_equal(attr(name_map, "removed_suffix"), "")
})

testthat::test_that("make_name_map still removes common prefixes", {
  labels <- c(
    "RUN_C3D9-Definitely_Not_Final_groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
    "RUN_C3D9-Definitely_Not_Final_groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
  )

  name_map <- util_tools$make_name_map(labels, delim = "[._\\s]")

  testthat::expect_equal(
    unname(as.character(name_map)),
    c(
      "groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
      "groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
    )
  )
  testthat::expect_equal(attr(name_map, "removed_prefix"), "RUN_C3D9-Definitely_Not_Final_")
  testthat::expect_equal(attr(name_map, "removed_suffix"), "")
})

testthat::test_that("make_name_map suffix stripping is explicit opt-in", {
  labels <- c("groupA_vs_control", "groupB_vs_control")

  default_map <- util_tools$make_name_map(labels, delim = "[._\\s]")
  suffix_map <- util_tools$make_name_map(labels, delim = "[._\\s]", strip_suffix = TRUE)

  testthat::expect_equal(unname(as.character(default_map)), labels)
  testthat::expect_equal(unname(as.character(suffix_map)), c("groupA", "groupB"))
  testthat::expect_equal(attr(suffix_map, "removed_suffix"), "_vs_control")
})

testthat::test_that("safe_filename normalizes troublesome label characters", {
  labels <- c(
    "A4B8-Sunk_My_Battleship",
    "A4B8-Sunk My Battleship",
    "A4B8-Sunk.My.Battleship",
    "A4B8-Sunk(My)Battleship",
    "A4B8-Sunk+My+Battleship",
    "RUN_C3D9-Definitely_Not_Final",
    "BATCH_X9Q2-Hotdog_Is_A_Sandwich",
    "GROUP_A4B8-Sunk__My__Battleship",
    "GROUP_A4B8-_Leading_Underscore-ish",
    "GROUP_A4B8-Trailing_Dash-",
    "GROUP_A4B8-Multiple---Dashes",
    "GROUP_A4B8-Bjork_Gudmundsdottir",
    "GROUP_A4B8-Cafe_Munchen",
    "GROUP_A4B8-Space At The Function",
    "GROUP_A4B8-Final_FINAL_v2_REALLY_FINAL"
  )

  sanitized <- stats::setNames(
    vapply(labels, util_tools$safe_filename, character(1), fallback = "label"),
    labels
  )

  testthat::expect_true(all(nzchar(sanitized)))
  testthat::expect_false(any(grepl("[[:space:]()+]", sanitized)))
  testthat::expect_equal(sanitized[["A4B8-Sunk My Battleship"]], "A4B8-Sunk_My_Battleship")
  testthat::expect_equal(sanitized[["A4B8-Sunk.My.Battleship"]], "A4B8-Sunk.My.Battleship")
  testthat::expect_equal(sanitized[["A4B8-Sunk(My)Battleship"]], "A4B8-Sunk_My_Battleship")
  testthat::expect_equal(sanitized[["A4B8-Sunk+My+Battleship"]], "A4B8-Sunk_My_Battleship")
  testthat::expect_equal(sanitized[["GROUP_A4B8-Space At The Function"]], "GROUP_A4B8-Space_At_The_Function")
  testthat::expect_equal(
    sanitized[["GROUP_A4B8-Final_FINAL_v2_REALLY_FINAL"]],
    "GROUP_A4B8-Final_FINAL_v2_REALLY_FINAL"
  )
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

testthat::test_that("plot selection variants keep default and suffix extra variants", {
  variants <- util_tools$normalize_plot_selection_variants(list(
    pstat_cutoff = 1,
    pstat_usetype = "fdr",
    sort_by = "nes",
    variants = list(
      list(name = "fdr25", pstat_cutoff = 0.25),
      list(sort_by = "pvalue")
    )
  ))

  testthat::expect_equal(length(variants), 3)
  testthat::expect_equal(variants[[1]]$name, "")
  testthat::expect_equal(variants[[1]]$pstat_usetype, "padj")
  testthat::expect_equal(variants[[1]]$sort_by, "NES")
  testthat::expect_equal(variants[[2]]$name, "fdr25")
  testthat::expect_equal(variants[[2]]$pstat_cutoff, 0.25)
  testthat::expect_equal(variants[[3]]$name, "by_pval")
  testthat::expect_equal(variants[[3]]$sort_by, "pval")
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

testthat::test_that("clean_args resolves savedir logfile sentinel", {
  root_dir <- tempdir()
  params <- list(
    savedir = "plots",
    advanced = list(logfile = "savedir"),
    genesets = list(list(category = "H", subcategory = "", collapse = FALSE))
  )

  cleaned <- util_tools$clean_args(params, root_dir = root_dir)

  testthat::expect_equal(
    cleaned$advanced$logfile,
    file.path(root_dir, "plots", "run.log")
  )
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
