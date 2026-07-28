suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))
plot_utils <- new.env()
source(file.path(src_dir, "./plot_utils.R"), local = plot_utils)


util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)
make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)


context("Utility functions")

test_that("get_args returns correct preset arguments", {
  f <- function() {}
  attr(f, "preset_args") <- list(a = 1, b = 2)

  expect_equal(get_args(f), list(a = 1, b = 2))
  expect_equal(get_args(function() {}), list()) # Should return an empty list for functions without preset_args
})

test_that("get_arg returns correct argument values", {
  f <- function() {}
  attr(f, "preset_args") <- list(a = 1, b = 2)

  expect_equal(get_arg(f, "a"), 1)
  expect_equal(get_arg(f, "b"), 2)
  expect_equal(get_arg(f, "c"), "") # Test for non-existent arg should return ""
})

test_that("pathway count summaries report rank and retained pathway statistics", {
  rank_a <- data.frame(
    pathway = c("A", "B", "B", "C"),
    padj = c(0.01, 0.20, 0.30, NA_real_),
    mainpathway = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  rank_b <- data.frame(
    pathway = c("A", "D"),
    padj = c(0.40, 0.04),
    mainpathway = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  summary <- plot_utils$make_pathway_count_summary(
    results_by_rank = list(rank_a = rank_a, rank_b = rank_b),
    rank_labels = c(rank_a = "Rank A", rank_b = "Rank B", rank_missing = "Missing"),
    expected_ranknames = c("rank_a", "rank_b", "rank_missing"),
    collapse = TRUE,
    combined = TRUE,
    main_pathway_ratio = 0.5,
    generated_at = as.POSIXct("2026-07-27 12:00:00", tz = "UTC")
  )

  expect_equal(summary$statistics$n_retained_pathways, 3L)
  expect_equal(summary$statistics$ranks$rank_a, list(
    label = "Rank A",
    n_pathways = 3L,
    n_main_pathways = 2L,
    n_main_padj_lt_0_25 = 2L,
    n_main_padj_lt_0_05 = 1L
  ))
  expect_equal(summary$statistics$ranks$rank_missing, list(
    label = "Missing",
    n_pathways = 0L,
    n_main_pathways = 0L,
    n_main_padj_lt_0_25 = 0L,
    n_main_padj_lt_0_05 = 0L
  ))
  expect_equal(summary$schema_version, 1L)
  expect_equal(summary$generated_at, "2026-07-27T12:00:00Z")
})

test_that("non-collapsed pathway summaries retain all tested pathways", {
  df <- data.frame(
    pathway = c("A", "B"),
    padj = c(0.01, 0.30),
    mainpathway = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  summary <- plot_utils$make_pathway_count_summary(
    results_by_rank = list(rank_a = df),
    expected_ranknames = "rank_a",
    collapse = FALSE,
    combined = FALSE
  )

  expect_equal(summary$statistics$n_retained_pathways, 2L)
  expect_equal(summary$statistics$ranks$rank_a$n_main_pathways, 0L)
})

test_that("pathway count writer is pretty and honors replace", {
  output_dir <- withr::local_tempdir()
  first <- list(
    statistics = list(n_retained_pathways = 1L, ranks = list()),
    schema_version = 1L,
    generated_at = "2026-07-27T12:00:00Z"
  )
  second <- list(
    statistics = list(n_retained_pathways = 2L, ranks = list()),
    schema_version = 1L,
    generated_at = "2026-07-27T13:00:00Z"
  )

  target <- plot_utils$write_pathway_count_summary(first, output_dir, replace = TRUE)
  original_text <- paste(readLines(target, warn = FALSE), collapse = "\n")
  expect_match(original_text, "\n  \"statistics\"", fixed = TRUE)

  plot_utils$write_pathway_count_summary(second, output_dir, replace = FALSE)
  expect_identical(
    paste(readLines(target, warn = FALSE), collapse = "\n"),
    original_text
  )

  plot_utils$write_pathway_count_summary(second, output_dir, replace = TRUE)
  rewritten <- jsonlite::read_json(target, simplifyVector = TRUE)
  expect_equal(rewritten$statistics$n_retained_pathways, 2L)
  expect_equal(rewritten$generated_at, "2026-07-27T13:00:00Z")
})
