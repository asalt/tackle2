suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))

plot_tools <- new.env()
source(file.path(src_dir, "./plot.R"), local = plot_tools)

bubble_tools <- new.env()
source(file.path(src_dir, "./plot_bubble.R"), local = bubble_tools)

make_faceted_plot_df <- function() {
  df <- data.frame(
    pathway = c(
      "HALLMARK_ALPHA_RESPONSE",
      "HALLMARK_BETA_SIGNALING",
      "HALLMARK_ALPHA_RESPONSE",
      "HALLMARK_BETA_SIGNALING"
    ),
    NES = c(1.4, -1.1, 1.6, -1.3),
    padj = c(0.01, 0.03, 0.02, 0.04),
    size = c(12, 10, 11, 9),
    rankname = c(
      "group_alpha_x-y",
      "group_alpha_x-y",
      "group_beta_x-y",
      "group_beta_x-y"
    ),
    stringsAsFactors = FALSE
  )
  df$leadingEdge <- I(list(
    c("A", "B", "C"),
    c("D", "E"),
    c("A", "F"),
    c("E", "G")
  ))
  df
}

make_many_faceted_bar_df <- function(n_pathways = 40, n_facets = 4) {
  pathways <- paste0(
    "HALLMARK_LONG_PATHWAY_NAME_WITH_MULTIPLE_TOKENS_",
    seq_len(n_pathways)
  )
  ranknames <- paste0(
    "group_alpha_long_context_",
    seq_len(n_facets),
    "_minus_group_beta_long_context_",
    seq_len(n_facets)
  )
  grid <- expand.grid(
    pathway = pathways,
    rankname = ranknames,
    stringsAsFactors = FALSE
  )
  grid$NES <- rep(seq(-2, 2, length.out = n_pathways), times = n_facets)
  grid$padj <- 0.01
  grid$size <- 10
  grid$leadingEdge <- I(rep(list(c("A", "B", "C")), nrow(grid)))
  grid
}

test_that("display labels preserve hyphenated contrast tokens", {
  expect_equal(
    plot_tools$format_display_label("group_alpha_x-y"),
    "group alpha x-y"
  )
  expect_match(
    plot_tools$format_display_label("group_alpha_x-y", wrap_width = 10),
    "x-y"
  )
})

test_that("display name map does not strip hyphenated suffixes", {
  labels <- c("group_alpha_x-y", "group_beta_x-y")
  name_map <- plot_tools$make_display_name_map(labels)
  name_values <- unname(as.character(name_map))

  expect_equal(unname(names(name_map)), labels)
  expect_true(all(grepl("x-y$", name_values)))
  expect_false(any(name_values %in% c("alpha", "beta")))
})

test_that("display name map preserves common cell-line suffixes", {
  labels <- c(
    "groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
    "groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
  )

  name_map <- plot_tools$make_display_name_map(labels)
  name_values <- unname(as.character(name_map))
  display_labels <- plot_tools$format_display_label(name_values)

  expect_equal(name_values, labels)
  expect_true(all(grepl("Control_X9Q2$", name_values)))
  expect_true(all(grepl("Control X9Q2", display_labels)))
  expect_false(any(grepl("minus group(ABC|DEF)$", display_labels)))
})

test_that("display name map preserves common control denominators", {
  labels <- c(
    "groupABC_minus_groupControl_X9Q2",
    "groupDEF_minus_groupControl_X9Q2"
  )

  name_map <- plot_tools$make_display_name_map(labels)
  name_values <- unname(as.character(name_map))

  expect_equal(name_values, labels)
  expect_true(all(grepl("groupControl_X9Q2$", name_values)))
})

test_that("display name map still shortens shared leading context", {
  labels <- c(
    "RUN_C3D9_Definitely_Not_Final_groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
    "RUN_C3D9_Definitely_Not_Final_groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
  )

  name_map <- plot_tools$make_display_name_map(labels)
  name_values <- unname(as.character(name_map))

  expect_equal(
    name_values,
    c(
      "groupABC_A4B8-Sunk_My_Battleship_minus_groupABC_Control_X9Q2",
      "groupDEF_A4B8-Sunk_My_Battleship_minus_groupDEF_Control_X9Q2"
    )
  )
  expect_true(all(grepl("Control_X9Q2$", name_values)))
})

test_that("faceted barplots disable strip clipping and keep strip padding", {
  plt <- plot_tools$barplot_with_numbers(make_faceted_plot_df(), title = "group_alpha_x-y")

  expect_equal(plt$theme$strip.clip, "off")
  expect_false(is.null(plt$theme$strip.text.x$margin))
  expect_true(all(grepl("x-y", plot_tools$format_display_label(unique(plt$data$rankname)))))
})

test_that("faceted barplot sizing grows more slowly for many pathways", {
  saved_dims <- NULL
  fake_save <- function(plot_code, width = NULL, height = NULL, ...) {
    saved_dims <<- c(width = width, height = height)
    TRUE
  }

  plt <- plot_tools$barplot_with_numbers(
    make_many_faceted_bar_df(n_pathways = 40, n_facets = 4),
    save_func = fake_save,
    title = "many faceted pathways"
  )

  expect_true(any(grepl("\n", as.character(plt$data$pathway), fixed = TRUE)))
  expect_lt(saved_dims[["height"]], 16)
  expect_gt(saved_dims[["height"]], 7)
})

test_that("faceted bubble plots disable strip clipping and keep strip padding", {
  plt <- bubble_tools$bubble_plot(
    make_faceted_plot_df(),
    title = "group_alpha_x-y",
    subtitle = "rank: group_alpha_x-y"
  )

  expect_equal(plt$theme$strip.clip, "off")
  expect_false(is.null(plt$theme$strip.text$margin))
  expect_true(all(grepl("x-y", plot_tools$format_display_label(unique(plt$data$rankname)))))
})

test_that("faceted bubble plot sizing grows more slowly for many pathways", {
  saved_dims <- NULL
  fake_save <- function(plot_code, width = NULL, height = NULL, ...) {
    saved_dims <<- c(width = width, height = height)
    TRUE
  }

  plt <- bubble_tools$bubble_plot(
    make_many_faceted_bar_df(n_pathways = 40, n_facets = 4),
    save_func = fake_save,
    title = "many faceted pathways"
  )

  expect_true(any(grepl("\n", as.character(plt$data$pathway), fixed = TRUE)))
  expect_lt(saved_dims[["height"]], 10)
  expect_gt(saved_dims[["height"]], 7)
})

test_that("bubble height scale compacts saved plot height", {
  saved_dims <- list()
  fake_save <- function(plot_code, width = NULL, height = NULL, ...) {
    saved_dims[[length(saved_dims) + 1]] <<- c(width = width, height = height)
    TRUE
  }

  bubble_tools$bubble_plot(
    make_many_faceted_bar_df(n_pathways = 20, n_facets = 4),
    save_func = fake_save,
    title = "unscaled bubble pathways",
    height_scale = 1
  )
  bubble_tools$bubble_plot(
    make_many_faceted_bar_df(n_pathways = 20, n_facets = 4),
    save_func = fake_save,
    title = "scaled bubble pathways",
    height_scale = 0.75
  )

  expect_equal(saved_dims[[2]][["height"]], saved_dims[[1]][["height"]] * 0.75)
})
