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

  expect_equal(unname(names(name_map)), labels)
  expect_true(all(grepl("x-y$", unname(name_map))))
  expect_false(any(unname(name_map) %in% c("alpha", "beta")))
})

test_that("faceted barplots disable strip clipping and keep strip padding", {
  plt <- plot_tools$barplot_with_numbers(make_faceted_plot_df(), title = "group_alpha_x-y")

  expect_equal(plt$theme$strip.clip, "off")
  expect_false(is.null(plt$theme$strip.text.x$margin))
  expect_true(all(grepl("x-y", plot_tools$format_display_label(unique(plt$data$rankname)))))
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
