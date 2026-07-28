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
    "HALLMARK_VERY_LONG_PATHWAY_NAME_WITH_MULTIPLE_TOKENS_",
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

bar_text_layer_size <- function(plot) {
  text_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomText"),
    plot$layers
  )
  stopifnot(length(text_layers) == 1)
  text_layers[[1]]$aes_params$size
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

test_that("adaptive facet labels use the tallest wrapped label as their size budget", {
  labels <- c(
    "one line",
    "two\nlines",
    "three\nline\nlabel",
    "four\nline\nlabel\nhere"
  )

  expect_equal(
    plot_tools$adaptive_facet_label_sizes(labels, base_size = 7),
    c(12, 12, 28 / 3, 7)
  )
})

test_that("adaptive facet label maps preserve text safely as rich labels", {
  values <- c(
    "short_&_direct",
    "a_much_longer_facet_label_that_wraps_across_multiple_lines"
  )
  formatted <- plot_tools$format_display_label(values, wrap_width = 18)
  expected_sizes <- plot_tools$adaptive_facet_label_sizes(formatted, base_size = 7)
  label_map <- plot_tools$make_adaptive_facet_label_map(
    values,
    wrap_width = 18,
    base_size = 7
  )

  expect_match(unname(label_map[[1]]), "&amp;", fixed = TRUE)
  expect_true(any(grepl("<br>", unname(label_map), fixed = TRUE)))
  expect_match(
    unname(label_map[[2]]),
    sprintf("font-size:%.2fpt", expected_sizes[[2]]),
    fixed = TRUE
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

test_that("rank display names extract contrasts unless labels are explicit", {
  ranknames <- c(
    "RUN_A4B8_20260714_nz2_dir_B_fna_T_group_groupIR_left_minus_groupnone_pre_imv_T_med",
    "RUN_A4B8_20260714_nz2_dir_B_fna_T_group_groupIR_right_minus_groupnone_pre_imv_T_med"
  )

  default_map <- plot_tools$make_rank_display_name_map(ranknames)
  expect_equal(names(default_map), ranknames)
  expect_equal(
    unname(default_map),
    c("groupIR_left_minus_groupnone_pre", "groupIR_right_minus_groupnone_pre")
  )

  rank_metadata <- data.frame(
    rankname = ranknames,
    rank_label = c("IR left vs pre", "IR right vs pre"),
    rank_order = seq_along(ranknames),
    rank_label_source = "names.txt",
    stringsAsFactors = FALSE
  )
  explicit_map <- plot_tools$make_rank_display_name_map(ranknames, rank_metadata)
  expect_equal(unname(explicit_map), rank_metadata$rank_label)
})

test_that("rankorder display frames use automatic contrast labels", {
  ranknames <- c(
    "RUN_A4B8_group_groupIR_left_minus_groupnone_pre_imv_T_med",
    "RUN_A4B8_group_groupIR_right_minus_groupnone_pre_imv_T_med"
  )
  df <- data.frame(rankname = ranknames, value = c(1, 2), stringsAsFactors = FALSE)

  labeled <- plot_tools$apply_rank_labels_to_df(df)

  expect_equal(
    labeled$rankname,
    c("groupIR_left_minus_groupnone_pre", "groupIR_right_minus_groupnone_pre")
  )
})

test_that("faceted barplots disable strip clipping and keep strip padding", {
  input <- make_faceted_plot_df()
  plt <- plot_tools$barplot_with_numbers(input, title = "group_alpha_x-y")

  expect_equal(plt$theme$strip.clip, "off")
  expect_false(is.null(plt$theme$strip.text.x$margin))
  expect_setequal(as.character(plt$data$rankname), c("alpha_x-y", "beta_x-y"))
  expect_true(all(grepl("x-y", plot_tools$format_display_label(unique(plt$data$rankname)))))
  if (requireNamespace("ggtext", quietly = TRUE)) {
    expect_true(inherits(plt$theme$strip.text.x, "element_markdown"))
  }
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
  expect_lte(plot_tools$max_wrapped_label_lines(unique(plt$data$pathway)), 2)
  expect_lt(plt$theme$axis.text.y$size, 5)
  expect_lt(bar_text_layer_size(plt), 2.2)
  expect_lt(saved_dims[["height"]], 16)
  expect_gt(saved_dims[["height"]], 7)
})

test_that("faceted barplot typography stays unchanged when rows have room", {
  plt <- plot_tools$barplot_with_numbers(
    make_many_faceted_bar_df(n_pathways = 12, n_facets = 4),
    title = "short faceted pathways"
  )

  expect_lte(plot_tools$max_wrapped_label_lines(unique(plt$data$pathway)), 2)
  expect_equal(plt$theme$axis.text.y$size, 6.2)
  expect_equal(bar_text_layer_size(plt), 2.2)
})

test_that("single barplot typography is not compacted", {
  df <- make_many_faceted_bar_df(n_pathways = 12, n_facets = 1)
  plt <- plot_tools$barplot_with_numbers(df, title = "single pathways")

  expect_equal(plt$theme$axis.text.y$size, 6.6)
  expect_equal(bar_text_layer_size(plt), 2.2)
  expect_false(inherits(plt$theme$strip.text.x, "element_markdown"))
})

test_that("faceted bubble plots disable strip clipping and keep strip padding", {
  input <- make_faceted_plot_df()
  plt <- bubble_tools$bubble_plot(
    input,
    title = "group_alpha_x-y",
    subtitle = "rank: group_alpha_x-y"
  )

  expect_equal(plt$theme$strip.clip, "off")
  expect_false(is.null(plt$theme$strip.text$margin))
  expect_setequal(as.character(plt$data$rankname), c("alpha_x-y", "beta_x-y"))
  expect_true(all(grepl("x-y", plot_tools$format_display_label(unique(plt$data$rankname)))))
  if (requireNamespace("ggtext", quietly = TRUE)) {
    expect_true(inherits(plt$theme$strip.text, "element_markdown"))
  }
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
  expect_lte(plot_tools$max_wrapped_label_lines(unique(plt$data$pathway)), 2)
  expect_lt(plt$theme$axis.text.y$size, 4.5)
  expect_lt(saved_dims[["height"]], 10)
  expect_gt(saved_dims[["height"]], 7)
})

test_that("single bubble typography is not compacted", {
  df <- make_many_faceted_bar_df(n_pathways = 12, n_facets = 1)
  plt <- bubble_tools$bubble_plot(df, title = "single pathways")

  expect_equal(plt$theme$axis.text.y$size, 6.6)
  expect_false(inherits(plt$theme$strip.text, "element_markdown"))
})

test_that("short bubble plots use a legend-safe total height", {
  saved_dims <- list()
  fake_save <- function(plot_code, width = NULL, height = NULL, ...) {
    saved_dims[[length(saved_dims) + 1]] <<- c(width = width, height = height)
    TRUE
  }

  bubble_tools$bubble_plot(
    make_many_faceted_bar_df(n_pathways = 12, n_facets = 1),
    save_func = fake_save,
    title = "short single pathways"
  )
  bubble_tools$bubble_plot(
    make_many_faceted_bar_df(n_pathways = 12, n_facets = 4),
    save_func = fake_save,
    title = "short faceted pathways"
  )

  expect_equal(saved_dims[[1]][["height"]], 5)
  expect_equal(saved_dims[[2]][["height"]], 5)
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
