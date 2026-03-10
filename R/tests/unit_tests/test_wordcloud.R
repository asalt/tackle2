suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(withr))

src_dir <- file.path(here("R"))

fgsea_tools <- new.env()
source(file.path(src_dir, "fgsea.R"), local = fgsea_tools)

sim_tools <- new.env()
source(file.path(src_dir, "simulate.R"), local = sim_tools)
generate_test_data <- sim_tools$generate_test_data

plot_utils <- new.env()
source(file.path(src_dir, "plot_utils.R"), local = plot_utils)

wordcloud_tools <- new.env()
source(file.path(src_dir, "plot_wordcloud.R"), local = wordcloud_tools)


context("Wordcloud summary plots")

TEST_DATA <- generate_test_data(collapse = FALSE, pathways = c("H"))


test_that("prepare_terms_for_wordcloud produces non-empty term table", {
  df <- TEST_DATA$H_ %>% fgsea_tools$concat_results_one_collection()
  terms <- wordcloud_tools$prepare_terms_for_wordcloud(
    df,
    padj_cutoff = 0.5,
    top_n_pathways = 40
  )

  expect_true(
    nrow(terms) > 0,
    info = "no terms returned from prepare_terms_for_wordcloud"
  )
  expect_true(
    all(c("word", "weight", "direction") %in% colnames(terms)),
    info = "term table is missing required columns"
  )
})


test_that("wordcloud_plot returns ggplot and saves example file", {
  df <- TEST_DATA$H_ %>% fgsea_tools$concat_results_one_collection()
  expect_true(
    nrow(df) > 0,
    info = "no rows available for wordcloud test"
  )

  withr::with_tempdir({
    output_dir <- file.path(getwd(), "wordcloud_examples")

    save_example <- function(plot_code, width, height, ...) {
      plot_utils$plot_and_save(
        plot_code = plot_code,
        filename = "wordcloud_H_example",
        path = output_dir,
        type = "pdf",
        width = width,
        height = height,
        replace = TRUE
      )
    }

    plt <- wordcloud_tools$wordcloud_plot(
      df,
      title = "H_wordcloud_example",
      padj_cutoff = 0.5,
      top_n_pathways = 40,
      max_words = 50,
      save_func = save_example
    )

    expect_true(
      "ggplot" %in% class(plt),
      info = "wordcloud_plot did not return a ggplot object"
    )

    expect_true(
      file.exists(file.path(output_dir, "wordcloud_H_example.pdf")),
      info = "expected wordcloud example pdf was not created"
    )
  })
}
)
