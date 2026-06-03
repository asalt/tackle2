# not fully fleshed out
# options(error = function() {
#             # print(dump.frames())
#             traceback(1)
#             cat("\n")
# })
options(error = NULL) # or recover

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(withr))

io_tools <- new.env()
source(file.path(here("R"), "io.R"), local = io_tools)

geneset_tools <- new.env()
source(file.path(here("R"), "./geneset_utils.R"), local = geneset_tools)

fgsea_tools <- new.env()
source(file.path(here("R"), "./fgsea.R"), local = fgsea_tools)

sim_tools <- new.env()
source(file.path(here("R"), "./simulate.R"), local = sim_tools)

run_env <- new.env()
source(file.path(here("R"), "./run.R"), local = run_env)

source(file.path(here("R"), "tests", "unit_tests", "helpers", "fgsea_test_helpers.R"))

# ==================================

trycatch <- function(expr, silent = TRUE) {
  tryCatch(expr, error = function(e) {
    if (!silent) {
      message("An error occurred: ", e$message)
    }
  })
}

testthat::test_that("enplot combine specs preserve facet labels while ordering by rankname", {
  metadata <- data.frame(
    id = c("rank_b", "rank_a", "rank_c"),
    group = c("G2", "G1", "G2"),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$id

  process_cut_by_fn <- function(cut_by, cdesc) {
    interaction(cdesc[, cut_by, drop = FALSE], drop = TRUE)
  }

  specs <- run_env$make_enplot_combine_specs(
    combine_by = "group",
    metadata = metadata,
    rankname_order = c("rank_a", "rank_b", "rank_c"),
    process_cut_by_fn = process_cut_by_fn
  )

  testthat::expect_equal(length(specs), 1)
  testthat::expect_equal(specs[[1]]$name, "group")
  testthat::expect_equal(specs[[1]]$combine_by$rankname, c("rank_a", "rank_b", "rank_c"))
  testthat::expect_equal(as.character(specs[[1]]$combine_by$facet), c("G1", "G2", "G2"))
  testthat::expect_false(any(is.na(specs[[1]]$combine_by$facet)))
  testthat::expect_equal(levels(specs[[1]]$combine_by$facet), c("G1", "G2"))
})

testthat::test_that("enplot combine specs return an all-route when combine_by is empty or invalid", {
  metadata <- data.frame(id = c("rank_a", "rank_b"), group = c("G1", "G2"))
  rownames(metadata) <- metadata$id
  null_process_cut_by <- function(cut_by, cdesc) NULL

  empty_spec <- run_env$make_enplot_combine_specs(
    combine_by = FALSE,
    metadata = metadata,
    rankname_order = rownames(metadata),
    process_cut_by_fn = null_process_cut_by
  )
  invalid_spec <- run_env$make_enplot_combine_specs(
    combine_by = "missing_column",
    metadata = metadata,
    rankname_order = rownames(metadata),
    process_cut_by_fn = null_process_cut_by
  )

  testthat::expect_equal(empty_spec[[1]]$name, "all")
  testthat::expect_null(empty_spec[[1]]$combine_by)
  testthat::expect_equal(invalid_spec[[1]]$name, "all")
  testthat::expect_null(invalid_spec[[1]]$combine_by)
})



setup <- function(data_dir, output_dir, gct_path) {
  if (!fs::dir_exists(data_dir)) fs::dir_create(data_dir, recurse = TRUE)
  if (!fs::dir_exists(output_dir)) fs::dir_create(output_dir, recurse = TRUE)

  datas1 <- purrr::map(1:3, ~ sim_tools$simulate_preranked_data(seed = 4321))
  datas1 %<>% purrr::map(~ .x %>% mutate(value = value + rnorm(nrow(.), sd = .1)))
  datas2 <- purrr::map(1:3, ~ sim_tools$simulate_preranked_data(seed = 1234))
  datas2 %<>% purrr::map(~ .x %>% mutate(value = value + rnorm(nrow(.), sd = .1)))
  datas3 <- purrr::map(1:3, ~ sim_tools$simulate_preranked_data(seed = 9999))
  datas3 %<>% purrr::map(~ .x %>% mutate(value = value + rnorm(nrow(.), sd = .1)))
  datas <-
    c(
      datas1,
      datas2,
      datas3
    )

  names(datas) <- c(
    paste0("group_A_vs_B_dirB_", seq_along(datas1)),
    paste0("group_C_vs_D_dirB_", seq_along(datas2)),
    paste0("group_E_vs_F_dirB_", seq_along(datas3))
  )

  ._ <- datas %>% purrr::imap(~ {
    # if (!file.exists(
    write_tsv(.x, file.path(data_dir, paste0(.y, ".tsv")))
  })

  .mat <- base::Reduce(
    accumulate = F,
    f = function(...) full_join(..., by = "id"),
    x = datas
  ) %>% as.data.frame()
  colnames(.mat) <- c("id", names(datas)) # lost the names

  rownames(.mat) <- .mat$id
  .mat$id <- NULL
  .meta <- data.frame(
    id = colnames(.mat),
    group = c(rep("A", 3), rep("B", 3), rep("C", 3))
  )
  rownames(.meta) <- .meta$id
  .rdesc <- data.frame(id = rownames(.mat), dummy = "X")
  rownames(.rdesc) <- rownames(.mat)
  gct <- new("GCT",
    mat = .mat %>% as.matrix(),
    cdesc = .meta,
    rdesc = .rdesc
  )
  # test_gct_path <- file.path(data_dir, "test.gct")
  cmapR::write_gct(gct, gct_path, appenddim = F)
}

testthat::test_that("test main function with valid parameters", {
  testthat::with_mocked_bindings({
    keep_artifacts <- identical(tolower(Sys.getenv("TACKLE2_KEEP_RUN_TEST_ARTIFACTS", "false")), "true")
    root_dir <- fs::path(tempdir(), sprintf("tackle2-run-%s-%04d",
      format(Sys.time(), "%Y%m%d%H%M%S"),
      sample.int(9999, 1)
    ))
    fs::dir_create(root_dir, recurse = TRUE)

    if (!keep_artifacts) {
      defer(fs::dir_delete(root_dir), envir = parent.frame())
    } else {
      message("Keeping tackle2 run-loop test artifacts in ", root_dir)
    }

    data_dir <- fs::path(root_dir, "data")
    output_dir <- fs::path(root_dir, "output")
    fs::dir_create(data_dir, recurse = TRUE)
    fs::dir_create(output_dir, recurse = TRUE)
    test_gct_path <- fs::path(data_dir, "test.gct")

    setup(data_dir, output_dir, test_gct_path)

    params <- list(
      rankfiledir = data_dir,
      volcanodir = data_dir,
      savedir = output_dir,
      gct_path = test_gct_path,
      ranks_from = "volcano",
      advanced = list(
        cache = FALSE,
        quiet = TRUE,
        parallel = FALSE,
        replace = TRUE,
        cachedir = fs::path(output_dir, "cache"),
        logfile = fs::path(output_dir, "run.log")
      ),
      barplot = list(do_individual = FALSE, do_combined = FALSE),
      heatmap_gsea = list(do = FALSE),
      heatmap_gene = list(do = FALSE),
      enplot = list(do_individual = FALSE, do_combined = FALSE),
      pca = list(do = FALSE),
      genesets = list(
        list(category = "H", subcategory = "", collapse = FALSE)
      )
    )

    testthat::expect_no_error({
      print("Function is being called")
      run_env$run(params)
    })
  }, fgsea = fgsea_test_fake_fgsea, .package = "fgsea")
})

testthat::test_that("sample_exclude is applied before table export", {
  testthat::with_mocked_bindings({
    root_dir <- fs::path(tempdir(), sprintf("tackle2-run-exclude-%s-%04d",
      format(Sys.time(), "%Y%m%d%H%M%S"),
      sample.int(9999, 1)
    ))
    fs::dir_create(root_dir, recurse = TRUE)
    defer(fs::dir_delete(root_dir), envir = parent.frame())

    data_dir <- fs::path(root_dir, "data")
    output_dir <- fs::path(root_dir, "output")
    fs::dir_create(data_dir, recurse = TRUE)
    fs::dir_create(output_dir, recurse = TRUE)
    test_gct_path <- fs::path(data_dir, "test.gct")
    setup(data_dir, output_dir, test_gct_path)

    params <- list(
      rankfiledir = data_dir,
      volcanodir = data_dir,
      savedir = output_dir,
      gct_path = test_gct_path,
      ranks_from = "volcano",
      sample_exclude = c("group_A_vs_B_dirB_1"),
      advanced = list(
        cache = FALSE,
        quiet = TRUE,
        parallel = FALSE,
        replace = TRUE,
        cachedir = fs::path(output_dir, "cache"),
        logfile = fs::path(output_dir, "run.log")
      ),
      barplot = list(do_individual = FALSE, do_combined = FALSE),
      heatmap_gsea = list(do = FALSE),
      heatmap_gene = list(do = FALSE),
      enplot = list(do_individual = FALSE, do_combined = FALSE),
      pca = list(do = FALSE),
      genesets = list(
        list(category = "H", subcategory = "", collapse = FALSE)
      )
    )

    testthat::expect_no_error(run_env$run(params))
    exported <- fs::dir_ls(fs::path(output_dir, "gsea_tables"), glob = "*.tsv")
    testthat::expect_equal(length(exported), 8)
  }, fgsea = fgsea_test_fake_fgsea, .package = "fgsea")
})


# testthat::test_that("test main function with valid parameters", {
#   params <- list(
#     rankfiledir = data_dir,
#     volcanodir = data_dir,
#     savedir = output_dir,
#     gct_path = test_gct_path,
#     ranks_from = "volcano",
#     genesets = list( list(category="H", subcategory="", collapse=FALSE) )
#   )
#   testthat::expect_no_error({
#     print("Function is being called")
#     run_env$run(params)
#   })
# })

# test_invalid_dir()
# test_one()

# teardown()

# test_invalid_dir <- function() {
#   rmarkdown::render("../../run.Rmd",
#     output_format = "html_document",
#     output_dir = output_dir,
#     output_file = "test_invalid_dir.html",
#     params = list(
#       rankfiledir = "does-not-exist",
#       savedir = output_dir,
#       genesets_json = '[{"category": "H", "subcategory": ""}]'
#     )
#   )
# }

# test_render <- function() {
#   rmarkdown::render("../../run.Rmd",
#     # output_format = "html_document",
#     output_dir = output_dir,
#     output_file = "test_defaults.html",
#     params = list(
#       rankfiledir = data_dir,
#       volcanodir = data_dir,
#       gct_path = NULL,
#       savedir = output_dir,
#       ranks_from = "volcano",
#       genesets_json = '[{"category": "H", "subcategory": ""}]'
#     )
#   )
# }


# test_one <- function() {
# geneset <- msigdbr::msigdbr(
#   species = "Homo sapiens",
#   category = "H",
#   subcategory = ""
# )


# genes <- geneset %>%
#   pull(entrez_gene) %>%
#   unique()

# for (f in c("groupA_vs_B_dirB.tsv", "group_A_vs_C_dirB.tsv")) {
#   .values <- rnorm(n = length(genes))
#   .data <- data.frame(
#     geneid = genes,
#     signedlogp = .values
#   )
#   .out <- file.path(data_dir, f)
#   write_tsv(.data, .out)
# }

#   rmarkdown::render("../../run.Rmd",
#     output_format = "html_document",
#     output_dir = output_dir,
#     output_file = "test_full.html",
#     params = list(
#       rankfiledir = data_dir,
#       savedir = file.path(output_dir, "test_full"),
#       genesets_json = '[{"category": "H", "subcategory": ""}]'
#     )
#   )
# }



# testthat::test_that(
#   "test base render",
#   {
#     testthat::expect_no_error(test_render())
#     # assertthat::file_ex
#     assertthat::are_equal(
#       fs::file_exists(
#         file.path(output_dir, "test_defaults.html")
#       ),
#       TRUE
#     )
#   }
# )

# test_render()
