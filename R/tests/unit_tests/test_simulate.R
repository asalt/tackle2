suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))



# src_dir <- file.path(here("R"))
# sim_tools <- new.env()
# source(file.path(src_dir, "simulate.R"), local = sim_tools)

source(file.path(here("R", "lazyloader.R")))
sim_tools <- get_tool_env("simulate")

# ==

testthat::test_that("test simulate", {
  data <- sim_tools$simulate_preranked_data()
  testthat::expect_true(
    "data.frame" %in% class(data)
  )

  testthat::expect_true(
    "id" %in% colnames(data)
  )

  testthat::expect_true(
    "value" %in% colnames(data)
  )
})


# unfortunately right now this depends on io.R and fgesa.R working
test_that("test generate testdata", {
  res <- sim_tools$generate_test_data()
  testthat::expect_true("C5_GO:BP" %in% names(res))
  testthat::expect_true(length(res) == 2) # we expect two geneset collections

  testthat::expect_true(
    res[[names(res)[1]]] %>% length() == 2 # we expect two rankobjs in each collection
  )

  testthat::expect_true(
    res[[names(res)[2]]] %>% length() == 2 # we expect two rankobjs in each collection
  )

  testthat::expect_true(
    "data.frame" %in% class(res[[names(res)[1]]][[1]]) # there are alot of brackets here
  )
})

test_that("test test data with collapse", {
  TEST_DATA_withCollapse <- sim_tools$generate_test_data(collapse = TRUE)
  TEST_DATA_withCollapse
  results <- TEST_DATA_withCollapse$`C5_GO:BP` %>%
    map(~ .x %>%
      pull(mainpathway) %>%
      table())

  # Iterate through each result element and perform checks
  walk(results, ~ {
    expect_true(all(c("TRUE", "FALSE") %in% names(.x)), info = "Both TRUE and FALSE should be present in mainpathway.")
    expect_true(.x["TRUE"] > 0, info = "There should be more than 0 TRUE values.")
    expect_true(.x["FALSE"] > 0, info = "There should be more than 0 FALSE values.")
  })
})


test_that("test test data without collapse", {
  TEST_DATA <- sim_tools$generate_test_data()
  results <- TEST_DATA$`C5_GO:BP` %>%
    map(~ .x %>%
      pull(mainpathway) %>%
      table())

  # Iterate through each result element and perform checks
  walk(results, ~ {
    expect_true("TRUE" %in% names(.x), info = "only TRUE should be present in mainpathway count.")
    expect_false("FALSE" %in% names(.x), info = "FALSE should not be present in mainpathway count.")
    expect_true(.x["TRUE"] > 0, info = "There should be more than 0 TRUE values.")
  })
  #
})

test_that("signed Hallmark simulator returns the expected structure", {
  sim <- sim_tools$simulate_signed_hallmark_scenarios(
    seed = 2026,
    n_replicates = 2
  )

  expect_true(is.list(sim))
  expect_true(all(c("rank_dfs", "gene_metadata", "membership_metadata", "scenario_metadata") %in% names(sim)))
  expect_true(length(sim$rank_dfs) == 10)
  expect_true(all(c("id", "value") %in% colnames(sim$rank_dfs[[1]])))
  expect_true(all(c("gene_id", "signal_mean", "value") %in% colnames(sim$gene_metadata)))
  expect_true(all(c("pathway", "weight", "activity") %in% colnames(sim$membership_metadata)))
  expect_true(all(c("scenario_id", "pathway", "activity") %in% colnames(sim$scenario_metadata)))
})

test_that("signed Hallmark simulator normalizes weights and overlap offsets correctly", {
  sim <- sim_tools$simulate_signed_hallmark_scenarios(
    scenarios = list(
      overlap_demo = c(
        HALLMARK_MYC_TARGETS_V1 = 1.2,
        HALLMARK_E2F_TARGETS = 0.8,
        HALLMARK_G2M_CHECKPOINT = 0.7
      )
    ),
    seed = 77,
    n_replicates = 2
  )

  weight_sums <- sim$membership_metadata %>%
    dplyr::group_by(scenario_id, gene_id) %>%
    dplyr::summarise(weight_sum = sum(weight), .groups = "drop")

  expect_true(all(abs(weight_sums$weight_sum - 1) < 1e-8))

  gene_static <- sim$gene_metadata %>%
    dplyr::distinct(
      scenario_id,
      gene_id,
      n_active_memberships,
      overlap_offset
    )

  expect_true(all(gene_static$overlap_offset[gene_static$n_active_memberships <= 1] == 0))
  expect_true(any(abs(gene_static$overlap_offset[gene_static$n_active_memberships > 1]) > 0))
})

test_that("signed Hallmark simulator preserves directionality for positive and negative themes", {
  sim <- sim_tools$simulate_signed_hallmark_scenarios(
    scenarios = list(
      polarity_demo = c(
        HALLMARK_OXIDATIVE_PHOSPHORYLATION = 1.1,
        HALLMARK_GLYCOLYSIS = -1.0
      )
    ),
    seed = 101,
    n_replicates = 1
  )

  positive_genes <- sim$membership_metadata %>%
    dplyr::filter(
      scenario_id == "polarity_demo",
      pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      n_active_memberships == 1
    ) %>%
    dplyr::pull(gene_id) %>%
    unique()

  negative_genes <- sim$membership_metadata %>%
    dplyr::filter(
      scenario_id == "polarity_demo",
      pathway == "HALLMARK_GLYCOLYSIS",
      n_active_memberships == 1
    ) %>%
    dplyr::pull(gene_id) %>%
    unique()

  gene_static <- sim$gene_metadata %>%
    dplyr::filter(rank_name == "polarity_demo_rep01") %>%
    dplyr::select(gene_id, weighted_activity)

  expect_true(length(positive_genes) > 0)
  expect_true(length(negative_genes) > 0)
  expect_true(all(gene_static$weighted_activity[gene_static$gene_id %in% positive_genes] > 0))
  expect_true(all(gene_static$weighted_activity[gene_static$gene_id %in% negative_genes] < 0))
})

test_that("signed Hallmark simulator is reproducible and replicate noise varies as expected", {
  sim1 <- sim_tools$simulate_signed_hallmark_scenarios(
    seed = 5150,
    n_replicates = 2
  )
  sim2 <- sim_tools$simulate_signed_hallmark_scenarios(
    seed = 5150,
    n_replicates = 2
  )

  expect_identical(sim1$rank_dfs, sim2$rank_dfs)
  expect_identical(sim1$membership_metadata, sim2$membership_metadata)
  expect_identical(sim1$gene_metadata, sim2$gene_metadata)

  rep1 <- sim1$rank_dfs[["myc_proliferation_rep01"]]
  rep2 <- sim1$rank_dfs[["myc_proliferation_rep02"]]
  expect_false(identical(rep1$value, rep2$value))
})

test_that("signed Hallmark expression simulator returns a log10 matrix and exportable GCT", {
  withr::with_tempdir({
    gct_path <- file.path(getwd(), "signed_hallmark_expression.gct")
    sim <- sim_tools$simulate_signed_hallmark_expression(
      seed = 9001,
      target_source = "signal_mean",
      n_samples_per_group = 4,
      export_gct_path = gct_path
    )

    expect_true(is.list(sim))
    expect_true(all(c(
      "expression_gct",
      "expression_matrix",
      "sample_metadata",
      "feature_metadata",
      "scenario_gene_metadata",
      "target_rank_dfs",
      "recovered_rank_dfs",
      "limma_tables",
      "recovery_metrics",
      "scenario_metadata",
      "score_simulation"
    ) %in% names(sim)))

    expect_true(file.exists(gct_path))
    parsed <- cmapR::parse_gctx(gct_path)
    expect_equal(dim(parsed@mat), dim(sim$expression_matrix))
    expect_true(all(c("group", "scenario_label", "is_baseline") %in% colnames(parsed@cdesc)))
    expect_true(all(c("gene_symbol", "baseline_log10", "sigma_gene") %in% colnames(parsed@rdesc)))

    expr_quantiles <- stats::quantile(as.numeric(sim$expression_matrix), c(0.05, 0.5, 0.95), na.rm = TRUE)
    expect_gt(expr_quantiles[[1]], 3.3)
    expect_gt(expr_quantiles[[2]], 4.4)
    expect_lt(expr_quantiles[[3]], 7.4)
  })
})

test_that("signed Hallmark expression simulator encodes abundance-dependent variance and limma-like recovery", {
  sim <- sim_tools$simulate_signed_hallmark_expression(
    seed = 4242,
    target_source = "signal_mean",
    n_samples_per_group = 5
  )

  abundance_cor <- stats::cor(
    sim$feature_metadata$baseline_log10,
    sim$feature_metadata$sigma_gene,
    method = "spearman"
  )

  expect_lt(abundance_cor, -0.95)
  expect_true(all(sim$feature_metadata$sigma_gene >= 0.06))
  expect_true(all(sim$feature_metadata$sigma_gene <= 0.4))

  recovery <- sim$recovery_metrics
  expect_true(all(recovery$spearman_t > 0.45))
  expect_true(all(recovery$spearman_signedlogP > 0.35))
  expect_true(all(recovery$pearson_logFC > 0.6))
})

test_that("teaching dataset catalog defines four A/B/C/D datasets with a shared mild vehicle group", {
  catalog <- sim_tools$get_hallmark_teaching_dataset_catalog()

  expect_equal(length(catalog), 4)

  group_names <- purrr::map(catalog, ~ names(.x$group_activities))
  purrr::walk(group_names, ~ expect_equal(.x, c("A", "B", "C", "D")))

  vehicle_profiles <- purrr::map(catalog, ~ .x$group_activities$B)
  expect_equal(vehicle_profiles[[1]], vehicle_profiles[[2]])
  expect_equal(vehicle_profiles[[1]], vehicle_profiles[[3]])
  expect_equal(vehicle_profiles[[1]], vehicle_profiles[[4]])
})

test_that("teaching dataset simulator returns balanced four-group matrices with inspectable metadata", {
  withr::with_tempdir({
    export_dir <- file.path(getwd(), "teaching_gcts")
    sim <- sim_tools$simulate_hallmark_teaching_datasets(
      seed = 3030,
      n_samples_per_group_batch = 2,
      export_dir = export_dir
    )

    expect_true(is.list(sim))
    expect_true(all(c("datasets", "feature_metadata", "group_activity_metadata") %in% names(sim)))
    expect_equal(length(sim$datasets), 4)

    for (dataset_id in names(sim$datasets)) {
      ds <- sim$datasets[[dataset_id]]
      expect_equal(dim(ds$expression_matrix), c(nrow(sim$feature_metadata), 24))
      expect_true(file.exists(ds$export_gct_path))
      expect_true(all(c(
        "sample",
        "group",
        "batch",
        "sample_biology_scale",
        "sample_batch_effect_mean",
        "sample_mean_shift_mean",
        "sample_latent_mean_log10_mean"
      ) %in% colnames(ds$sample_metadata)))
      expect_true(all(c("baseline_log10", "sigma_gene", "batch_module") %in% colnames(ds$feature_metadata)))
      expect_true(all(c("group", "signal_mean", "mean_shift", "group_mean_log10") %in% colnames(ds$gene_group_metadata)))
      expect_true(all(c("pathway", "weight", "n_union_memberships") %in% colnames(ds$membership_metadata)))
      expect_true(all(c("batch_offset", "batch_module_offset_sd") %in% colnames(ds$batch_metadata)))
      expect_true(all(ds$sample_metadata$sample_biology_scale >= 0.75))
      expect_true(all(ds$sample_metadata$sample_biology_scale <= 1.25))

      group_batch_counts <- table(ds$sample_metadata$group, ds$sample_metadata$batch)
      expect_true(all(group_batch_counts == 2))

      expr_quantiles <- stats::quantile(as.numeric(ds$expression_matrix), c(0.05, 0.5, 0.95), na.rm = TRUE)
      expect_gt(expr_quantiles[[1]], 3.4)
      expect_gt(expr_quantiles[[2]], 4.5)
      expect_lt(expr_quantiles[[3]], 6.9)

      activity_meta <- ds$group_activity_metadata %>%
        dplyr::filter(!is.na(pathway), !is_vehicle)

      d_vs_c <- activity_meta %>%
        dplyr::filter(group %in% c("C", "D")) %>%
        dplyr::select(group, pathway, activity) %>%
        tidyr::pivot_wider(names_from = group, values_from = activity) %>%
        dplyr::mutate(
          C = dplyr::coalesce(C, 0),
          D = dplyr::coalesce(D, 0)
        )

      expect_true(all(abs(d_vs_c$D) >= abs(d_vs_c$C)))
    }

    abundance_cor <- stats::cor(
      sim$feature_metadata$baseline_log10,
      sim$feature_metadata$sigma_gene,
      method = "spearman"
    )
    expect_lt(abundance_cor, -0.95)
  })
})

test_that("teaching dataset simulator adds mild sample-level biology jitter while leaving baseline means centered", {
  sim <- sim_tools$simulate_hallmark_teaching_datasets(
    seed = 3031,
    dataset_ids = "dataset2_inflammatory_interferon",
    n_samples_per_group_batch = 2
  )

  ds <- sim$datasets[[1]]
  sample_meta <- ds$sample_metadata

  expect_true(all(sample_meta$sample_biology_scale >= 0.75))
  expect_true(all(sample_meta$sample_biology_scale <= 1.25))

  non_baseline <- sample_meta %>%
    dplyr::filter(group %in% c("B", "C", "D"))
  expect_gt(stats::sd(non_baseline$sample_mean_shift_mean), 0)

  group_c_shift_values <- sample_meta %>%
    dplyr::filter(group == "C") %>%
    dplyr::pull(sample_mean_shift_mean)
  expect_gt(dplyr::n_distinct(round(group_c_shift_values, 8)), 1)

  baseline_shift_values <- sample_meta %>%
    dplyr::filter(group == "A") %>%
    dplyr::pull(sample_mean_shift_mean)
  expect_true(all(abs(baseline_shift_values) < 1e-12))
})

test_that("teaching dataset batch spike creates gene-specific crude batch structure", {
  sim <- sim_tools$simulate_hallmark_teaching_datasets(
    seed = 3032,
    dataset_ids = "dataset1_proliferation_suppression",
    n_samples_per_group_batch = 2,
    batch_effect_sd = 0.08
  )

  ds <- sim$datasets[[1]]
  expect_gt(dplyr::n_distinct(ds$feature_metadata$batch_module), 1)

  meta <- ds$sample_metadata
  baseline_idx <- meta$group == "A"
  batch_ids <- sort(unique(meta$batch))
  batch_gene_means <- vapply(
    batch_ids,
    function(batch_id) {
      rowMeans(ds$expression_matrix[, baseline_idx & meta$batch == batch_id, drop = FALSE])
    },
    numeric(nrow(ds$expression_matrix))
  )

  batch_delta <- batch_gene_means[, 1] - batch_gene_means[, 2]
  expect_gt(stats::sd(batch_delta), 0.02)
  expect_gt(mean(ds$batch_metadata$batch_module_offset_sd), 0.01)
})

test_that("grouped Hallmark simulator supports arbitrary two-group dataset specs with labels", {
  custom_datasets <- list(
    dataset1_epithelial_innate_activation = list(
      label = "Dataset 1: Epithelial Innate Inflammation",
      description = paste(
        "Balanced two-group practice dataset with strong inflammatory and TNF/NFkB activity,",
        "moderate IL6/JAK/STAT3, and mild ROS plus IFN-gamma seasoning."
      ),
      groups = list(
        control = list(
          label = "Sham",
          activities = numeric(),
          treatment_level = 0
        ),
        case = list(
          label = "Innate Activated",
          activities = c(
            HALLMARK_INFLAMMATORY_RESPONSE = 1.2,
            HALLMARK_TNFA_SIGNALING_VIA_NFKB = 1.0,
            HALLMARK_IL6_JAK_STAT3_SIGNALING = 0.6,
            HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = 0.3,
            HALLMARK_INTERFERON_GAMMA_RESPONSE = 0.25
          ),
          treatment_level = 1
        )
      ),
      expected_up_pathways = c(
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE"
      )
    )
  )

  sim <- sim_tools$simulate_hallmark_grouped_datasets(
    datasets = custom_datasets,
    seed = 3033,
    n_samples_per_group_batch = 3,
    batch_effect_sd = 0.12
  )

  expect_true(all(c("datasets", "feature_metadata", "group_metadata", "group_activity_metadata") %in% names(sim)))
  expect_equal(length(sim$datasets), 1)

  ds <- sim$datasets[[1]]
  expect_equal(dim(ds$expression_matrix), c(nrow(sim$feature_metadata), 18))
  expect_equal(ds$group_metadata$group, c("control", "case"))
  expect_equal(ds$group_metadata$group_label, c("Sham", "Innate Activated"))
  expect_true(all(ds$sample_metadata$is_vehicle == FALSE))
  expect_equal(unique(ds$sample_metadata$treatment_level[ds$sample_metadata$group == "control"]), 0)
  expect_equal(unique(ds$sample_metadata$treatment_level[ds$sample_metadata$group == "case"]), 1)

  group_batch_counts <- table(ds$sample_metadata$group, ds$sample_metadata$batch)
  expect_true(all(group_batch_counts == 3))

  case_pathways <- ds$group_activity_metadata %>%
    dplyr::filter(group == "case", !is.na(pathway)) %>%
    dplyr::pull(pathway)
  expect_setequal(
    case_pathways,
    c(
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
      "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE"
    )
  )

  expect_gt(mean(ds$batch_metadata$batch_module_offset_sd), 0.01)
})

test_that("grouped Hallmark simulator applies named gene spikes alongside pathway programs", {
  custom_datasets <- list(
    dataset_gene_spike_demo = list(
      label = "Dataset Gene Spike Demo",
      description = "Two-group custom dataset with pathway and gene-anchor signal.",
      groups = list(
        control = list(
          label = "Control",
          activities = numeric(),
          treatment_level = 0
        ),
        case = list(
          label = "Glycolytic Shift",
          activities = c(
            HALLMARK_GLYCOLYSIS = 1.2,
            HALLMARK_HYPOXIA = 0.6
          ),
          gene_spikes = c(
            HK2 = 1.4,
            LDHA = 1.2
          ),
          treatment_level = 1
        )
      )
    )
  )

  sim <- sim_tools$simulate_hallmark_grouped_datasets(
    datasets = custom_datasets,
    seed = 3034,
    n_samples_per_group_batch = 3,
    batch_effect_sd = 0.12
  )

  expect_true("group_gene_spike_metadata" %in% names(sim))
  expect_equal(nrow(sim$group_gene_spike_metadata), 2)

  ds <- sim$datasets[[1]]
  expect_true("group_gene_spike_metadata" %in% names(ds))
  expect_true(all(c("gene_spike_value", "gene_spike_mean_shift", "pathway_mean_shift") %in% colnames(ds$gene_group_metadata)))

  hk2_case <- ds$gene_group_metadata %>%
    dplyr::filter(group == "case", gene_symbol == "HK2") %>%
    dplyr::slice(1)
  hk2_control <- ds$gene_group_metadata %>%
    dplyr::filter(group == "control", gene_symbol == "HK2") %>%
    dplyr::slice(1)

  expect_equal(hk2_case$gene_spike_target_name, "HK2")
  expect_gt(hk2_case$gene_spike_value, 1)
  expect_gt(hk2_case$gene_spike_mean_shift, 0)
  expect_equal(hk2_control$gene_spike_mean_shift, 0)
  expect_gt(hk2_case$mean_shift, hk2_control$mean_shift)
})


test_that("grouped Hallmark simulator supports null datasets with no pathway or gene-spike signal", {
  custom_datasets <- list(
    dataset_null_demo = list(
      label = "Dataset Null Demo",
      description = "Two-group null dataset with no intentional biological contrast.",
      groups = list(
        control = list(
          label = "Control",
          activities = numeric(),
          treatment_level = 0
        ),
        case = list(
          label = "Null Contrast",
          activities = numeric(),
          treatment_level = 0
        )
      )
    )
  )

  sim <- sim_tools$simulate_hallmark_grouped_datasets(
    datasets = custom_datasets,
    seed = 3035,
    n_samples_per_group_batch = 3,
    batch_effect_sd = 0.12
  )

  ds <- sim$datasets[[1]]
  expect_true(all(is.na(ds$group_activity_metadata$pathway)))
  expect_true(all(ds$group_activity_metadata$activity == 0))
  expect_equal(nrow(ds$group_gene_spike_metadata), 0)
  expect_true(all(abs(ds$gene_group_metadata$mean_shift) < 1e-12))
  expect_true(all(abs(ds$sample_metadata$sample_mean_shift_mean) < 1e-12))
})
