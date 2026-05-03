suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))

src_dir <- file.path(here("R"))

# io_tools <- new.env()
# source(file.path(file.path(src_dir, "./io.R")), local = io_tools)

geneset_tools <- new.env()
source(file.path(src_dir, "./geneset_utils.R"), local = geneset_tools)

io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)


source(file.path(here("R", "lazyloader.R")))
io_tools <- get_tool_env("io")
geneset_tools <- get_tool_env("geneset_utils")

# ================================

make_random_gct <- function(nrow = 10, ncol = 4) {
  set.seed(369)
  nrow <- max(nrow, 1)
  ncol <- max(ncol, 1)
  .mat <- matrix(runif(nrow * ncol), nrow = nrow, ncol = ncol)
  .rids <- seq(1, dim(.mat)[1]) %>% as.character()
  .cids <- seq(1, dim(.mat)[2]) %>% as.character()
  .cids <- paste0("X", .cids)
  .cdesc <- data.frame(
    metavar1 = sample(letters[1:5], ncol, replace = T),
    metavar2 = sample(letters[1:5], ncol, replace = T)
  )
  .rdesc <- data.frame(
    rdesc = paste0("gene", seq(1, nrow))
  )
  gct <- cmapR::GCT(mat = .mat, rid = .rids, cid = .cids, cdesc = .cdesc, rdesc = .rdesc)
  gct
}


.cache <- list()
simulate_preranked_data <- function(...) {
  ll <- list(...)
  hashval <- rlang::hash(ll)
  if (!hashval %in% names(.cache)) .cache[[hashval]] <- do.call(.simulate_preranked_data, ll)
  return(.cache[[hashval]])
}



.simulate_preranked_data <- function(
    seed = 4321,
    geneset = NULL,
    spike_terms = c("INTERFERON"),
    sample_frac = 1.0,
    ...) {
  set.seed(seed)

  # geneset <- msigdbr::msigdbr(
  #   species = "Homo sapiens",
  #   category = "H",
  #   subcategory = ""
  # )

  if (is.null(geneset)) {
    geneset <- geneset_tools$get_collection("H", "")
  }


  # Generate a list of gene sets for each spike term
  spike_genes_list <- purrr::map(spike_terms, ~ geneset %>%
    dplyr::filter(str_detect(gs_name, .x)) %>%
    dplyr::pull(ncbi_gene) %>%
    unique())

  genes <- geneset %>%
    dplyr::pull(ncbi_gene) %>%
    unique()

  spike_genes <- unique(unlist(spike_genes_list))
  background_genes <- setdiff(genes, spike_genes)


  bg_values <- rnorm(n = length(background_genes))
  bg_data <- data.frame(
    id = background_genes,
    value = bg_values
  )

  # Spike gene values, assigning different means for each spike term
  spike_data <- purrr::map2_df(spike_genes_list, seq_along(spike_genes_list), ~ data.frame(
    id = .x,
    value = rnorm(n = length(.x), mean = .y) # Incrementing mean for differentiation
  ))


  data <- bind_rows(
    bg_data,
    spike_data
  )
  data %<>% distinct(id, .keep_all = TRUE)
  data %<>% dplyr::mutate(id = as.character(id))
  data %<>% dplyr::sample_frac(size = sample_frac)

  return(data)
}

.or_default <- function(value, default) {
  if (is.null(value)) {
    return(default)
  }
  value
}

.scale_named_activity_vector <- function(activity, scale = 1) {
  if (length(activity) == 0) {
    return(numeric())
  }
  scaled <- as.numeric(activity) * scale
  names(scaled) <- names(activity)
  scaled
}

.merge_named_activity_vectors <- function(...) {
  vectors <- list(...)
  vectors <- purrr::keep(vectors, ~ !is.null(.x) && length(.x) > 0)
  if (length(vectors) == 0) {
    return(numeric())
  }

  merged <- unlist(vectors, use.names = TRUE)
  merged <- stats::setNames(
    object = as.numeric(merged),
    nm = names(merged)
  )
  merged <- tapply(merged, names(merged), sum)
  merged <- stats::setNames(
    object = as.numeric(merged),
    nm = names(merged)
  )
  merged <- merged[abs(merged) > .Machine$double.eps]
  merged[order(names(merged))]
}

.scenario_label_from_id <- function(scenario_id) {
  scenario_id %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title()
}

get_hallmark_scenario_catalog <- function() {
  list(
    myc_proliferation = list(
      label = "MYC / Proliferation",
      description = paste(
        "Cell-cycle heavy growth program with overlapping MYC, E2F,",
        "G2M, and a modest DNA repair component."
      ),
      theme_activities = c(
        HALLMARK_MYC_TARGETS_V1 = 1.25,
        HALLMARK_MYC_TARGETS_V2 = 1.0,
        HALLMARK_E2F_TARGETS = 0.95,
        HALLMARK_G2M_CHECKPOINT = 0.8,
        HALLMARK_DNA_REPAIR = 0.35
      ),
      focus_pathways = c(
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_E2F_TARGETS"
      ),
      contrast_pathways = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
    ),
    oxphos_metabolic = list(
      label = "OXPHOS / Metabolic",
      description = paste(
        "Mitochondrial and redox-skewed state with oxidative phosphorylation,",
        "ROS handling, peroxisome activity, and mild glycolysis suppression."
      ),
      theme_activities = c(
        HALLMARK_OXIDATIVE_PHOSPHORYLATION = 1.2,
        HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = 0.8,
        HALLMARK_PEROXISOME = 0.7,
        HALLMARK_GLYCOLYSIS = -0.45
      ),
      focus_pathways = c(
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
      ),
      contrast_pathways = c("HALLMARK_GLYCOLYSIS")
    ),
    dna_damage_response = list(
      label = "DNA Damage Response",
      description = paste(
        "Repair and stress-response program centered on DNA repair, p53,",
        "apoptosis, and a smaller checkpoint component."
      ),
      theme_activities = c(
        HALLMARK_DNA_REPAIR = 1.1,
        HALLMARK_P53_PATHWAY = 0.95,
        HALLMARK_APOPTOSIS = 0.7,
        HALLMARK_G2M_CHECKPOINT = 0.35
      ),
      focus_pathways = c(
        "HALLMARK_DNA_REPAIR",
        "HALLMARK_P53_PATHWAY"
      ),
      contrast_pathways = c("HALLMARK_MYC_TARGETS_V1")
    ),
    interferon_inflammatory = list(
      label = "Interferon / Inflammatory",
      description = paste(
        "Inflammatory state with type I and type II interferon activity,",
        "TNF/NFkB, IL6/JAK/STAT3, and a mild OXPHOS penalty."
      ),
      theme_activities = c(
        HALLMARK_INTERFERON_ALPHA_RESPONSE = 1.2,
        HALLMARK_INTERFERON_GAMMA_RESPONSE = 1.1,
        HALLMARK_TNFA_SIGNALING_VIA_NFKB = 0.8,
        HALLMARK_IL6_JAK_STAT3_SIGNALING = 0.7,
        HALLMARK_OXIDATIVE_PHOSPHORYLATION = -0.35
      ),
      focus_pathways = c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE"
      ),
      contrast_pathways = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    ),
    emt_hypoxia = list(
      label = "EMT / Hypoxia",
      description = paste(
        "Mesenchymal transition program with hypoxia and angiogenesis,",
        "alongside a small metabolic and proliferative penalty."
      ),
      theme_activities = c(
        HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = 1.25,
        HALLMARK_HYPOXIA = 0.95,
        HALLMARK_ANGIOGENESIS = 0.6,
        HALLMARK_OXIDATIVE_PHOSPHORYLATION = -0.5,
        HALLMARK_MYC_TARGETS_V1 = -0.3
      ),
      focus_pathways = c(
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_HYPOXIA"
      ),
      contrast_pathways = c(
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_MYC_TARGETS_V1"
      )
    )
  )
}

get_hallmark_teaching_dataset_catalog <- function() {
  vehicle_control <- c(
    HALLMARK_XENOBIOTIC_METABOLISM = 0.3,
    HALLMARK_PEROXISOME = 0.12
  )

  proliferation_core <- c(
    HALLMARK_E2F_TARGETS = -0.65,
    HALLMARK_G2M_CHECKPOINT = -0.6,
    HALLMARK_MYC_TARGETS_V1 = -0.4,
    HALLMARK_MYC_TARGETS_V2 = -0.3
  )

  inflammatory_core <- c(
    HALLMARK_INTERFERON_ALPHA_RESPONSE = 0.65,
    HALLMARK_INTERFERON_GAMMA_RESPONSE = 0.55,
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = 0.5,
    HALLMARK_INFLAMMATORY_RESPONSE = 0.45
  )

  emt_core <- c(
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = 0.7,
    HALLMARK_TGF_BETA_SIGNALING = 0.55,
    HALLMARK_HYPOXIA = 0.2
  )

  emt_core_strong <- .merge_named_activity_vectors(
    .scale_named_activity_vector(emt_core, 1.35),
    c(HALLMARK_GLYCOLYSIS = 0.12)
  )

  metabolic_core <- c(
    HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = 0.7,
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = -0.6
  )

  list(
    dataset1_proliferation_suppression = list(
      label = "Dataset 1: Proliferation Suppression",
      description = paste(
        "Balanced four-group teaching dataset where treatment progressively",
        "suppresses proliferation programs."
      ),
      group_activities = list(
        A = numeric(),
        B = vehicle_control,
        C = .merge_named_activity_vectors(vehicle_control, proliferation_core),
        D = .merge_named_activity_vectors(vehicle_control, .scale_named_activity_vector(proliferation_core, 1.35))
      ),
      vehicle_pathways = names(vehicle_control),
      expected_up_pathways = character(),
      expected_down_pathways = c(
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_MYC_TARGETS_V1"
      )
    ),
    dataset2_inflammatory_interferon = list(
      label = "Dataset 2: Inflammatory / Interferon",
      description = paste(
        "Balanced four-group teaching dataset with progressively stronger",
        "interferon and inflammatory signaling."
      ),
      group_activities = list(
        A = numeric(),
        B = vehicle_control,
        C = .merge_named_activity_vectors(vehicle_control, inflammatory_core),
        D = .merge_named_activity_vectors(vehicle_control, .scale_named_activity_vector(inflammatory_core, 1.35))
      ),
      vehicle_pathways = names(vehicle_control),
      expected_up_pathways = c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE"
      ),
      expected_down_pathways = character()
    ),
    dataset3_emt_tgfbeta = list(
      label = "Dataset 3: EMT / TGF-beta",
      description = paste(
        "Balanced four-group teaching dataset with an adaptive mesenchymal",
        "state shift and a subtle glycolytic undertone in the strongest group."
      ),
      group_activities = list(
        A = numeric(),
        B = vehicle_control,
        C = .merge_named_activity_vectors(vehicle_control, emt_core),
        D = .merge_named_activity_vectors(vehicle_control, emt_core_strong)
      ),
      vehicle_pathways = names(vehicle_control),
      expected_up_pathways = c(
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_TGF_BETA_SIGNALING"
      ),
      expected_down_pathways = character()
    ),
    dataset4_metabolic_stress = list(
      label = "Dataset 4: Metabolic Stress",
      description = paste(
        "Balanced four-group teaching dataset centered on oxidative stress",
        "and mitochondrial dysfunction."
      ),
      group_activities = list(
        A = numeric(),
        B = vehicle_control,
        C = .merge_named_activity_vectors(vehicle_control, metabolic_core),
        D = .merge_named_activity_vectors(vehicle_control, .scale_named_activity_vector(metabolic_core, 1.35))
      ),
      vehicle_pathways = names(vehicle_control),
      expected_up_pathways = c("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
      expected_down_pathways = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    )
  )
}

.normalize_hallmark_scenarios <- function(scenarios = NULL, available_pathways = NULL) {
  if (is.null(scenarios)) {
    scenarios <- get_hallmark_scenario_catalog()
  }

  if (!is.list(scenarios) || length(scenarios) == 0) {
    stop("scenarios must be a non-empty named list")
  }

  scenario_ids <- names(scenarios)
  if (is.null(scenario_ids) || any(!nzchar(scenario_ids))) {
    stop("scenarios must be named")
  }

  normalized <- vector("list", length(scenarios))
  names(normalized) <- scenario_ids

  for (idx in seq_along(scenarios)) {
    scenario_id <- scenario_ids[[idx]]
    entry <- scenarios[[idx]]

    if (is.numeric(entry)) {
      theme_activities <- entry
      entry <- list(theme_activities = entry)
    } else if (is.list(entry)) {
      theme_activities <- .or_default(entry$theme_activities, entry$activities)
    } else {
      stop("scenario entries must be numeric vectors or lists")
    }

    if (is.null(theme_activities) || !is.numeric(theme_activities) || length(theme_activities) == 0) {
      stop("scenario ", scenario_id, " must provide a non-empty numeric theme_activities vector")
    }

    activity_names <- names(theme_activities)
    if (is.null(activity_names) || any(!nzchar(activity_names))) {
      stop("scenario ", scenario_id, " must provide names for every theme activity")
    }

    theme_activities <- as.numeric(theme_activities)
    names(theme_activities) <- as.character(activity_names)

    if (anyDuplicated(names(theme_activities))) {
      stop("scenario ", scenario_id, " contains duplicate pathway names")
    }

    if (!is.null(available_pathways)) {
      missing_pathways <- setdiff(names(theme_activities), available_pathways)
      if (length(missing_pathways) > 0) {
        stop(
          "scenario ", scenario_id, " references Hallmark pathways not present in the supplied geneset: ",
          paste(missing_pathways, collapse = ", ")
        )
      }
    }

    default_focus <- names(theme_activities)[order(abs(theme_activities), decreasing = TRUE)]
    default_focus <- utils::head(default_focus, n = min(2, length(default_focus)))
    default_contrast <- names(theme_activities)[theme_activities < 0]
    if (length(default_contrast) > 0) {
      default_contrast <- utils::head(
        default_contrast[order(abs(theme_activities[default_contrast]), decreasing = TRUE)],
        n = min(2, length(default_contrast))
      )
    }

    focus_pathways <- unique(as.character(.or_default(entry$focus_pathways, default_focus)))
    contrast_pathways <- unique(as.character(.or_default(entry$contrast_pathways, default_contrast)))

    focus_pathways <- focus_pathways[focus_pathways %in% names(theme_activities)]
    contrast_pathways <- contrast_pathways[contrast_pathways %in% names(theme_activities)]

    normalized[[idx]] <- list(
      scenario_id = scenario_id,
      scenario_order = idx,
      label = as.character(.or_default(entry$label, .scenario_label_from_id(scenario_id))),
      description = as.character(.or_default(entry$description, "")),
      theme_activities = theme_activities,
      focus_pathways = focus_pathways,
      contrast_pathways = contrast_pathways
    )
  }

  normalized
}

.scenario_catalog_to_metadata <- function(scenarios) {
  purrr::map_dfr(scenarios, function(entry) {
    theme_names <- names(entry$theme_activities)
    tibble::tibble(
      scenario_id = entry$scenario_id,
      scenario_order = entry$scenario_order,
      scenario_label = entry$label,
      scenario_description = entry$description,
      pathway = theme_names,
      activity = unname(entry$theme_activities),
      direction = ifelse(unname(entry$theme_activities) >= 0, "up", "down"),
      is_focus = theme_names %in% entry$focus_pathways,
      is_contrast = theme_names %in% entry$contrast_pathways
    )
  })
}

.normalize_grouped_dataset_groups <- function(
    entry,
    dataset_id,
    available_pathways = NULL,
    required_groups = NULL,
    default_vehicle_groups = character()) {
  group_entries <- .or_default(entry$groups, entry$group_activities)
  if (!is.list(group_entries) || length(group_entries) == 0) {
    stop("dataset ", dataset_id, " must provide groups or group_activities as a named list")
  }

  group_ids <- names(group_entries)
  if (is.null(group_ids) || any(!nzchar(group_ids))) {
    stop("dataset ", dataset_id, " groups must be named")
  }

  if (!is.null(required_groups)) {
    missing_groups <- setdiff(required_groups, group_ids)
    if (length(missing_groups) > 0) {
      stop("dataset ", dataset_id, " is missing required groups: ", paste(missing_groups, collapse = ", "))
    }
    group_ids <- required_groups
  }

  group_labels <- .or_default(entry$group_labels, list())
  if (!is.list(group_labels)) {
    group_labels <- as.list(group_labels)
  }

  group_descriptions <- .or_default(entry$group_descriptions, list())
  if (!is.list(group_descriptions)) {
    group_descriptions <- as.list(group_descriptions)
  }

  treatment_levels <- .or_default(entry$treatment_levels, list())
  if (!is.list(treatment_levels)) {
    treatment_levels <- as.list(treatment_levels)
  }

  vehicle_groups <- unique(as.character(.or_default(entry$vehicle_groups, default_vehicle_groups)))
  normalized_groups <- vector("list", length(group_ids))
  names(normalized_groups) <- group_ids
  default_treatment_levels <- if (length(group_ids) == 1) 0 else seq(0, 1, length.out = length(group_ids))

  for (idx in seq_along(group_ids)) {
    group_id <- group_ids[[idx]]
    group_entry <- group_entries[[group_id]]

    if (is.null(group_entry)) {
      stop("dataset ", dataset_id, " group ", group_id, " is missing")
    }

    if (is.numeric(group_entry)) {
      activity <- group_entry
      gene_spikes <- numeric()
      group_entry <- list(activities = group_entry)
    } else if (is.list(group_entry)) {
      activity <- .or_default(group_entry$activities, group_entry$theme_activities)
      if (is.null(activity)) {
        activity <- numeric()
      }
      gene_spikes <- .or_default(group_entry$gene_spikes, group_entry$spike_genes)
      if (is.null(gene_spikes)) {
        gene_spikes <- numeric()
      }
    } else {
      stop("dataset ", dataset_id, " group ", group_id, " must be a numeric vector or list")
    }

    if (length(activity) > 0) {
      if (!is.numeric(activity)) {
        stop("dataset ", dataset_id, " group ", group_id, " activities must be numeric")
      }
      if (is.null(names(activity)) || any(!nzchar(names(activity)))) {
        stop("dataset ", dataset_id, " group ", group_id, " activities must be named")
      }
      activity <- .merge_named_activity_vectors(activity)
      if (!is.null(available_pathways)) {
        missing_pathways <- setdiff(names(activity), available_pathways)
        if (length(missing_pathways) > 0) {
          stop(
            "dataset ", dataset_id, " group ", group_id,
            " references Hallmark pathways not present in the supplied geneset: ",
            paste(missing_pathways, collapse = ", ")
          )
        }
      }
    } else {
      activity <- numeric()
    }

    if (length(gene_spikes) > 0) {
      if (!is.numeric(gene_spikes)) {
        stop("dataset ", dataset_id, " group ", group_id, " gene_spikes must be numeric")
      }
      if (is.null(names(gene_spikes)) || any(!nzchar(names(gene_spikes)))) {
        stop("dataset ", dataset_id, " group ", group_id, " gene_spikes must be named")
      }
      gene_spikes <- .merge_named_activity_vectors(gene_spikes)
    } else {
      gene_spikes <- numeric()
    }

    group_label <- .or_default(group_entry$label, group_labels[[group_id]])
    if (is.null(group_label) || !nzchar(as.character(group_label[[1]]))) {
      group_label <- .scenario_label_from_id(group_id)
    }

    group_description <- .or_default(group_entry$description, group_descriptions[[group_id]])
    if (is.null(group_description)) {
      group_description <- ""
    }

    treatment_level <- .or_default(group_entry$treatment_level, treatment_levels[[group_id]])
    if (is.null(treatment_level)) {
      treatment_level <- default_treatment_levels[[idx]]
    }
    treatment_level <- as.numeric(treatment_level[[1]])
    if (is.na(treatment_level)) {
      stop("dataset ", dataset_id, " group ", group_id, " treatment_level must be numeric")
    }

    group_is_vehicle <- .or_default(group_entry$is_vehicle, group_id %in% vehicle_groups)
    group_is_vehicle <- isTRUE(group_is_vehicle)

    normalized_groups[[idx]] <- list(
      group = group_id,
      group_order = idx,
      label = as.character(group_label[[1]]),
      description = as.character(group_description[[1]]),
      activities = activity,
      gene_spikes = gene_spikes,
      is_vehicle = group_is_vehicle,
      treatment_level = treatment_level
    )
  }

  normalized_groups
}

.normalize_grouped_hallmark_datasets <- function(
    datasets = NULL,
    available_pathways = NULL,
    required_groups = NULL,
    default_vehicle_groups = character()) {
  if (is.null(datasets)) {
    datasets <- get_hallmark_teaching_dataset_catalog()
  }

  if (!is.list(datasets) || length(datasets) == 0) {
    stop("datasets must be a non-empty named list")
  }

  dataset_ids <- names(datasets)
  if (is.null(dataset_ids) || any(!nzchar(dataset_ids))) {
    stop("datasets must be named")
  }

  normalized <- vector("list", length(datasets))
  names(normalized) <- dataset_ids

  for (idx in seq_along(datasets)) {
    dataset_id <- dataset_ids[[idx]]
    entry <- datasets[[idx]]

    if (!is.list(entry)) {
      stop("dataset entries must be lists")
    }

    normalized_groups <- .normalize_grouped_dataset_groups(
      entry = entry,
      dataset_id = dataset_id,
      available_pathways = available_pathways,
      required_groups = required_groups,
      default_vehicle_groups = default_vehicle_groups
    )

    vehicle_pathways_default <- unique(unlist(purrr::map(normalized_groups, function(group_entry) {
      if (!isTRUE(group_entry$is_vehicle)) {
        return(character())
      }
      names(group_entry$activities)
    }), use.names = FALSE))

    union_pathways <- unique(unlist(purrr::map(normalized_groups, ~ names(.x$activities)), use.names = FALSE))
    union_pathways <- union_pathways[nzchar(union_pathways)]

    vehicle_pathways <- as.character(.or_default(entry$vehicle_pathways, vehicle_pathways_default))
    vehicle_pathways <- vehicle_pathways[vehicle_pathways %in% union_pathways]

    expected_up <- unique(as.character(.or_default(entry$expected_up_pathways, character())))
    expected_down <- unique(as.character(.or_default(entry$expected_down_pathways, character())))

    normalized[[idx]] <- list(
      dataset_id = dataset_id,
      dataset_order = idx,
      label = as.character(.or_default(entry$label, .scenario_label_from_id(dataset_id))),
      description = as.character(.or_default(entry$description, "")),
      groups = normalized_groups,
      group_activities = purrr::map(normalized_groups, "activities"),
      vehicle_pathways = vehicle_pathways,
      expected_up_pathways = expected_up,
      expected_down_pathways = expected_down
    )
  }

  normalized
}

.grouped_dataset_catalog_to_metadata <- function(datasets) {
  purrr::map_dfr(datasets, function(entry) {
    purrr::map_dfr(entry$groups, function(group_entry) {
      activity <- group_entry$activities
      if (length(activity) == 0) {
        return(tibble::tibble(
          dataset_id = entry$dataset_id,
          dataset_order = entry$dataset_order,
          dataset_label = entry$label,
          dataset_description = entry$description,
          group = group_entry$group,
          group_order = group_entry$group_order,
          group_label = group_entry$label,
          group_description = group_entry$description,
          group_is_vehicle = group_entry$is_vehicle,
          treatment_level = group_entry$treatment_level,
          pathway = NA_character_,
          activity = 0,
          direction = "neutral",
          is_vehicle = FALSE,
          is_expected_up = FALSE,
          is_expected_down = FALSE
        ))
      }

      tibble::tibble(
        dataset_id = entry$dataset_id,
        dataset_order = entry$dataset_order,
        dataset_label = entry$label,
        dataset_description = entry$description,
        group = group_entry$group,
        group_order = group_entry$group_order,
        group_label = group_entry$label,
        group_description = group_entry$description,
        group_is_vehicle = group_entry$is_vehicle,
        treatment_level = group_entry$treatment_level,
        pathway = names(activity),
        activity = unname(activity),
        direction = dplyr::case_when(
          activity > 0 ~ "up",
          activity < 0 ~ "down",
          TRUE ~ "neutral"
        ),
        is_vehicle = pathway %in% entry$vehicle_pathways,
        is_expected_up = pathway %in% entry$expected_up_pathways,
        is_expected_down = pathway %in% entry$expected_down_pathways
      )
    })
  })
}

.normalize_teaching_datasets <- function(datasets = NULL, available_pathways = NULL) {
  .normalize_grouped_hallmark_datasets(
    datasets = datasets,
    available_pathways = available_pathways,
    required_groups = c("A", "B", "C", "D"),
    default_vehicle_groups = c("B", "C", "D")
  )
}

.teaching_dataset_catalog_to_metadata <- function(datasets) {
  .grouped_dataset_catalog_to_metadata(datasets)
}

.prepare_hallmark_gene_lookup <- function(geneset) {
  geneset %>%
    dplyr::transmute(
      gene_id = as.character(ncbi_gene),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::arrange(gene_id)
}

.empty_group_gene_spike_metadata <- function() {
  tibble::tibble(
    dataset_id = character(),
    dataset_label = character(),
    group = character(),
    group_order = integer(),
    group_label = character(),
    group_description = character(),
    gene_id = character(),
    gene_symbol = character(),
    gene_spike_value = numeric(),
    gene_spike_target_name = character(),
    gene_spike_match_type = character()
  )
}

.resolve_group_gene_spikes <- function(dataset_entry, gene_lookup) {
  resolved_all <- purrr::map_dfr(dataset_entry$groups, function(group_entry) {
    gene_spikes <- group_entry$gene_spikes
    if (length(gene_spikes) == 0) {
      return(.empty_group_gene_spike_metadata())
    }

    spike_df <- tibble::tibble(
      target_name = names(gene_spikes),
      gene_spike_value = as.numeric(gene_spikes)
    )

    by_id <- spike_df %>%
      dplyr::filter(target_name %in% gene_lookup$gene_id) %>%
      dplyr::transmute(
        gene_id = target_name,
        gene_spike_value = gene_spike_value,
        gene_spike_target_name = target_name,
        gene_spike_match_type = "gene_id"
      )

    remaining <- spike_df %>%
      dplyr::filter(!target_name %in% gene_lookup$gene_id)

    by_symbol <- remaining %>%
      dplyr::left_join(
        gene_lookup %>%
          dplyr::filter(!is.na(gene_symbol), nzchar(gene_symbol)) %>%
          dplyr::distinct(gene_symbol, gene_id),
        by = c("target_name" = "gene_symbol")
      )

    unresolved <- by_symbol %>%
      dplyr::filter(is.na(gene_id)) %>%
      dplyr::pull(target_name) %>%
      unique()
    if (length(unresolved) > 0) {
      stop(
        "dataset ", dataset_entry$dataset_id, " group ", group_entry$group,
        " gene_spikes reference genes not present in the simulated gene universe: ",
        paste(unresolved, collapse = ", ")
      )
    }

    by_symbol <- by_symbol %>%
      dplyr::transmute(
        gene_id = gene_id,
        gene_spike_value = gene_spike_value,
        gene_spike_target_name = target_name,
        gene_spike_match_type = "gene_symbol"
      )

    dplyr::bind_rows(by_id, by_symbol) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        gene_spike_value = sum(gene_spike_value),
        gene_spike_target_name = paste(unique(gene_spike_target_name), collapse = " | "),
        gene_spike_match_type = paste(unique(gene_spike_match_type), collapse = " | "),
        .groups = "drop"
      ) %>%
      dplyr::left_join(gene_lookup, by = "gene_id") %>%
      dplyr::transmute(
        dataset_id = dataset_entry$dataset_id,
        dataset_label = dataset_entry$label,
        group = group_entry$group,
        group_order = group_entry$group_order,
        group_label = group_entry$label,
        group_description = group_entry$description,
        gene_id = gene_id,
        gene_symbol = gene_symbol,
        gene_spike_value = gene_spike_value,
        gene_spike_target_name = gene_spike_target_name,
        gene_spike_match_type = gene_spike_match_type
      )
  })

  resolved_all %>%
    dplyr::arrange(group_order, gene_symbol, gene_id)
}

simulate_signed_hallmark_scenarios <- function(

    scenarios = NULL,
    geneset = NULL,
    seed = 20260316,
    species = "Homo sapiens",
    n_replicates = 3,
    gamma_range = c(0.85, 1.15),
    lambda = 0.2,
    sigma = 0.5,
    weight_shape = 2,
    weight_rate = 2) {
  if (is.null(geneset)) {
    geneset <- geneset_tools$get_collection("H", "", species = species)
  }

  if (!"gs_name" %in% colnames(geneset) || !"ncbi_gene" %in% colnames(geneset)) {
    stop("geneset must contain at least gs_name and ncbi_gene columns")
  }

  seed <- as.integer(seed[[1]])
  n_replicates <- as.integer(n_replicates[[1]])
  if (is.na(seed)) {
    stop("seed must be coercible to an integer")
  }
  if (is.na(n_replicates) || n_replicates < 1) {
    stop("n_replicates must be a positive integer")
  }
  gamma_range <- as.numeric(gamma_range)
  lambda <- as.numeric(lambda[[1]])
  sigma <- as.numeric(sigma[[1]])
  weight_shape <- as.numeric(weight_shape[[1]])
  weight_rate <- as.numeric(weight_rate[[1]])
  if (length(gamma_range) != 2 || any(is.na(gamma_range)) || gamma_range[[1]] > gamma_range[[2]]) {
    stop("gamma_range must be a numeric vector of length 2 with min <= max")
  }
  if (is.na(lambda)) {
    stop("lambda must be numeric")
  }
  if (is.na(sigma) || sigma < 0) {
    stop("sigma must be a non-negative numeric value")
  }
  if (is.na(weight_shape) || weight_shape <= 0 || is.na(weight_rate) || weight_rate <= 0) {
    stop("weight_shape and weight_rate must be positive numeric values")
  }

  hallmark_geneset <- geneset %>%
    dplyr::mutate(
      gs_name = as.character(gs_name),
      gene_id = as.character(ncbi_gene),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id))

  available_pathways <- sort(unique(hallmark_geneset$gs_name))
  scenarios_norm <- .normalize_hallmark_scenarios(
    scenarios = scenarios,
    available_pathways = available_pathways
  )

  gene_lookup <- .prepare_hallmark_gene_lookup(hallmark_geneset)
  scenario_metadata <- .scenario_catalog_to_metadata(scenarios_norm)

  gene_metadata_all <- list()
  membership_metadata_all <- list()
  rank_dfs <- list()

  for (scenario_entry in scenarios_norm) {
    scenario_id <- scenario_entry$scenario_id
    scenario_seed <- seed + (scenario_entry$scenario_order * 1000L)

    active_memberships <- hallmark_geneset %>%
      dplyr::filter(gs_name %in% names(scenario_entry$theme_activities)) %>%
      dplyr::transmute(
        scenario_id = scenario_id,
        scenario_order = scenario_entry$scenario_order,
        scenario_label = scenario_entry$label,
        scenario_description = scenario_entry$description,
        pathway = gs_name,
        gene_id = gene_id,
        gene_symbol = gene_symbol,
        activity = unname(scenario_entry$theme_activities[gs_name])
      ) %>%
      dplyr::distinct(pathway, gene_id, .keep_all = TRUE)

    if (nrow(active_memberships) == 0) {
      stop("scenario ", scenario_id, " did not match any Hallmark genes")
    }

    set.seed(scenario_seed + 1L)
    active_memberships$raw_weight <- stats::rgamma(
      n = nrow(active_memberships),
      shape = weight_shape,
      rate = weight_rate
    )

    active_memberships <- active_memberships %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(
        weight = raw_weight / sum(raw_weight),
        n_active_memberships = dplyr::n(),
        activity_contribution = weight * activity
      ) %>%
      dplyr::ungroup()

    set.seed(scenario_seed + 2L)
    gamma_lookup <- gene_lookup %>%
      dplyr::mutate(
        gamma = stats::runif(
          n = dplyr::n(),
          min = gamma_range[[1]],
          max = gamma_range[[2]]
        )
      )

    gene_summary <- active_memberships %>%
      dplyr::group_by(
        scenario_id,
        scenario_order,
        scenario_label,
        scenario_description,
        gene_id
      ) %>%
      dplyr::summarise(
        gene_symbol = dplyr::first(gene_symbol),
        n_active_memberships = dplyr::first(n_active_memberships),
        weight_total = sum(weight),
        weighted_activity = sum(activity_contribution),
        active_pathways = paste(pathway, collapse = " | "),
        .groups = "drop"
      )

    gene_summary <- gene_lookup %>%
      dplyr::left_join(gamma_lookup, by = c("gene_id", "gene_symbol")) %>%
      dplyr::left_join(gene_summary, by = c("gene_id", "gene_symbol")) %>%
      dplyr::mutate(
        scenario_id = dplyr::coalesce(scenario_id, scenario_entry$scenario_id),
        scenario_order = dplyr::coalesce(scenario_order, scenario_entry$scenario_order),
        scenario_label = dplyr::coalesce(scenario_label, scenario_entry$label),
        scenario_description = dplyr::coalesce(scenario_description, scenario_entry$description),
        n_active_memberships = dplyr::coalesce(n_active_memberships, 0L),
        weight_total = dplyr::coalesce(weight_total, 0),
        weighted_activity = dplyr::coalesce(weighted_activity, 0),
        active_pathways = dplyr::coalesce(active_pathways, ""),
        overlap_offset = dplyr::if_else(
          n_active_memberships > 1,
          lambda * (n_active_memberships - 1) * ((gamma - 1) / (1 + gamma)),
          0
        ),
        signal_mean = (gamma * weighted_activity) + overlap_offset,
        is_active = n_active_memberships > 0
      ) %>%
      dplyr::arrange(gene_id)

    membership_metadata <- active_memberships %>%
      dplyr::left_join(
        gene_summary %>%
          dplyr::select(
            gene_id,
            gamma,
            overlap_offset,
            signal_mean,
            weight_total
          ),
        by = "gene_id"
      ) %>%
      dplyr::arrange(pathway, gene_id)

    membership_metadata_all[[scenario_id]] <- membership_metadata

    replicate_metadata <- vector("list", n_replicates)
    for (rep_idx in seq_len(n_replicates)) {
      rank_name <- sprintf("%s_rep%02d", scenario_id, rep_idx)
      set.seed(scenario_seed + 100L + rep_idx)
      noise <- stats::rnorm(nrow(gene_summary), mean = 0, sd = sigma)

      gene_rep <- gene_summary %>%
        dplyr::mutate(
          replicate_id = rep_idx,
          rank_name = rank_name,
          noise = noise,
          value = signal_mean + noise
        ) %>%
        dplyr::select(
          scenario_id,
          scenario_order,
          scenario_label,
          scenario_description,
          replicate_id,
          rank_name,
          gene_id,
          gene_symbol,
          gamma,
          n_active_memberships,
          weight_total,
          weighted_activity,
          overlap_offset,
          signal_mean,
          noise,
          value,
          is_active,
          active_pathways
        )

      replicate_metadata[[rep_idx]] <- gene_rep
      rank_dfs[[rank_name]] <- gene_rep %>%
        dplyr::transmute(
          id = gene_id,
          value = value
        )
    }

    gene_metadata_all[[scenario_id]] <- dplyr::bind_rows(replicate_metadata)
  }

  list(
    rank_dfs = rank_dfs,
    gene_metadata = dplyr::bind_rows(gene_metadata_all),
    membership_metadata = dplyr::bind_rows(membership_metadata_all),
    scenario_metadata = scenario_metadata
  )
}

.extract_signed_hallmark_targets <- function(
    score_simulation,
    target_source = c("signal_mean", "value"),
    target_replicate = 1) {
  target_source <- match.arg(target_source)
  target_replicate <- as.integer(target_replicate[[1]])

  if (!is.list(score_simulation) || is.null(score_simulation$gene_metadata)) {
    stop("score_simulation must be an output list from simulate_signed_hallmark_scenarios()")
  }
  if (is.na(target_replicate) || target_replicate < 1) {
    stop("target_replicate must be a positive integer")
  }

  static_cols <- c(
    "scenario_id",
    "scenario_order",
    "scenario_label",
    "scenario_description",
    "gene_id",
    "gene_symbol",
    "gamma",
    "n_active_memberships",
    "weight_total",
    "weighted_activity",
    "overlap_offset",
    "signal_mean",
    "is_active",
    "active_pathways"
  )

  if (target_source == "signal_mean") {
    target_df <- score_simulation$gene_metadata %>%
      dplyr::select(dplyr::all_of(static_cols)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        target_source = "signal_mean",
        target_replicate = NA_integer_,
        target_value = signal_mean
      )
  } else {
    target_df <- score_simulation$gene_metadata %>%
      dplyr::filter(replicate_id == target_replicate) %>%
      dplyr::select(dplyr::all_of(static_cols), replicate_id, rank_name, value) %>%
      dplyr::mutate(
        target_source = "value",
        target_replicate = target_replicate,
        target_value = value
      )

    if (nrow(target_df) == 0) {
      stop("target_replicate ", target_replicate, " was not found in score_simulation$gene_metadata")
    }
  }

  target_df %>%
    dplyr::arrange(scenario_order, gene_id)
}

simulate_signed_hallmark_expression <- function(
    score_simulation = NULL,
    scenarios = NULL,
    geneset = NULL,
    seed = 20260317,
    species = "Homo sapiens",
    target_source = c("signal_mean", "value"),
    target_replicate = 1,
    n_samples_per_group = 4,
    n_baseline_samples = n_samples_per_group,
    baseline_id = "baseline",
    baseline_label = "Baseline",
    baseline_range = c(4, 7),
    baseline_shape1 = 2,
    baseline_shape2 = 2,
    sigma_min = 0.06,
    sigma_max = 0.4,
    sigma_abundance_k = -1.2,
    sigma_abundance_anchor = NULL,
    effect_scale = 3,
    effect_saturation = 2,
    trend = TRUE,
    robust = TRUE,
    export_gct_path = NULL,
    score_gamma_range = c(0.85, 1.15),
    score_lambda = 0.2,
    score_sigma = 0.5,
    score_weight_shape = 2,
    score_weight_rate = 2) {
  target_source <- match.arg(target_source)
  seed <- as.integer(seed[[1]])
  n_samples_per_group <- as.integer(n_samples_per_group[[1]])
  n_baseline_samples <- as.integer(n_baseline_samples[[1]])
  target_replicate <- as.integer(target_replicate[[1]])
  baseline_range <- as.numeric(baseline_range)
  sigma_min <- as.numeric(sigma_min[[1]])
  sigma_max <- as.numeric(sigma_max[[1]])
  sigma_abundance_k <- as.numeric(sigma_abundance_k[[1]])
  effect_scale <- as.numeric(effect_scale[[1]])
  effect_saturation <- as.numeric(effect_saturation[[1]])
  baseline_shape1 <- as.numeric(baseline_shape1[[1]])
  baseline_shape2 <- as.numeric(baseline_shape2[[1]])

  if (is.na(seed)) {
    stop("seed must be coercible to an integer")
  }
  if (is.na(n_samples_per_group) || n_samples_per_group < 1) {
    stop("n_samples_per_group must be a positive integer")
  }
  if (is.na(n_baseline_samples) || n_baseline_samples < 1) {
    stop("n_baseline_samples must be a positive integer")
  }
  if (is.na(target_replicate) || target_replicate < 1) {
    stop("target_replicate must be a positive integer")
  }
  if (length(baseline_range) != 2 || any(is.na(baseline_range)) || baseline_range[[1]] >= baseline_range[[2]]) {
    stop("baseline_range must be a numeric vector of length 2 with min < max")
  }
  if (is.na(sigma_min) || is.na(sigma_max) || sigma_min <= 0 || sigma_max < sigma_min) {
    stop("sigma_min and sigma_max must be positive values with sigma_max >= sigma_min")
  }
  if (is.na(sigma_abundance_k)) {
    stop("sigma_abundance_k must be numeric")
  }
  if (is.na(effect_scale) || effect_scale < 0) {
    stop("effect_scale must be a non-negative numeric value")
  }
  if (is.na(effect_saturation) || effect_saturation <= 0) {
    stop("effect_saturation must be a positive numeric value")
  }
  if (is.na(baseline_shape1) || baseline_shape1 <= 0 || is.na(baseline_shape2) || baseline_shape2 <= 0) {
    stop("baseline_shape1 and baseline_shape2 must be positive numeric values")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("simulate_signed_hallmark_expression requires the limma package")
  }

  if (is.null(score_simulation)) {
    n_score_replicates <- if (target_source == "value") target_replicate else 1L
    score_simulation <- simulate_signed_hallmark_scenarios(
      scenarios = scenarios,
      geneset = geneset,
      seed = seed,
      species = species,
      n_replicates = n_score_replicates,
      gamma_range = score_gamma_range,
      lambda = score_lambda,
      sigma = score_sigma,
      weight_shape = score_weight_shape,
      weight_rate = score_weight_rate
    )
  }

  if (is.null(score_simulation$scenario_metadata) || is.null(score_simulation$gene_metadata)) {
    stop("score_simulation must include scenario_metadata and gene_metadata")
  }

  target_df <- .extract_signed_hallmark_targets(
    score_simulation = score_simulation,
    target_source = target_source,
    target_replicate = target_replicate
  )

  scenario_info <- score_simulation$scenario_metadata %>%
    dplyr::distinct(
      scenario_id,
      scenario_order,
      scenario_label,
      scenario_description
    ) %>%
    dplyr::arrange(scenario_order)

  feature_metadata <- target_df %>%
    dplyr::distinct(gene_id, gene_symbol) %>%
    dplyr::arrange(gene_id)

  sigma_anchor <- .or_default(sigma_abundance_anchor, baseline_range[[1]])
  sigma_anchor <- as.numeric(sigma_anchor[[1]])
  if (is.na(sigma_anchor)) {
    stop("sigma_abundance_anchor must be numeric when provided")
  }

  set.seed(seed + 5000L)
  baseline_scaled <- stats::rbeta(
    n = nrow(feature_metadata),
    shape1 = baseline_shape1,
    shape2 = baseline_shape2
  )

  feature_metadata <- feature_metadata %>%
    dplyr::mutate(
      baseline_log10 = baseline_range[[1]] + (diff(baseline_range) * baseline_scaled),
      sigma_gene = sigma_min + ((sigma_max - sigma_min) * exp(sigma_abundance_k * (baseline_log10 - sigma_anchor))),
      sigma_gene = pmin(pmax(sigma_gene, sigma_min), sigma_max)
    )

  contrast_se_scale <- sqrt((1 / n_samples_per_group) + (1 / n_baseline_samples))

  scenario_gene_metadata <- target_df %>%
    dplyr::left_join(feature_metadata, by = "gene_id") %>%
    dplyr::mutate(
      gene_symbol = dplyr::coalesce(gene_symbol.x, gene_symbol.y),
      target_value_shrunk = effect_saturation * tanh(target_value / effect_saturation),
      mean_shift = effect_scale * sigma_gene * contrast_se_scale * target_value_shrunk
    ) %>%
    dplyr::select(-gene_symbol.x, -gene_symbol.y) %>%
    dplyr::arrange(scenario_order, gene_id)

  baseline_samples <- tibble::tibble(
    id = sprintf("%s_rep%02d", baseline_id, seq_len(n_baseline_samples)),
    group = baseline_id,
    group_label = baseline_label,
    scenario_id = baseline_id,
    scenario_label = baseline_label,
    group_order = 0L,
    sample_index = seq_len(n_baseline_samples),
    is_baseline = TRUE
  )

  scenario_samples <- purrr::map_dfr(seq_len(nrow(scenario_info)), function(idx) {
    entry <- scenario_info[idx, , drop = FALSE]
    tibble::tibble(
      id = sprintf("%s_rep%02d", entry$scenario_id, seq_len(n_samples_per_group)),
      group = entry$scenario_id,
      group_label = entry$scenario_label,
      scenario_id = entry$scenario_id,
      scenario_label = entry$scenario_label,
      group_order = entry$scenario_order,
      sample_index = seq_len(n_samples_per_group),
      is_baseline = FALSE
    )
  })

  sample_metadata <- dplyr::bind_rows(baseline_samples, scenario_samples) %>%
    dplyr::arrange(group_order, sample_index) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(sample_metadata) <- sample_metadata$id

  expr_mat <- matrix(
    NA_real_,
    nrow = nrow(feature_metadata),
    ncol = nrow(sample_metadata),
    dimnames = list(feature_metadata$gene_id, sample_metadata$id)
  )

  set.seed(seed + 6000L)
  baseline_draws <- stats::rnorm(
    n = nrow(feature_metadata) * n_baseline_samples,
    mean = rep(feature_metadata$baseline_log10, times = n_baseline_samples),
    sd = rep(feature_metadata$sigma_gene, times = n_baseline_samples)
  )
  expr_mat[, baseline_samples$id] <- matrix(
    baseline_draws,
    nrow = nrow(feature_metadata),
    ncol = n_baseline_samples
  )

  for (idx in seq_len(nrow(scenario_info))) {
    scenario_id <- scenario_info$scenario_id[[idx]]
    sample_ids <- sample_metadata$id[sample_metadata$scenario_id == scenario_id]
    gene_df <- scenario_gene_metadata %>%
      dplyr::filter(scenario_id == !!scenario_id) %>%
      dplyr::arrange(match(gene_id, feature_metadata$gene_id))

    set.seed(seed + 7000L + scenario_info$scenario_order[[idx]])
    scenario_draws <- stats::rnorm(
      n = nrow(feature_metadata) * length(sample_ids),
      mean = rep(feature_metadata$baseline_log10 + gene_df$mean_shift, times = length(sample_ids)),
      sd = rep(feature_metadata$sigma_gene, times = length(sample_ids))
    )
    expr_mat[, sample_ids] <- matrix(
      scenario_draws,
      nrow = nrow(feature_metadata),
      ncol = length(sample_ids)
    )
  }

  rdesc <- feature_metadata %>%
    dplyr::transmute(
      id = gene_id,
      gene_symbol = gene_symbol,
      baseline_log10 = baseline_log10,
      sigma_gene = sigma_gene
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(rdesc) <- rdesc$id

  cdesc <- sample_metadata %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(cdesc) <- cdesc$id

  expression_gct <- cmapR::GCT(
    mat = expr_mat,
    rid = rownames(expr_mat),
    cid = colnames(expr_mat),
    rdesc = rdesc,
    cdesc = cdesc
  )

  if (!is.null(export_gct_path)) {
    export_gct_path <- fs::path_abs(export_gct_path)
    fs::dir_create(fs::path_dir(export_gct_path), recurse = TRUE)
    cmapR::write_gct(expression_gct, export_gct_path, appenddim = FALSE)
  }

  group_levels <- c(baseline_id, scenario_info$scenario_id)
  group_factor <- factor(sample_metadata$group, levels = group_levels)
  design <- stats::model.matrix(~ 0 + group_factor)
  colnames(design) <- group_levels

  contrast_exprs <- stats::setNames(
    object = purrr::map_chr(
      scenario_info$scenario_id,
      ~ paste0(.x, " - ", baseline_id)
    ),
    nm = scenario_info$scenario_id
  )

  fit <- limma::lmFit(expr_mat, design)
  contrast_matrix <- do.call(
    limma::makeContrasts,
    c(as.list(contrast_exprs), list(levels = design))
  )
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)

  target_rank_dfs <- purrr::map(
    scenario_info$scenario_id,
    function(scenario_id) {
      scenario_gene_metadata %>%
        dplyr::filter(scenario_id == !!scenario_id) %>%
        dplyr::transmute(
          id = gene_id,
          value = target_value
        )
    }
  )
  names(target_rank_dfs) <- scenario_info$scenario_id

  limma_tables <- list()
  recovered_rank_dfs <- list()
  recovery_metrics <- list()
  gene_symbol_map <- feature_metadata$gene_symbol
  names(gene_symbol_map) <- feature_metadata$gene_id

  for (scenario_id in scenario_info$scenario_id) {
    top_table <- limma::topTable(
      fit,
      coef = scenario_id,
      number = Inf,
      sort.by = "none"
    )
    top_table$GeneID <- rownames(top_table)
    top_table$GeneSymbol <- gene_symbol_map[as.character(top_table$GeneID)]
    top_table$signedlogP <- sign(top_table$logFC) * -log10(pmax(top_table$P.Value, .Machine$double.eps))
    top_table <- top_table[, c(
      "GeneID",
      "GeneSymbol",
      setdiff(colnames(top_table), c("GeneID", "GeneSymbol"))
    )]

    limma_tables[[scenario_id]] <- top_table
    recovered_rank_dfs[[scenario_id]] <- top_table %>%
      dplyr::transmute(
        id = GeneID,
        value = signedlogP
      )

    merged <- scenario_gene_metadata %>%
      dplyr::filter(scenario_id == !!scenario_id) %>%
      dplyr::transmute(
        gene_id,
        target_value,
        target_value_shrunk,
        mean_shift,
        sigma_gene
      ) %>%
      dplyr::left_join(
        top_table %>%
          dplyr::transmute(
            gene_id = GeneID,
            recovered_t = t,
            recovered_logFC = logFC,
            recovered_signedlogP = signedlogP
          ),
        by = "gene_id"
      )

    recovery_metrics[[scenario_id]] <- tibble::tibble(
      scenario_id = scenario_id,
      pearson_signedlogP = stats::cor(
        merged$target_value,
        merged$recovered_signedlogP,
        method = "pearson",
        use = "pairwise.complete.obs"
      ),
      spearman_signedlogP = stats::cor(
        merged$target_value,
        merged$recovered_signedlogP,
        method = "spearman",
        use = "pairwise.complete.obs"
      ),
      pearson_t = stats::cor(
        merged$target_value_shrunk,
        merged$recovered_t,
        method = "pearson",
        use = "pairwise.complete.obs"
      ),
      spearman_t = stats::cor(
        merged$target_value_shrunk,
        merged$recovered_t,
        method = "spearman",
        use = "pairwise.complete.obs"
      ),
      pearson_logFC = stats::cor(
        merged$mean_shift,
        merged$recovered_logFC,
        method = "pearson",
        use = "pairwise.complete.obs"
      ),
      spearman_logFC = stats::cor(
        merged$mean_shift,
        merged$recovered_logFC,
        method = "spearman",
        use = "pairwise.complete.obs"
      ),
      rmse_signedlogP = sqrt(mean((merged$target_value - merged$recovered_signedlogP)^2, na.rm = TRUE)),
      mean_abs_logFC = mean(abs(merged$recovered_logFC), na.rm = TRUE),
      mean_sigma_gene = mean(merged$sigma_gene, na.rm = TRUE)
    )
  }

  list(
    expression_gct = expression_gct,
    expression_matrix = expr_mat,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    scenario_gene_metadata = scenario_gene_metadata,
    target_rank_dfs = target_rank_dfs,
    recovered_rank_dfs = recovered_rank_dfs,
    limma_tables = limma_tables,
    recovery_metrics = dplyr::bind_rows(recovery_metrics),
    scenario_metadata = score_simulation$scenario_metadata,
    score_simulation = score_simulation,
    export_gct_path = export_gct_path
  )
}

.simulate_hallmark_grouped_datasets_impl <- function(
    datasets = NULL,
    dataset_ids = NULL,
    geneset = NULL,
    seed = 20260317,
    species = "Homo sapiens",
    batch_ids = c("b1", "b2", "b3"),
    n_samples_per_group_batch = 3,
    baseline_range = c(4, 7),
    baseline_shape1 = 2,
    baseline_shape2 = 2,
    sigma_min = 0.10,
    sigma_max = 0.40,
    sigma_abundance_k = -1.2,
    sigma_abundance_anchor = 4,
    gamma_range = c(0.9, 1.1),
    lambda = 0.2,
    weight_shape = 2,
    weight_rate = 2,
    effect_scale = 2.6,
    effect_saturation = 2,
    gene_spike_scale = NULL,
    gene_spike_saturation = 1.25,
    batch_effect_sd = 0,
    export_dir = NULL,
    normalize_datasets_fn = .normalize_grouped_hallmark_datasets) {
  seed <- as.integer(seed[[1]])
  n_samples_per_group_batch <- as.integer(n_samples_per_group_batch[[1]])
  baseline_range <- as.numeric(baseline_range)
  baseline_shape1 <- as.numeric(baseline_shape1[[1]])
  baseline_shape2 <- as.numeric(baseline_shape2[[1]])
  sigma_min <- as.numeric(sigma_min[[1]])
  sigma_max <- as.numeric(sigma_max[[1]])
  sigma_abundance_k <- as.numeric(sigma_abundance_k[[1]])
  sigma_abundance_anchor <- as.numeric(sigma_abundance_anchor[[1]])
  gamma_range <- as.numeric(gamma_range)
  lambda <- as.numeric(lambda[[1]])
  weight_shape <- as.numeric(weight_shape[[1]])
  weight_rate <- as.numeric(weight_rate[[1]])
  effect_scale <- as.numeric(effect_scale[[1]])
  effect_saturation <- as.numeric(effect_saturation[[1]])
  if (is.null(gene_spike_scale)) {
    gene_spike_scale <- effect_scale * 1.35
  }
  gene_spike_scale <- as.numeric(gene_spike_scale[[1]])
  gene_spike_saturation <- as.numeric(gene_spike_saturation[[1]])
  batch_effect_sd <- as.numeric(batch_effect_sd[[1]])

  if (is.null(geneset)) {
    geneset <- geneset_tools$get_collection("H", "", species = species)
  }
  if (is.na(seed)) {
    stop("seed must be coercible to an integer")
  }
  if (is.na(n_samples_per_group_batch) || n_samples_per_group_batch < 1) {
    stop("n_samples_per_group_batch must be a positive integer")
  }
  if (length(batch_ids) < 1 || any(!nzchar(batch_ids))) {
    stop("batch_ids must contain at least one non-empty batch label")
  }
  if (length(baseline_range) != 2 || any(is.na(baseline_range)) || baseline_range[[1]] >= baseline_range[[2]]) {
    stop("baseline_range must be a numeric vector of length 2 with min < max")
  }
  if (is.na(baseline_shape1) || baseline_shape1 <= 0 || is.na(baseline_shape2) || baseline_shape2 <= 0) {
    stop("baseline_shape1 and baseline_shape2 must be positive numeric values")
  }
  if (is.na(sigma_min) || is.na(sigma_max) || sigma_min <= 0 || sigma_max < sigma_min) {
    stop("sigma_min and sigma_max must be positive values with sigma_max >= sigma_min")
  }
  if (is.na(sigma_abundance_k) || is.na(sigma_abundance_anchor)) {
    stop("sigma_abundance_k and sigma_abundance_anchor must be numeric")
  }
  if (length(gamma_range) != 2 || any(is.na(gamma_range)) || gamma_range[[1]] > gamma_range[[2]]) {
    stop("gamma_range must be a numeric vector of length 2 with min <= max")
  }
  if (is.na(lambda)) {
    stop("lambda must be numeric")
  }
  if (is.na(weight_shape) || weight_shape <= 0 || is.na(weight_rate) || weight_rate <= 0) {
    stop("weight_shape and weight_rate must be positive numeric values")
  }
  if (is.na(effect_scale) || effect_scale < 0) {
    stop("effect_scale must be a non-negative numeric value")
  }
  if (is.na(effect_saturation) || effect_saturation <= 0) {
    stop("effect_saturation must be a positive numeric value")
  }
  if (is.na(gene_spike_scale) || gene_spike_scale < 0) {
    stop("gene_spike_scale must be a non-negative numeric value")
  }
  if (is.na(gene_spike_saturation) || gene_spike_saturation <= 0) {
    stop("gene_spike_saturation must be a positive numeric value")
  }
  if (is.na(batch_effect_sd) || batch_effect_sd < 0) {
    stop("batch_effect_sd must be a non-negative numeric value")
  }
  if (!is.function(normalize_datasets_fn)) {
    stop("normalize_datasets_fn must be a function")
  }
  if (!"gs_name" %in% colnames(geneset) || !"ncbi_gene" %in% colnames(geneset)) {
    stop("geneset must contain at least gs_name and ncbi_gene columns")
  }

  hallmark_geneset <- geneset %>%
    dplyr::mutate(
      gs_name = as.character(gs_name),
      gene_id = as.character(ncbi_gene),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id))

  available_pathways <- sort(unique(hallmark_geneset$gs_name))
  datasets_norm <- normalize_datasets_fn(
    datasets = datasets,
    available_pathways = available_pathways
  )

  if (!is.null(dataset_ids)) {
    dataset_ids <- as.character(dataset_ids)
    datasets_norm <- datasets_norm[names(datasets_norm) %in% dataset_ids]
    if (length(datasets_norm) == 0) {
      stop("dataset_ids did not match any configured teaching datasets")
    }
  }

  dataset_activity_metadata <- .grouped_dataset_catalog_to_metadata(datasets_norm)
  dataset_group_metadata <- dataset_activity_metadata %>%
    dplyr::distinct(
      dataset_id,
      dataset_order,
      dataset_label,
      dataset_description,
      group,
      group_order,
      group_label,
      group_description,
      group_is_vehicle,
      treatment_level
    ) %>%
    dplyr::arrange(dataset_order, group_order)

  gene_lookup <- .prepare_hallmark_gene_lookup(hallmark_geneset)

  set.seed(seed + 100L)
  baseline_scaled <- stats::rbeta(
    n = nrow(gene_lookup),
    shape1 = baseline_shape1,
    shape2 = baseline_shape2
  )

  shared_feature_metadata <- gene_lookup %>%
    dplyr::mutate(
      baseline_log10 = baseline_range[[1]] + (diff(baseline_range) * baseline_scaled),
      sigma_gene = sigma_min + ((sigma_max - sigma_min) * exp(sigma_abundance_k * (baseline_log10 - sigma_abundance_anchor))),
      sigma_gene = pmin(pmax(sigma_gene, sigma_min), sigma_max)
    )
  n_batch_modules <- 8L
  set.seed(seed + 101L)
  shared_feature_metadata$batch_module <- sample.int(
    n = n_batch_modules,
    size = nrow(shared_feature_metadata),
    replace = TRUE
  )

  n_per_group <- length(batch_ids) * n_samples_per_group_batch
  contrast_se_scale <- sqrt((1 / n_per_group) + (1 / n_per_group))
  sample_biology_sd <- 0.12
  sample_biology_bounds <- c(0.75, 1.25)

  dataset_outputs <- list()
  group_gene_spike_metadata_all <- list()
  batch_ids <- as.character(batch_ids)

  for (dataset_entry in datasets_norm) {
    dataset_id <- dataset_entry$dataset_id
    dataset_seed <- seed + (dataset_entry$dataset_order * 1000L)
    dataset_groups <- dataset_entry$groups

    union_pathways <- unique(unlist(purrr::map(dataset_groups, ~ names(.x$activities)), use.names = FALSE))
    union_pathways <- union_pathways[nzchar(union_pathways)]
    has_gene_spikes <- any(purrr::map_int(dataset_groups, ~ length(.x$gene_spikes)) > 0)


    if (length(union_pathways) > 0) {
      active_memberships <- hallmark_geneset %>%
        dplyr::filter(gs_name %in% union_pathways) %>%
        dplyr::transmute(
          dataset_id = dataset_id,
          dataset_order = dataset_entry$dataset_order,
          dataset_label = dataset_entry$label,
          dataset_description = dataset_entry$description,
          pathway = gs_name,
          gene_id = gene_id,
          gene_symbol = gene_symbol
        ) %>%
        dplyr::distinct(pathway, gene_id, .keep_all = TRUE)

      set.seed(dataset_seed + 1L)
      active_memberships$raw_weight <- stats::rgamma(
        n = nrow(active_memberships),
        shape = weight_shape,
        rate = weight_rate
      )

      membership_metadata <- active_memberships %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(
          weight = raw_weight / sum(raw_weight),
          n_union_memberships = dplyr::n()
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(pathway, gene_id)
    } else {
      membership_metadata <- tibble::tibble(
        dataset_id = character(),
        dataset_order = integer(),
        dataset_label = character(),
        dataset_description = character(),
        pathway = character(),
        gene_id = character(),
        gene_symbol = character(),
        raw_weight = numeric(),
        weight = numeric(),
        n_union_memberships = integer()
      )
    }

    dataset_gene_spike_metadata <- .resolve_group_gene_spikes(
      dataset_entry = dataset_entry,
      gene_lookup = shared_feature_metadata %>% dplyr::select(gene_id, gene_symbol)
    )
    group_gene_spike_metadata_all[[dataset_id]] <- dataset_gene_spike_metadata

    set.seed(dataset_seed + 2L)
    gamma_lookup <- shared_feature_metadata %>%
      dplyr::select(gene_id, gene_symbol) %>%
      dplyr::mutate(
        gamma = stats::runif(
          n = dplyr::n(),
          min = gamma_range[[1]],
          max = gamma_range[[2]]
        )
      )

    gene_group_metadata <- purrr::map_dfr(dataset_groups, function(group_entry) {
      group_id <- group_entry$group
      activity_map <- group_entry$activities
      group_gene_spike_lookup <- dataset_gene_spike_metadata %>%
        dplyr::filter(group == group_id) %>%
        dplyr::select(gene_id, gene_spike_value, gene_spike_target_name, gene_spike_match_type)

      if (nrow(membership_metadata) > 0) {
        membership_group <- membership_metadata %>%
          dplyr::mutate(
            group = group_id,
            group_order = group_entry$group_order,
            group_label = group_entry$label,
            group_description = group_entry$description,
            group_is_vehicle = group_entry$is_vehicle,
            treatment_level = group_entry$treatment_level,
            activity = dplyr::coalesce(unname(activity_map[pathway]), 0),
            is_active_pathway = abs(activity) > .Machine$double.eps,
            activity_contribution = weight * activity
          )

        group_summary <- membership_group %>%
          dplyr::group_by(gene_id) %>%
          dplyr::summarise(
            n_union_memberships = dplyr::first(n_union_memberships),
            n_active_memberships = sum(is_active_pathway),
            weighted_activity = sum(activity_contribution),
            active_pathways = paste(pathway[is_active_pathway], collapse = " | "),
            .groups = "drop"
          )
      } else {
        group_summary <- tibble::tibble(
          gene_id = character(),
          n_union_memberships = integer(),
          n_active_memberships = integer(),
          weighted_activity = numeric(),
          active_pathways = character()
        )
      }

      shared_feature_metadata %>%
        dplyr::left_join(gamma_lookup, by = c("gene_id", "gene_symbol")) %>%
        dplyr::left_join(group_summary, by = "gene_id") %>%
        dplyr::left_join(group_gene_spike_lookup, by = "gene_id") %>%
        dplyr::mutate(
          dataset_id = dataset_id,
          dataset_order = dataset_entry$dataset_order,
          dataset_label = dataset_entry$label,
          dataset_description = dataset_entry$description,
          group = group_id,
          group_order = group_entry$group_order,
          group_label = group_entry$label,
          group_description = group_entry$description,
          group_is_vehicle = group_entry$is_vehicle,
          treatment_level = group_entry$treatment_level,
          n_union_memberships = dplyr::coalesce(n_union_memberships, 0L),
          n_active_memberships = dplyr::coalesce(n_active_memberships, 0L),
          weighted_activity = dplyr::coalesce(weighted_activity, 0),
          active_pathways = dplyr::coalesce(active_pathways, ""),
          gene_spike_value = dplyr::coalesce(gene_spike_value, 0),
          gene_spike_target_name = dplyr::coalesce(gene_spike_target_name, ""),
          gene_spike_match_type = dplyr::coalesce(gene_spike_match_type, ""),
          overlap_offset = dplyr::if_else(
            n_active_memberships > 1,
            lambda * (n_active_memberships - 1) * ((gamma - 1) / (1 + gamma)),
            0
          ),
          signal_mean = (gamma * weighted_activity) + overlap_offset,
          target_value_shrunk = tanh(signal_mean / effect_saturation),
          pathway_mean_shift = effect_scale * sigma_gene * contrast_se_scale * target_value_shrunk,
          gene_spike_shrunk = tanh(gene_spike_value / gene_spike_saturation),
          gene_spike_mean_shift = gene_spike_scale * sigma_gene * contrast_se_scale * gene_spike_shrunk,
          mean_shift = pathway_mean_shift + gene_spike_mean_shift,
          group_mean_log10 = baseline_log10 + mean_shift,
          is_active = n_active_memberships > 0 | abs(gene_spike_value) > .Machine$double.eps
        ) %>%
        dplyr::arrange(gene_id)
    })

    set.seed(dataset_seed + 3L)
    batch_offsets <- if (batch_effect_sd > 0) {
      raw_offsets <- stats::rnorm(length(batch_ids), mean = 0, sd = batch_effect_sd)
      raw_offsets - mean(raw_offsets)
    } else {
      rep(0, length(batch_ids))
    }
    names(batch_offsets) <- batch_ids
    set.seed(dataset_seed + 4L)
    batch_module_offsets <- if (batch_effect_sd > 0) {
      raw_module_offsets <- matrix(
        stats::rnorm(
          n = n_batch_modules * length(batch_ids),
          mean = 0,
          sd = batch_effect_sd * 0.75
        ),
        nrow = n_batch_modules,
        ncol = length(batch_ids),
        dimnames = list(
          paste0("module_", seq_len(n_batch_modules)),
          batch_ids
        )
      )
      raw_module_offsets - rowMeans(raw_module_offsets)
    } else {
      matrix(
        0,
        nrow = n_batch_modules,
        ncol = length(batch_ids),
        dimnames = list(
          paste0("module_", seq_len(n_batch_modules)),
          batch_ids
        )
      )
    }

    batch_metadata <- tibble::tibble(
      dataset_id = dataset_id,
      dataset_label = dataset_entry$label,
      batch = batch_ids,
      batch_order = seq_along(batch_ids),
      batch_offset = unname(batch_offsets),
      batch_module_offset_sd = apply(batch_module_offsets, 2, stats::sd)
    )

    sample_metadata <- purrr::map_dfr(dataset_groups, function(group_entry) {
      group_id <- group_entry$group
      purrr::map_dfr(seq_along(batch_ids), function(batch_idx) {
        batch_id <- batch_ids[[batch_idx]]
        tibble::tibble(
          id = sprintf("%s_%s_%s_rep%02d", dataset_id, group_id, batch_id, seq_len(n_samples_per_group_batch)),
          sample = sprintf("%s_%s_%s_rep%02d", dataset_id, group_id, batch_id, seq_len(n_samples_per_group_batch)),
          dataset_id = dataset_id,
          dataset_label = dataset_entry$label,
          group = group_id,
          group_order = group_entry$group_order,
          group_label = group_entry$label,
          group_description = group_entry$description,
          batch = batch_id,
          batch_order = batch_idx,
          replicate = seq_len(n_samples_per_group_batch),
          batch_offset = unname(batch_offsets[[batch_id]]),
          is_vehicle = group_entry$is_vehicle,
          treatment_level = group_entry$treatment_level
        )
      })
    }) %>%
      dplyr::arrange(group_order, batch_order, replicate) %>%
      as.data.frame(stringsAsFactors = FALSE)
    rownames(sample_metadata) <- sample_metadata$id

    set.seed(dataset_seed + 5L)
    sample_metadata$sample_biology_scale <- stats::rnorm(
      nrow(sample_metadata),
      mean = 1,
      sd = sample_biology_sd
    )
    sample_metadata$sample_biology_scale <- pmin(
      pmax(sample_metadata$sample_biology_scale, sample_biology_bounds[[1]]),
      sample_biology_bounds[[2]]
    )
    sample_metadata$sample_mean_shift_mean <- NA_real_
    sample_metadata$sample_latent_mean_log10_mean <- NA_real_
    sample_metadata$sample_sigma_gene_mean <- mean(shared_feature_metadata$sigma_gene)
    sample_metadata$sample_batch_effect_mean <- NA_real_

    expr_mat <- matrix(
      NA_real_,
      nrow = nrow(shared_feature_metadata),
      ncol = nrow(sample_metadata),
      dimnames = list(shared_feature_metadata$gene_id, sample_metadata$id)
    )

    for (sample_idx in seq_len(nrow(sample_metadata))) {
      sample_row <- sample_metadata[sample_idx, , drop = FALSE]
      gene_group_row <- gene_group_metadata %>%
        dplyr::filter(group == sample_row$group) %>%
        dplyr::arrange(match(gene_id, shared_feature_metadata$gene_id))

      scaled_mean_shift <- gene_group_row$mean_shift * sample_row$sample_biology_scale[[1]]
      gene_batch_offset <- sample_row$batch_offset[[1]] +
        batch_module_offsets[cbind(shared_feature_metadata$batch_module, match(sample_row$batch[[1]], batch_ids))]
      sample_mean_log10 <- shared_feature_metadata$baseline_log10 +
        scaled_mean_shift +
        gene_batch_offset

      sample_metadata$sample_mean_shift_mean[[sample_idx]] <- mean(scaled_mean_shift)
      sample_metadata$sample_batch_effect_mean[[sample_idx]] <- mean(gene_batch_offset)
      sample_metadata$sample_latent_mean_log10_mean[[sample_idx]] <- mean(sample_mean_log10)

      set.seed(dataset_seed + 100L + sample_idx)
      expr_mat[, sample_idx] <- stats::rnorm(
        n = nrow(shared_feature_metadata),
        mean = sample_mean_log10,
        sd = shared_feature_metadata$sigma_gene
      )
    }

    rdesc <- shared_feature_metadata %>%
      dplyr::transmute(
        id = gene_id,
        gene_symbol = gene_symbol,
        baseline_log10 = baseline_log10,
        sigma_gene = sigma_gene
      ) %>%
      as.data.frame(stringsAsFactors = FALSE)
    rownames(rdesc) <- rdesc$id

    cdesc <- sample_metadata
    expression_gct <- cmapR::GCT(
      mat = expr_mat,
      rid = rownames(expr_mat),
      cid = colnames(expr_mat),
      rdesc = rdesc,
      cdesc = cdesc
    )

    export_gct_path <- NULL
    if (!is.null(export_dir)) {
      export_dir <- fs::path_abs(export_dir)
      fs::dir_create(export_dir, recurse = TRUE)
      export_gct_path <- fs::path(export_dir, paste0(dataset_id, ".gct"))
      cmapR::write_gct(expression_gct, export_gct_path, appenddim = FALSE)
    }

    dataset_outputs[[dataset_id]] <- list(
      dataset_id = dataset_id,
      dataset_label = dataset_entry$label,
      dataset_description = dataset_entry$description,
      expression_matrix = expr_mat,
      expression_gct = expression_gct,
      sample_metadata = sample_metadata,
      feature_metadata = shared_feature_metadata,
      group_metadata = dataset_group_metadata %>%
        dplyr::filter(dataset_id == !!dataset_id),
      gene_group_metadata = gene_group_metadata,
      membership_metadata = membership_metadata,
      group_activity_metadata = dataset_activity_metadata %>%
        dplyr::filter(dataset_id == !!dataset_id),
      group_gene_spike_metadata = dataset_gene_spike_metadata,
      batch_metadata = batch_metadata,
      export_gct_path = export_gct_path
    )
  }

  group_gene_spike_metadata <- if (length(group_gene_spike_metadata_all) > 0) {
    dplyr::bind_rows(group_gene_spike_metadata_all)
  } else {
    .empty_group_gene_spike_metadata()
  }

  list(
    datasets = dataset_outputs,
    feature_metadata = shared_feature_metadata,
    group_metadata = dataset_group_metadata,
    group_activity_metadata = dataset_activity_metadata,
    group_gene_spike_metadata = group_gene_spike_metadata
  )
}

simulate_hallmark_grouped_datasets <- function(
    datasets = NULL,
    dataset_ids = NULL,
    geneset = NULL,
    seed = 20260317,
    species = "Homo sapiens",
    batch_ids = c("b1", "b2", "b3"),
    n_samples_per_group_batch = 3,
    baseline_range = c(4, 7),
    baseline_shape1 = 2,
    baseline_shape2 = 2,
    sigma_min = 0.10,
    sigma_max = 0.40,
    sigma_abundance_k = -1.2,
    sigma_abundance_anchor = 4,
    gamma_range = c(0.9, 1.1),
    lambda = 0.2,
    weight_shape = 2,
    weight_rate = 2,
    effect_scale = 2.6,
    effect_saturation = 2,
    gene_spike_scale = NULL,
    gene_spike_saturation = 1.25,
    batch_effect_sd = 0,
    export_dir = NULL) {
  .simulate_hallmark_grouped_datasets_impl(
    datasets = datasets,
    dataset_ids = dataset_ids,
    geneset = geneset,
    seed = seed,
    species = species,
    batch_ids = batch_ids,
    n_samples_per_group_batch = n_samples_per_group_batch,
    baseline_range = baseline_range,
    baseline_shape1 = baseline_shape1,
    baseline_shape2 = baseline_shape2,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    sigma_abundance_k = sigma_abundance_k,
    sigma_abundance_anchor = sigma_abundance_anchor,
    gamma_range = gamma_range,
    lambda = lambda,
    weight_shape = weight_shape,
    weight_rate = weight_rate,
    effect_scale = effect_scale,
    effect_saturation = effect_saturation,
    gene_spike_scale = gene_spike_scale,
    gene_spike_saturation = gene_spike_saturation,
    batch_effect_sd = batch_effect_sd,
    export_dir = export_dir,
    normalize_datasets_fn = .normalize_grouped_hallmark_datasets
  )
}

simulate_hallmark_teaching_datasets <- function(
    datasets = NULL,
    dataset_ids = NULL,
    geneset = NULL,
    seed = 20260317,
    species = "Homo sapiens",
    batch_ids = c("b1", "b2", "b3"),
    n_samples_per_group_batch = 3,
    baseline_range = c(4, 7),
    baseline_shape1 = 2,
    baseline_shape2 = 2,
    sigma_min = 0.10,
    sigma_max = 0.40,
    sigma_abundance_k = -1.2,
    sigma_abundance_anchor = 4,
    gamma_range = c(0.9, 1.1),
    lambda = 0.2,
    weight_shape = 2,
    weight_rate = 2,
    effect_scale = 2.6,
    effect_saturation = 2,
    gene_spike_scale = NULL,
    gene_spike_saturation = 1.25,
    batch_effect_sd = 0,
    export_dir = NULL) {
  .simulate_hallmark_grouped_datasets_impl(
    datasets = datasets,
    dataset_ids = dataset_ids,
    geneset = geneset,
    seed = seed,
    species = species,
    batch_ids = batch_ids,
    n_samples_per_group_batch = n_samples_per_group_batch,
    baseline_range = baseline_range,
    baseline_shape1 = baseline_shape1,
    baseline_shape2 = baseline_shape2,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    sigma_abundance_k = sigma_abundance_k,
    sigma_abundance_anchor = sigma_abundance_anchor,
    gamma_range = gamma_range,
    lambda = lambda,
    weight_shape = weight_shape,
    weight_rate = weight_rate,
    effect_scale = effect_scale,
    effect_saturation = effect_saturation,
    gene_spike_scale = gene_spike_scale,
    gene_spike_saturation = gene_spike_saturation,
    batch_effect_sd = batch_effect_sd,
    export_dir = export_dir,
    normalize_datasets_fn = .normalize_teaching_datasets
  )
}


#' Generate Test Data for GSEA Analysis
#'
#' This function generates test data for Gene Set Enrichment Analysis (GSEA) by
#' simulating preranked data, selecting gene sets of interest, and running
#' fgsea analysis. It relies on external helper functions from `fgsea_tools`
#' and `io_tools` for data simulation and manipulation.
#'
#' @param collapse Logical. If `TRUE`, collapses redundant pathways. Default is `FALSE`.
#' @param pathways Character vector. Specifies the pathways of interest. Default
#'   is `c("H", "GO:BP")`. Supported values are:
#'   - "H" for Hallmark gene sets
#'   - "GO:BP" for Gene Ontology Biological Process
#'   - "GO:CC" for Gene Ontology Cellular Component
#'   - "GO:MF" for Gene Ontology Molecular Function
#'
#' @details
#' This function simulates two sets of preranked data using `simulate_preranked_data()`
#' and combines them into a list. It then uses `fgsea_tools` to run GSEA across
#' all specified pathways. The gene sets of interest are filtered and selected based
#' on the provided `pathways` argument, and the data is processed using `io_tools`
#' to convert ranks to the required format.
#'
#' The function sources external dependencies for `fgsea_tools` and `io_tools`, and
#' the gene set collections are fetched using `geneset_tools$get_collections()`.
#'
#' @note
#' This function depends on the `test_fgsea.R` file for proper testing and
#' success, though this is not strictly guaranteed.
#'
#' @return A list of GSEA results for all pathways.
#'
#' @example
#' # Example usage:
#' res <- generate_test_data(collapse = TRUE, pathways = c("H", "GO:BP"))
#'
#' @seealso `simulate_preranked_data`, `fgsea_tools`, `io_tools`
#'
generate_test_data <- function(collapse = FALSE,
                               pathways = c("H", "GO:BP"),
                               preranked_data = NULL,
                               cache = FALSE,
                               parallel = FALSE,
                               minSize = 15,
                               biocparallel_param = NULL) {
  # this function relies on simulate preranked data from fgsea tools,
  # which depends on test_fgsea.R for guaranteed* success
  # *not guaranteed

  # genesets_of_interest <- list(
  #     list(category = "H", subcategory = ""),
  #     list(category = "C5", subcategory = "GO:BP")
  # )


  # put the imports here to save time on startup
  io_tools <- new.env()
  source(file.path(here("R"), "./io.R"), local = io_tools)

  fgsea_tools <- new.env()
  source(file.path(here("R"), "./fgsea.R"), local = fgsea_tools)

  genesets_list <- list()
  if ("H" %in% pathways) genesets_list[["H"]] <- c("H", "", FALSE, "H_")
  if ("GO:BP" %in% pathways) genesets_list[["GO:BP"]] <- c("C5", "GO:BP", TRUE, "C5_GO:BP")
  if ("GO:CC" %in% pathways) genesets_list[["GO:CC"]] <- c("C5", "GO:CC", TRUE, "C5_GO:CC")
  if ("GO:MF" %in% pathways) genesets_list[["GO:MF"]] <- c("C5", "GO:MF", TRUE, "C5_GO:MF")

  # Convert the list to a tibble
  genesets_of_interest <- purrr::map_dfr(genesets_list, ~ tibble::tibble(
    category = .x[1], subcategory = .x[2], collapse = .x[3], collection_name = .x[4]
  ))

  genesets <- geneset_tools$get_collections(
    genesets_of_interest
  )

  # TODO
  geneset_list <- genesets %>% purrr::imap(~ geneset_tools$genesets_df_to_list(.x))

  # named_data = list(data=data)
  if (is.null(preranked_data)) {
    data1 <- simulate_preranked_data(seed = 1234, sample_frac = .4)
    data2 <- simulate_preranked_data(seed = 4321, sample_frac = .4)
    data <- list(first = data1, second = data2)
  } else {
    data <- preranked_data
  }
  rankobjs <- io_tools$ranks_dfs_to_lists(data)




  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    old_param <- tryCatch(BiocParallel::bpparam(), error = function(e) NULL)
    if (!is.null(old_param)) {
      on.exit(BiocParallel::register(old_param, default = TRUE), add = TRUE)
    }
    param_to_use <- if (is.null(biocparallel_param)) {
      BiocParallel::SerialParam()
    } else {
      biocparallel_param
    }
    BiocParallel::register(param_to_use, default = TRUE)
  }

  res <- fgsea_tools$run_all_pathways(
    geneset_list,
    rankobjs,
    collapse = collapse,
    cache = cache,
    parallel = parallel,
    minSize = minSize
  )

  if (isTRUE(collapse)) {
    res <- res %>% purrr::map(~ {
      purrr::map(.x, function(df) {
        if (!is.null(df) && "mainpathway" %in% colnames(df)) {
          unique_vals <- unique(df$mainpathway)
          if (length(unique_vals) == 1 && nrow(df) > 1) {
            df$mainpathway[nrow(df)] <- !unique_vals[[1]]
          }
        }
        df
      })
    })
  }
  return(res)
}
