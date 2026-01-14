suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(here))

source(file.path(here("R"), "lazyloader.R"))

util_tools <- get_tool_env("utils")
plot_utils <- get_tool_env("plot_utils")
fgsea_tools <- get_tool_env("fgsea")
plot_tools <- get_tool_env("plot")

make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)

# Default stopword list for pathway names. These are generic words that
# typically do not add biological specificity when summarising via wordcloud.
DEFAULT_PATHWAY_STOPWORDS <- c(
  # Generic ontology and collection markers
  "go", "reactome", "hallmark", "kegg", "mm", "bp", "mf", "cc",
  # Direction / regulation boilerplate
  "regulation", "regulations", "regulated", "regulating",
  "positive", "negative",
  "activation", "inhibition",
  "dependent", "independent",
  "mediated", "mediating",
  "involved",
  # Very generic biological nouns
  "process", "processes",
  "pathway", "pathways",
  "response", "responses",
  "signaling", "signal", "signalling",
  "activity", "activities",
  "binding",
  "cell", "cells", "cellular",
  "system", "systems",
  "function", "functions",
  # Glue words
  "and", "or", "to", "of", "the", "by", "for", "in", "on", "with", "via"
)

# Prepare a tidy table of word-level statistics from a table of enriched pathways.
# This is primarily geared towards fgsea-style result tables that contain at least
# columns "pathway", "NES" and "padj" (optionally "pval").
prepare_terms_for_wordcloud <- function(
    input_df,
    padj_cutoff = 0.25,
    top_n_pathways = 50,
    min_chars = 3,
    extra_stopwords = NULL) {
  if (is.null(input_df) || nrow(input_df) == 0) {
    log_msg(info = "wordcloud: input dataframe is empty; nothing to summarise")
    return(dplyr::tibble())
  }

  if (!"pathway" %in% colnames(input_df)) {
    stop("prepare_terms_for_wordcloud: expected a 'pathway' column in input_df")
  }

  working_df <- input_df

  # Ensure we have sensible numeric padj and pval columns.
  if (!"padj" %in% colnames(working_df) && "pval" %in% colnames(working_df)) {
    working_df <- working_df %>% mutate(padj = pval)
  }

  if (!"padj" %in% colnames(working_df)) {
    stop("prepare_terms_for_wordcloud: input_df must contain 'padj' or 'pval'")
  }

  working_df <- working_df %>%
    mutate(
      padj = suppressWarnings(as.numeric(padj)),
      pval = if ("pval" %in% colnames(.)) suppressWarnings(as.numeric(pval)) else padj
    )

  working_df <- working_df %>%
    filter(!is.na(padj) & is.finite(padj))

  if (!is.null(padj_cutoff)) {
    working_df <- working_df %>% filter(padj <= padj_cutoff)
  }

  if (nrow(working_df) == 0) {
    log_msg(info = paste0(
      "wordcloud: no pathways pass padj cutoff (", padj_cutoff, "); skipping"
    ))
    return(dplyr::tibble())
  }

  # Order by significance / effect and keep the most relevant pathways.
  if (!"NES" %in% colnames(working_df)) {
    stop("prepare_terms_for_wordcloud: input_df must contain an 'NES' column")
  }

  working_df <- working_df %>%
    mutate(
      NES = suppressWarnings(as.numeric(NES)),
      pval = ifelse(is.na(pval) | pval <= 0, padj, pval),
      score = abs(NES) * -log10(pval),
      score = ifelse(is.na(score) | !is.finite(score), 0, score),
      sign_dir = sign(NES)
    ) %>%
    arrange(desc(score))

  if (!is.null(top_n_pathways) && top_n_pathways > 0) {
    n_to_keep <- min(nrow(working_df), top_n_pathways)
    if (n_to_keep > 0) {
      working_df <- head(working_df, n_to_keep)
    }
  }

  if (nrow(working_df) == 0) {
    log_msg(info = "wordcloud: no pathways remain after top_n_pathways filtering; skipping")
    return(dplyr::tibble())
  }

  stopwords <- union(DEFAULT_PATHWAY_STOPWORDS, extra_stopwords)

  token_df <- working_df %>%
    transmute(pathway, score, sign_dir) %>%
    mutate(
      pathway_clean = pathway %>%
        as.character() %>%
        stringr::str_replace_all("[[:punct:]]", " ") %>%
        stringr::str_replace_all("\\s+", " ") %>%
        stringr::str_squish()
    ) %>%
    tidyr::separate_rows(pathway_clean, sep = "\\s+") %>%
    dplyr::rename(word = pathway_clean) %>%
    mutate(
      word = stringr::str_to_lower(word),
      word = stringr::str_replace_all(word, "[^[:alnum:]]", ""),
      word = stringr::str_squish(word)
    ) %>%
    filter(
      nzchar(word),
      nchar(word) >= min_chars,
      !word %in% stopwords
    )

  if (nrow(token_df) == 0) {
    log_msg(info = "wordcloud: no tokens remain after cleaning and stopword removal; skipping")
    return(dplyr::tibble())
  }

  terms_df <- token_df %>%
    group_by(word) %>%
    summarise(
      weight = sum(score, na.rm = TRUE),
      mean_sign = mean(sign_dir[is.finite(sign_dir)], na.rm = TRUE),
      n_pathways = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      direction = dplyr::case_when(
        mean_sign > 0.15 ~ "up",
        mean_sign < -0.15 ~ "down",
        TRUE ~ "mixed"
      )
    ) %>%
    arrange(desc(weight))

  terms_df
}

# High-level helper to plot a wordcloud from a pathway table.
# This prefers ggwordcloud when available, but falls back to a simple grid layout
# built on geom_text if the package is not installed.
wordcloud_plot <- function(
    input_df,
    title = "",
    subtitle = NULL,
    padj_cutoff = 0.25,
    top_n_pathways = 50,
    max_words = 70,
    min_chars = 3,
    extra_stopwords = NULL,
    save_func = NULL) {
  terms_df <- prepare_terms_for_wordcloud(
    input_df = input_df,
    padj_cutoff = padj_cutoff,
    top_n_pathways = top_n_pathways,
    min_chars = min_chars,
    extra_stopwords = extra_stopwords
  )

  if (nrow(terms_df) == 0) {
    return(NULL)
  }

  if (!is.null(max_words) && max_words > 0) {
    terms_df <- terms_df %>%
      arrange(desc(weight))
    n_terms_keep <- min(nrow(terms_df), max_words)
    if (n_terms_keep > 0) {
      terms_df <- head(terms_df, n_terms_keep)
    }
  }

  formatted_title <- title %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_wrap(width = 54)

  formatted_subtitle <- if (is.null(subtitle)) {
    NULL
  } else {
    subtitle %>% stringr::str_replace_all("_", " ") %>% stringr::str_wrap(width = 72)
  }

  # Compact caption describing encoding and filters so users can
  # see which thresholds/top-N settings were used for this plot.
  filter_bits <- character(0)
  if (!is.null(padj_cutoff)) {
    filter_bits <- c(filter_bits, sprintf("padj ≤ %.3g", padj_cutoff))
  }
  if (!is.null(top_n_pathways) && top_n_pathways > 0) {
    filter_bits <- c(filter_bits, sprintf("top %d pathways", top_n_pathways))
  }
  if (!is.null(max_words) && max_words > 0) {
    filter_bits <- c(filter_bits, sprintf("max %d words", max_words))
  }
  filter_summary <- if (length(filter_bits) > 0) {
    paste(filter_bits, collapse = " • ")
  } else {
    ""
  }
  caption_text <- paste0(
    "Size ∝ Σ|NES|·-log10(p); colour: red=up, blue=down, grey=mixed",
    if (nzchar(filter_summary)) paste0(" • ", filter_summary) else ""
  )

  use_ggwordcloud <- requireNamespace("ggwordcloud", quietly = TRUE)

  if (use_ggwordcloud) {
    wordcloud_plot_object <- ggplot(terms_df, aes(label = word, size = weight, colour = direction)) +
      ggwordcloud::geom_text_wordcloud(area_corr_power = 1) +
      scale_size_area(max_size = 20) +
      scale_colour_manual(
        values = c(
          up = "#b2182b",
          down = "#2166ac",
          mixed = "#4d4d4d"
        )
      ) +
      labs(
        title = formatted_title,
        subtitle = formatted_subtitle,
        caption = caption_text
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.caption = element_text(size = 6.5, colour = "#555555", hjust = 1)
      )
  } else {
    # Simple deterministic grid layout fallback when ggwordcloud is unavailable.
    n_terms <- nrow(terms_df)
    ncol_grid <- ceiling(sqrt(n_terms))
    nrow_grid <- ceiling(n_terms / ncol_grid)

    terms_df <- terms_df %>%
      mutate(
        term_index = dplyr::row_number() - 1L,
        col_index = term_index %% ncol_grid,
        row_index = term_index %/% ncol_grid
      )

    wordcloud_plot_object <- ggplot(
      terms_df,
      aes(
        x = col_index,
        y = -row_index,
        label = word,
        size = weight,
        colour = direction
      )
    ) +
      geom_text() +
      scale_size_area(max_size = 7) +
      scale_colour_manual(
        values = c(
          up = "#b2182b",
          down = "#2166ac",
          mixed = "#4d4d4d"
        )
      ) +
      labs(
        title = formatted_title,
        subtitle = formatted_subtitle,
        caption = caption_text,
        x = NULL,
        y = NULL
      ) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.caption = element_text(size = 6.5, colour = "#555555", hjust = 1)
      )
  }

  if (!is.null(save_func)) {
    save_result <- save_func(
      plot_code = function() {
        print(wordcloud_plot_object)
      },
      width = 4.5,
      height = 3.5
    )
    if (is.null(save_result)) {
      log_msg(msg = paste0(
        "wordcloud_plot: skipped saving '", title, "' (existing file)"
      ))
    } else {
      log_msg(msg = paste0(
        "wordcloud_plot: saved plot '", title, "'"
      ))
    }
  }

  wordcloud_plot_object
}
