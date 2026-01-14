suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(scales))

source(file.path(here("R"), "lazyloader.R"))

util_tools <- get_tool_env("utils")
plot_utils <- get_tool_env("plot_utils")
fgsea_tools <- get_tool_env("fgsea" )
plot_tools <- get_tool_env("plot" )

make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)

DEFAULT_BUBBLE_GLYPH <- "⁕"

# Format a readable source/collection label for subtitles
format_source_label <- function(name) {
  name %>%
    stringr::str_replace_all("[:._]", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_squish()
}

prepare_data_for_bubble <- function(df, glyph = DEFAULT_BUBBLE_GLYPH) {
  if (!"pval" %in% colnames(df)) {
    df <- df %>% mutate(pval = padj)
  }

  plot_tools$prepare_data_for_barplot(df) %>%
    mutate(
      padj = {
        tmp <- padj
        if (is.list(tmp)) {
          tmp <- vapply(tmp, function(x) as.numeric(x)[1], numeric(1), USE.NAMES = FALSE)
        }
        tmp <- suppressWarnings(as.numeric(tmp))
        ifelse(is.na(tmp), 1, tmp)
      },
      pval = {
        tmp <- pval
        if (is.list(tmp)) {
          tmp <- vapply(tmp, function(x) as.numeric(x)[1], numeric(1), USE.NAMES = FALSE)
        }
        tmp <- suppressWarnings(as.numeric(tmp))
        ifelse(is.na(tmp), 1, tmp)
      },
      plot_leading_edge = pmax(leadingEdgeNum, 1),
      sig_category = dplyr::case_when(
        padj < 0.05 ~ "<0.05",
        padj < 0.25 ~ "<0.25",
        TRUE ~ ">=0.25"
      ),
      fill_value = sign(NES) * (1 - pmin(pval, 1)),
      sig_label = ifelse(sig_category == "<0.05", glyph, ""),
      text_color = ifelse(abs(fill_value) > 0.55, "#FFFFFF", "#111111")
    )
}

bubble_plot <- function(
    df,
    title = "",
    subtitle = NULL,
    save_func = NULL,
    facet_order = NULL,
    nes_range = NULL,
    size_range = c(3.0, 9.0),
    glyph = DEFAULT_BUBBLE_GLYPH,
    ...) {

  sel <- prepare_data_for_bubble(df, glyph = glyph)
  log_msg(msg = paste0("bubble_plot: received ", nrow(df), " rows, plotting ", nrow(sel), " after prep"))

  if (nrow(sel) == 0) {
    log_msg(warning = paste0("bubble_plot: no rows to plot for title '", title, "'"))
    return(NULL)
  }

  formatted_title <- title %>%
    str_replace_all("_", " ") %>%
    str_wrap(width = 54)

  formatted_subtitle <- if (is.null(subtitle)) {
    NULL
  } else {
    subtitle %>% str_replace_all("_", " ") %>% str_wrap(width = 72)
  }

  custom_labeller <- function(value) {
    vapply(value, function(label) {
      label %>%
        str_replace_all("_", " ") %>%
        str_wrap(width = 54)
    }, character(1))
  }

  if (!is.null(facet_order) && "rankname" %in% colnames(sel)) {
    sel <- sel %>%
      mutate(rankname = factor(rankname, levels = facet_order, ordered = TRUE)) %>%
      arrange(rankname)
  }

  if ("rankname" %in% names(sel)) {
    # right here we make a mapping and rename for display
    # for facetted plots
    name_map <- util_tools$make_name_map(unique(sel$rankname))
    sel <- sel %>% dplyr::mutate(
                          rankname=dplyr::recode(
                                               rankname,
                                               !!!name_map  # named vector
                                               )
                          )
    }



  sel <- sel %>% mutate(plot_leading_edge = pmax(1, plot_leading_edge))

  sel_text <- sel %>% filter(sig_category == "<0.05")
  sel_text_dark <- sel_text %>% filter(text_color == "#FFFFFF")
  sel_text_light <- sel_text %>% filter(text_color != "#FFFFFF")

  # Dynamically adjust y-axis text size based on label length
  # Remove line breaks when measuring to approximate visual width
  pathway_chars <- sel$pathway %>% as.character() %>%
    stringr::str_replace_all("\n", " ") %>%
    stringr::str_squish()
  max_label_chars <- if (length(pathway_chars) > 0) max(nchar(pathway_chars)) else 0
  # Align size scaling with heatmap/barplot conventions
  axis_text_y_size <- dplyr::case_when(
    max_label_chars < 36 ~ 7.6,
    max_label_chars < 64 ~ 6.6,
    max_label_chars < 84 ~ 6.2,
    TRUE ~ 5.2
  )
  strip_label_chars <- if ("rankname" %in% names(sel)) {
    sel$rankname %>% as.character() %>%
      stringr::str_replace_all("\n", " ") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_squish()
  } else {
    character(0)
  }
  max_strip_chars <- if (length(strip_label_chars) > 0) max(nchar(strip_label_chars)) else 0
  strip_text_size <- dplyr::case_when(
    max_strip_chars < 28 ~ 10,
    max_strip_chars < 48 ~ 9,
    max_strip_chars < 72 ~ 8,
    TRUE ~ 7
  )

  p <- ggplot(sel, aes(x = NES, y = pathway)) +
    # Reference line at x=0, layered behind points
    geom_vline(xintercept = 0, colour = scales::alpha("#555555", 0.6), size = 0.4, show.legend = FALSE) +
    geom_point(
      aes(
        size = plot_leading_edge,
        fill = fill_value,
        colour = sig_category
      ),
      shape = 21,
      stroke = 0.8,
      alpha = 0.9
    ) +
    scale_fill_gradient2(
      low = "#084594",
      mid = "#ffffff",
      high = "#b30000",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "#f7f7f7",
      guide = guide_colourbar(title = "1 - pval")
    ) +
    scale_size(
      range = size_range,
      guide = guide_legend(title = "Leading edge genes")
    ) +
    scale_colour_manual(
      values = c(
        "<0.25" = scales::alpha("#3f3f3f", 0.45),
        "<0.05" = scales::alpha("#000000", 0.55),
        ">=0.25" = scales::alpha("#000000", 0)
      ),
      breaks = c("<0.25", "<0.05"),
      labels = c("padj < 0.25", paste0("padj < 0.05 (", glyph, ")")),
      guide = guide_legend(
        title = "Significance",
        override.aes = list(
          shape = 21,
          fill = "grey70",
          size = 4.5,
          stroke = 0.8,
          alpha = 1,
          colour = scales::alpha("#000000", 0.55)
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    scale_y_discrete(position = "right") +
    labs(
      title = formatted_title,
      subtitle = formatted_subtitle,
      x = "NES",
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = axis_text_y_size, face = "bold"),
      axis.text.x = element_text(size = 7.0),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      strip.text = element_text(size = strip_text_size, face = "bold", hjust = 0.5),
      legend.position = "right"
    )

  if (nrow(sel_text_dark) > 0) {
    p <- p + geom_text(
      data = sel_text_dark,
      inherit.aes = FALSE,
      aes(x = NES, y = pathway, label = sig_label),
      colour = "#FFFFFF",
      size = 3.0,
      vjust = 0.5,
      show.legend = FALSE
    )
  }

  if (nrow(sel_text_light) > 0) {
    p <- p + geom_text(
      data = sel_text_light,
      inherit.aes = FALSE,
      aes(x = NES, y = pathway, label = sig_label),
      colour = "#111111",
      size = 3.0,
      vjust = 0.5,
      show.legend = FALSE
    )
  }

  if (is.null(nes_range)) {
    nes_max <- suppressWarnings(max(sel$NES, na.rm = TRUE))
    nes_min <- suppressWarnings(min(sel$NES, na.rm = TRUE))
    if (!is.finite(nes_max) || !is.finite(nes_min)) {
      nes_range <- NULL
    } else if (nes_min < 0 && nes_max > 0) {
      max_abs <- max(abs(nes_min), abs(nes_max))
      nes_range <- c(-max_abs, max_abs)
    } else {
      nes_range <- c(min(nes_min, 0, na.rm = TRUE), max(nes_max, 0, na.rm = TRUE))
    }
  } else {
    if (length(nes_range) != 2) {
      stop("nes_range should be a numeric vector of length 2")
    }
  }

  if (!is.null(nes_range) && all(is.finite(nes_range))) {
    p <- p + coord_cartesian(xlim = c(nes_range[1], nes_range[2]))
  }

  if ("rankname" %in% colnames(sel) && length(unique(sel$rankname)) > 1) {
    p <- p + facet_wrap(~rankname, labeller = as_labeller(custom_labeller))
    num_panels <- length(unique(sel$rankname))
    ncol <- ceiling(sqrt(num_panels))
    nrow <- ceiling(num_panels / ncol)
  } else {
    ncol <- 1
    nrow <- 1
  }

  panel_width_in <- 4.0
  # Compute height so each pathway row has reasonable space
  n_pathways <- dplyr::n_distinct(sel$pathway)
  text_scale <- axis_text_y_size / 6.6
  per_row_in <- dplyr::case_when(
    n_pathways <= 20 ~ 0.28,
    n_pathways <= 60 ~ 0.25,
    TRUE ~ 0.26
  ) * text_scale
  per_row_in <- max(min(per_row_in, 0.34), 0.18)
  min_panel_height_in <- 3.2
  facet_strip_pad_in <- 0.35
  panel_height_in_calc <- max(min_panel_height_in, per_row_in * n_pathways + facet_strip_pad_in)
  total_width_in <- 2.2 + (panel_width_in * ncol)
  total_height_in <- panel_height_in_calc * nrow

  if (!is.null(save_func)) {
    save_result <- save_func(
      plot_code = function() {
        print(p)
      },
      width = total_width_in,
      height = total_height_in
    )
    if (is.null(save_result)) {
      log_msg(msg = paste0(
        "bubble_plot: skipped saving '", title, "' (existing file)"
      ))
    } else {
      log_msg(msg = paste0(
        "bubble_plot: saved plot '", title, "' with dimensions ",
        round(total_width_in, 2), "x", round(total_height_in, 2)
      ))
    }
  }

  p
}

all_bubble_plots <- function(
    results_list,
    save_func = NULL,
    facet_order = NULL,
    limit = 20,
    size_range = c(3.0, 9.0),
    glyph = DEFAULT_BUBBLE_GLYPH,
    ...) {
  if (!is.null(save_func)) {
    existing_filename <- get_arg(save_func, "filename")
    if (!nzchar(existing_filename)) {
      base_filename <- util_tools$safe_filename("bubble", fallback = "bubble")
      save_func <- make_partial(save_func, filename = base_filename)
    }
  }

  results_list %>% purrr::imap(
    ~ {
      collection_name <- .y
      collection_label <- format_source_label(collection_name)
      list_of_comparisons <- .x
      # Build a friendly mapping to strip shared affixes from comparison names
      comparison_names <- names(list_of_comparisons)
      name_map <- util_tools$make_name_map(comparison_names)
      list_of_comparisons %>% purrr::imap(
        ~ {
          dataframe <- .x
          comparison_name <- .y
          comparison_label <- name_map[[comparison_name]] %||% comparison_name

          purrr::map(limit, function(.limit) {
            sel <- fgsea_tools$select_topn(dataframe, limit = .limit, pstat_cutoff = 1)
            log_msg(msg = paste0(
              "bubble all: collection=", collection_name,
              " comparison=", comparison_name,
              " limit=", .limit,
              " selected_rows=", nrow(sel)
            ))
            if (nrow(sel) == 0) {
              log_msg(warning = paste0(
                "no pathways available for bubble plot (requested top ", .limit, ") in ",
                collection_name, " / ", comparison_name
              ))
              return(NULL)
            }
            n_pathways <- dplyr::n_distinct(sel$pathway)
            effective_limit <- min(.limit, n_pathways)
            nes_max <- suppressWarnings(max(abs(dataframe$NES), na.rm = TRUE))
            nes_range <- if (is.finite(nes_max)) c(-nes_max, nes_max) else NULL

            local_save_func <- save_func
            save_path <- NULL
            if (!is.null(local_save_func)) {
              collection_dir <- util_tools$safe_path_component(collection_name)
              comparison_dir <- util_tools$safe_path_component(comparison_label)
              filename <- util_tools$safe_filename(
                "bubble",
                paste0("top", effective_limit),
                paste0("n", n_pathways),
                fallback = "bubble"
              )
              save_path <- util_tools$safe_subdir(
                get_arg(local_save_func, "path"),
                collection_dir,
                "bubble",
                comparison_dir
              )
              local_save_func <- make_partial(local_save_func, filename = filename, path = save_path)
              log_msg(msg = paste0(
                "bubble all: saving to ", file.path(save_path, paste0(filename, ".pdf"))
              ))
            }

            rank_label <- sel$rankname %>% na.omit() %>% unique()
            rank_label <- if (length(rank_label) == 0) {
              NULL
            } else {
              rank_label %>%
                str_replace_all("_", " ") %>%
                paste(collapse = ", ")
            }
            base_subtitle <- if (!is.null(rank_label)) {
              paste0("rank: ", rank_label, " • top ", effective_limit)
            } else {
              paste0("top ", effective_limit, " pathways")
            }
            subtitle_text <- paste0(base_subtitle, " • source: ", collection_label)

            bubble_plot(
              sel,
              # title = comparison_name,
              title = comparison_label,
              subtitle = subtitle_text,
              save_func = local_save_func,
              facet_order = facet_order,
              nes_range = nes_range,
              size_range = size_range,
              glyph = glyph,
              ...
            )
          })
        }
      )
    }
  )
}

do_combined_bubble_plots <- function(
    results_list,
    save_func = NULL,
    facet_order = NULL,
    limit = 20,
    size_range = c(3.0, 9.0),
    glyph = DEFAULT_BUBBLE_GLYPH,
    ...) {
  genesets <- names(results_list)

  purrr::map(genesets, function(geneset_name) {
    fgsea_res_list <- results_list[[geneset_name]]
    collection_label <- format_source_label(geneset_name)

    purrr::map(limit, function(.limit) {
      res <- fgsea_res_list %>% bind_rows(.id = "rankname")
      res <- fgsea_tools$select_topn(res, limit = .limit, pstat_cutoff = 1)
      n_sel <- res %>% distinct(pathway) %>% nrow()
      effective_limit <- min(.limit, n_sel)
      if (n_sel == 0) {
        log_msg(warning = paste0(
          "no pathways available for combined bubble plot (requested top ", .limit, ") for ",
          geneset_name
        ))
        return(NULL)
      }
      log_msg(msg = paste0(
        "bubble combined: geneset=", geneset_name,
        " limit=", .limit,
        " selected_rows=", nrow(res),
        " distinct_pathways=", n_sel
      ))
      nes_max <- suppressWarnings(max(abs(res$NES), na.rm = TRUE))
      nes_range <- if (is.finite(nes_max)) c(-nes_max, nes_max) else NULL

      local_save_func <- save_func
      if (!is.null(local_save_func)) {
        geneset_dir <- util_tools$safe_path_component(geneset_name)
        filename <- util_tools$safe_filename(
          "bubble",
          paste0("top", effective_limit),
          paste0("n", n_sel),
          "all",
          fallback = "bubble_all"
        )
        path <- util_tools$safe_subdir(get_arg(local_save_func, "path"), geneset_dir, "bubble")
        local_save_func <- make_partial(local_save_func, filename = filename, path = path)
        log_msg(msg = paste0(
          "bubble combined: saving to ", file.path(path, paste0(filename, ".pdf"))
        ))
      }

      bubble_plot(
        res,
        title = geneset_name,
        subtitle = paste0("top ", effective_limit, " pathways • source: ", collection_label),
        save_func = local_save_func,
        facet_order = facet_order,
        nes_range = nes_range,
        size_range = size_range,
        glyph = glyph,
        ...
      )
    })
  })
}
