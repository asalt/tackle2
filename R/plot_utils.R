suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(colorspace)

suppressPackageStartupMessages(library(here))

basedir <- file.path(here())
src_dir <- file.path(here("R"))

util_tools <- new.env()
source(file.path(src_dir, "./utils.R"), local = util_tools)
make_partial <- util_tools$make_partial
get_args <- util_tools$get_args
get_arg <- util_tools$get_arg
log_msg <- util_tools$make_partial(util_tools$log_msg)

fgsea_tools <- new.env()
source(file.path(src_dir, "./fgsea.R"), local = fgsea_tools)

# moved to utils.R
# # Helper function to get current preset arguments or an empty list if none are set
# get_args <- function(f) {
#   if (!is.null(attr(f, "preset_args"))) {
#     return(attr(f, "preset_args"))
#   } else {
#     return(list()) # Return an empty list for easier handling
#   }
# }
# get_arg <- function(f, arg) {
#   args <- get_args(f)
#   val <- args[[arg]]
#   if (is.null(val)) {
#     return("")
#   }
#   return(val)
# }

# # Custom partial function with dynamic argument handling
# make_partial <- function(.f, ...) {
#   # Retrieve current preset arguments, if any
#   current_args <- get_args(.f)

#   # New fixed arguments
#   args_fixed <- list(...)

#   # Ensure that named arguments are handled properly
#   if (!is.null(names(args_fixed))) {
#     # Overwrite or add new arguments
#     current_args[names(args_fixed)] <- args_fixed
#   }

#   # Inner function to call .f with the correct arguments
#   inner <- function(...) {
#     # Combine fixed arguments with any new ones provided at call time
#     args <- modifyList(current_args, list(...))
#     do.call(.f, args)
#   }

#   # Attach updated preset arguments to the inner function for later inspection
#   attr(inner, "preset_args") <- current_args

#   return(inner)
# }


plot_and_save_unsafe <- function(
    plot_code,
    filename,
    path = file.path(
      basedir,
      "plots/"
    ),
    type = "pdf",
    width = 8,
    height = 6,
    replace = T,
    ...) {
  # Setup: Open the appropriate graphics device

  #log_msg(msg = "plot_and_save")
  if (is.null(filename)){
      log_msg(warning = "filename is null")
  }

  filename <- util_tools$safe_path_component(filename, fallback = "plot")

  full_path <- file.path(path, paste0(filename, ".", type))
  log_msg(msg = paste0("plot_and_save: target ", full_path))

  max_total <- 240
  if (nchar(full_path) > max_total) {
    base_len <- nchar(file.path(path, ""))
    extension_len <- nchar(paste0(".", type))
    allowed <- max_total - base_len - extension_len
    allowed <- max(16, allowed)
    filename <- util_tools$safe_path_component(filename, fallback = "plot", max_chars = allowed)
    full_path <- file.path(path, paste0(filename, ".", type))
  }
  # log_msg(msg = paste0("full_path: ", full_path))

  if (!fs::dir_exists(path)) fs::dir_create(path)

  if (file.exists(full_path) && replace == FALSE) {
    log_msg(msg = paste0("plot_and_save: skipping existing file ", full_path))
    graphics.off() # turn off anything that opened
    return()
  }
  # ??
  on.exit(dev.off(), add = TRUE)

  if (type == "pdf") {
    # pdf(full_path, width = width, height = height)
    cairo_pdf(full_path, width = width, height = height)
  } else if (type == "png") {
    png(full_path, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file type")
  }


  # Execute the plot code (passed as a function)
  h <- plot_code()
  # if (inherits(h, "ggplot")) {  # this may be included already in the plot_code closure and unnecessary and not worth pytting here
  #   print(h)  # Required for ggplot rendering inside functions
  # }

  # Teardown: Close the graphics device
  # dev.off()
  # Ensure device closes even if an error occurs by putting the command here
  # is this the proper way to do this to ensure avoiding error of having too many open devices?
  # check all open devices with dev.list()
  # on.exit(dev.off(), add = TRUE)  moved to above

  # log_msg(msg = paste0("done"))

  log_msg(msg = paste0("plot_and_save: wrote ", full_path))
  return(h)
}

safe_plot_and_save <- function(...) {
  tryCatch(
    plot_and_save_unsafe(...),
    error = function(e) {
      message("\n❗️ Error caught: ", conditionMessage(e))
      message("\nLast traceback:")
      print(rlang::last_trace(drop = FALSE))

      message("\nOpen graphics devices:")
      print(dev.list())

      # Force close them
      if (length(dev.list()) > 0) {
        message("\nClosing all devices...")
        graphics.off()
      }

      # Optional: Write to log file
      log_msg(warning = paste("Plotting failed:", conditionMessage(e)))

      # Optionally re-throw if you want upstream failure
      # stop(e)
    }
  )
}
plot_and_save <- safe_plot_and_save #


# Matplotlib default colors
#matplotlib_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
#                       "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
#                       "#bcbd22", "#17becf")


is_numericish <- function(x) {
  if (is.numeric(x)) return(TRUE)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]                       # drop empty strings
  if (length(x) == 0) return(FALSE)
  suppressWarnings(all(!is.na(as.numeric(x))))
}

create_named_color_list <- function(df, columns, c=80) {

  # Initialize an empty list to store the result
  color_list <- list()

  # Iterate over each column
  user_map <- getOption(util_tools$pkg_option_name("user_colormap"), NULL)
  for (col_name in columns) {
    col_data <- df[[col_name]]

    # If numeric-like, build a continuous color mapping;
    # otherwise create a discrete palette with explicit NA color.
    if (is_numericish(col_data)) {
      suppressWarnings({ vals <- as.numeric(as.character(col_data)) })
      vals <- vals[is.finite(vals)]
      n_unique <- length(unique(vals))
      tol <- .Machine$double.eps^0.5
      is_integerish <- length(vals) > 0 && all(abs(vals - round(vals)) < tol)

      # Treat small-cardinality integer-like as categorical; otherwise continuous
      if (n_unique <= 5 && is_integerish) {
        # Discrete palette for small integer sets (e.g., 0/1 groups)
        unique_vals <- sort(unique(df[[col_name]]))
        n_vals <- length(unique_vals)
        colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
        assigned <- setNames(colors_assigned, as.character(unique_vals))
        # Apply user overrides (column-specific first, then global)
        if (!is.null(user_map)) {
          col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
          if (!is.null(col_override) && length(col_override) > 0) {
            matches <- intersect(names(assigned), names(col_override))
            assigned[matches] <- col_override[matches]
          }
          global_override <- tryCatch(user_map$global, error = function(e) NULL)
          if (!is.null(global_override) && length(global_override) > 0) {
            # apply only for values still not overridden
            remaining <- setdiff(names(assigned), names(col_override))
            matches2 <- intersect(remaining, names(global_override))
            assigned[matches2] <- global_override[matches2]
          }
        }
        if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
        if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
        color_list[[col_name]] <- assigned
      } else if (length(vals) >= 2) {
        minv <- min(vals, na.rm = TRUE)
        maxv <- max(vals, na.rm = TRUE)

        if (is.finite(minv) && is.finite(maxv) && minv != maxv) {
          # Diverging if range spans zero; otherwise sequential
          if (minv < 0 && maxv > 0) {
            brks <- c(minv, 0, maxv)
            cols <- c("#2166ac", "#f7f7f7", "#b2182b")
          } else {
            brks <- c(minv, maxv)
            cols <- c("#eff3ff", "#08519c")
          }
          color_list[[col_name]] <- circlize::colorRamp2(brks, cols)
          next
        } else {
          # Fallback for constant range
          unique_vals <- sort(unique(df[[col_name]]))
          n_vals <- length(unique_vals)
          colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
          assigned <- setNames(colors_assigned, as.character(unique_vals))
          if (!is.null(user_map)) {
            col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
            if (!is.null(col_override) && length(col_override) > 0) {
              matches <- intersect(names(assigned), names(col_override))
              assigned[matches] <- col_override[matches]
            }
            global_override <- tryCatch(user_map$global, error = function(e) NULL)
            if (!is.null(global_override) && length(global_override) > 0) {
              remaining <- setdiff(names(assigned), names(col_override))
              matches2 <- intersect(remaining, names(global_override))
              assigned[matches2] <- global_override[matches2]
            }
          }
          if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
          if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
          color_list[[col_name]] <- assigned
        }
      } else {
        # Not enough values to infer; treat as discrete
        unique_vals <- sort(unique(df[[col_name]]))
        n_vals <- length(unique_vals)
        colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
        assigned <- setNames(colors_assigned, as.character(unique_vals))
        if (!is.null(user_map)) {
          col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
          if (!is.null(col_override) && length(col_override) > 0) {
            matches <- intersect(names(assigned), names(col_override))
            assigned[matches] <- col_override[matches]
          }
          global_override <- tryCatch(user_map$global, error = function(e) NULL)
          if (!is.null(global_override) && length(global_override) > 0) {
            remaining <- setdiff(names(assigned), names(col_override))
            matches2 <- intersect(remaining, names(global_override))
            assigned[matches2] <- global_override[matches2]
          }
        }
        if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
        if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
        color_list[[col_name]] <- assigned
      }
    } else {
      # Discrete categorical mapping
      col_data[is.na(col_data)] <- "NA"
      unique_vals <- sort(unique(col_data))
      n_vals <- length(unique_vals)
      colors_assigned <- colorspace::qualitative_hcl(palette='Dynamic', n=n_vals, c=c)
      assigned <- setNames(colors_assigned, as.character(unique_vals))
      # Apply user overrides (column-specific first, then global)
      if (!is.null(user_map)) {
        col_override <- tryCatch(user_map$by_column[[col_name]], error = function(e) NULL)
        if (!is.null(col_override) && length(col_override) > 0) {
          matches <- intersect(names(assigned), names(col_override))
          assigned[matches] <- col_override[matches]
        }
        global_override <- tryCatch(user_map$global, error = function(e) NULL)
        if (!is.null(global_override) && length(global_override) > 0) {
          remaining <- setdiff(names(assigned), names(col_override))
          matches2 <- intersect(remaining, names(global_override))
          assigned[matches2] <- global_override[matches2]
        }
      }
      if ("NA" %in% names(assigned)) assigned[["NA"]] <- "grey40"
      if ("na" %in% names(assigned)) assigned[["na"]] <- "grey40"
      color_list[[col_name]] <- assigned
    }
  }

  return(color_list)
}



process_cut_by <- function(cut_by, cdesc) {
  # Return NULL immediately if cut_by is NULL
  if (is.null(cut_by) || length(cut_by) == 0) {
    return(NULL)
  }

  # Treat NA and FALSE as "no cut"
  if (is.logical(cut_by) && length(cut_by) == 1 && isFALSE(cut_by)) {
    return(NULL)
  }
  if (all(is.na(cut_by))) {
    return(NULL)
  }

  # If cut_by is a single string containing ':', split it into a vector
  if (is.character(cut_by) && length(cut_by) == 1 && grepl(":", cut_by)) {
    cut_by <- strsplit(cut_by, ":")[[1]]
  }

  # Ensure cut_by is now a character vector
  if (!is.character(cut_by)) {
    # warning("cut_by should be a character string or vector.")
    # return(NULL)
    # this is fine
    cut_by <- as.character(cut_by)
  }

  # Drop empty/NA values after coercion
  cut_by <- cut_by[!is.na(cut_by)]
  cut_by <- cut_by[nzchar(trimws(cut_by))]
  if (length(cut_by) == 0) {
    return(NULL)
  }

  # Check if all elements in cut_by are valid column names
  invalid_cols <- setdiff(cut_by, colnames(cdesc))
  if (length(invalid_cols) > 0) {
    warning(
      "The following cut_by elements are not column names in cdesc: ",
      paste(invalid_cols, collapse = ", ")
    )
    return(NULL)
  }

  # Subset the relevant columns and create the interaction factor
  cut_by_factor <- interaction(cdesc[, cut_by, drop = FALSE], drop = TRUE)

  return(cut_by_factor)
}
