# lazyloader.R
# TODO have every tool in the R/ directory be loaded in a lazy manner with this script
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))

# Only define get_tool_env if it doesn't already exist
if (!exists("get_tool_env", envir = .GlobalEnv)) {
  get_tool_env <- local({
    # all_tools <- file.path(here("R"), list.files(here("R"), pattern = "\\.R$", full.names = FALSE))
    all_tools <- fs::dir_ls(path = here("R"), glob = "*.R") %>%
      map(basename) %>%
      purrr::list_c() %>%
      str_remove(".R")

    if (interactive()) {
      print(getwd())
      print(all_tools)
    }

    # all_tools <- c("heatmap", "io", "utils")
    #heatmap.R  io.R  lazyloader.R  reduce.R  summarize_rmd.R  utils.R

    # Initialize cache with empty environments
    tools_cache <- setNames(vector("list", length(all_tools)), all_tools)
    for (tool_name in all_tools) {
      tools_cache[[tool_name]] <- new.env()
    }

    # Keep track of tools currently being loaded
    tools_loading <- list()

    function(tool_name) {
      if (!tool_name %in% names(tools_cache)) {
        stop("Unknown tool: ", tool_name)
      }

      # Return environment if already loaded
      if (exists(".__loaded__", envir = tools_cache[[tool_name]])) {
        return(tools_cache[[tool_name]])
      }

      # Return environment if already loading (circular dependency)
      if (tool_name %in% tools_loading) {
        return(tools_cache[[tool_name]])
      }

      # Mark as currently loading
      tools_loading <<- c(tools_loading, tool_name)

      # Source the tool file into its environment
      src_dir <- file.path(here("R"))
      source_file <- paste0(tool_name, ".R")
      full_path <- file.path(src_dir, source_file)
      if (!file.exists(full_path)) {
        tools_loading <<- setdiff(tools_loading, tool_name)
        stop("Source file does not exist for tool: ", tool_name)
      }

      # Provide access to get_tool_env
      env <- tools_cache[[tool_name]]
      env$get_tool_env <- get_tool_env
      env$tools_cache <- tools_cache

      # Source the file
      source(full_path, local = env)

      # Mark as loaded
      assign(".__loaded__", TRUE, envir = env)

      # Remove from loading state
      tools_loading <<- setdiff(tools_loading, tool_name)

      return(env)
    }
  })

  # Assign get_tool_env to the global environment
  assign("get_tool_env", get_tool_env, envir = .GlobalEnv)
}
