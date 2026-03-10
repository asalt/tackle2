suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(RcppTOML))
suppressPackageStartupMessages(library(here))
options(error = rlang::entrace)


source(file.path(here("R"), "run.R")) #
source(file.path(here("R"), "lazyloader.R")) #
util_tools <- get_tool_env("utils")

setClass("filetype")
setMethod(
  "coerce", c(from = "ANY", to = "filetype"),
  function(from, to) {
    if (!file.exists(from)) {
      stop(paste0(from, " does not exist, exiting.."))
    }
    return(from)
  }
)

get_parser <- function() {
  parser <- arg_parser("tackle2")
  parser <- add_argument(parser, "--quiet",
    help = "suppress console output (still writes run.log)",
    flag = TRUE,
    short = "-q"
  )
  parser <- add_argument(parser, "--verbose",
    help = "enable verbose console output (overrides --quiet)",
    flag = TRUE,
    short = "-v"
  )
  parser <- add_argument(parser, "--loglevel",
    help = "logging threshold: DEBUG, INFO, WARNING, ERROR",
    default = NULL,
    short = "-l"
  )
  parser <- add_argument(parser, "config", help = "toml config file", type = "filetype")
  return(parser)
}



main <- function() {
  parser <- get_parser()
  argv <- parse_args(parser)

  params <- RcppTOML::parseTOML(argv$config)
  params_in <- params$params %||% list()
  params_in$advanced <- params_in$advanced %||% list()

  if (isTRUE(argv$quiet) && isTRUE(argv$verbose)) {
    stop("Cannot specify both --quiet and --verbose")
  }

  if (isTRUE(argv$verbose)) {
    params_in$advanced$verbose <- TRUE
    params_in$advanced$quiet <- FALSE
  } else if (isTRUE(argv$quiet)) {
    params_in$advanced$quiet <- TRUE
    params_in$advanced$verbose <- FALSE
  }

  if (!is.null(argv$loglevel) && nzchar(argv$loglevel)) {
    params_in$advanced$loglevel <- toupper(argv$loglevel)
  }

  run(params_in) # named list with first order [params] and nested subsections
}

if (sys.nframe() == 0) { # if ran directly, not sourced, equivalent to python if __name__ == "__main__"
  main()
}
