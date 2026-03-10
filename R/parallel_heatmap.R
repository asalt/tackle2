suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(purrr))

suppressPackageStartupMessages(library(here))
source(file.path(here("R"), "lazyloader.R"))

resolve_multisession_workers <- function(job_count, workers = NULL, reserve = 1L) {
  if (is.null(job_count) || job_count <= 0) {
    return(0L)
  }

  if (!is.null(workers)) {
    workers <- as.integer(workers)
    workers <- workers[!is.na(workers)]
    if (length(workers) == 0) return(0L)
    return(max(0L, min(job_count, max(workers))))
  }

  available <- future::availableCores()
  if (!is.numeric(available) || length(available) == 0) {
    available <- 1L
  }

  available <- max(1L, as.integer(available[1]))
  target <- max(1L, available - reserve)
  min(job_count, target)
}

with_multisession_plan <- function(workers, expr) {
  workers <- as.integer(workers)
  if (length(workers) == 0 || workers <= 1L) {
    return(eval.parent(substitute(expr)))
  }

  prev_plan <- future::plan()
  on.exit(future::plan(prev_plan), add = TRUE)

  future::plan(future::multisession, workers = workers)
  eval.parent(substitute(expr))
}

run_multisession_jobs <- function(
    jobs,
    job_fn,
    workers = NULL,
    packages = NULL,
    globals = TRUE,
    stop_on_error = FALSE,
    reserve = 1L,
    seed = TRUE) {

  if (!is.function(job_fn)) {
    stop("job_fn must be a function")
  }

  job_count <- length(jobs)
  resolved_workers <- resolve_multisession_workers(job_count, workers, reserve = reserve)

  results <- vector("list", job_count)
  errors <- vector("list", job_count)

  if (resolved_workers <= 1L || job_count <= 1L) {
    for (idx in seq_len(job_count)) {
      job <- jobs[[idx]]
      results[[idx]] <- tryCatch(
        job_fn(job),
        error = function(e) {
          errors[[idx]] <<- e
          if (isTRUE(stop_on_error)) stop(e)
          NULL
        }
      )
    }
    return(list(results = results, errors = errors, workers = 1L))
  }

  run_expr <- quote({
    futures <- purrr::imap(jobs, function(job, idx) {
      future::future(
        expr = job_fn(job),
        packages = packages,
        globals = globals,
        seed = seed
      )
    })

    purrr::iwalk(futures, function(fut, idx) {
      results[[idx]] <<- tryCatch(
        future::value(fut),
        error = function(e) {
          errors[[idx]] <<- e
          if (isTRUE(stop_on_error)) stop(e)
          NULL
        }
      )
    })
  })

  with_multisession_plan(resolved_workers, eval(run_expr))

  list(results = results, errors = errors, workers = resolved_workers)
}
