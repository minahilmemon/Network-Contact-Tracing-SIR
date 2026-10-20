# Run with: Rscript tests/check_comparison.R
# Uses synthetic CSVs and a stub ABM runner; no simulations are rerun.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(root, "compare.R"))

check_comparison <- function() {
  directory <- tempfile("comparison-check-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  abm_file <- file.path(directory, "It_p030_all_neighbors-2.csv")
  pde_file <- file.path(directory, "PDE_contact_tracing_dt0.5.csv")
  write.csv(data.frame(t = 0:2, run_1 = c(20, 30, 10), run_2 = c(20, 50, 0)),
            abm_file, row.names = FALSE)
  write.csv(data.frame(t = seq(0, 2, 0.5), I_PDE = c(20, 25, 35, 25, 5),
                       S_PDE = 99900), pde_file, row.names = FALSE)
  write.csv(data.frame(t = 0:2, mean_I = 999, sd_I = 500),
            file.path(directory, "It_p030_all_neighbors-2-summary.csv"), row.names = FALSE)
  write.csv(data.frame(t = 0:2, irrelevant = 999),
            file.path(directory, "unrelated.csv"), row.names = FALSE)
  stopifnot(length(find_result_csvs(c(directory, abm_file))) == 2L)
  device <- grDevices::dev.cur()
  png_file <- file.path(directory, "comparison.png")
  series <- compare_results(c(directory, abm_file), png_file)
  abm <- Filter(function(x) x$kind == "ABM", series)[[1]]
  pde <- Filter(function(x) x$kind == "PDE", series)[[1]]
  stopifnot(length(series) == 2L, identical(abm$I, c(20, 40, 5)),
            isTRUE(all.equal(pde$I, c(20, 25, 35, 25, 5))),
            length(abm$t) == 3L, length(pde$t) == 5L,
            file.info(png_file)$size > 0,
            identical(grDevices::dev.cur(), device))
  pdf_file <- file.path(directory, "pde-only.pdf")
  stopifnot(length(compare_results(pde_file, pdf_file)) == 1L,
            file.info(pdf_file)$size > 0)

  # Extinction points stay zero in the data but become gaps on a log axis.
  zero_file <- file.path(directory, "PDE_with_zeros.csv")
  write.csv(data.frame(t = 0:3, I_PDE = c(20, 0, 5, 0)), zero_file, row.names = FALSE)
  log_file <- file.path(directory, "comparison-log.png")
  logged <- withCallingHandlers(compare_results(c(abm_file, zero_file), log_file, log_y = TRUE),
                               warning = function(w) stop(conditionMessage(w)))
  stopifnot(file.info(log_file)$size > 0,
            identical(logged[[2]]$I, c(20L, 0L, 5L, 0L)),
            identical(grDevices::dev.cur(), device))
  empty_file <- file.path(directory, "PDE_extinct.csv")
  write.csv(data.frame(t = 0:2, I_PDE = 0), empty_file, row.names = FALSE)
  message <- tryCatch(compare_results(empty_file, log_file, log_y = TRUE), error = conditionMessage)
  stopifnot(is.character(message), grepl("positive infectious count", message))

  invalid <- file.path(directory, "PDE_invalid.csv")
  write.csv(data.frame(t = 0:2, I_PDE = c(20, NA, 0)), invalid, row.names = FALSE)
  warned <- FALSE
  result <- withCallingHandlers(read_result_series(invalid), warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  stopifnot(is.null(result), warned)

  # The ABM entry point must write/cache its data without opening a device.
  # Evaluate just this function definition so plotting tests do not need ABM.
  expressions <- parse(file.path(root, "network-contact-tracing.R"))
  definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                         identical(x[[2]], as.name("run_contact_tracing_comparison")),
                       as.list(expressions))
  env <- new.env(parent = globalenv())
  env$run_contact_tracing_replicates <- function(parameters, I0, times, runs, workers, seeds) {
    lapply(seq_len(runs), function(run) rep(run, length(times)))
  }
  eval(definition[[1]], env)
  data_dir <- file.path(directory, "simulation-only")
  dir.create(data_dir)
  data_file <- file.path(data_dir, "It_p030_all_neighbors-2.csv")
  result <- env$run_contact_tracing_comparison(runs = 2L, workers = 1L, file = data_file)
  cached <- env$run_contact_tracing_comparison(runs = 2L, workers = 1L, file = data_file)
  stopifnot(nrow(result) == 201L, ncol(result) == 3L,
            isTRUE(all.equal(result, cached)),
            identical(grDevices::dev.cur(), device),
            !any(grepl("\\.(png|pdf)$", list.files(data_dir))))
}

check_comparison()
cat("PASS: discovery, duplicate/summary exclusion, ABM means, separate time grids,\n",
    "      PNG/PDF/log output, zero handling, and headless ABM data output.\n")
