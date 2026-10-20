# Run with: Rscript tests/check_abm_fresh_sessions.R
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))

check_fresh_sessions <- function() {
  source(file.path(root, "network-contact-tracing.R"), local = TRUE)
  model <- network_contact_tracing
  original_plan <- future::plan()
  # Observe the worker process directly, without creating a large simulation.
  network_contact_tracing <- function(parameters, I0, times) {
    draw <- stats::runif(1)
    if (I0 < 0 && draw < 0.5) stop("expected simulation failure")
    list(pid = Sys.getpid(), draw = draw)
  }
  seeds <- as.list(1:4)
  first <- run_contact_tracing_replicates(
    list(), 0L, 0, runs = 4L, workers = 1L, seeds = seeds)
  second <- run_contact_tracing_replicates(
    list(), 0L, 0, runs = 4L, workers = 2L, seeds = seeds)
  stopifnot(length(unique(vapply(first, `[[`, integer(1), "pid"))) == 4L,
            length(unique(vapply(second, `[[`, integer(1), "pid"))) == 4L,
            all(vapply(first, `[[`, integer(1), "pid") != Sys.getpid()),
            identical(lapply(first, `[[`, "draw"), lapply(second, `[[`, "draw")),
            identical(class(future::plan()), class(original_plan)))

  # With no seed list, fresh workers initialize distinct RNG streams.
  unseeded <- run_contact_tracing_replicates(
    list(), 0L, 0, runs = 4L, workers = 2L)
  stopifnot(length(unique(vapply(unseeded, `[[`, integer(1), "pid"))) == 4L,
            length(unique(vapply(unseeded, `[[`, numeric(1), "draw"))) == 4L)

  results <- run_contact_tracing_replicates(
    list(), -1L, 0, runs = 4L, workers = 2L, seeds = seeds)
  stopifnot(length(results) == 4L,
            inherits(results[[1]], "error"), results[[1]]$run == 1L,
            conditionMessage(results[[1]]) == "expected simulation failure",
            !inherits(results[[4]], "error"),
            identical(results[[4]]$draw, first[[4]]$draw))

  # Check the actual ABM model against an identically seeded direct run.
  network_contact_tracing <- model
  parameters <- list(N = 1000L, mean_degree = 5, beta = 0.1,
                     gamma = 0.1, tau = 0.15, theta = 1, p = 0.3)
  isolated <- run_contact_tracing_replicates(
    parameters, 20L, 0:5, runs = 2L, workers = 2L,
    seeds = list(2026L, 2027L))
  for (i in 1:2) {
    set.seed(2025L + i)
    expected <- network_contact_tracing(parameters, 20L, 0:5)
    stopifnot(identical(isolated[[i]], expected))
  }
}

check_fresh_sessions()
cat("PASS: fresh sessions at one/two workers, independent unseeded RNG, seed reproducibility, isolated\n",
    "      errors, restored future plan, and seeded ABM equivalence.\n")
