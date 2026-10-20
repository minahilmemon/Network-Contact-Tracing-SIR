if (!requireNamespace("ABM", quietly = TRUE)) {
  stop("Install the ABM package before running this simulation.")
}

# Contact tracing for an SIR epidemic on a static contact network.
#
# This uses the current ABM main-branch API: the per-neighbour contact rate is
# attached to the Contact object, and the contact transition refers to its type.

# Register the production transitions separately so small fixed-network
# examples can test exactly the same transmission and tracing callbacks.
add_contact_tracing_transitions <- function(sim, network, gamma, tau, theta, p) {
  if (!"getNeighbors" %in% getNamespaceExports("ABM")) {
    stop("Update ABM to a build exporting getNeighbors() before running this model.")
  }
  sim$addTransition(
    "I" -> "R", gamma, logging = list(ABM::dec("I"))
  )
  sim$addTransition(
    "I" -> "T", tau, logging = list(ABM::dec("I"), ABM::inc("T"))
  )
  sim$addTransition(
    "I" + "S" -> "I" + "I" ~ "network",
    logging = list(ABM::inc("I"))
  )
  # A process-wide callback key must distinguish independent simulations.
  # Use a small explicit environment so the registry can release the graph
  # after a replicate, rather than retaining the whole simulation forever.
  context <- list2env(list(sim = sim, network = network, p = p), parent = baseenv())
  trace_neighbors <- function(time, agent) {
    # The source is already X. Trace all infected network neighbours,
    # regardless of which node infected them or any prior contact history.
    for (contact in ABM::getNeighbors(network, agent)) {
      if (ABM::matchState(contact, "I") && stats::runif(1) < p) {
        ABM::setState(contact, "T")
        state <- sim$state
        state[["I"]] <- state[["I"]] - 1
        state[["T"]] <- state[["T"]] + 1
        sim$state <- state
      }
    }
  }
  environment(trace_neighbors) <- context
  trace_neighbors <- ABM::registerRFunction(
    trace_neighbors, name = tempfile("network-contact-tracing-"))
  sim$addTransition(
    "T" -> "X", theta, changed_callback = trace_neighbors,
    logging = list(ABM::dec("T"))
  )
  # Call only after the simulation is finished and will not be resumed.
  function() rm(list = ls(context, all.names = TRUE), envir = context)
}

#' Simulate contact tracing on a contact-network SIR model
#'
#' @param parameters A named list with elements `N`, `mean_degree`, `beta`,
#'   `gamma`, `tau`, `theta`, and `p`.
#'   `beta` is the per-neighbour contact rate; `gamma` is
#'   the recovery rate; `tau` is the diagnosis rate; `theta` is the tracing
#'   completion rate; and `p` is the probability that an infected network
#'   neighbour is successfully traced during a tracing episode.
#' @param I0 The number infectious at time zero.
#' @param times A nondecreasing numeric vector of observation times.
#'
#' @return A named numeric vector containing I(t) at the requested time points.
network_contact_tracing <- function(parameters, I0, times) {
  with (parameters, {
    N <- as.integer(N)

    if (!is.numeric(I0) || length(I0) != 1L || !is.finite(I0) ||
        I0 != floor(I0) || I0 < 0 || I0 > N) {
      stop("I0 must be an integer between 0 and N", call. = FALSE)
    }
    I0 <- as.integer(I0)

    if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times)) ||
        any(diff(times) < 0)) {
      stop("times must be a nonempty, nondecreasing finite numeric vector",
           call. = FALSE)
    }

    # first I0 individuals are I, the remaining ones are S
    initial_states = c(
      rep(list(list("I")), I0),
      rep(list(list("S")), N-I0))
    # IDs are stable because this model has no births, deaths, or migration.
    sim <- ABM::Simulation$new(initial_states)

    sim$state <- list(
      I = I0,
      T = 0
    )
    for (state_name in names(sim$state)) {
      sim$addLogger(state_name)
    }

    # The same stored network supplies transmission and tracing neighbours.
    network <- ABM::newConfigurationModel(
      rng = function(n) stats::rpois(n, mean_degree),
      rate = beta,
      type = "network"
    )
    sim$addContact(network)
    release_tracing <- add_contact_tracing_transitions(sim, network, gamma, tau, theta, p)
    on.exit(release_tracing(), add = TRUE)

    result = sim$run(times)
    setNames(as.numeric(result$I), as.character(result$times))
  })
}

# Each replicate runs in a fresh R session so ABM's native allocations are
# reclaimed when that process exits. NULL seeds retain fresh-session RNG setup;
# a list of one integer seed per replicate permits reproducible reruns.
run_contact_tracing_replicates <- function(parameters, I0, times, runs = 1000L,
                                         workers = 8L, seeds = NULL) {
  for (package in c("future", "future.apply", "future.callr")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the ", package, " package before running replicates.")
    }
  }
  is_count <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x <= .Machine$integer.max && x == floor(x)
  }
  if (!is_count(runs) || !is_count(workers)) {
    stop("runs and workers must be positive integers.", call. = FALSE)
  }
  if (!is.null(seeds) &&
      (!is.list(seeds) || length(seeds) != runs ||
       !all(vapply(seeds, function(seed) {
         is.numeric(seed) && length(seed) == 1L && is.finite(seed) &&
           abs(seed) <= .Machine$integer.max && seed == floor(seed)
       }, logical(1))))) {
    stop("seeds must be NULL or a list of one integer seed per run.",
         call. = FALSE)
  }

  old_plan <- future::plan(future.callr::callr, workers = min(workers, runs))
  on.exit(future::plan(old_plan), add = TRUE)
  future.apply::future_lapply(
    seq_len(runs),
    function(run, parameters, I0, times, seeds) {
      result <- tryCatch({
        if (!is.null(seeds)) set.seed(seeds[[run]])
        network_contact_tracing(parameters, I0, times)
      }, error = function(error) {
        error$run <- run
        error
      })
      cat("run:", run, if (inherits(result, "error")) "failed" else "finished", "\n")
      result
    },
    parameters = parameters, I0 = I0, times = times, seeds = seeds,
    future.globals = list(
      network_contact_tracing = network_contact_tracing,
      add_contact_tracing_transitions = add_contact_tracing_transitions
    ),
    future.chunk.size = 1L,
    future.seed = NULL
  )
}

# Source this file to load the functions without starting 1000 replicates.
# Rscript network-contact-tracing.R runs the default comparison below.
run_contact_tracing_comparison <- function(runs = 1000L, workers = 8L,
                                         seeds = NULL,
                                         file = sprintf("It_p030_all_neighbors-%d.csv", runs)) {
  parameters <- list(
    N = 100000L, mean_degree = 5, beta = 0.1,
    gamma = 0.1, tau = 0.15, theta = 1, p = 0.3
  )

  # Do not reuse simulations produced by the former transmission-partner model.
  times <- 0:200
  if (file.exists(file)) {
    data <- read.csv(file)
    if (ncol(data) != runs + 1L || !identical(data$t, times)) {
      stop("Existing results do not match the requested runs and observation times.")
    }
  } else {
    results <- run_contact_tracing_replicates(
      parameters, I0 = 20L, times = times, runs = runs,
      workers = workers, seeds = seeds
    )
    # Preserve every result, including error conditions, before making the CSV.
    results_file <- paste0(tools::file_path_sans_ext(file), ".rds")
    saveRDS(list(parameters = parameters, I0 = 20L, times = times,
                 seeds = seeds, results = results,
                 session_info = sessionInfo()), results_file)
    failed <- vapply(results, inherits, logical(1), what = "error")
    if (any(failed)) {
      stop(sum(failed), " runs failed; all results are saved in ", results_file,
           ". Failed run IDs: ", paste(which(failed), collapse = ", "))
    }
    trajectories <- do.call(cbind, results)
    colnames(trajectories) <- paste0("run_", seq_len(runs))
    data <- data.frame(t = times, trajectories)
    write.csv(data, file = file, row.names = FALSE)
  }
  # each column of data contains I(t) for one simulation.
  # remove the time column then calculate an average at each time point
  I = rowMeans(data[,-1, drop=FALSE])
  plot(data$t, I, type="l", xlab="Time", ylab="Mean infectious individuals")
  invisible(data)
}

if (sys.nframe() == 0L) run_contact_tracing_comparison()
