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

# Source this file to load the functions without starting 100 replicates.
# Rscript network-contact-tracing.R runs the default comparison below.
run_contact_tracing_comparison <- function() {
  parameters <- list(
    N = 100000L, mean_degree = 5, beta = 0.1,
    gamma = 0.1, tau = 0.15, theta = 1, p = 0.3
  )

  # Do not reuse simulations produced by the former transmission-partner model.
  file="It_p030_all_neighbors-1000.csv"
  if (file.exists(file)) {
    data = read.csv(file)
  } else {
    times = 0:200
    data = data.frame(t=times)
    # Generate 100 independent replicates.
    for (run in 1:1000) {
      cat("run:", run, "\n")
      I <- network_contact_tracing(parameters, I0 = 20L, times=times)
      data = cbind(data, I)
    }
    write.csv(data, file=file, row.names = FALSE)
  }
  # each column of data contains I(t) for one simulation.
  # remove the time column then calculate an average at each time point
  I = rowMeans(data[,-1])
  plot(data$t, I, type="l")
}

if (sys.nframe() == 0L) run_contact_tracing_comparison()
