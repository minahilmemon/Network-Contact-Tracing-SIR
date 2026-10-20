# Run: Rscript tests/check_abm_dynamics.R
# Test the production transitions against independent probability calculations.
# Uses small diagnostic populations; does not overwrite epidemic CSVs.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else normalizePath(".")
source(file.path(root, "network-contact-tracing.R"))

make_case <- function(states, degree_rng, beta, gamma, tau, theta, p) {
  sim <- ABM::Simulation$new(lapply(states, function(x) list(x)))
  sim$state <- list(I = sum(states == "I"), T = sum(states == "T"))
  sim$addLogger("I"); sim$addLogger("T")
  network <- ABM::newConfigurationModel(degree_rng, rate = beta, type = "network")
  sim$addContact(network)
  cleanup <- add_contact_tracing_transitions(sim, network, gamma, tau, theta, p)
  invisible(sim$run(0))
  list(sim = sim, network = network, cleanup = cleanup, N = length(states))
}
agent_states <- function(case) {
  vapply(seq_len(case$N), function(i) ABM::getState(case$sim$agent(i))[[1]],
         character(1))
}
observe <- function(case, time) {
  result <- case$sim$resume(time)
  states <- agent_states(case)
  stopifnot(result$I == sum(states == "I"), result$T == sum(states == "T"),
            all(states %in% c("S", "I", "T", "R", "X")))
  c(I = result$I, T = result$T)
}
within_sampling_error <- function(actual, expected, variance, label) {
  # Fixed seeds make these checks reproducible; six standard deviations
  # allow ordinary Monte Carlo variation while detecting rate mistakes.
  limit <- 6*sqrt(variance) + 2
  if (abs(actual-expected) > limit)
    stop(sprintf("%s: observed %.3f, expected %.3f, tolerance %.3f",
                 label, actual, expected, limit))
}
adjacency <- function(case) {
  lapply(seq_len(case$N), function(i) {
    vapply(ABM::getNeighbors(case$network, case$sim$agent(i)),
           ABM::getID, numeric(1))
  })
}

check_rates <- function() {
  # Every susceptible has degree one. Infection can only arrive from an
  # initially infected hub; there is no onward transmission to another S.
  # Each I--S edge therefore infects independently with probability 1-e^-bt,
  # including when its source has multiple neighbours.
  set.seed(9301)
  hubs <- 2000L
  states <- c(rep("I", hubs), rep("S", 3*hubs))
  degrees <- c(rep(3L, hubs), rep(1L, 3*hubs))
  beta <- .4
  case <- make_case(states, function(n) degrees, beta, 0, 0, 0, 0)
  on.exit(case$cleanup())
  neighbours <- adjacency(case)
  susceptible <- seq.int(hubs+1L, length(states))
  exposed <- sum(vapply(neighbours[susceptible], function(x) any(x <= hubs), logical(1)))
  stopifnot(exposed > 1000)
  for (time in c(.5, 1, 2, 4)) {
    counts <- observe(case, time)
    probability <- -expm1(-beta*time)
    within_sampling_error(counts["I"]-hubs, exposed*probability,
                          exposed*probability*(1-probability), "per-edge transmission")
  }
  cat("PASS transmission rate beta per I--S edge, including multi-neighbour sources.\n")
}

check_tracing_probability <- function() {
  set.seed(9302)
  states <- c(rep("T", 1000), rep("I", 3000), rep("S", 400), rep("R", 300), rep("X", 300))
  p <- .3
  case <- make_case(states, function(n) rpois(n, 5), 0, 0, 0, function(time) 1, p)
  on.exit(case$cleanup())
  neighbours <- adjacency(case)
  infected <- which(states == "I")
  # Each initial T fires at t=1. A target with m such neighbours escapes
  # with probability (1-p)^m. Newly traced targets cannot fire before t=2.
  m <- vapply(neighbours[infected], function(x) sum(states[x] == "T"), integer(1))
  traced_probability <- 1-(1-p)^m
  stopifnot(any(m == 0), any(m >= 2), observe(case, .5)["I"] == length(infected))
  counts <- observe(case, 1.5)
  actual <- length(infected)-counts["I"]
  within_sampling_error(actual, sum(traced_probability),
                        sum(traced_probability*(1-traced_probability)),
                        "Bernoulli tracing of infected neighbours")
  final <- agent_states(case)
  stopifnot(counts["T"] == actual, all(final[states == "T"] == "X"),
            identical(final[states %in% c("S","R","X")],
                      states[states %in% c("S","R","X")]))
  cat("PASS p=.3 tracing, shared targets, state eligibility and one episode per source.\n")
}

check_pair_markov_model <- function() {
  # Degree-one networks are disjoint pairs. Enumerate all 25 pair states
  # and evolve their exact continuous-time Markov generator.
  labels <- c("S","I","T","R","X")
  pairs <- expand.grid(left=labels, right=labels, stringsAsFactors=FALSE)
  keys <- paste(pairs$left, pairs$right, sep=":")
  index <- function(state) match(paste(state, collapse=":"), keys)
  beta <- .4; gamma <- .2; tau <- .3; theta <- .7; p <- .3
  Q <- matrix(0, 25, 25)
  for (i in seq_len(25)) {
    state <- unlist(pairs[i, ], use.names=FALSE)
    add <- function(next_state, rate) {
      j <- index(next_state)
      Q[i,j] <<- Q[i,j]+rate
    }
    for (focal in 1:2) {
      other <- 3-focal
      if (state[focal] == "I") {
        next_state <- state; next_state[focal] <- "R"; add(next_state, gamma)
        next_state <- state; next_state[focal] <- "T"; add(next_state, tau)
        if (state[other] == "S") {
          next_state <- state; next_state[other] <- "I"; add(next_state, beta)
        }
      }
      if (state[focal] == "T") {
        next_state <- state; next_state[focal] <- "X"
        if (state[other] == "I") {
          add(next_state, theta*(1-p))
          next_state[other] <- "T"; add(next_state, theta*p)
        } else add(next_state, theta)
      }
    }
  }
  diag(Q) <- -rowSums(Q)
  nI <- rowSums(pairs == "I"); nT <- rowSums(pairs == "T")
  # Uniformization evaluates exp(tQ) through a Poisson-weighted series.
  transition_matrix <- function(time) {
    rate <- max(-diag(Q))
    P <- diag(25)+Q/rate
    power <- diag(25)
    transition <- dpois(0, rate*time)*power
    for (i in seq_len(qpois(1-1e-14, rate*time))) {
      power <- power %*% P
      transition <- transition+dpois(i, rate*time)*power
    }
    stopifnot(max(abs(rowSums(transition)-1)) < 1e-12)
    transition
  }
  set.seed(9303)
  states <- rep(labels, each=2000)
  case <- make_case(states, function(n) rep(1L,n), beta, gamma, tau, theta, p)
  on.exit(case$cleanup())
  neighbour <- vapply(adjacency(case), function(x) x[[1]], numeric(1))
  first <- which(seq_along(states) < neighbour)
  stopifnot(length(first) == length(states)/2)
  initial <- tabulate(match(paste(states[first], states[neighbour[first]], sep=":"), keys),
                      nbins=25)
  for (time in c(.5, 1, 2, 4)) {
    counts <- observe(case, time)
    transition <- transition_matrix(time)
    for (label in c("I","T")) {
      values <- if (label == "I") nI else nT
      means <- as.vector(transition %*% values)
      variances <- as.vector(transition %*% (values^2))-means^2
      within_sampling_error(counts[label], sum(initial*means),
                            sum(initial*variances), paste("exact pair CTMC", label))
    }
  }
  cat("PASS exact two-node Markov model: transmission, recovery, diagnosis and tracing clocks.\n")
}

check_observation_times <- function() {
  parameters <- list(N=1000L, mean_degree=0, beta=0, gamma=1, tau=0, theta=0, p=0)
  set.seed(31)
  full <- network_contact_tracing(parameters, 1000L, 0:3)
  set.seed(31)
  later <- network_contact_tracing(parameters, 1000L, c(2,2,3))
  stopifnot(identical(later, full[c("2","2","3")]))
  for (times in list(-1, c(0, -1), NA_real_, Inf, numeric())) {
    error <- tryCatch(network_contact_tracing(parameters, 1000L, times), error=identity)
    stopifnot(inherits(error, "error"), grepl("times must", conditionMessage(error)))
  }
  cat("PASS observations start from the same time-zero epidemic, including duplicate times.\n")
}

check_rates()
check_tracing_probability()
check_pair_markov_model()
check_observation_times()
cat("All ABM dynamics checks passed; logged I/T equal actual agent counts at every check.\n")
