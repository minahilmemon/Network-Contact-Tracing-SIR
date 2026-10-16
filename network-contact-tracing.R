# install the ABM library
if (!require(ABM)) devtools::install_github("junlingm/ABM")

# Contact tracing for an SIR epidemic on a static contact network.
#
# This uses the current ABM main-branch API: the per-neighbour contact rate is
# attached to the Contact object, and the contact transition refers to its type.

#' Simulate contact tracing on a contact-network SIR model
#'
#' @param parameters A named list with elements `N`, `mean_degree`, `beta`,
#'   `gamma`, `tau`, `theta`, `p`, and
#'   ``beta` is the per-neighbour contact rate; `gamma` is
#'   the recovery rate; `tau` is the diagnosis rate; `theta` is the tracing
#'   completion rate; and `p` is the probability that a recorded contact is
#'   successfully traced.
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

    # Static Poisson-degree configuration-model contact network.
    sim$addContact(ABM::newConfigurationModel(
      rng = function(n) stats::rpois(n, mean_degree),
      rate = beta,
      type = "network"
    ))
  
    sim$addTransition(
      "I" -> "R",
      gamma,
      logging = list(ABM::dec("I"))
    )
    sim$addTransition(
      "I" -> "T",
      tau,
      logging = list(ABM::dec("I"), ABM::inc("T"))
    )

    # Record both endpoints of each I-S contact, including non-transmitting ones.
    sim$addTransition(
      "I" + "S" -> "I" + "I" ~ "network",
      to_change_callback = function(time, agent, contact) {
        agent_state <- ABM::getState(agent)
        agent_state$contact_ids <- unique(c(
          agent_state$contact_ids,
          ABM::getID(contact)
        ))
        ABM::setState(agent, agent_state)
  
        contact_state <- ABM::getState(contact)
        contact_state$contact_ids <- unique(c(
          contact_state$contact_ids,
          ABM::getID(agent)
        ))
        ABM::setState(contact, contact_state)

        TRUE
      },
      logging = list(ABM::inc("I"))
    )
  
    sim$addTransition(
      "T" -> "X",
      theta,
      changed_callback = function(time, agent) {
        for (id in ABM::getState(agent)$contact_ids) {
          contact <- sim$agent(id)
          if (stats::runif(1) >= p) {
            next
          }
  
          if (ABM::matchState(contact, "I")) {
            ABM::setState(contact, "T")
            state = sim$state
            state[["I"]] = state[["I"]] - 1
            state[["T"]] = state[["T"]] + 1
            sim$state = state
          }
        }
      },
      logging = list(ABM::dec("T"))
    )
  
    result = sim$run(times)
    setNames(as.numeric(result$I), as.character(result$times))
  })
}

parameters <- list(
  N = 100000L, mean_degree = 5, beta = 0.1,
  gamma = 0.1, tau = 0.15, theta = 1, p = 0.3
)

file="It_p030.csv"
if (file.exists(file)) {
  data = read.csv(file)
} else {
  times = 0:200
  data = data.frame(t=times)
  # make 10 runs #make more runs
  for (run in 1:100) {
    I <- network_contact_tracing(parameters, I0 = 20L, times=times)
    data = cbind(data, I)
  } 
  write.csv(data, file=file, row.names = FALSE)
}
# each column of data contains I(t) for one simulation.
# remove the time column then calculate an average at each time point
I = rowMeans(data[,-1])
plot(data$t, I, type="l")
