# Run with: Rscript tests/check_abm_tracing.R
# Requires an ABM build exporting getNeighbors(). No epidemic CSVs are written.
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1]))))
} else normalizePath('.')
source(file.path(root, 'network-contact-tracing.R'))

make_case <- function(states, degrees, p, rate = 0, seed = 1) {
  set.seed(seed)
  sim <- ABM::Simulation$new(lapply(states, function(x) list(x)))
  sim$state <- list(I = sum(states == 'I'), T = sum(states == 'T'))
  sim$addLogger('I')
  sim$addLogger('T')
  network <- ABM::newConfigurationModel(function(n) degrees,
                                       rate = rate, type = 'network')
  sim$addContact(network)
  cleanup <- add_contact_tracing_transitions(sim, network, gamma = 0, tau = 0,
                                 theta = function(time) 1, p = p)
  invisible(sim$run(0))
  list(sim = sim, network = network, cleanup = cleanup)
}

# Neither endpoint has any contact history. With p=1 the I neighbour must
# nevertheless be traced; with p=0 only the original T becomes X.
for (p in c(0, 1)) {
  case <- make_case(c('T', 'I'), c(1L, 1L), p)
  result <- case$sim$resume(c(0.5, 1.5, 2.5))
  stopifnot(identical(as.numeric(result$I), c(1, 1-p, 1-p)),
            identical(as.numeric(result$T), c(1, p, 0)),
            ABM::matchState(case$sim$agent(1), 'X'),
            ABM::matchState(case$sim$agent(2), if (p == 1) 'X' else 'I'))
  case$cleanup()
}

# S, R, T and X neighbours are never changed to T by another node's episode.
for (state in c('S', 'R', 'T', 'X')) {
  case <- make_case(c('T', state), c(1L, 1L), p = 1)
  result <- case$sim$resume(c(1.5, 3.5))
  stopifnot(all(result$I == 0), all(result$T == 0),
            ABM::matchState(case$sim$agent(2), if (state == 'T') 'X' else state))
  case$cleanup()
}

# Construct the simple path A--B--C from configuration-model pairing.
# A starts T, B starts S, C starts I. B is infected by C at t=.25,
# then traced by A at t=1 despite never transmitting to/from A.
case <- NULL
for (seed in 1:20) {
  candidate <- make_case(c('T', 'S', 'I'), c(1L, 2L, 1L), p = 1,
                         rate = function(time) 0.25, seed = seed)
  if (length(ABM::getNeighbors(candidate$network, candidate$sim$agent(2))) == 2) {
    case <- candidate
    break
  }
  candidate$cleanup()
}
stopifnot(!is.null(case))
first <- case$sim$resume(0.5)
stopifnot(first$I == 2, ABM::matchState(case$sim$agent(2), 'I'))
second <- case$sim$resume(1.5)
stopifnot(second$I == 1, second$T == 1,
          ABM::matchState(case$sim$agent(1), 'X'),
          ABM::matchState(case$sim$agent(2), 'T'),
          ABM::matchState(case$sim$agent(3), 'I'))
last <- case$sim$resume(c(2.5, 3.5, 5.5))
stopifnot(all(last$I == 0), identical(as.numeric(last$T), c(1, 0, 0)))
for (i in 1:3) stopifnot(ABM::matchState(case$sim$agent(i), 'X'))
case$cleanup()

# Exercise the public entry point on a stochastic network as well.
for (seed in 2026:2027) {
  set.seed(seed)
  result <- network_contact_tracing(
    list(N = 1000L, mean_degree = 5, beta = 0.1, gamma = 0.1,
         tau = 0.15, theta = 1, p = 0.3), I0 = 20L, times = 0:5)
  stopifnot(length(result) == 6, result[1] == 20,
            all(is.finite(result)), all(result >= 0), all(result <= 1000))
}
cat('PASS: all-neighbour tracing, p=0/1, state eligibility, outside infection,\n',
    '      tracing chains, one episode per T, and stochastic model smoke test.\n')
