# Run with: Rscript tests/check_pde_consistency.R
# Evaluate only initialization and the time loop, avoiding plots, CSV writes,
# package installs, and full-size runs. No solver source is rewritten by tests.
script_arg <- grep('^--file=', commandArgs(), value = TRUE)
root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1]))))
} else {
  normalizePath('.')
}

initialize <- function(overrides = list()) {
  config <- modifyList(list(dt = 0.4, da = 0.4, t_max = 4,
                            a_max = 8, Kmax = 6, mean_degree = 2), overrides)
  exprs <- parse(file.path(root, 'pde.R'))
  loop <- which(vapply(exprs, function(x) {
    is.call(x) && identical(x[[1]], as.name('for')) &&
      identical(x[[2]], as.name('n'))
  }, logical(1)))
  stopifnot(length(loop) == 1)
  env <- new.env(parent = globalenv())
  invisible(capture.output(for (expr in exprs[seq_len(loop - 1L)]) {
    if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names(config)) {
      expr[[3]] <- config[[as.character(expr[[2]])]]
    }
    eval(expr, env)
  }))
  list(env = env, loop = exprs[[loop]])
}
run <- function(case) {
  invisible(capture.output(eval(case$loop, case$env)))
  case$env
}
close_to <- function(actual, expected, label, tolerance = 1e-10) {
  error <- max(abs(actual - expected))
  if (!is.finite(error) || error > tolerance * max(1, abs(expected))) {
    stop(sprintf('%s: error %.12g', label, error))
  }
}

case <- initialize()
e <- case$env
mass <- sum(e$I) * e$da + e$I0 * e$p0
quadrature <- function(y) e$da * sum(y)
close_to(e$Lambda_history[1, ], e$eps, 'initial Lambda')
close_to(mass, e$I0, 'initial infected mass')
close_to(e$I_total[1], mass, 'reported initial mass')
cohort <- numeric(e$A)
cohort[1] <- 1 / e$da
for (step in seq_len(e$A)) {
  close_to(quadrature(cohort), 1, 'transported unit mass')
  cohort <- c(0, head(cohort, -1))
}

# Infection-free S-T decay must be exp(-delta*t); also test delta=0.
for (decay_rate in c(1, 0)) {
  case <- initialize(list(dt = 1, da = 1, t_max = 2, delta = decay_rate))
  e <- case$env
  e$I[] <- 0; e$SI[] <- 0; e$II[] <- 0; e$IT[] <- 0
  e$I0 <- 0
  e$ST <- outer(e$ks * e$S, e$q * 0.05)
  e$SS <- outer(e$ks * e$S, e$q * 0.95)
  initial <- e$ST
  e <- run(case)
  close_to(e$ST, initial * exp(-decay_rate * tail(e$times, 1)),
           'infection-free ST decay')
}

# With exposure removed and only recovery active, the initial cohort
# must follow I0*exp(-gamma*t), without age-grid oscillations.
case <- initialize(list(alpha = 0, gamma = 0.25, p = 0))
e <- case$env
e$SI[] <- 0; e$ST[] <- 0; e$II[] <- 0; e$IT[] <- 0
e <- run(case)
close_to(e$I_total, e$I0 * exp(-e$gamma * e$times), 'recovering seed cohort')

# A short nonlinear refinement check (not a full epidemic validation).
endpoints <- numeric(3)
for (i in seq_along(endpoints)) {
  h <- c(0.4, 0.2, 0.1)[i]
  e <- run(initialize(list(dt = h, da = h)))
  fields <- c('S', 'I', 'SS', 'ST', 'SI', 'IT', 'II')
  stopifnot(all(vapply(fields, function(name) {
    x <- e[[name]]
    all(is.finite(x)) && all(x >= 0)
  }, logical(1))))
  endpoints[i] <- tail(e$I_total, 1)
}
stopifnot(abs(endpoints[3] - endpoints[2]) < abs(endpoints[2] - endpoints[1]))
cat(sprintf('PASS pde.R: mass, normalization, decay, recovery; I(4) at h=.4/.2/.1 = %s\n',
            paste(sprintf('%.8f', endpoints), collapse = ', ')))
cat('All PDE consistency checks passed.\n')
