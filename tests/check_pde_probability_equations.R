# Run: Rscript tests/check_pde_probability_equations.R
# Independent checks of manuscript equations (13)-(18), their inflow
# boundaries, and two limiting models. No CSVs or full epidemics are produced.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else normalizePath(".")

initialize <- function(config = list()) {
  overrides <- modifyList(list(engine = "R", dt = .2, da = .2,
    t_max = 4, a_max = 8, Kmax = 6, mean_degree = 2), config)
  expressions <- parse(file.path(root, "pde.R"))
  loop <- which(vapply(expressions, function(x) {
    is.call(x) && identical(x[[1]], as.name("for")) &&
      identical(x[[2]], as.name("n"))
  }, logical(1)))
  stopifnot(length(loop) == 1)
  e <- new.env(parent = globalenv())
  invisible(capture.output(for (expr in expressions[seq_len(loop - 1L)]) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names(overrides)) {
      expr[[3]] <- overrides[[as.character(expr[[2]])]]
    }
    eval(expr, e)
  }))
  list(e = e, loop = expressions[[loop]])
}
run <- function(case) {
  invisible(capture.output(eval(case$loop, case$e)))
  case$e
}

# Reconstruct every probability from its defining ratio, including IS
# from the reverse orientation of SI. These are not solver state variables.
probabilities <- function(e) {
  out <- list(
    SS = sweep(e$SS, 1, e$ks * e$S, "/"),
    ST = sweep(e$ST, 1, e$ks * e$S, "/"),
    SI = sweep(e$SI, 1, e$ks * e$S, "/"),
    IS = array(0, dim(e$SI)),
    IT = array(0, dim(e$IT)),
    II = array(0, dim(e$II)))
  for (k in e$ks) for (j in e$ks) {
    out$IS[k, j, ] <- e$SI[j, k, ] / (k * e$I[k, ])
    out$IT[k, j, ] <- e$IT[k, j, ] / (k * e$I[k, ])
    out$II[k, j, , ] <- sweep(e$II[k, j, , ], 1, k * e$I[k, ], "/")
  }
  out$T <- apply(out$IT, c(1, 3), sum)
  out$Lambda <- apply(out$SI, 1, sum) * e$da
  out
}

# Dense, nonconstant manufactured states exercise every degree factor and
# infection-age argument, instead of checking only a sparse seed cohort.
fixture <- function(h) {
  case <- initialize(list(dt = h, da = h, t_max = h, a_max = .4, Kmax = 3))
  e <- case$e
  e$S <- c(12, 14, 16)
  e$I <- outer(e$ks, e$ages, function(k, a) 2 + .2*k + .1*a)
  e$SS <- .3 * outer(e$ks, e$ks)
  e$ST <- .02 * outer(e$ks, e$ks)
  for (k in e$ks) for (j in e$ks) for (a in seq_len(e$A)) {
    e$SI[k, j, a] <- .1 * (k + j + e$ages[a])
    e$IT[k, j, a] <- k * e$I[k, a] * .003 * (k + j + e$ages[a])
    for (b in seq_len(e$A)) e$II[k, j, a, b] <-
      .02 * (k + j + e$ages[a] + e$ages[b])
  }
  case
}

generator_error <- function(h) {
  case <- fixture(h)
  e <- case$e
  z <- probabilities(e)
  beta <- e$beta; mu <- e$alpha + e$gamma; dp <- e$delta_p
  interior <- setNames(rep(0, 6), c("13:SI", "14:ST", "15:IS",
                                    "16:SS", "17:II", "18:IT"))
  boundary <- setNames(rep(0, 5), c("SI", "IS", "IT", "II_focal", "II_neighbour"))
  e <- run(case)
  zz <- probabilities(e)
  old <- seq_len(e$A - 1L)
  # Relative to a unit scale for small probabilities, to the reference
  # density for the singular 1/Lambda inflow.
  error <- function(actual, expected) max(abs(actual - expected)) / max(1, abs(expected))
  for (k in e$ks) for (j in e$ks) {
    rhs <- (-(beta + mu) - dp*(j-1)*z$T[j, ] + beta*z$Lambda[k]) * z$SI[k,j, ]
    interior[1] <- max(interior[1],
      error((zz$SI[k,j,old+1] - z$SI[k,j,old])/h, rhs[old]))
    rhs <- (beta*z$Lambda[k] - e$delta)*z$ST[k,j] +
      e$alpha*sum(z$SI[k,j, ])*h +
      dp*(j-1)*sum(z$SI[k,j, ]*z$T[j, ])*h
    interior[2] <- max(interior[2], error((zz$ST[k,j] - z$ST[k,j])/h, rhs))
    rhs <- (-beta - beta*(j-1)*z$Lambda[j] + dp*z$T[k, ]) * z$IS[k,j, ]
    interior[3] <- max(interior[3],
      error((zz$IS[k,j,old+1] - z$IS[k,j,old])/h, rhs[old]))
    rhs <- (-beta*(j-1)*z$Lambda[j] + beta*z$Lambda[k]) * z$SS[k,j]
    interior[4] <- max(interior[4], error((zz$SS[k,j] - z$SS[k,j])/h, rhs))
    rhs <- z$II[k,j, , ] * outer(z$T[k, ], z$T[j, ],
      function(tk, tj) -mu - dp*(j-1)*tj + dp*tk)
    interior[5] <- max(interior[5],
      error((zz$II[k,j,old+1,old+1] - z$II[k,j,old,old])/h, rhs[old,old]))
    W <- rowSums(z$II[k,j, , ])*h
    Q <- as.vector(z$II[k,j, , ] %*% z$T[j, ])*h
    rhs <- (-e$delta + dp*z$T[k, ])*z$IT[k,j, ] + e$alpha*W + dp*(j-1)*Q
    interior[6] <- max(interior[6],
      error((zz$IT[k,j,old+1] - z$IT[k,j,old])/h, rhs[old]))

    boundary[1] <- max(boundary[1],
      error(zz$SI[k,j,1], beta*(j-1)*z$SS[k,j]*z$Lambda[j]))
    boundary[2] <- max(boundary[2], error(zz$IS[k,j,1], (k-1)/k*z$SS[k,j]))
    boundary[3] <- max(boundary[3], error(zz$IT[k,j,1], (k-1)/k*z$ST[k,j]))
    boundary[4] <- max(boundary[4], error(zz$II[k,j,1,old+1],
      z$SI[k,j,old] * ((k-1)/k + 1/(k*z$Lambda[k]))))
    boundary[5] <- max(boundary[5], error(zz$II[k,j,old+1,1],
      beta*(1+(j-1)*z$Lambda[j])*z$IS[k,j,old]))
  }
  c(interior, boundary)
}
h <- c(.02, .01, .005)
errors <- vapply(h, generator_error, numeric(11))
stopifnot(all(is.finite(errors)), max(errors[,3]) < .01)
# A first-order method's characteristic derivative and boundary residuals
# should decrease at least proportionately to h, allowing rounding error.
stopifnot(all(errors[,2] <= .65*errors[,1] + 1e-11),
          all(errors[,3] <= .65*errors[,2] + 1e-11))
cat("PASS all six probability equations and five inflow boundaries under refinement.\n")
print(signif(errors, 4))

rk4 <- function(rhs, y, duration, h = .001) {
  steps <- ceiling(duration/h); h <- duration/steps
  for (step in seq_len(steps)) {
    a <- rhs(y); b <- rhs(y + h*a/2)
    c <- rhs(y + h*b/2); d <- rhs(y + h*c)
    y <- y + h*(a + 2*b + 2*c + d)/6
  }
  y
}

# With p=0, integrate the age equations analytically to get a finite ODE.
# This reference has no age grid, no two-age transport, and no boundary
# quadrature. Variables are node S/I and integrated SS/SI edge counts.
sir_reference <- function(e, duration, h = .001) {
  K <- e$K
  unpack <- function(y) list(S=y[seq_len(K)], I=y[K+seq_len(K)],
    SS=matrix(y[2*K+seq_len(K*K)], K, K),
    SI=matrix(y[2*K+K*K+seq_len(K*K)], K, K))
  initial <- c(e$S, rowSums(e$I)*e$da, e$SS, apply(e$SI,c(1,2),sum)*e$da)
  rhs <- function(y) {
    z <- unpack(y)
    exposure <- rowSums(z$SI)
    Lambda <- exposure/(e$ks*z$S)
    external <- e$beta*(e$ks-1)*Lambda
    dSS <- -outer(external, external, "+")*z$SS
    inflow <- sweep(z$SS, 2, external, "*")
    dSI <- inflow - sweep(z$SI, 1, e$beta+e$alpha+e$gamma+external, "*")
    c(-e$beta*exposure, e$beta*exposure-(e$alpha+e$gamma)*z$I, dSS, dSI)
  }
  z <- unpack(rk4(rhs, initial, duration, h))
  c(S=sum(z$S)+e$N*e$p0*(1-e$eps),
    I=sum(z$I)+e$I0*e$p0*exp(-(e$alpha+e$gamma)*duration))
}

# A degree-one graph consists of isolated vertices and disjoint edges.
# Tracing has no external triples, so these integrated equations are an
# exact Markov expectation for independently seeded edges, even for p>0.
degree_one_reference <- function(e, duration, h = .001) {
  initial <- c(S=sum(e$S), I=sum(e$I)*e$da, SI=sum(e$SI)*e$da,
    II=sum(e$II)*e$da^2, IT=sum(e$IT)*e$da)
  mu <- e$alpha+e$gamma
  rhs <- function(y) c(
    -e$beta*y["SI"],
    e$beta*y["SI"] - mu*y["I"] - e$delta_p*y["IT"],
    -(e$beta+mu)*y["SI"],
    2*e$beta*y["SI"] - 2*mu*y["II"],
    e$alpha*y["II"] - (e$delta+mu)*y["IT"])
  y <- rk4(rhs, initial, duration, h)
  c(S=unname(y["S"])+e$N*e$p0*(1-e$eps),
    I=unname(y["I"])+e$I0*e$p0*exp(-mu*duration))
}

for (reference_name in c("SIR without tracing", "degree-one graph with tracing")) {
  config <- if (reference_name == "SIR without tracing") list(p=0, I0=10000) else
    list(Kmax=1, p=.8, I0=10000)
  reference <- if (reference_name == "SIR without tracing") sir_reference else
    degree_one_reference
  case <- initialize(config)
  expected <- reference(case$e, 4)
  finer_reference <- reference(case$e, 4, h=.0005)
  stopifnot(max(abs(expected-finer_reference)) < 1e-6)
  residual <- vapply(c(.2,.1,.05), function(h) {
    e <- run(initialize(modifyList(config,list(dt=h,da=h))))
    actual <- c(S=tail(e$S_total,1), I=tail(e$I_total,1))
    abs(actual-expected)
  }, numeric(2))
  stopifnot(all(residual[,2] < .6*residual[,1]),
            all(residual[,3] < .6*residual[,2]))
  cat("PASS", reference_name, ": independent ODE reference and first-order convergence.\n")
  print(signif(residual,6))
}
cat("All probability-equation checks passed.\n")
