# Run: Rscript tests/check_pde_engines.R [--benchmark]
# Check the compiled transport kernel against the retained R implementation.
script_arg <- grep('^--file=', commandArgs(), value=TRUE)
root <- if (length(script_arg)) dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1])))) else normalizePath('.')
Rcpp::sourceCpp(file.path(root, 'pde_ii.cpp'), showOutput=FALSE)
close_to <- function(x, y, label) {
  if (!identical(dim(x),dim(y)) || any(!is.finite(x)) ||
      max(abs(x-y))/max(1,max(abs(y))) > 1e-11) stop(label)
}

# Dense, asymmetric, nonconstant states exercise all ages and degree factors,
# independently of the special symmetry and empty cells of epidemic seeds.
reference_pair <- function(II,SI,PT,Lambda,SI_in,survival,common,beta,alpha,gamma,delta_p,dt,da) {
  K <- dim(II)[1]; A <- dim(II)[3]; old <- seq_len(A-1); next_age <- old+1
  out <- array(0,dim(II)); W <- Q <- array(0,c(K,K,A))
  for(k in seq_len(K)) for(j in seq_len(K)) {
    pair <- II[k,j,,]
    W[k,j,] <- rowSums(pair)*da
    Q[k,j,] <- as.vector(pair %*% PT[j,])*da
    out[k,j,next_age,next_age] <- pair[old,old]*(common^2)*outer(survival[k,old]^2,survival[j,old]^2)
    loss <- beta+alpha+gamma+delta_p*(j-1)*PT[j,]+beta*(k-1)*Lambda[k]
    reverse_loss <- beta+alpha+gamma+delta_p*(k-1)*PT[k,]+beta*(j-1)*Lambda[j]
    fk <- beta*(1+(k-1)*Lambda[k]); fj <- beta*(1+(j-1)*Lambda[j])
    out[k,j,1,next_age] <- fk*SI[k,j,old]*exp(-loss[old]*dt/2)*common*survival[k,1]*survival[j,old]
    out[k,j,next_age,1] <- fj*SI[j,k,old]*exp(-reverse_loss[old]*dt/2)*common*survival[k,old]*survival[j,1]
    out[k,j,1,1] <- .5*(fk*SI_in[k,j]+fj*SI_in[j,k])*common*survival[k,1]*survival[j,1]
  }
  list(II=out,W=W,Q=Q)
}
set.seed(2026)
K <- 4L; A <- 7L; beta <- .1; alpha <- .15; gamma <- .1; dp <- .3; h <- .2
II <- array(runif(K*K*A*A),c(K,K,A,A)); SI <- array(runif(K*K*A),c(K,K,A))
PT <- matrix(runif(K*A,0,.2),K,A); Lambda <- runif(K,0,.1)
SI_in <- matrix(runif(K*K),K,K)
survival <- exp(-dp*((seq_len(K)-1)*PT)*h/2); common <- exp(-(alpha+gamma)*h)
args <- list(II,SI,PT,Lambda,SI_in,survival,common,beta,alpha,gamma,dp,h,h)
r <- do.call(reference_pair,args); c <- do.call(pde_ii_step,c(args,list(A)))
for(name in names(r)) close_to(c[[name]],r[[name]],paste('dense kernel',name))
close_to(pde_ii_symmetry(II),max(abs(II-aperm(II,c(2,1,4,3))))/max(1,max(II)), 'symmetry diagnostic')
cat('PASS dense kernel, both incoming faces, corner, weighted integrals and symmetry.\n')

initialize <- function(engine, config=list()) {
  overrides <- modifyList(list(engine=engine,dt=.2,da=.2,t_max=2,a_max=2,Kmax=4,mean_degree=2),config)
  exprs <- parse(file.path(root,'pde.R'))
  loop <- which(vapply(exprs,function(x) is.call(x) && identical(x[[1]],as.name('for')) && identical(x[[2]],as.name('n')),logical(1)))
  stopifnot(length(loop)==1)
  e <- new.env(parent=globalenv())
  invisible(capture.output(for(expr in exprs[seq_len(loop-1)]) {
    if(is.call(expr) && identical(expr[[1]],as.name('<-')) && is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names(overrides)) expr[[3]] <- overrides[[as.character(expr[[2]])]]
    eval(expr,e)
  }))
  list(env=e,loop=exprs[[loop]])
}
run <- function(case) {
  invisible(capture.output(eval(case$loop,case$env)))
  case$env
}
fields <- c('S','I','SS','ST','SI','IT','II','I_total','S_total','Lambda_history','age_outflow')
for(config in list(list(),list(p=0),list(p=1),list(delta=0),list(alpha=0))) {
  r <- run(initialize('R',config)); c <- run(initialize('cpp',config))
  for(name in fields) close_to(c[[name]],r[[name]],paste('trajectory',name))
}
cat('PASS complete trajectories for default rates, p=0/1, delta=0 and alpha=0, including age-cutoff outflow.\n')

if ('--benchmark' %in% commandArgs(trailingOnly=TRUE)) {
  # Timings exclude compilation and initialization. Keep the benchmark small
  # enough to coexist with a user's fine-grid simulation.
  config <- list(Kmax=20,mean_degree=5,dt=.1,da=.1,t_max=4,a_max=15)
  elapsed <- allocated <- numeric(2)
  for(i in 1:2) {
    case <- initialize(c('R','cpp')[i],config)
    profile <- tempfile()
    gc()
    Rprofmem(profile)
    elapsed[i] <- system.time(answer <- run(case))['elapsed']
    Rprofmem(NULL)
    sizes <- suppressWarnings(as.numeric(sub(' .*','',readLines(profile))))
    allocated[i] <- sum(sizes,na.rm=TRUE)/1e9
    unlink(profile)
    if(i==1) reference <- answer$I_total else close_to(answer$I_total,reference,'benchmark trajectory')
    rm(answer,case); gc()
  }
  cat(sprintf('40-step benchmark, K=20, A=150: R %.3fs, cpp %.3fs, speedup %.2fx.\n',elapsed[1],elapsed[2],elapsed[1]/elapsed[2]))
  cat(sprintf('Cumulative R allocations: R %.3f GB, cpp %.3f GB (not peak memory).\n',allocated[1],allocated[2]))
  # Also benchmark a fully occupied two-age array, rather than relying on
  # the empty high-age cells present early in a simulation.
  K <- 20; A <- 150
  II <- array(runif(K*K*A*A),c(K,K,A,A)); SI <- array(runif(K*K*A),c(K,K,A))
  PT <- matrix(runif(K*A,0,.2),K,A); Lambda <- runif(K,0,.1); SI_in <- matrix(runif(K*K),K,K)
  survival <- exp(-dp*((seq_len(K)-1)*PT)*h/2)
  args <- list(II,SI,PT,Lambda,SI_in,survival,common,beta,alpha,gamma,dp,h,h)
  a <- system.time(r <- do.call(reference_pair,args))['elapsed']
  b <- system.time(c <- do.call(pde_ii_step,c(args,list(A)) ))['elapsed']
  for(name in names(r)) close_to(c[[name]],r[[name]],paste('dense benchmark',name))
  cat(sprintf('Fully occupied II kernel: R %.3fs, cpp %.3fs, speedup %.2fx.\n',a,b,a/b))
}
