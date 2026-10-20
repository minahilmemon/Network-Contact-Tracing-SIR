# Age-structured network contact-tracing model: quick exploratory solver.
# Uses the corrected edge equations and derives conditional probabilities
# from their defining ratios. No independent Lambda/R evolution, floors,
# probability rescaling, or Simpson quadrature is used.
#
# Cell-centred, first-order characteristic/inflow discretization; dt = da.
# Cell integrals are da * sum(density). Boundary inflows are evaluated
# at frozen midpoint estimates. Newly infected cohorts receive half a
# step of reaction; existing cohorts receive a full step.
# Refine dt=da before using numerical values in the manuscript.

beta <- 0.10
alpha <- 0.15
gamma <- 0.10
delta <- 1.0
p <- 0.3
delta_p <- delta * p
N <- 100000
I0 <- 20
mean_degree <- 5
Kmax <- 20

dt <- 0.1
da <- 0.1
t_max <- 200
a_max <- 60

# The compiled engine avoids copying two-age slices and large temporary arrays.
# Use --engine=R for the original numerical implementation as a reference.
engine <- "auto"
# Example: Rscript pde_new.R --dt=0.1 --tmax=3 --amax=12 --engine=cpp
args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (startsWith(arg, '--tmax=')) t_max <- as.numeric(sub('--tmax=', '', arg))
  if (startsWith(arg, '--amax=')) a_max <- as.numeric(sub('--amax=', '', arg))
  if (startsWith(arg, '--dt=')) dt <- da <- as.numeric(sub('--dt=', '', arg))
  if (startsWith(arg, '--engine=')) engine <- sub('--engine=', '', arg)
}
stopifnot(dt > 0, abs(dt-da) < 1e-12, a_max >= 2*da,
          t_max > 0, beta > 0, alpha >= 0, gamma >= 0,
          delta >= 0, p >= 0, p <= 1, I0 > 0, I0 < N)
stopifnot(engine %in% c("auto", "R", "cpp"))
use_cpp <- engine == "cpp" || (engine == "auto" && requireNamespace("Rcpp", quietly=TRUE))
if (use_cpp) {
  if (!requireNamespace("Rcpp", quietly=TRUE)) stop("The cpp engine requires Rcpp.")
  # Resolve the companion source for Rscript, source(), or the test harness.
  script_files <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
  source_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  candidate_dirs <- unique(c(dirname(c(script_files, source_files)),
                             dirname(dirname(script_files)), getwd()))
  candidates <- file.path(candidate_dirs, "pde_ii.cpp")
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Cannot find pde_ii.cpp beside the solver.")
  cpp <- new.env()
  Rcpp::sourceCpp(normalizePath(candidates[1]), env=cpp, showOutput=FALSE)
}
A <- ceiling(a_max/da)
a_max <- A*da
ages <- ((seq_len(A))-0.5)*da
steps <- ceiling(t_max/dt)
times <- (0:steps)*dt
K <- Kmax
ks <- seq_len(K)
old <- seq_len(A-1L)
next_age <- old+1L
pk_all <- dpois(0:K, mean_degree)
pk_all <- pk_all/sum(pk_all)
p0 <- pk_all[1]
pk <- pk_all[-1]
q <- ks*pk/sum(ks*pk)
eps <- I0/N
S <- N*pk*(1-eps)
I <- matrix(0, K, A)
# Represent the seed cohort by its mass in the first age cell [0, da].
# This is a finite-volume approximation, not a nodal delta integrated
# with Simpson weights. Its initial integral is exactly I0 * sum(pk).
I[,1] <- I0*pk/da
# Edge densities: SI[k,j,a], IT[k,j,a], II[k,j,a1,a2].
SS <- outer(ks*S, q*(1-eps))
ST <- matrix(0,K,K)
SI <- array(0,c(K,K,A))
IT <- array(0,c(K,K,A))
II <- array(0,c(K,K,A,A))
SI[,,1] <- outer(ks*S, q*eps/da)
II[,,1,1] <- outer(ks*I[,1], q*eps/da)

# Ratios are only evaluated on occupied conditioning cells. The zero
# value outside their support is a computational extension, not a DFE
# assignment; products involving those empty edge cells must also vanish.
conditional_PT <- function(IT, I) {
  counts <- apply(IT,c(1,3),sum)
  denom <- ks*I
  out <- matrix(0,K,A)
  occupied <- denom > 0
  out[occupied] <- counts[occupied]/denom[occupied]
  if (any(counts[!occupied] > 1e-9))
    stop('Traced edges exist without an infected conditioning cohort.')
  out
}
exponential_source <- function(x, loss, source, h) {
  z <- loss*h
  phi <- rep(h,length(z))
  nz <- abs(z) > 1e-8
  phi[nz] <- -expm1(-z[nz])/loss[nz]
  x*exp(-z)+source*phi
}

I_total <- numeric(steps+1)
S_total <- numeric(steps+1)
I_total[1] <- sum(I)*da+I0*p0
S_total[1] <- sum(S)+N*p0*(1-eps)
Lambda_history <- matrix(0,steps+1,K)
Lambda_history[1,] <- apply(SI,1,sum)*da/(ks*S)
max_PT <- 0
max_S_neighbours <- 0
symmetry_SS <- 0
symmetry_II <- 0
age_outflow <- 0
started <- proc.time()[3]
cat(sprintf('Running dt=da=%.2f, %d age cells, %d degree classes, %d steps (%s engine).\n',
            dt,A,K,steps,if (use_cpp) 'cpp' else 'R'))

for (n in seq_len(steps)) {
  PT <- conditional_PT(IT,I)
  max_PT <- max(max_PT,PT)
  if (max(PT)>1+1e-6)
    stop(sprintf('Conditional traced-neighbour probability exceeds 1 at t=%.3f; refine the grid or investigate closure consistency.',times[n]))
  exposure <- apply(SI,1,sum)*da
  Lambda <- exposure/(ks*S)
  if (any(!is.finite(Lambda)) || any(Lambda < 0) || any(Lambda>1+1e-6))
    stop('Invalid susceptible infected-neighbour probability.')
  susceptible_neighbours <- (rowSums(SS)+rowSums(ST)+exposure)/(ks*S)
  max_S_neighbours <- max(max_S_neighbours,susceptible_neighbours)
  # Missing probability corresponds to recovered/removed neighbours.
  if (max(susceptible_neighbours)>1+0.02)
    stop('Susceptible conditional neighbour mass exceeds 1 by over 0.02; refine grid or investigate closure.')

  hazard_S <- beta*ks*Lambda
  S_new <- S*exp(-hazard_S*dt)
  # The incoming node mass equals the susceptible loss exactly.
  incidence_avg <- -S*expm1(-hazard_S*dt)/dt
  node_loss <- alpha+gamma+delta_p*(ks*PT)
  age_outflow <- age_outflow+sum(I[,A])*da
  I_new <- matrix(0,K,A)
  I_new[,next_age] <- I[,old,drop=FALSE]*exp(-node_loss[,old,drop=FALSE]*dt)
  I_new[,1] <- incidence_avg*exp(-node_loss[,1]*dt/2)

  SS_loss <- beta*(outer((ks-1)*Lambda,rep(1,K))+
                    outer(rep(1,K),(ks-1)*Lambda))
  SS_mid <- SS*exp(-SS_loss*dt/2)
  SS_new <- SS*exp(-SS_loss*dt)
  SI_new <- array(0,c(K,K,A))
  IT_new <- array(0,c(K,K,A))
  if (!use_cpp) II_new <- array(0,c(K,K,A,A))
  ST_new <- matrix(0,K,K)
  # SI inflow: infection of the neighbour by another neighbour.
  SI_in <- sweep(SS_mid,2,beta*(ks-1)*Lambda,'*')

  ii_surv <- exp(-delta_p*((ks-1)*PT)*dt/2)
  ii_common <- exp(-(alpha+gamma)*dt)
  if (use_cpp) {
    # Initially only age cell 1 is occupied; unit-speed transport grows the
    # occupied range by one cell per step. No nonzero density is discarded.
    pair_update <- cpp$pde_ii_step(II,SI,PT,Lambda,SI_in,ii_surv,
                                  ii_common,beta,alpha,gamma,delta_p,dt,da,min(n,A))
    II_new <- pair_update$II
  }
  for (k in ks) {
    for (j in ks) {
      si_loss <- beta+alpha+gamma+delta_p*(j-1)*PT[j,]+
                 beta*(k-1)*Lambda[k]
      SI_new[k,j,next_age] <- SI[k,j,old]*exp(-si_loss[old]*dt)
      SI_new[k,j,1] <- SI_in[k,j]*exp(-si_loss[1]*dt/2)

      # W and Q are integrals of edge densities here; division by
      # k I is implicit when reconstructing probability equations.
      if (use_cpp) {
        W_edges <- pair_update$W[k,j,]
        Q_edges <- pair_update$Q[k,j,]
      } else {
        pair <- II[k,j,,]
        W_edges <- rowSums(pair)*da
        Q_edges <- as.vector(pair %*% PT[j,])*da
      }
      it_source <- alpha*W_edges+delta_p*(j-1)*Q_edges
      it_loss <- delta+alpha+gamma+delta_p*(k-1)*PT[k,]
      IT_new[k,j,next_age] <- exponential_source(
        IT[k,j,old],it_loss[old],it_source[old],dt)
      # ST source integrated over the neighbour's infection age.
      st_source <- alpha*sum(SI[k,j,])*da+
                   delta_p*(j-1)*sum(SI[k,j,]*PT[j,])*da
      # The traced neighbour triggers its episode and enters X at rate delta.
      st_loss <- delta+beta*(k-1)*Lambda[k]
      ST_new[k,j] <- exponential_source(ST[k,j],st_loss,st_source,dt)

      if (!use_cpp) {
        II_new[k,j,next_age,next_age] <- pair[old,old]*
          (ii_common^2)*outer(ii_surv[k,old]^2,ii_surv[j,old]^2)
        # Midpoint estimates for the two incoming face strips. EIS[k,j]
        # equals SI[j,k] by undirected edge symmetry.
        si_reverse_loss <- beta+alpha+gamma+delta_p*(k-1)*PT[k,]+
                           beta*(j-1)*Lambda[j]
        factor_k <- beta*(1+(k-1)*Lambda[k])
        factor_j <- beta*(1+(j-1)*Lambda[j])
        II_new[k,j,1,next_age] <- factor_k*SI[k,j,old]*
          exp(-si_loss[old]*dt/2)*ii_common*ii_surv[k,1]*ii_surv[j,old]
        II_new[k,j,next_age,1] <- factor_j*SI[j,k,old]*
          exp(-si_reverse_loss[old]*dt/2)*ii_common*ii_surv[k,old]*ii_surv[j,1]
        # The corner is a cell average, not a shared boundary node.
        # Each face occupies half the corner cell under unit-speed
        # diagonal transport; both contributions are included once.
        II_new[k,j,1,1] <- 0.5*(factor_k*SI_in[k,j]+
                                 factor_j*SI_in[j,k])*
                            ii_common*ii_surv[k,1]*ii_surv[j,1]
      }
    }
  }
  # Boundary IT comes from infection of the S endpoint of ST.
  # Node incidence is used consistently to avoid dividing by Lambda.
  ST_mid <- (ST+ST_new)/2
  for (k in ks) {
    IT_new[k,,1] <- ((k-1)/k)*(ST_mid[k,]/S[k])*
                   incidence_avg[k]*exp(-node_loss[k,1]*dt/2)
  }

  state_arrays <- list(S_new,I_new,SS_new,ST_new,SI_new,IT_new)
  # The cpp kernel validates II while updating it, without allocating huge
  # temporary logical arrays for a second scan of the four-dimensional state.
  if (!use_cpp) state_arrays <- c(state_arrays,list(II_new))
  if (any(vapply(state_arrays,function(x) any(!is.finite(x)) || any(x<0),logical(1))))
    stop(sprintf('Nonfinite or negative state at t=%.3f.',times[n+1]))
  S <- S_new; I <- I_new; SS <- SS_new; ST <- ST_new
  SI <- SI_new; IT <- IT_new; II <- II_new
  I_total[n+1] <- sum(I)*da+I0*p0*exp(-(alpha+gamma)*times[n+1])
  S_total[n+1] <- sum(S)+N*p0*(1-eps)
  Lambda_history[n+1,] <- apply(SI,1,sum)*da/(ks*S)
  if (n %% 25 == 0 || n == steps) {
    cat(sprintf('t=%6.1f  I=%10.3f  elapsed=%5.1fs\n',
                times[n+1],I_total[n+1],proc.time()[3]-started))
    flush.console()
  }
}

# Final diagnostics; no post-hoc probability projection is performed.
max_PT <- max(max_PT,conditional_PT(IT,I))
symmetry_SS <- max(abs(SS-t(SS)))/max(1,max(SS))
symmetry_II <- if (use_cpp) cpp$pde_ii_symmetry(II) else
  max(abs(II-aperm(II,c(2,1,4,3))))/max(1,max(II))
results <- data.frame(t=times,I_PDE=I_total,S_PDE=S_total)
# Outputs are saved next to the script when run with Rscript, otherwise
# in the current working directory (for example when sourced in RStudio).
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
out_dir <- if(length(script_arg)) dirname(normalizePath(sub('^--file=','',script_arg[1]))) else getwd()
prefix <- file.path(out_dir,paste0('PDE_contact_tracing_dt',format(dt,trim=TRUE,scientific=FALSE)))
write.csv(results,paste0(prefix,'.csv'),row.names=FALSE)
plot_result <- function() {
  plot(times,I_total,type='l',lwd=2,col='#166678',
       xlab='Time',ylab='Number infectious',
       main=sprintf('Contact-tracing PDE: exploratory grid dt = da = %.1f',dt))
  abm_file <- file.path(out_dir,'It_p030_all_neighbors.csv')
  if(file.exists(abm_file)) {
    abm <- read.csv(abm_file)
    if('t' %in% names(abm) && ncol(abm)>1) {
      values <- abm[,setdiff(names(abm),'t'),drop=FALSE]
      if(all(vapply(values,is.numeric,logical(1)))) {
        abm_mean <- rowMeans(values)
        plot(range(c(times,abm$t)),range(c(0,I_total,abm_mean),finite=TRUE),
             type='n',xlab='Time',ylab='Number infectious',main='PDE and supplied ABM data')
        lines(times,I_total,col='#166678',lwd=2)
        lines(abm$t,abm_mean,col='#D35D35',lwd=2,lty=2)
        legend('topright',c('PDE (coarse grid)','ABM mean'),
               col=c('#166678','#D35D35'),lty=c(1,2),lwd=2,bty='n')
      }
    }
  }
}
png(paste0(prefix,'.png'),width=1400,height=850,res=150)
plot_result()
dev.off()
if(interactive()) plot_result()
cat(sprintf('\nFinished in %.1f seconds. Initial I=%.6f; peak I=%.3f at t=%.1f.\n',
            proc.time()[3]-started,I_total[1],max(I_total),times[which.max(I_total)]))
cat(sprintf('Max traced-neighbour probability: %.6g\n',max_PT))
cat(sprintf('Max tracked susceptible neighbour mass: %.6g\n',max_S_neighbours))
cat(sprintf('Relative SS / II symmetry errors: %.3g / %.3g\n',symmetry_SS,symmetry_II))
cat(sprintf('Cumulative infected mass exiting age cutoff: %.6g\n',age_outflow))
cat('Saved:',paste0(prefix,'.csv'),'and',paste0(prefix,'.png'),'\n')
if(age_outflow>1) warning('Non-negligible infected age outflow: increase a_max.')
