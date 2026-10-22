# Run: Rscript tests/check_theta_full_pde.R --dt=0.1
# Reuses pde.R's unmodified initialization and full nonlinear time loop.
# A tiny deterministic seed separates invasion growth from susceptible loss.
# Writes diagnostic trajectories under output/theta-profiles/full-pde;
# existing epidemic CSVs are never overwritten.
script_arg <- grep('^--file=', commandArgs(), value=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1]))))
args <- commandArgs(trailingOnly=TRUE)
value <- function(name, default) {
  option <- args[startsWith(args,paste0('--',name,'='))]
  if (length(option)) as.numeric(sub('^[^=]*=', '', tail(option,1))) else default
}
h <- value('dt', .1)
duration <- value('tmax', 30)
seed_fraction <- value('seed-fraction', 1e-11)
stopifnot(h > 0, duration >= 30, seed_fraction > 0, seed_fraction < .01)
config <- list(engine='cpp', dt=h, da=h, t_max=duration,
               a_max=duration+2*h, Kmax=20, mean_degree=5,
               N=100000, I0=100000*seed_fraction)
expressions <- parse(file.path(root,'pde.R'))
loop <- which(vapply(expressions,function(x) is.call(x) &&
  identical(x[[1]],as.name('for')) && identical(x[[2]],as.name('n')),logical(1)))
stopifnot(length(loop)==1)
e <- new.env(parent=globalenv())
for(expr in expressions[seq_len(loop-1L)]) {
  if(is.call(expr) && identical(expr[[1]],as.name('<-')) &&
     is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% names(config))
    expr[[3]] <- config[[as.character(expr[[2]])]]
  eval(expr,e)
}
eval(expressions[[loop]],e)

trajectory <- data.frame(t=e$times, I_PDE=e$I_total, S_PDE=e$S_total,
                        Lambda=rowMeans(e$Lambda_history),
                        Lambda_relative_degree_spread=apply(e$Lambda_history,1,
                          function(x) (max(x)-min(x))/mean(x)))
fits <- do.call(rbind,lapply(list(c(10,15),c(15,20),c(20,25),c(25,30),c(20,30)),
  function(window) {
    z <- trajectory[trajectory$t >= window[1]-1e-9 & trajectory$t <= window[2]+1e-9,]
    I_fit <- lm(log(I_PDE)~t,data=z)
    L_fit <- lm(log(Lambda)~t,data=z)
    data.frame(dt=h, start=window[1], end=window[2],
               growth_I=unname(coef(I_fit)[2]), growth_Lambda=unname(coef(L_fit)[2]),
               max_log_I_fit_residual=max(abs(resid(I_fit))),
               susceptible_loss=1-tail(z$S_PDE,1)/e$N)
  }))
stopifnot(e$age_outflow==0, all(is.finite(as.matrix(trajectory))),
          max(fits$susceptible_loss) < 1e-5)
PT <- e$conditional_PT(e$IT,e$I)
profiles <- expand.grid(degree=e$ks,infection_age=e$ages)
profiles$theta <- as.vector(PT)
profiles <- profiles[profiles$infection_age <= duration-h,]
directory <- file.path(root,'output','theta-profiles','full-pde')
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
suffix <- paste0('dt',format(h,trim=TRUE,scientific=FALSE))
if (seed_fraction != 1e-11) suffix <- paste0(suffix,'-seed',format(seed_fraction,scientific=TRUE))
write.csv(trajectory,file.path(directory,paste0('trajectory-',suffix,'.csv')),row.names=FALSE)
write.csv(fits,file.path(directory,paste0('growth-',suffix,'.csv')),row.names=FALSE)
write.csv(profiles,file.path(directory,paste0('theta-',suffix,'.csv')),row.names=FALSE)
cat('\nFull nonlinear PDE invasion-growth estimates:\n')
print(fits,digits=10,row.names=FALSE)
cat('Initial fraction:',seed_fraction,'; age-cutoff loss:',e$age_outflow,'\n')
