# Compare completed ABM batches and PDE solutions using infectious fractions.
# No simulations are run by this script. Rscript compare-population-sizes.R
# Reads committed CSVs only; no local RDS files or ABM installation required.
# --output-dir=DIR selects where the PNG, PDF and comparison tables are saved.
args <- commandArgs(trailingOnly=TRUE)
if("--help" %in% args) {
  cat("Usage: Rscript compare-population-sizes.R [--output-dir=DIR]\n")
  quit(status=0)
}
stopifnot(length(args)<=1L,all(startsWith(args,"--output-dir=")))
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
workspace <- dirname(normalizePath(sub("^--file=","",script_arg[1])))
output_dir <- if(length(args)) sub("^--output-dir=","",args[1]) else workspace
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
read_batch <- function(stem,metadata) {
  data <- read.csv(file.path(workspace, paste0(stem, ".csv")))
  trajectories <- as.matrix(data[,-1])
  stopifnot(all(is.finite(trajectories)), all(trajectories>=0),
            identical(data$t, metadata$times),
            all(trajectories[1,]==metadata$I0))
  N <- metadata$parameters$N
  means <- rowMeans(trajectories)
  se <- apply(trajectories,1,sd)/sqrt(ncol(trajectories))
  version <- metadata$ABM_version
  if (is.null(version)) version <- metadata$session_info$loadedOnly$ABM$Version
  list(t=data$t, I=means, se=se, N=N, runs=ncol(trajectories),
       trajectories=trajectories, metadata=metadata, version=version)
}
batch_info <- read.csv(file.path(workspace,"population-size-batches.csv"))
batch_info <- batch_info[batch_info$N==1000000 & batch_info$I0==200,,drop=FALSE]
stopifnot(nrow(batch_info)==1L)
large_metadata <- list(
  parameters=as.list(batch_info[1,c("N","mean_degree","beta","gamma","tau","theta","p")]),
  I0=batch_info$I0,times=seq(batch_info$first_time,batch_info$last_time,
                           by=batch_info$observation_interval),
  ABM_version=batch_info$ABM_version)
# Historical 1000-run batch: same rates and degree distribution, N=100000,
# I0=20, ABM 0.6.99 (recorded in the original batch metadata).
baseline_metadata <- large_metadata
baseline_metadata$parameters$N <- 100000L
baseline_metadata$I0 <- 20L
baseline_metadata$ABM_version <- "0.6.99"
baseline <- read_batch("It_p030_all_neighbors-1000",baseline_metadata)
large <- read_batch(tools::file_path_sans_ext(batch_info$file),large_metadata)
stopifnot(identical(baseline$t,large$t),
          isTRUE(all.equal(baseline$metadata$I0/baseline$N,
                           large$metadata$I0/large$N)),
          isTRUE(all.equal(baseline$metadata$parameters[-1],large$metadata$parameters[-1])))
batches <- list(baseline,large)
# The saved PDE solutions were generated at N=100000 and I0=20 by pde.R.
pde <- lapply(c("0.1","0.4"),function(dt) {
  data <- read.csv(file.path(workspace,paste0("PDE_contact_tracing_dt",dt,".csv")))
  stopifnot(all(c("t","I_PDE") %in% names(data)),
            all(is.finite(data$t)),all(diff(data$t)>0),
            all(is.finite(data$I_PDE)),all(data$I_PDE>=0),
            data$t[1]==0,abs(data$I_PDE[1]-20)<1e-8,
            isTRUE(all.equal(20/100000,large$metadata$I0/large$N)))
  list(t=data$t,I=data$I_PDE,N=100000L,I0=20L,dt=dt)
})
colors <- c("#0072B2","#D55E00","#009E73","#666666")
labels <- c(vapply(batches,function(x) sprintf("ABM: N = %s, I0 = %d (%s runs)",
  format(x$N,big.mark=",",scientific=FALSE),x$metadata$I0,
  format(x$runs,big.mark=",",scientific=FALSE)),character(1)),
  vapply(pde,function(x) sprintf("PDE: dt = %s days",x$dt),character(1)))
styles <- c(1,1,2,3)
percent <- lapply(batches,function(x) 100*x$I/x$N)
pde_percent <- lapply(pde,function(x) 100*x$I/x$N)
upper <- lapply(batches,function(x) 100*(x$I+1.96*x$se)/x$N)
lower <- lapply(batches,function(x) 100*pmax(0,x$I-1.96*x$se)/x$N)
draw_comparison <- function() {
  par(mfrow=c(2,1),mar=c(3.8,4.5,2,1),oma=c(4.6,0,1.8,0),las=1)
  all_percent <- c(percent,pde_percent)
  all_times <- c(lapply(batches,`[[`,"t"),lapply(pde,`[[`,"t"))
  for(log_y in c(FALSE,TRUE)) {
    positive <- unlist(all_percent)
    positive <- positive[positive>0]
    limits <- if(log_y) range(positive)*c(.8,1.2) else
      c(0,max(unlist(c(upper,pde_percent)))*1.13)
    plot(range(large$t),limits,type="n",log=if(log_y) "y" else "",
         xlab="Time (days)",ylab="Infectious population (%)",
         main=if(log_y) "Logarithmic scale" else "Linear scale")
    # Bands are pointwise intervals for the ABM ensemble mean.
    for(i in seq_along(batches)) {
      lo <- lower[[i]]
      if(log_y) lo <- pmax(lo,limits[1])
      polygon(c(batches[[i]]$t,rev(batches[[i]]$t)),c(lo,rev(upper[[i]])),
              col=adjustcolor(colors[i],alpha.f=.15),border=NA)
    }
    for(i in seq_along(all_percent)) {
      y <- all_percent[[i]]
      if(log_y) y[y<=0] <- NA_real_
      lines(all_times[[i]],y,col=colors[i],lty=styles[i],lwd=2)
    }
    legend("topright",legend=labels,col=colors,lty=styles,lwd=2,bty="n",cex=.78)
  }
  mtext("Contact tracing: ABM population size and PDE comparison",
        side=3,outer=TRUE,line=.15,font=2,cex=1.1)
  mtext(sprintf("All curves start with %.3f%% infectious. Shading: approximate pointwise 95%% intervals for ABM means.",
                100*large$metadata$I0/large$N),side=1,outer=TRUE,line=1,cex=.73)
  mtext("Mean degree = 5; beta = 0.10, gamma = 0.10, tau = 0.15, theta = 1, p = 0.30.",
        side=1,outer=TRUE,line=2.1,cex=.73)
  mtext(sprintf("ABM versions: N = 100,000 used %s; N = 1,000,000 used %s.",
                baseline$version,large$version),side=1,outer=TRUE,line=3.2,cex=.73)
}
for(extension in c("png","pdf")) {
  outfile <- file.path(output_dir,paste0("population-size-comparison.",extension))
  if(extension=="png") png(outfile,width=2200,height=1800,res=220)
  else pdf(outfile,width=10,height=1800/220,useDingbats=FALSE)
  draw_comparison()
  invisible(dev.off())
  cat("Saved:",outfile,"\n")
}
comparison <- data.frame(t=large$t,
  mean_I_N100000=baseline$I,se_I_N100000=baseline$se,
  mean_I_N1000000=large$I,se_I_N1000000=large$se,
  fraction_N100000=baseline$I/baseline$N,
  fraction_N1000000=large$I/large$N)
write.csv(comparison,file.path(output_dir,"population-size-comparison.csv"),row.names=FALSE)
metrics <- do.call(rbind,lapply(batches,function(x) data.frame(
  N=x$N,I0=x$metadata$I0,runs=x$runs,ABM_version=x$version,
  mean_curve_peak=max(x$I),peak_fraction=max(x$I)/x$N,
  peak_time=x$t[which.max(x$I)],final_mean=tail(x$I,1),
  extinct_at_200=sum(x$trajectories[nrow(x$trajectories),]==0))))
print(metrics,row.names=FALSE)
pde_metrics <- do.call(rbind,lapply(pde,function(x) data.frame(
  dt=x$dt,N=x$N,I0=x$I0,peak_fraction=max(x$I)/x$N,
  peak_time=x$t[which.max(x$I)],final_mean=tail(x$I,1))))
print(pde_metrics,row.names=FALSE)
write.csv(metrics,file.path(output_dir,"population-size-comparison-ABM-metrics.csv"),row.names=FALSE)
write.csv(pde_metrics,file.path(output_dir,"population-size-comparison-PDE-metrics.csv"),row.names=FALSE)
