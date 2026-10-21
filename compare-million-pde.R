# Rscript compare-million-pde.R [--output=figure.pdf]
# Compare saved fractions; shading is the central 95% of stochastic runs,
# computed separately at each observation time, not a confidence interval.
args <- commandArgs(trailingOnly=TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript compare-million-pde.R [--output=figure.pdf]\n")
  quit(status=0)
}
stopifnot(length(args)<=1L,all(startsWith(args,"--output=")))
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
root <- dirname(normalizePath(sub("^--file=","",script_arg[1])))
output <- if(length(args)) sub("^--output=","",args[1]) else
  file.path(root,"output","pdf","abm-million-pde-dt01.pdf")
stopifnot(tolower(tools::file_ext(output))=="pdf")

metadata <- read.csv(file.path(root,"population-size-batches.csv"))
metadata <- metadata[metadata$N==1000000 & metadata$I0==200,,drop=FALSE]
stopifnot(nrow(metadata)==1L,metadata$runs==100L)
abm <- read.csv(file.path(root,metadata$file))
trajectories <- as.matrix(abm[,-1])
pde <- read.csv(file.path(root,"PDE_contact_tracing_dt0.1.csv"))
pde_N <- 100000L # Population size used to generate this saved PDE solution.
stopifnot(identical(abm$t,0:200),ncol(trajectories)==metadata$runs,
          all(is.finite(trajectories)),all(trajectories>=0),
          all(trajectories<=metadata$N),all(trajectories==floor(trajectories)),
          all(trajectories[1,]==metadata$I0),
          all(is.finite(pde$t)),all(is.finite(pde$I_PDE)),all(pde$I_PDE>=0),
          pde$t[1]==0,max(abs(diff(pde$t)-.1))<1e-10,
          isTRUE(all.equal(pde$I_PDE[1]/pde_N,metadata$I0/metadata$N)))

# R's type-7 empirical quantiles interpolate adjacent order statistics.
quantiles <- apply(trajectories,1,quantile,probs=c(.025,.975),type=7,names=FALSE)
lower <- 100*quantiles[1,]/metadata$N
upper <- 100*quantiles[2,]/metadata$N
mean_abm <- 100*rowMeans(trajectories)/metadata$N
pde_fraction <- 100*pde$I_PDE/pde_N
stopifnot(all(lower<=upper))

dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
pdf(output,width=9,height=8,pointsize=10,useDingbats=FALSE,
    title="Million-agent ABM ensemble and PDE comparison")
par(mfrow=c(2,1),mar=c(3.6,4.6,1.8,.9),oma=c(4.3,0,2.0,0),
    mgp=c(2.9,.7,0),las=1,bty="l",family="Helvetica")
abm_color <- "#D55E00"
pde_color <- "#0072B2"
band_color <- "#D3D3D3"
for(log_y in c(FALSE,TRUE)) {
  top <- max(upper,pde_fraction)*1.12
  values <- c(lower,mean_abm,pde_fraction)
  bottom <- if(log_y) 10^floor(log10(min(values[values>0]))) else 0
  plot(range(abm$t),c(bottom,top),type="n",log=if(log_y) "y" else "",
       xlab="Time (days)",ylab="Infectious fraction (%)",
       main=if(log_y) "Logarithmic scale" else "Linear scale")
  # Zero quantiles are clipped at the lower plotting boundary on the log
  # panel; their actual value remains zero. No positive pseudocount is added.
  band_lower <- if(log_y) pmax(lower,bottom) else lower
  polygon(c(abm$t,rev(abm$t)),c(band_lower,rev(upper)),
          col=band_color,border=NA)
  lines(abm$t,mean_abm,col=abm_color,lwd=2)
  lines(pde$t,pde_fraction,col=pde_color,lwd=2,lty=2)
  legend("topright",bty="n",cex=.92,
         legend=c("ABM ensemble mean (100 runs)","PDE (dt = 0.1 days)",
                  "ABM 2.5th-97.5th percentiles"),
         col=c(abm_color,pde_color,band_color),lty=c(1,2,NA),lwd=c(2,2,NA),
         pch=c(NA,NA,15),pt.cex=1.6)
}
mtext("Contact tracing: million-agent ensemble and PDE",side=3,outer=TRUE,
      line=.35,font=2,cex=1.2)
mtext("ABM: N = 1,000,000; I0 = 200. Both models start with 0.020% infectious.",
      side=1,outer=TRUE,line=.9,cex=.83)
mtext("Grey band: pointwise empirical 95% quantile interval across the stochastic runs.",
      side=1,outer=TRUE,line=2.0,cex=.83)
mtext("Zero lower quantiles extend below the logarithmic plotting range (from day 172).",
      side=1,outer=TRUE,line=3.1,cex=.78)
invisible(dev.off())
cat("Saved:",normalizePath(output),"\n")
