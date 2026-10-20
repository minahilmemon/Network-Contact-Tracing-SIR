# Run from any directory:
# Rscript run-abm-1million.R --workers=4
# --pilot runs only replicate 1; subsequent invocations resume checkpoints.
# --output-dir=DIR keeps outputs/checkpoints outside the source directory.
# --dry-run prints the setup without starting simulations or writing files.
args <- commandArgs(trailingOnly=TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript run-abm-1million.R [--workers=4] [--output-dir=DIR] [--pilot] [--dry-run]\n")
  quit(status=0)
}
stopifnot(all(grepl("^--(workers=|output-dir=|pilot$|dry-run$)",args)))
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=","",script_arg[1])))
output_arg <- grep("^--output-dir=",args,value=TRUE)
stopifnot(length(output_arg)<=1L)
workspace <- if(length(output_arg)) sub("^--output-dir=","",output_arg) else script_dir
source_file <- file.path(script_dir,"network-contact-tracing.R")
source(source_file)
workers_arg <- grep("^--workers=", args, value=TRUE)
workers <- if (length(workers_arg)) as.integer(sub("^--workers=", "", workers_arg[1])) else 4L
stopifnot(length(workers)==1L, is.finite(workers), workers>=1L)
pilot <- "--pilot" %in% args
parameters <- list(N=1000000L, mean_degree=5, beta=.1, gamma=.1,
                   tau=.15, theta=1, p=.3)
I0 <- 200L
times <- 0:200
runs <- 100L
seeds <- as.list(20261020L + seq_len(runs))
checkpoint_dir <- file.path(workspace, sprintf("abm-1million-I0%d-100-checkpoints", I0))
metadata_file <- file.path(checkpoint_dir, "metadata.rds")
metadata <- list(parameters=parameters, I0=I0, times=times, runs=runs, seeds=seeds,
                 ABM_version=as.character(utils::packageVersion("ABM")),
                 ABM_commit=utils::packageDescription("ABM")$RemoteSha,
                 source_md5=unname(tools::md5sum(source_file)))
if ("--dry-run" %in% args) {
  print(metadata)
  cat("Output directory:",workspace,"\nWorkers:",workers,"\n")
  quit(status=0)
}
dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
stem <- sprintf("It_p030_all_neighbors-N%d-I0_%d-%d",parameters$N,I0,runs)
if (file.exists(metadata_file)) {
  stopifnot(identical(readRDS(metadata_file), metadata))
} else {
  saveRDS(metadata, metadata_file)
  stopifnot(file.copy(source_file, file.path(checkpoint_dir, "model-source.R")))
}
checkpoint_file <- function(id) file.path(checkpoint_dir, sprintf("run-%03d.rds", id))
valid_result <- function(x) {
  is.numeric(x) && length(x)==length(times) && identical(names(x),as.character(times)) &&
    all(is.finite(x)) && all(x>=0) && all(x<=parameters$N) &&
    all(x==floor(x)) && x[1]==I0
}
completed <- vapply(seq_len(runs), function(id) {
  path <- checkpoint_file(id)
  if (!file.exists(path)) return(FALSE)
  saved <- readRDS(path)
  identical(saved$id,id) && identical(saved$seed,seeds[[id]]) && valid_result(saved$result)
}, logical(1))
pending <- which(!completed)
if (pilot) pending <- intersect(pending, 1L)
cat(sprintf("N=%d; I0=%d; %d/%d checkpoints complete; workers=%d; ABM=%s (%s)\n",
            parameters$N,I0,sum(completed),runs,workers,
            metadata$ABM_version,metadata$ABM_commit))
flush.console()
started <- proc.time()[3]
while (length(pending)) {
  ids <- head(pending, workers)
  batch_started <- proc.time()[3]
  results <- run_contact_tracing_replicates(parameters, I0, times,
    runs=length(ids), workers=workers, seeds=seeds[ids])
  for (i in seq_along(ids)) {
    id <- ids[i]
    result <- results[[i]]
    if (inherits(result,"error")) result$run <- id
    saved <- list(id=id, seed=seeds[[id]], result=result, finished=Sys.time())
    temp <- tempfile("completed-",tmpdir=checkpoint_dir)
    saveRDS(saved,temp)
    stopifnot(file.rename(temp,checkpoint_file(id)))
    if (valid_result(result)) completed[id] <- TRUE
  }
  failed <- ids[!completed[ids]]
  cat(sprintf("[%s] Completed %d/%d; batch %.1fs; elapsed %.1fs\n",
              format(Sys.time(),"%H:%M:%S"),sum(completed),runs,
              proc.time()[3]-batch_started,proc.time()[3]-started))
  flush.console()
  if (length(failed)) stop("Failed or invalid runs saved in checkpoints: ",paste(failed,collapse=", "))
  pending <- setdiff(pending,ids)
}
if (pilot) {
  cat("Pilot complete. Its checkpoint will be reused by the full batch.\n")
} else {
  stopifnot(all(completed))
  saved <- lapply(seq_len(runs), function(id) readRDS(checkpoint_file(id)))
  results <- lapply(saved, `[[`, "result")
  stopifnot(all(vapply(results,valid_result,logical(1))))
  trajectories <- do.call(cbind,results)
  colnames(trajectories) <- paste0("run_",seq_len(runs))
  data <- data.frame(t=times,trajectories)
  output <- file.path(workspace,paste0(stem,".csv"))
  temporary <- tempfile("validated-abm-",tmpdir=workspace,fileext=".csv")
  write.csv(data,temporary,row.names=FALSE)
  roundtrip <- read.csv(temporary)
  stopifnot(identical(roundtrip$t,times),ncol(roundtrip)==runs+1L,
            all(as.matrix(roundtrip[,-1])==trajectories))
  # Preserve the prior 100-run file outside compare.R's nonrecursive search.
  if (file.exists(output)) {
    archive <- file.path(workspace,"archives",paste0("before-1million-I0",I0,"-",format(Sys.time(),"%Y%m%d-%H%M%S")))
    dir.create(archive,recursive=TRUE)
    prior <- file.path(workspace,c(paste0(stem,c(".csv",".rds","-summary.csv")),
                                  "population-size-batches.csv"))
    prior <- prior[file.exists(prior)]
    stopifnot(all(file.copy(prior,file.path(archive,basename(prior)))))
    cat("Previous CSV and metadata saved in:",archive,"\n")
  }
  saveRDS(c(metadata,list(results=results,finished=Sys.time(),
                         session_info=sessionInfo())),
          file.path(workspace,paste0(stem,".rds")))
  summary <- data.frame(t=times,mean_I=rowMeans(trajectories),
                        sd_I=apply(trajectories,1,sd))
  summary$se_I <- summary$sd_I/sqrt(runs)
  write.csv(summary,file.path(workspace,paste0(stem,"-summary.csv")),row.names=FALSE)
  stopifnot(file.rename(temporary,output))
  batch_info <- data.frame(file=paste0(stem,".csv"),summary_file=paste0(stem,"-summary.csv"),
    N=parameters$N,I0=I0,initial_infectious_fraction=I0/parameters$N,runs=runs,
    mean_degree=parameters$mean_degree,beta=parameters$beta,gamma=parameters$gamma,
    tau=parameters$tau,theta=parameters$theta,p=parameters$p,
    first_time=min(times),last_time=max(times),observation_interval=1,
    first_seed=seeds[[1]],seed_increment=1,ABM_version=metadata$ABM_version,
    ABM_commit=if(is.null(metadata$ABM_commit)) NA_character_ else metadata$ABM_commit,
    model_source_md5=metadata$source_md5)
  manifest_file <- file.path(workspace,"population-size-batches.csv")
  if(file.exists(manifest_file)) {
    prior_info <- read.csv(manifest_file)
    batch_info <- rbind(prior_info[prior_info$file!=batch_info$file,,drop=FALSE],batch_info)
  }
  write.csv(batch_info,manifest_file,row.names=FALSE)
  cat("Replaced:",output,"\n")
  cat(sprintf("Mean peak: %.3f at t=%g; final mean: %.3f; extinct by t=200: %d/%d\n",
              max(summary$mean_I),times[which.max(summary$mean_I)],
              tail(summary$mean_I,1),sum(trajectories[nrow(trajectories),]==0),runs))
}
