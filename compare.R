# Compare saved simulations without loading ABM or rerunning either model.
# Rscript compare.R [--log] [--output=comparison.png] [CSV_FILE_OR_DIRECTORY ...]
# --log uses a logarithmic y axis and defaults to comparison-log.png.
# With no paths, search the working directory, this script's directory, and
# its parent. Pass paths explicitly to compare a selected set of results.

find_result_csvs <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Cannot find: ", paste(missing, collapse = ", "))
  files <- unlist(lapply(paths, function(path) {
    if (dir.exists(path)) list.files(path, pattern = "\\.csv$", full.names = TRUE)
    else path
  }), use.names = FALSE)
  files <- files[grepl("^(It_|PDE_).*\\.csv$", basename(files)) &
                   !grepl("[-_]summary\\.csv$", basename(files), ignore.case = TRUE)]
  sort(unique(normalizePath(files, mustWork = TRUE)))
}

read_result_series <- function(file) {
  tryCatch({
    data <- read.csv(file)
    if (!"t" %in% names(data) || !is.numeric(data$t) || !nrow(data) ||
        any(!is.finite(data$t)) || any(diff(data$t) <= 0)) {
      stop("expected a finite, increasing numeric t column")
    }
    stem <- tools::file_path_sans_ext(basename(file))
    if (startsWith(stem, "PDE_")) {
      if (!"I_PDE" %in% names(data) || !is.numeric(data$I_PDE)) {
        stop("expected a numeric I_PDE column")
      }
      values <- data$I_PDE
      label <- if (startsWith(stem, "PDE_contact_tracing_dt")) {
        paste0("PDE (dt = ", sub("^PDE_contact_tracing_dt", "", stem), ")")
      } else gsub("_", " ", stem)
      kind <- "PDE"
    } else {
      # Summary statistics are not independent simulation trajectories.
      if ("mean_I" %in% names(data)) return(NULL)
      replicates <- data[, names(data) != "t", drop = FALSE]
      if (!ncol(replicates) || !all(vapply(replicates, is.numeric, logical(1))) ||
          any(!is.finite(as.matrix(replicates))) || any(as.matrix(replicates) < 0)) {
        stop("expected finite, nonnegative numeric ABM replicate columns")
      }
      values <- rowMeans(replicates)
      name <- sub("-[0-9]+$", "", sub("^It_", "", stem))
      label <- sprintf("ABM mean (%d runs): %s", ncol(replicates), gsub("_", " ", name))
      kind <- "ABM"
    }
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("infectious counts must be finite and nonnegative")
    }
    list(file = file, kind = kind, label = label, t = data$t, I = values)
  }, error = function(error) {
    warning("Skipping ", basename(file), ": ", conditionMessage(error), call. = FALSE)
    NULL
  })
}

compare_results <- function(paths = ".", output = NULL, log_y = FALSE) {
  if (is.null(output)) output <- if (log_y) "comparison-log.png" else "comparison.png"
  files <- find_result_csvs(paths)
  series <- Filter(Negate(is.null), lapply(files, read_result_series))
  if (!length(series)) stop("No usable ABM or PDE result CSVs found.", call. = FALSE)
  labels <- vapply(series, `[[`, character(1), "label")
  repeated <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[repeated] <- paste0(labels[repeated], " [",
                            vapply(series[repeated], function(x) dirname(x$file), character(1)), "]")

  counts <- unlist(lapply(series, `[[`, "I"))
  if (log_y) {
    positive <- counts[counts > 0]
    if (!length(positive)) stop("Log scale requires at least one positive infectious count.",
                                call. = FALSE)
    limits <- range(positive)
    limits <- limits * if (diff(limits) == 0) c(0.5, 2) else c(1 / 1.06, 1.06)
  } else limits <- c(0, max(1, max(counts) * 1.06))

  extension <- tolower(tools::file_ext(output))
  if (extension == "png") {
    grDevices::png(output, width = 1800, height = 1100, res = 180)
  } else if (extension == "pdf") {
    grDevices::pdf(output, width = 10, height = 6.1)
  } else stop("output must have a .png or .pdf extension.", call. = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  colors <- grDevices::hcl.colors(length(series), "Dark 3")
  styles <- ifelse(vapply(series, `[[`, character(1), "kind") == "ABM", 1L, 2L)
  time_range <- range(unlist(lapply(series, `[[`, "t")))
  if (diff(time_range) == 0) time_range <- time_range + c(-0.5, 0.5)
  graphics::plot(time_range, limits, type = "n", log = if (log_y) "y" else "",
                 xlab = "Time",
                 ylab = if (log_y) "Number infectious (log scale)" else "Number infectious",
                 main = "Contact tracing: ABM and PDE simulations")
  for (i in seq_along(series)) {
    values <- series[[i]]$I
    # Keep gaps at zeros rather than inventing a positive pseudocount.
    if (log_y) values[values == 0] <- NA_real_
    graphics::lines(series[[i]]$t, values,
                    col = colors[i], lty = styles[i], lwd = 2)
  }
  if (log_y && any(counts == 0)) {
    graphics::mtext("Zero counts omitted on logarithmic axis", side = 3, line = 0.3, cex = 0.7)
  }
  graphics::legend("topright", legend = labels, col = colors, lty = styles,
                   lwd = 2, bty = "n", cex = 0.8)
  cat("Included CSVs:\n", paste0("  ", vapply(series, `[[`, character(1), "file"),
                                  collapse = "\n"), "\n", sep = "")
  cat("Saved:", normalizePath(output, mustWork = FALSE), "\n")
  invisible(series)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat("Usage: Rscript compare.R [--log] [--output=FILE.png|FILE.pdf] [CSV_FILE_OR_DIRECTORY ...]\n")
  } else {
    output_arg <- grep("^--output=", args, value = TRUE)
    if (length(output_arg) > 1L) stop("Specify --output only once.")
    log_y <- "--log" %in% args
    output <- if (length(output_arg)) sub("^--output=", "", output_arg) else NULL
    paths <- args[!startsWith(args, "--output=") & args != "--log"]
    if (any(startsWith(paths, "--"))) stop("Unknown option; use --help for usage.")
    if (!length(paths)) {
      script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
      script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
      paths <- unique(c(getwd(), script_dir, dirname(script_dir)))
    }
    compare_results(paths, output, log_y = log_y)
  }
}
