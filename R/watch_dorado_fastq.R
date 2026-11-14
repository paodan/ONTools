#' Demultiplex one classified BAM directly to FASTQ
#'
#' `demux_one_bam_fastq()` runs `dorado demux --no-classify --emit-fastq` on a
#' BAM whose reads were classified during basecalling. Each invocation should
#' use a directory dedicated to one input BAM.
#'
#' @param bam_file Classified BAM created from one POD5 file.
#' @param demux_dir Dedicated Dorado output directory for this BAM.
#' @param threads Optional positive integer passed to `dorado demux --threads`.
#' @param dorado Command name or executable path for Dorado.
#' @param conda_env Optional Conda environment name.
#' @param conda Conda command name or executable path.
#' @param dry_run Logical. Return the planned command without running it.
#' @param echo Logical. Print the planned command.
#' @param stderr Passed to [system2()] for Dorado's standard error stream.
#'
#' @return Invisibly returns command details and discovered FASTQ paths.
#' @export
demux_one_bam_fastq <- function(bam_file,
                                demux_dir,
                                threads = NULL,
                                dorado = "dorado",
                                conda_env = NULL,
                                conda = "conda",
                                dry_run = FALSE,
                                echo = TRUE,
                                stderr = "") {
  check_scalar_character(bam_file, "bam_file")
  check_scalar_character(demux_dir, "demux_dir")
  check_scalar_character(dorado, "dorado")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  threads <- validate_optional_positive_integer(threads, "threads")
  validate_system2_stderr(stderr)

  if (isTRUE(dry_run)) {
    bam_file <- normalizePath(bam_file, mustWork = FALSE)
    demux_dir <- normalizePath(demux_dir, mustWork = FALSE)
  } else {
    check_file_arg(bam_file, "bam_file")
    bam_file <- normalizePath(bam_file, mustWork = TRUE)
    if (dir.exists(demux_dir) && length(list.files(demux_dir, all.files = TRUE, no.. = TRUE))) {
      stop("`demux_dir` must be empty: ", demux_dir, call. = FALSE)
    }
    dir.create(demux_dir, recursive = TRUE, showWarnings = FALSE)
    demux_dir <- normalizePath(demux_dir, mustWork = TRUE)
  }

  args <- c(
    "demux", bam_file,
    "--output-dir", demux_dir,
    "--no-classify", "--emit-fastq"
  )
  if (!is.null(threads)) args <- c(args, "--threads", as.character(threads))

  call <- realtime_dorado_external_call(dorado, args, conda_env, conda)
  command_string <- paste(c(shQuote(call$command), shQuote(call$args)), collapse = " ")
  if (isTRUE(echo)) message(command_string)

  result <- list(
    status = NA_integer_,
    command = call$command,
    args = call$args,
    command_string = command_string,
    paths = list(bam_file = bam_file, demux_dir = demux_dir, fastq_files = character()),
    conda_env = conda_env
  )
  if (isTRUE(dry_run)) return(invisible(result))

  if (is.null(conda_env)) require_external_command(dorado) else require_external_command(conda)
  status <- suppressWarnings(system2(
    command = call$command,
    args = shQuote(call$args),
    stdout = "",
    stderr = stderr
  ))
  if (!identical(status, 0L)) {
    stop("dorado demux failed with exit status: ", status, call. = FALSE)
  }

  fastq_files <- realtime_fastq_files(demux_dir)
  if (!length(fastq_files)) {
    stop("dorado demux completed but produced no FASTQ files in: ", demux_dir,
         call. = FALSE)
  }

  result$status <- status
  result$paths$fastq_files <- fastq_files
  invisible(result)
}

#' Publish one POD5 batch of demultiplexed FASTQ files
#'
#' FASTQ chunks produced by Dorado are grouped by alias and concatenated into
#' one immutable file per POD5 and alias, named
#' `<pod5_stem>__<alias>.fastq[.gz]`.
#'
#' @param demux_dir Dorado demultiplexing output directory.
#' @param pod5_stem POD5 basename without the `.pod5` extension.
#' @param fastq_dir Root directory in which alias directories are created.
#' @param include_unclassified Include the `unclassified` group.
#' @param overwrite Replace an existing file for the same POD5 and alias.
#'
#' @return A data frame describing the published files.
#' @export
publish_demux_fastqs <- function(demux_dir,
                                 pod5_stem,
                                 fastq_dir,
                                 include_unclassified = TRUE,
                                 overwrite = FALSE) {
  check_dir_arg(demux_dir, "demux_dir")
  check_scalar_character(pod5_stem, "pod5_stem")
  check_scalar_character(fastq_dir, "fastq_dir")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(overwrite, "overwrite")
  if (!nzchar(pod5_stem) || grepl("[/\\\\]", pod5_stem)) {
    stop("`pod5_stem` must be a non-empty file stem.", call. = FALSE)
  }

  inputs <- realtime_fastq_files(demux_dir)
  if (!length(inputs)) {
    stop("No FASTQ files found under: ", demux_dir, call. = FALSE)
  }
  aliases <- basename(dirname(inputs))
  keep <- include_unclassified | aliases != "unclassified"
  inputs <- inputs[keep]
  aliases <- aliases[keep]
  if (!length(inputs)) {
    stop("No publishable FASTQ files found under: ", demux_dir, call. = FALSE)
  }

  dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
  fastq_dir <- normalizePath(fastq_dir, mustWork = TRUE)
  groups <- split(inputs, aliases)

  plan <- lapply(names(groups), function(alias) {
    files <- sort(groups[[alias]])
    compressed <- grepl("[.]gz$", files)
    if (length(unique(compressed)) != 1L) {
      stop("Mixed compressed and uncompressed FASTQ chunks for alias: ", alias,
           call. = FALSE)
    }
    alias_dir <- file.path(fastq_dir, alias)
    output <- file.path(
      alias_dir,
      paste0(pod5_stem, "__", alias, if (compressed[[1L]]) ".fastq.gz" else ".fastq")
    )
    list(alias = alias, inputs = files, alias_dir = alias_dir, output = output)
  })

  existing <- vapply(plan, function(x) file.exists(x$output), logical(1))
  if (any(existing) && !isTRUE(overwrite)) {
    stop("Published FASTQ already exists: ", plan[[which(existing)[[1L]]]]$output,
         call. = FALSE)
  }

  outputs <- character(length(plan))
  for (i in seq_along(plan)) {
    item <- plan[[i]]
    dir.create(item$alias_dir, recursive = TRUE, showWarnings = FALSE)
    realtime_concatenate_files(item$inputs, item$output, overwrite = overwrite)
    outputs[[i]] <- normalizePath(item$output, mustWork = TRUE)
  }

  data.frame(
    alias = names(groups),
    output_fastq = outputs,
    n_chunks = vapply(plan, function(x) length(x$inputs), integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Write a completion marker for one POD5 batch
#'
#' @param done_file Path to the `.done` marker.
#' @param pod5_file Source POD5 file.
#' @param bam_file Completed BAM file.
#' @param fastq_files Published FASTQ files.
#' @param overwrite Replace an existing marker.
#'
#' @return Invisibly returns the normalized marker path.
#' @export
write_done_marker <- function(done_file,
                              pod5_file,
                              bam_file,
                              fastq_files,
                              overwrite = FALSE) {
  check_scalar_character(done_file, "done_file")
  check_file_arg(pod5_file, "pod5_file")
  check_file_arg(bam_file, "bam_file")
  check_logical_scalar(overwrite, "overwrite")
  if (!is.character(fastq_files) || !length(fastq_files) || anyNA(fastq_files)) {
    stop("`fastq_files` must contain at least one path.", call. = FALSE)
  }
  missing <- fastq_files[!file.exists(fastq_files)]
  if (length(missing)) stop("FASTQ file does not exist: ", missing[[1L]], call. = FALSE)

  dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
  done_file <- normalizePath(done_file, mustWork = FALSE)
  if (file.exists(done_file) && !isTRUE(overwrite)) {
    stop("Done marker already exists: ", done_file, call. = FALSE)
  }

  lines <- c(
    paste("completed_at", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
    paste("pod5", normalizePath(pod5_file, mustWork = TRUE), sep = "\t"),
    paste("bam", normalizePath(bam_file, mustWork = TRUE), sep = "\t"),
    paste("fastq", normalizePath(fastq_files, mustWork = TRUE), sep = "\t")
  )
  temporary <- tempfile(pattern = paste0(basename(done_file), "."), tmpdir = dirname(done_file))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (file.exists(done_file)) unlink(done_file)
  if (!file.rename(temporary, done_file)) {
    stop("Unable to publish done marker: ", done_file, call. = FALSE)
  }
  invisible(normalizePath(done_file, mustWork = TRUE))
}

#' Process one completed POD5 file through FASTQ publication
#'
#' @param pod5_file Input POD5 file.
#' @param bam_dir Directory for the corresponding BAM file.
#' @param fastq_dir Root directory for per-alias FASTQ fragments.
#' @param state_dir Directory for `.done` markers and temporary demux output.
#' @param model Dorado basecalling model or alias.
#' @param kit_name Dorado barcode kit name.
#' @param sample_sheet MinKNOW-compatible sample sheet containing aliases.
#' @param threads Optional Dorado demux thread count.
#' @param include_unclassified Publish unclassified reads.
#' @param dorado,conda_env,conda,dry_run,echo,stderr Passed to the component
#'   Dorado functions.
#'
#' @return Invisibly returns paths and component results.
#' @export
process_one_pod5 <- function(pod5_file,
                             bam_dir,
                             fastq_dir,
                             state_dir,
                             model = "hac",
                             kit_name,
                             sample_sheet,
                             threads = NULL,
                             include_unclassified = TRUE,
                             dorado = "dorado",
                             conda_env = NULL,
                             conda = "conda",
                             dry_run = FALSE,
                             echo = TRUE,
                             stderr = "") {
  check_file_arg(pod5_file, "pod5_file")
  check_scalar_character(bam_dir, "bam_dir")
  check_scalar_character(fastq_dir, "fastq_dir")
  check_scalar_character(state_dir, "state_dir")
  check_scalar_character(model, "model")
  check_scalar_character(kit_name, "kit_name")
  check_file_arg(sample_sheet, "sample_sheet")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(dry_run, "dry_run")

  pod5_file <- normalizePath(pod5_file, mustWork = TRUE)
  sample_sheet <- normalizePath(sample_sheet, mustWork = TRUE)
  stem <- sub("[.]pod5$", "", basename(pod5_file), ignore.case = TRUE)
  if (identical(stem, basename(pod5_file))) {
    stop("`pod5_file` must end in .pod5.", call. = FALSE)
  }

  dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  bam_dir <- normalizePath(bam_dir, mustWork = TRUE)
  fastq_dir <- normalizePath(fastq_dir, mustWork = TRUE)
  state_dir <- normalizePath(state_dir, mustWork = TRUE)
  bam_file <- file.path(bam_dir, paste0(stem, ".bam"))
  done_file <- file.path(state_dir, paste0(stem, ".done"))

  if (file.exists(done_file) && !isTRUE(dry_run)) {
    return(invisible(list(
      status = "skipped",
      paths = list(pod5_file = pod5_file, bam_file = bam_file,
                   fastq_files = character(), done_file = done_file)
    )))
  }

  basecall_result <- NULL
  if (!file.exists(bam_file) || isTRUE(dry_run)) {
    basecall_result <- basecall_one(
      pod5_file = pod5_file,
      bam_file = bam_file,
      model = model,
      kit_name = kit_name,
      sample_sheet = sample_sheet,
      dorado = dorado,
      conda_env = conda_env,
      conda = conda,
      dry_run = dry_run,
      echo = echo,
      stderr = stderr
    )
  } else if (is.na(file.info(bam_file)$size) || file.info(bam_file)$size == 0) {
    stop("Existing BAM is empty: ", bam_file, call. = FALSE)
  }

  demux_work_root <- file.path(state_dir, ".demux_work")
  dir.create(demux_work_root, recursive = TRUE, showWarnings = FALSE)
  demux_dir <- tempfile(pattern = paste0(stem, "."), tmpdir = demux_work_root)
  on.exit(unlink(demux_dir, recursive = TRUE, force = TRUE), add = TRUE)

  demux_result <- demux_one_bam_fastq(
    bam_file = bam_file,
    demux_dir = demux_dir,
    threads = threads,
    dorado = dorado,
    conda_env = conda_env,
    conda = conda,
    dry_run = dry_run,
    echo = echo,
    stderr = stderr
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = "planned",
      basecall = basecall_result,
      demux = demux_result,
      paths = list(pod5_file = pod5_file, bam_file = bam_file,
                   fastq_files = character(), done_file = done_file)
    )))
  }

  publication <- publish_demux_fastqs(
    demux_dir = demux_dir,
    pod5_stem = stem,
    fastq_dir = fastq_dir,
    include_unclassified = include_unclassified,
    overwrite = TRUE
  )
  done_file <- write_done_marker(
    done_file = done_file,
    pod5_file = pod5_file,
    bam_file = bam_file,
    fastq_files = publication$output_fastq,
    overwrite = TRUE
  )

  invisible(list(
    status = "completed",
    basecall = basecall_result,
    demux = demux_result,
    publication = publication,
    paths = list(
      pod5_file = pod5_file,
      bam_file = normalizePath(bam_file, mustWork = TRUE),
      fastq_files = publication$output_fastq,
      done_file = done_file
    )
  ))
}

#' Watch a POD5 directory and publish per-alias FASTQ fragments
#'
#' The watcher preserves the existing ONTools stability rule: a POD5 is
#' processed after its size is unchanged across one checking interval.
#'
#' @param pod5_path Directory containing POD5 files.
#' @param bam_path Directory for one BAM per POD5.
#' @param fastq_path Root directory for per-alias FASTQ fragments.
#' @param state_path Directory for completion markers.
#' @param max_run_time Maximum watcher runtime in seconds.
#' @param check_interval Seconds between POD5 size checks.
#' @param idle_interval Maximum seconds to sleep when no work is available.
#' @inheritParams process_one_pod5
#'
#' @return Invisibly returns elapsed time and per-POD5 results.
#' @export
watch_dorado_fastq <- function(pod5_path,
                               bam_path,
                               fastq_path,
                               state_path,
                               kit_name,
                               sample_sheet,
                               max_run_time = 10,
                               check_interval = 3,
                               idle_interval = 30,
                               model = "hac",
                               threads = NULL,
                               include_unclassified = TRUE,
                               dorado = "dorado",
                               conda_env = NULL,
                               conda = "conda",
                               dry_run = FALSE,
                               echo = TRUE,
                               stderr = "") {
  check_dir_arg(pod5_path, "pod5_path")
  check_scalar_character(bam_path, "bam_path")
  check_scalar_character(fastq_path, "fastq_path")
  check_scalar_character(state_path, "state_path")
  check_scalar_character(kit_name, "kit_name")
  check_file_arg(sample_sheet, "sample_sheet")
  max_run_time <- validate_watch_duration(max_run_time, "max_run_time")
  check_interval <- validate_watch_duration(check_interval, "check_interval")
  idle_interval <- validate_watch_duration(idle_interval, "idle_interval")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  pod5_path <- normalizePath(pod5_path, mustWork = TRUE)
  dir.create(bam_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(fastq_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(state_path, recursive = TRUE, showWarnings = FALSE)
  bam_path <- normalizePath(bam_path, mustWork = TRUE)
  fastq_path <- normalizePath(fastq_path, mustWork = TRUE)
  state_path <- normalizePath(state_path, mustWork = TRUE)

  start_time <- Sys.time()
  elapsed_time <- function() as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  results <- list()

  repeat {
    elapsed <- elapsed_time()
    if (elapsed >= max_run_time && !isTRUE(dry_run)) break
    if (isTRUE(echo)) message("Watcher elapsed time: ", round(elapsed, 1), " seconds")

    pod5_files <- watch_file(pod5_path, ".pod5")
    done_files <- watch_file(state_path, ".done")
    new_names <- setdiff(pod5_files$name, done_files$name)

    if (isTRUE(dry_run)) {
      for (file_name in new_names) {
        input <- pod5_files$path[pod5_files$name == file_name][[1L]]
        results[[file_name]] <- process_one_pod5(
          input, bam_path, fastq_path, state_path,
          model = model, kit_name = kit_name, sample_sheet = sample_sheet,
          threads = threads, include_unclassified = include_unclassified,
          dorado = dorado, conda_env = conda_env, conda = conda,
          dry_run = TRUE, echo = echo, stderr = stderr
        )
      }
      break
    }

    if (!length(new_names)) {
      remaining <- max_run_time - elapsed_time()
      wait_time <- min(idle_interval, max(0, remaining))
      if (wait_time > 0) Sys.sleep(wait_time)
      next
    }

    for (file_name in new_names) {
      if (elapsed_time() >= max_run_time) break
      input <- pod5_files$path[pod5_files$name == file_name][[1L]]
      last_size <- file.info(input)$size
      stable <- FALSE

      while (elapsed_time() < max_run_time) {
        remaining <- max_run_time - elapsed_time()
        Sys.sleep(min(check_interval, max(0, remaining)))
        current_size <- file.info(input)$size
        if (!length(current_size) || is.na(current_size)) break
        if (isTRUE(current_size == last_size)) {
          stable <- TRUE
          break
        }
        if (isTRUE(echo)) message("POD5 file is still growing: ", input)
        last_size <- current_size
      }

      if (!stable || elapsed_time() >= max_run_time) next
      results[[file_name]] <- process_one_pod5(
        input, bam_path, fastq_path, state_path,
        model = model, kit_name = kit_name, sample_sheet = sample_sheet,
        threads = threads, include_unclassified = include_unclassified,
        dorado = dorado, conda_env = conda_env, conda = conda,
        dry_run = FALSE, echo = echo, stderr = stderr
      )
    }
  }

  invisible(list(
    elapsed_time = elapsed_time(),
    batches = results,
    paths = list(pod5_path = pod5_path, bam_path = bam_path,
                 fastq_path = fastq_path, state_path = state_path)
  ))
}

#' Finalize a run by merging per-POD5 FASTQ fragments by alias
#'
#' This function is intentionally separate from [watch_dorado_fastq()] and
#' should only be called by the user after sequencing and real-time processing
#' are complete.
#'
#' @param fastq_path Root directory created by `publish_demux_fastqs()`.
#' @param output_path Directory for one final FASTQ per alias.
#' @param include_unclassified Include the `unclassified` alias.
#' @param overwrite Replace existing final FASTQ files.
#'
#' @return A data frame describing the final FASTQ files.
#' @export
finalize_run_fastqs <- function(fastq_path,
                                output_path,
                                include_unclassified = TRUE,
                                overwrite = FALSE) {
  check_dir_arg(fastq_path, "fastq_path")
  check_scalar_character(output_path, "output_path")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(overwrite, "overwrite")

  fastq_path <- normalizePath(fastq_path, mustWork = TRUE)
  alias_dirs <- list.dirs(fastq_path, recursive = FALSE, full.names = TRUE)
  if (!isTRUE(include_unclassified)) {
    alias_dirs <- alias_dirs[basename(alias_dirs) != "unclassified"]
  }
  alias_dirs <- alias_dirs[vapply(alias_dirs, function(x) length(realtime_fastq_files(x)) > 0L,
                                  logical(1))]
  if (!length(alias_dirs)) stop("No per-alias FASTQ fragments found in: ", fastq_path,
                                call. = FALSE)

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path, mustWork = TRUE)
  plan <- lapply(alias_dirs, function(alias_dir) {
    alias <- basename(alias_dir)
    inputs <- realtime_fastq_files(alias_dir)
    compressed <- grepl("[.]gz$", inputs)
    if (length(unique(compressed)) != 1L) {
      stop("Mixed compressed and uncompressed FASTQ fragments for alias: ", alias,
           call. = FALSE)
    }
    output <- file.path(output_path, paste0(alias, if (compressed[[1L]]) ".fastq.gz" else ".fastq"))
    list(alias = alias, inputs = sort(inputs), output = output)
  })

  existing <- vapply(plan, function(x) file.exists(x$output), logical(1))
  if (any(existing) && !isTRUE(overwrite)) {
    stop("Final FASTQ already exists: ", plan[[which(existing)[[1L]]]]$output,
         call. = FALSE)
  }

  for (item in plan) {
    realtime_concatenate_files(item$inputs, item$output, overwrite = overwrite)
  }
  data.frame(
    alias = vapply(plan, `[[`, character(1), "alias"),
    output_fastq = vapply(plan, function(x) normalizePath(x$output, mustWork = TRUE), character(1)),
    n_fragments = vapply(plan, function(x) length(x$inputs), integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

realtime_dorado_external_call <- function(dorado, args, conda_env, conda) {
  if (is.null(conda_env)) {
    list(command = dorado, args = args)
  } else {
    list(
      command = conda,
      args = c("run", "--no-capture-output", "-n", conda_env, "--", dorado, args)
    )
  }
}

realtime_fastq_files <- function(path) {
  files <- list.files(
    path,
    pattern = "[.]fastq([.]gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  files <- files[file.exists(files) & !dir.exists(files)]
  sort(normalizePath(files, mustWork = TRUE))
}

realtime_concatenate_files <- function(inputs, output, overwrite = FALSE) {
  if (file.exists(output) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", output, call. = FALSE)
  }
  temporary <- tempfile(pattern = paste0(basename(output), "."), tmpdir = dirname(output))
  on.exit(unlink(temporary), add = TRUE)
  destination <- file(temporary, open = "wb")
  on.exit(close(destination), add = TRUE)

  for (input in inputs) {
    source <- file(input, open = "rb")
    repeat {
      chunk <- readBin(source, what = "raw", n = 1024L * 1024L)
      if (!length(chunk)) break
      writeBin(chunk, destination)
    }
    close(source)
  }
  close(destination)
  on.exit(NULL, add = FALSE)

  if (file.exists(output)) unlink(output)
  if (!file.rename(temporary, output)) {
    stop("Unable to publish output file: ", output, call. = FALSE)
  }
  invisible(output)
}
