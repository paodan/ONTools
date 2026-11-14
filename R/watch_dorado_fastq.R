#' Create a standardized Dorado run directory structure
#'
#' `create_dorado_run_structure()` creates all directories used by the
#' real-time POD5-to-FASTQ workflow and returns their resolved paths.
#'
#' @param run_path Root directory for one sequencing run.
#' @param create Logical. Create missing directories when `TRUE`; when `FALSE`,
#'   only return the paths that would be used.
#'
#' @return A named list containing the run root and standardized output
#'   directories.
#'
#' @details
#' The standardized layout is:
#'
#' ```
#' <run_path>/
#' |-- bam/
#' |-- fastq_live/
#' |-- fastq_final/
#' |-- qc_live/
#' |-- qc_cache/
#' |-- state/
#' `-- logs/
#' ```
#'
#' POD5 files and the sample sheet are inputs and are deliberately not copied,
#' moved, or created under `run_path`. Pass their existing paths directly to
#' [watch_dorado_fastq()] or [process_one_pod5()].
#'
#' @examples
#' layout <- create_dorado_run_structure(tempfile("dorado-run-"))
#' names(layout)
#' dir.exists(layout$fastq_path)
#' @export
create_dorado_run_structure <- function(run_path, create = TRUE) {
  check_scalar_character(run_path, "run_path")
  check_logical_scalar(create, "create")
  if (!nzchar(run_path)) stop("`run_path` must not be empty.", call. = FALSE)

  if (isTRUE(create)) {
    dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(run_path)) {
      stop("Unable to create `run_path`: ", run_path, call. = FALSE)
    }
  }
  root <- normalizePath(run_path, mustWork = create)
  layout <- list(
    run_path = root,
    bam_path = file.path(root, "bam"),
    fastq_path = file.path(root, "fastq_live"),
    final_fastq_path = file.path(root, "fastq_final"),
    qc_path = file.path(root, "qc_live"),
    qc_cache_path = file.path(root, "qc_cache"),
    state_path = file.path(root, "state"),
    log_path = file.path(root, "logs")
  )
  if (isTRUE(create)) {
    directories <- unlist(layout[c(
      "bam_path", "fastq_path", "final_fastq_path",
      "qc_path", "qc_cache_path", "state_path", "log_path"
    )], use.names = FALSE)
    for (directory in directories) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(directory)) {
        stop("Unable to create run directory: ", directory, call. = FALSE)
      }
    }
    layout <- lapply(layout, normalizePath, mustWork = TRUE)
  }
  class(layout) <- c("dorado_run_structure", "list")
  layout
}

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
#'
#' @details
#' This function must receive a BAM that was already classified during
#' basecalling. It deliberately uses `--no-classify`, so `kit_name` and
#' `sample_sheet` are not passed again. Use a unique `demux_dir` for every BAM
#' because Dorado may reuse batch file names such as `_0.fastq`.
#'
#' @examples
#' planned <- demux_one_bam_fastq(
#'   bam_file = "/data/project/bam/batch001.bam",
#'   demux_dir = "/data/project/state/demux/batch001",
#'   threads = 8,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' planned$command_string
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
#' `<pod5_stem>__<alias>.fastq.gz` by default.
#'
#' @param demux_dir Dorado demultiplexing output directory.
#' @param pod5_stem POD5 basename without the `.pod5` extension.
#' @param fastq_dir Root directory in which alias directories are created.
#' @param include_unclassified Include the `unclassified` group.
#' @param compress Logical. If `TRUE` (the default), gzip-compress each
#'   published per-POD5 FASTQ. If `FALSE`, write uncompressed `.fastq` files.
#' @param overwrite Replace an existing file for the same POD5 and alias.
#'
#' @return A data frame describing the published files.
#'
#' @details
#' Dorado may create several FASTQ chunks for one alias. They are concatenated
#' without parsing the records. With `compress = TRUE`, uncompressed Dorado
#' output is streamed directly through gzip into the published file. The
#' published files are organized as
#' `fastq_dir/<alias>/<pod5_stem>__<alias>.fastq.gz`.
#'
#' @examples
#' demux <- tempfile("demux-")
#' dir.create(file.path(demux, "run", "fastq_pass", "sampleA"), recursive = TRUE)
#' writeLines(
#'   c("@read1", "ACGT", "+", "IIII"),
#'   file.path(demux, "run", "fastq_pass", "sampleA", "chunk1.fastq")
#' )
#' writeLines(
#'   c("@read2", "TGCA", "+", "IIII"),
#'   file.path(demux, "run", "fastq_pass", "sampleA", "chunk2.fastq")
#' )
#' published <- publish_demux_fastqs(
#'   demux_dir = demux,
#'   pod5_stem = "batch001",
#'   fastq_dir = tempfile("fastq-live-"),
#'   compress = TRUE
#' )
#' published
#' con <- gzfile(published$output_fastq, open = "rt")
#' readLines(con)
#' close(con)
#' @export
publish_demux_fastqs <- function(demux_dir,
                                 pod5_stem,
                                 fastq_dir,
                                 include_unclassified = TRUE,
                                 compress = TRUE,
                                 overwrite = FALSE) {
  check_dir_arg(demux_dir, "demux_dir")
  check_scalar_character(pod5_stem, "pod5_stem")
  check_scalar_character(fastq_dir, "fastq_dir")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(compress, "compress")
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
    output <- file.path(alias_dir, paste0(
      pod5_stem, "__", alias,
      if (isTRUE(compress)) ".fastq.gz" else ".fastq"
    ))
    list(
      alias = alias,
      inputs = files,
      input_compressed = compressed[[1L]],
      alias_dir = alias_dir,
      output = output
    )
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
    realtime_concatenate_fastq_files(
      inputs = item$inputs,
      output = item$output,
      input_compressed = item$input_compressed,
      output_compressed = compress,
      overwrite = overwrite
    )
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

#' Update one cumulative QC plot for an entire FASTQ directory
#'
#' `update_fastq_directory_qc()` treats every immediate alias directory as one
#' sample. It temporarily merges all per-POD5 FASTQ fragments belonging to an
#' alias, then passes the complete set of one-file-per-alias inputs to
#' [plot_fastq_statistics()] in a single call.
#'
#' @param fastq_path Root directory containing one subdirectory per alias.
#'   May be `NULL` when `run_path` is supplied.
#' @param qc_path Directory in which the run-level QC plot is written. May be
#'   `NULL` when `run_path` is supplied.
#' @param run_path Optional standardized run root created by
#'   [create_dorado_run_structure()]. It supplies missing `fastq_path` and
#'   `qc_path` values.
#' @param include_unclassified Include the `unclassified` alias.
#' @param output_name Base name for the QC plot, without an extension.
#' @param format Plot format, either `"png"` or `"pdf"`.
#' @param width,height Positive plot dimensions in inches.
#' @param dpi Positive PNG resolution in dots per inch.
#' @param plot_args Named list of additional arguments passed to
#'   [plot_fastq_statistics()]. Do not include `fastqs` or `path`.
#' @param error_action Action when plotting fails: `"warn"` records an
#'   error log and emits a warning, `"silent"` only records the log, and
#'   `"stop"` records the log and stops.
#'
#' @return A one-row data frame describing the status, alias and fragment
#'   counts, plot path, error-log path, and any error message.
#'
#' @details
#' Merged alias FASTQs exist only in a temporary working directory and are
#' deleted when the function returns. `plot_fastq_statistics()` is called once
#' with these merged files, `barcode_column = "source_file"`, and explicit
#' sample labels equal to the alias directory names. Consequently, each alias
#' contributes one sample to both the read-depth and peak-statistics panels,
#' regardless of barcode tags stored in the FASTQ headers.
#'
#' The final plot is published atomically. If plotting fails, the previous
#' successful plot is preserved and `<qc_path>/<output_name>.error.log` is
#' updated. A later call automatically retries with all fragments then present.
#'
#' @examples
#' root <- tempfile("directory-qc-")
#' layout <- create_dorado_run_structure(root)
#' dir.create(file.path(layout$fastq_path, "sampleA"))
#' dir.create(file.path(layout$fastq_path, "sampleB"))
#' write_example <- function(path, prefix, lengths) {
#'   records <- unlist(lapply(seq_along(lengths), function(i) {
#'     sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'     quality <- paste(rep("I", lengths[[i]]), collapse = "")
#'     c(
#'       paste0("@", prefix, i, " qs:f:15.5 ch:i:1 SM:Z:barcode01"),
#'       sequence, "+", quality
#'     )
#'   }), use.names = FALSE)
#'   con <- gzfile(path, "wt")
#'   writeLines(records, con)
#'   close(con)
#' }
#' write_example(
#'   file.path(layout$fastq_path, "sampleA", "batch001__sampleA.fastq.gz"),
#'   "a", c(95:105, 295:305)
#' )
#' write_example(
#'   file.path(layout$fastq_path, "sampleB", "batch001__sampleB.fastq.gz"),
#'   "b", c(145:155, 395:405)
#' )
#'
#' qc <- update_fastq_directory_qc(
#'   run_path = root,
#'   include_unclassified = FALSE,
#'   plot_args = list(min_distance = 100, density_n = 1024)
#' )
#' qc[, c("status", "n_aliases", "n_fragments")]
#' file.exists(qc$plot_file)
#' @export
update_fastq_directory_qc <- function(fastq_path = NULL,
                                      qc_path = NULL,
                                      include_unclassified = TRUE,
                                      output_name = "fastq_statistics",
                                      format = c("png", "pdf"),
                                      width = 12,
                                      height = 10,
                                      dpi = 150,
                                      plot_args = list(),
                                      error_action = c("warn", "silent", "stop"),
                                      run_path = NULL) {
  if (!is.null(run_path)) {
    layout <- create_dorado_run_structure(run_path)
    if (is.null(fastq_path)) fastq_path <- layout$fastq_path
    if (is.null(qc_path)) qc_path <- layout$qc_path
  }
  if (is.null(fastq_path)) {
    stop("Supply `fastq_path` or `run_path`.", call. = FALSE)
  }
  if (is.null(qc_path)) {
    stop("Supply `qc_path` or `run_path`.", call. = FALSE)
  }
  check_dir_arg(fastq_path, "fastq_path")
  check_scalar_character(qc_path, "qc_path")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_scalar_character(output_name, "output_name")
  if (!nzchar(output_name) || grepl("[/\\\\]", output_name)) {
    stop("`output_name` must be a non-empty file name without a path.", call. = FALSE)
  }
  format <- match.arg(format)
  error_action <- match.arg(error_action)
  width <- realtime_positive_number(width, "width")
  height <- realtime_positive_number(height, "height")
  dpi <- realtime_positive_number(dpi, "dpi")
  if (!is.list(plot_args)) {
    stop("`plot_args` must be a list.", call. = FALSE)
  }
  if (length(plot_args)) {
    if (is.null(names(plot_args)) || any(!nzchar(names(plot_args))) || anyDuplicated(names(plot_args))) {
      stop("`plot_args` must be a uniquely named list.", call. = FALSE)
    }
    reserved <- intersect(names(plot_args), c(
      "fastqs", "path", "barcode_column", "sample_names"
    ))
    if (length(reserved)) {
      stop("`plot_args` must not contain: ", paste(reserved, collapse = ", "),
           call. = FALSE)
    }
  }

  fastq_path <- normalizePath(fastq_path, mustWork = TRUE)
  dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)
  qc_path <- normalizePath(qc_path, mustWork = TRUE)
  alias_dirs <- list.dirs(fastq_path, recursive = FALSE, full.names = TRUE)
  if (!isTRUE(include_unclassified)) {
    alias_dirs <- alias_dirs[basename(alias_dirs) != "unclassified"]
  }
  fragment_lists <- lapply(alias_dirs, realtime_fastq_files)
  keep <- lengths(fragment_lists) > 0L
  alias_dirs <- alias_dirs[keep]
  fragment_lists <- fragment_lists[keep]
  aliases <- basename(alias_dirs)
  if (!length(aliases)) {
    stop("No per-alias FASTQ fragments found under: ", fastq_path, call. = FALSE)
  }

  plot_file <- file.path(qc_path, paste0(output_name, ".", format))
  error_file <- file.path(qc_path, paste0(output_name, ".error.log"))
  work_dir <- tempfile("merged-alias-fastq-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  error_message <- tryCatch({
    merged_fastqs <- vapply(seq_along(aliases), function(i) {
      fragments <- fragment_lists[[i]]
      compressed <- grepl("[.]gz$", fragments, ignore.case = TRUE)
      if (length(unique(compressed)) != 1L) {
        stop("Mixed compressed and uncompressed FASTQ fragments for alias: ",
             aliases[[i]], call. = FALSE)
      }
      output <- file.path(work_dir, sprintf("sample_%06d.fastq.gz", i))
      realtime_concatenate_fastq_files(
        inputs = fragments,
        output = output,
        input_compressed = compressed[[1L]],
        output_compressed = TRUE,
        overwrite = FALSE
      )
      output
    }, character(1))

    args <- plot_args
    args$fastqs <- merged_fastqs
    args$barcode_column <- "source_file"
    args$sample_names <- aliases
    if (!"title" %in% names(args)) args$title <- "Cumulative FASTQ statistics"
    plot <- do.call(plot_fastq_statistics, args)
    temporary <- tempfile(
      pattern = paste0(output_name, "."),
      tmpdir = qc_path,
      fileext = paste0(".", format)
    )
    on.exit(unlink(temporary), add = TRUE)
    ggplot2::ggsave(
      filename = temporary,
      plot = plot,
      device = format,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE
    )
    if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
        file.info(temporary)$size == 0) {
      stop("QC renderer produced an empty plot file.", call. = FALSE)
    }
    realtime_publish_file(temporary, plot_file)
    if (file.exists(error_file)) unlink(error_file)
    NULL
  }, error = function(e) conditionMessage(e))

  status <- if (is.null(error_message)) "updated" else "error"
  if (!is.null(error_message)) {
    realtime_write_lines_atomic(
      c(
        paste("updated_at", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
        paste("n_aliases", length(aliases), sep = "\t"),
        paste("n_fragments", sum(lengths(fragment_lists)), sep = "\t"),
        paste("error", error_message, sep = "\t")
      ),
      error_file
    )
    message_text <- paste0("FASTQ directory QC failed: ", error_message)
    if (identical(error_action, "warn")) warning(message_text, call. = FALSE)
    if (identical(error_action, "stop")) stop(message_text, call. = FALSE)
  }

  data.frame(
    status = status,
    n_aliases = length(aliases),
    n_fragments = sum(lengths(fragment_lists)),
    plot_file = if (file.exists(plot_file)) normalizePath(plot_file) else NA_character_,
    error_file = if (file.exists(error_file)) normalizePath(error_file) else NA_character_,
    error_message = if (is.null(error_message)) NA_character_ else error_message,
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
#'
#' @details
#' The marker is written through a temporary file and renamed only after all
#' referenced files exist. A watcher therefore treats a POD5 as complete only
#' after BAM generation, demultiplexing, and FASTQ publication all succeed.
#'
#' @examples
#' root <- tempfile("done-example-")
#' dir.create(root)
#' pod5 <- file.path(root, "batch001.pod5")
#' bam <- file.path(root, "batch001.bam")
#' fastq <- file.path(root, "batch001__sampleA.fastq")
#' writeLines("pod5", pod5)
#' writeLines("bam", bam)
#' writeLines(c("@read1", "AC", "+", "II"), fastq)
#' marker <- write_done_marker(
#'   done_file = file.path(root, "batch001.done"),
#'   pod5_file = pod5,
#'   bam_file = bam,
#'   fastq_files = fastq
#' )
#' readLines(marker)
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
#' @param bam_dir Directory for the corresponding BAM file. May be `NULL` when
#'   `run_path` is supplied.
#' @param fastq_dir Root directory for per-alias FASTQ fragments. May be `NULL`
#'   when `run_path` is supplied.
#' @param state_dir Directory for `.done` markers and temporary demux output.
#'   May be `NULL` when `run_path` is supplied.
#' @param run_path Optional standardized run root created by
#'   [create_dorado_run_structure()]. It supplies missing output directories;
#'   it does not alter `pod5_file` or `sample_sheet`.
#' @param model Dorado basecalling model or alias.
#' @param kit_name Dorado barcode kit name.
#' @param sample_sheet Existing MinKNOW-compatible sample sheet containing
#'   aliases. The file is read in place and is not copied into `run_path`.
#' @param barcode_arrangement Optional custom Dorado barcode arrangement TOML
#'   file passed to [basecall_one()]. Requires `kit_name`.
#' @param barcode_sequences Optional custom barcode sequence FASTA file passed
#'   to [basecall_one()]. Requires `barcode_arrangement`.
#' @param extra_basecaller_params Optional character vector of additional
#'   `dorado basecaller` arguments passed to [basecall_one()]. Dedicated
#'   barcode and sample-sheet options must use their named arguments instead.
#' @param threads Optional Dorado demux thread count.
#' @param include_unclassified Publish unclassified reads.
#' @param compress_fastq Gzip-compress each published per-POD5 FASTQ fragment.
#' @param update_qc Logical. After the `.done` marker is written, update one
#'   cumulative FASTQ statistics plot for the entire `fastq_dir`.
#' @param qc_include_unclassified Include `unclassified` in the live QC plot.
#'   This is independent of `include_unclassified`, which controls whether
#'   unclassified FASTQ fragments are published.
#' @param qc_path Root directory for live QC plots. When `NULL`, use
#'   `qc_live` beside `fastq_dir`.
#' @param qc_every_n_pod5 Positive integer controlling run-level QC frequency.
#'   For example, `1` updates after every completed POD5 and `5` after every
#'   fifth completed POD5. A failed attempt is retried by a later eligible batch.
#' @param qc_format,qc_width,qc_height,qc_dpi Plot settings passed to
#'   [update_fastq_directory_qc()].
#' @param qc_plot_args Named list of additional arguments passed from
#'   [update_fastq_directory_qc()] to [plot_fastq_statistics()].
#' @param dorado,conda_env,conda,dry_run,echo,stderr Passed to the component
#'   Dorado functions.
#'
#' @return Invisibly returns paths and component results.
#'
#' @details
#' The function performs one transactional batch:
#'
#' 1. Create or reuse `<bam_dir>/<pod5_stem>.bam`.
#' 2. Run Dorado demultiplexing in a private temporary directory.
#' 3. Publish one FASTQ fragment per alias under `fastq_dir`.
#' 4. Write `<state_dir>/<pod5_stem>.done`.
#' 5. Optionally update one cumulative run-level QC plot, treating every alias
#'    directory as one merged sample.
#'
#' If the BAM exists but the marker does not, basecalling is not repeated; the
#' function retries demultiplexing and publication. If the marker exists, the
#' batch is returned with status `"skipped"`. The completion marker is written
#' before QC, so a plotting failure never causes basecalling or demultiplexing
#' to be repeated. QC errors are recorded under `qc_path` and retried later.
#'
#' The FASTQ fragments created here are per-POD5 files. They are not the final
#' run-level FASTQ files; call [finalize_run_fastqs()] only after sequencing is
#' complete.
#'
#' @examples
#' input_dir <- tempfile("one-pod5-input-")
#' dir.create(input_dir)
#' run_dir <- tempfile("one-pod5-output-")
#' pod5 <- file.path(input_dir, "batch001.pod5")
#' sheet <- file.path(input_dir, "sample_sheet.csv")
#' writeLines("pod5", pod5)
#' writeLines(
#'   c(
#'     "experiment_id,flow_cell_id,kit,barcode,alias",
#'     "run1,FLOWCELL01,SQK-NBD114-96,barcode01,sampleA"
#'   ),
#'   sheet
#' )
#'
#' # dry_run validates paths and returns both planned Dorado commands without
#' # requiring a Dorado installation or a real POD5 file.
#' planned <- process_one_pod5(
#'   pod5_file = pod5,
#'   sample_sheet = sheet,
#'   run_path = run_dir,
#'   kit_name = "SQK-NBD114-96",
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' planned$status
#' planned$basecall$command_string
#' planned$demux$command_string
#'
#' \dontrun{
#' # With Dorado installed and a real POD5, omit dry_run.
#' completed <- process_one_pod5(
#'   pod5_file = "/data/minknow/run01/pod5/batch001.pod5",
#'   sample_sheet = "/data/config/sample_sheet.csv",
#'   run_path = "/data/project/run01",
#'   model = "hac",
#'   kit_name = "SQK-NBD114-96",
#'   threads = 8,
#'   update_qc = TRUE,
#'   qc_every_n_pod5 = 1,
#'   qc_plot_args = list(min_distance = 100, density_n = 1024)
#' )
#'
#' # Custom barcodes: `kit_name` must match the name in the arrangement TOML.
#' process_one_pod5(
#'   pod5_file = "/data/minknow/run01/pod5/batch001.pod5",
#'   sample_sheet = "/data/config/sample_sheet.csv",
#'   run_path = "/data/project/run01",
#'   kit_name = "MY_CUSTOM_KIT",
#'   barcode_arrangement = "/data/config/custom_barcodes.toml",
#'   barcode_sequences = "/data/config/custom_barcodes.fasta",
#'   extra_basecaller_params = c("--device", "cuda:all")
#' )
#' completed$paths
#' }
#' @export
process_one_pod5 <- function(pod5_file,
                             bam_dir = NULL,
                             fastq_dir = NULL,
                             state_dir = NULL,
                             model = "hac",
                             kit_name,
                             sample_sheet,
                             threads = NULL,
                             include_unclassified = TRUE,
                             compress_fastq = TRUE,
                             update_qc = TRUE,
                             qc_include_unclassified = include_unclassified,
                             qc_path = NULL,
                             qc_every_n_pod5 = 1L,
                             qc_format = "png",
                             qc_width = 12,
                             qc_height = 10,
                             qc_dpi = 150,
                             qc_plot_args = list(),
                             dorado = "dorado",
                             conda_env = NULL,
                             conda = "conda",
                             dry_run = FALSE,
                             echo = TRUE,
                             stderr = "",
                             run_path = NULL,
                             barcode_arrangement = NULL,
                             barcode_sequences = NULL,
                             extra_basecaller_params = NULL) {
  layout <- NULL
  if (!is.null(run_path)) {
    layout <- create_dorado_run_structure(run_path)
    if (is.null(bam_dir)) bam_dir <- layout$bam_path
    if (is.null(fastq_dir)) fastq_dir <- layout$fastq_path
    if (is.null(state_dir)) state_dir <- layout$state_path
    if (is.null(qc_path)) qc_path <- layout$qc_path
  }
  if (is.null(bam_dir) || is.null(fastq_dir) || is.null(state_dir)) {
    stop("Supply `bam_dir`, `fastq_dir`, and `state_dir`, or supply `run_path`.",
         call. = FALSE)
  }
  if (is.null(sample_sheet)) {
    stop("Supply the existing `sample_sheet` path.", call. = FALSE)
  }
  check_file_arg(pod5_file, "pod5_file")
  check_scalar_character(bam_dir, "bam_dir")
  check_scalar_character(fastq_dir, "fastq_dir")
  check_scalar_character(state_dir, "state_dir")
  check_scalar_character(model, "model")
  check_scalar_character(kit_name, "kit_name")
  check_file_arg(sample_sheet, "sample_sheet")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(compress_fastq, "compress_fastq")
  check_logical_scalar(update_qc, "update_qc")
  check_logical_scalar(qc_include_unclassified, "qc_include_unclassified")
  check_logical_scalar(dry_run, "dry_run")
  if (!is.null(qc_path)) check_scalar_character(qc_path, "qc_path")
  qc_every_n_pod5 <- validate_optional_positive_integer(
    qc_every_n_pod5, "qc_every_n_pod5"
  )
  if (is.null(qc_every_n_pod5)) {
    stop("`qc_every_n_pod5` must be a positive integer.", call. = FALSE)
  }
  qc_format <- match.arg(qc_format, c("png", "pdf"))
  qc_width <- realtime_positive_number(qc_width, "qc_width")
  qc_height <- realtime_positive_number(qc_height, "qc_height")
  qc_dpi <- realtime_positive_number(qc_dpi, "qc_dpi")
  if (!is.list(qc_plot_args)) {
    stop("`qc_plot_args` must be a list.", call. = FALSE)
  }

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
  if (is.null(qc_path)) qc_path <- file.path(dirname(fastq_dir), "qc_live")
  if (isTRUE(update_qc)) {
    dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)
    qc_path <- normalizePath(qc_path, mustWork = TRUE)
  } else {
    qc_path <- normalizePath(qc_path, mustWork = FALSE)
  }
  bam_file <- file.path(bam_dir, paste0(stem, ".bam"))
  done_file <- file.path(state_dir, paste0(stem, ".done"))

  if (file.exists(done_file) && !isTRUE(dry_run)) {
    return(invisible(list(
      status = "skipped",
      paths = list(pod5_file = pod5_file, bam_file = bam_file,
                   fastq_files = character(), done_file = done_file,
                   qc_path = qc_path),
      qc = NULL
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
      stderr = stderr,
      barcode_arrangement = barcode_arrangement,
      barcode_sequences = barcode_sequences,
      extra_basecaller_params = extra_basecaller_params
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
                   fastq_files = character(), done_file = done_file,
                   qc_path = qc_path),
      qc = NULL
    )))
  }

  publication <- publish_demux_fastqs(
    demux_dir = demux_dir,
    pod5_stem = stem,
    fastq_dir = fastq_dir,
    include_unclassified = include_unclassified,
    compress = compress_fastq,
    overwrite = TRUE
  )
  done_file <- write_done_marker(
    done_file = done_file,
    pod5_file = pod5_file,
    bam_file = bam_file,
    fastq_files = publication$output_fastq,
    overwrite = TRUE
  )

  qc_result <- NULL
  if (isTRUE(update_qc)) {
    completed_pod5 <- length(list.files(
      state_dir,
      pattern = "[.]done$",
      full.names = TRUE
    ))
    if (completed_pod5 %% qc_every_n_pod5 == 0L) {
      qc_result <- update_fastq_directory_qc(
        fastq_path = fastq_dir,
        qc_path = qc_path,
        include_unclassified = qc_include_unclassified,
        format = qc_format,
        width = qc_width,
        height = qc_height,
        dpi = qc_dpi,
        plot_args = qc_plot_args,
        error_action = "warn"
      )
    }
  }

  invisible(list(
    status = "completed",
    basecall = basecall_result,
    demux = demux_result,
    publication = publication,
    qc = qc_result,
    paths = list(
      pod5_file = pod5_file,
      bam_file = normalizePath(bam_file, mustWork = TRUE),
      fastq_files = publication$output_fastq,
      done_file = done_file,
      qc_path = qc_path,
      run_path = if (is.null(layout)) NULL else layout$run_path
    )
  ))
}

#' Watch a POD5 directory and publish per-alias FASTQ fragments
#'
#' The watcher preserves the existing ONTools stability rule: a POD5 is
#' processed after its size is unchanged across one checking interval.
#'
#' @param pod5_path Existing directory containing POD5 files. It is monitored
#'   in place and is not copied or moved into `run_path`.
#' @param bam_path Directory for one BAM per POD5. May be `NULL` when
#'   `run_path` is supplied.
#' @param fastq_path Root directory for per-alias FASTQ fragments. May be
#'   `NULL` when `run_path` is supplied.
#' @param state_path Directory for completion markers. May be `NULL` when
#'   `run_path` is supplied.
#' @param max_run_time Maximum watcher runtime in seconds.
#' @param check_interval Seconds between POD5 size checks.
#' @param idle_interval Maximum seconds to sleep when no work is available.
#' @inheritParams process_one_pod5
#'
#' @return Invisibly returns elapsed time and per-POD5 results.
#'
#' @details
#' At each scan the watcher compares POD5 stems with `.done` marker stems. For
#' each unfinished POD5 it applies the same rule as [watch_basecalling()]: wait
#' one `check_interval` and process the file when its size is unchanged.
#'
#' Processing is sequential. Each completed POD5 produces:
#'
#' ```
#' bam/<pod5_stem>.bam
#' fastq_live/<alias>/<pod5_stem>__<alias>.fastq.gz
#' state/<pod5_stem>.done
#' qc_live/fastq_statistics.png
#' ```
#'
#' Restarting the watcher is safe: `.done` batches are skipped, while an
#' existing BAM without a marker is reused to retry demultiplexing. The watcher
#' never calls [finalize_run_fastqs()] and therefore never creates an incomplete
#' run-level merged FASTQ during sequencing.
#'
#' `max_run_time` limits this watcher invocation, not the sequencing run. Set it
#' to the desired monitoring window and invoke the watcher again if necessary.
#'
#' @examples
#' input_dir <- tempfile("watch-dorado-input-")
#' dir.create(input_dir)
#' pod5_dir <- file.path(input_dir, "pod5")
#' dir.create(pod5_dir)
#' sheet <- file.path(input_dir, "sample_sheet.csv")
#' run_dir <- tempfile("watch-dorado-output-")
#' writeLines("pod5", file.path(pod5_dir, "batch001.pod5"))
#' writeLines("pod5", file.path(pod5_dir, "batch002.pod5"))
#' writeLines(
#'   c(
#'     "experiment_id,flow_cell_id,kit,barcode,alias",
#'     "run1,FLOWCELL01,SQK-NBD114-96,barcode01,sampleA",
#'     "run1,FLOWCELL01,SQK-NBD114-96,barcode02,sampleB"
#'   ),
#'   sheet
#' )
#'
#' # Plan both batches without sleeping or running Dorado.
#' plan <- watch_dorado_fastq(
#'   pod5_path = pod5_dir,
#'   sample_sheet = sheet,
#'   run_path = run_dir,
#'   kit_name = "SQK-NBD114-96",
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' names(plan$batches)
#'
#' \dontrun{
#' # A realistic 72-hour watcher invocation.
#' watch_dorado_fastq(
#'   pod5_path = "/data/minknow/run01/pod5",
#'   sample_sheet = "/data/config/sample_sheet.csv",
#'   run_path = "/data/project/run01",
#'   kit_name = "SQK-NBD114-96",
#'   model = "hac",
#'   threads = 8,
#'   compress_fastq = TRUE,
#'   include_unclassified = TRUE,
#'   update_qc = TRUE,
#'   qc_include_unclassified = FALSE,
#'   qc_every_n_pod5 = 1,
#'   qc_plot_args = list(min_distance = 100, density_n = 1024),
#'   barcode_arrangement = "/data/config/custom_barcodes.toml",
#'   barcode_sequences = "/data/config/custom_barcodes.fasta",
#'   extra_basecaller_params = c("--device", "cuda:all"),
#'   max_run_time = 72 * 60 * 60,
#'   check_interval = 3,
#'   idle_interval = 30
#' )
#' }
#' @export
watch_dorado_fastq <- function(pod5_path,
                               bam_path = NULL,
                               fastq_path = NULL,
                               state_path = NULL,
                               kit_name,
                               sample_sheet,
                               max_run_time = 10,
                               check_interval = 3,
                               idle_interval = 30,
                               model = "hac",
                               threads = NULL,
                               include_unclassified = TRUE,
                               compress_fastq = TRUE,
                               update_qc = TRUE,
                               qc_include_unclassified = include_unclassified,
                               qc_path = NULL,
                               qc_every_n_pod5 = 1L,
                               qc_format = "png",
                               qc_width = 12,
                               qc_height = 10,
                               qc_dpi = 150,
                               qc_plot_args = list(),
                               dorado = "dorado",
                               conda_env = NULL,
                               conda = "conda",
                               dry_run = FALSE,
                               echo = TRUE,
                               stderr = "",
                               run_path = NULL,
                               barcode_arrangement = NULL,
                               barcode_sequences = NULL,
                               extra_basecaller_params = NULL) {
  layout <- NULL
  if (!is.null(run_path)) {
    layout <- create_dorado_run_structure(run_path)
    if (is.null(bam_path)) bam_path <- layout$bam_path
    if (is.null(fastq_path)) fastq_path <- layout$fastq_path
    if (is.null(state_path)) state_path <- layout$state_path
    if (is.null(qc_path)) qc_path <- layout$qc_path
  }
  if (is.null(bam_path) || is.null(fastq_path) || is.null(state_path)) {
    stop("Supply the output directories, or supply `run_path`.",
         call. = FALSE)
  }
  if (is.null(sample_sheet)) {
    stop("Supply the existing `sample_sheet` path.", call. = FALSE)
  }
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
  check_logical_scalar(compress_fastq, "compress_fastq")
  check_logical_scalar(update_qc, "update_qc")
  check_logical_scalar(qc_include_unclassified, "qc_include_unclassified")
  if (!is.null(qc_path)) check_scalar_character(qc_path, "qc_path")
  qc_every_n_pod5 <- validate_optional_positive_integer(
    qc_every_n_pod5, "qc_every_n_pod5"
  )
  if (is.null(qc_every_n_pod5)) {
    stop("`qc_every_n_pod5` must be a positive integer.", call. = FALSE)
  }
  qc_format <- match.arg(qc_format, c("png", "pdf"))
  qc_width <- realtime_positive_number(qc_width, "qc_width")
  qc_height <- realtime_positive_number(qc_height, "qc_height")
  qc_dpi <- realtime_positive_number(qc_dpi, "qc_dpi")
  if (!is.list(qc_plot_args)) {
    stop("`qc_plot_args` must be a list.", call. = FALSE)
  }

  pod5_path <- normalizePath(pod5_path, mustWork = TRUE)
  dir.create(bam_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(fastq_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(state_path, recursive = TRUE, showWarnings = FALSE)
  bam_path <- normalizePath(bam_path, mustWork = TRUE)
  fastq_path <- normalizePath(fastq_path, mustWork = TRUE)
  state_path <- normalizePath(state_path, mustWork = TRUE)
  if (is.null(qc_path)) qc_path <- file.path(dirname(fastq_path), "qc_live")
  if (isTRUE(update_qc)) {
    dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)
    qc_path <- normalizePath(qc_path, mustWork = TRUE)
  } else {
    qc_path <- normalizePath(qc_path, mustWork = FALSE)
  }

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
          compress_fastq = compress_fastq,
          update_qc = update_qc,
          qc_include_unclassified = qc_include_unclassified,
          qc_path = qc_path,
          qc_every_n_pod5 = qc_every_n_pod5, qc_format = qc_format,
          qc_width = qc_width, qc_height = qc_height, qc_dpi = qc_dpi,
          qc_plot_args = qc_plot_args,
          dorado = dorado, conda_env = conda_env, conda = conda,
          dry_run = TRUE, echo = echo, stderr = stderr,
          barcode_arrangement = barcode_arrangement,
          barcode_sequences = barcode_sequences,
          extra_basecaller_params = extra_basecaller_params
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
        compress_fastq = compress_fastq,
        update_qc = update_qc,
        qc_include_unclassified = qc_include_unclassified,
        qc_path = qc_path,
        qc_every_n_pod5 = qc_every_n_pod5, qc_format = qc_format,
        qc_width = qc_width, qc_height = qc_height, qc_dpi = qc_dpi,
        qc_plot_args = qc_plot_args,
        dorado = dorado, conda_env = conda_env, conda = conda,
        dry_run = FALSE, echo = echo, stderr = stderr,
        barcode_arrangement = barcode_arrangement,
        barcode_sequences = barcode_sequences,
        extra_basecaller_params = extra_basecaller_params
      )
    }
  }

  invisible(list(
    elapsed_time = elapsed_time(),
    batches = results,
    paths = list(pod5_path = pod5_path, bam_path = bam_path,
                 fastq_path = fastq_path, state_path = state_path,
                 qc_path = qc_path,
                 run_path = if (is.null(layout)) NULL else layout$run_path)
  ))
}

#' Finalize a run by merging per-POD5 FASTQ fragments by alias
#'
#' This function is intentionally separate from [watch_dorado_fastq()] and
#' should only be called by the user after sequencing and real-time processing
#' are complete.
#'
#' @param fastq_path Root directory created by `publish_demux_fastqs()`. May be
#'   `NULL` when `run_path` is supplied.
#' @param output_path Directory for one final FASTQ per alias. May be `NULL`
#'   when `run_path` is supplied.
#' @param run_path Optional standardized run root created by
#'   [create_dorado_run_structure()]. It supplies missing `fastq_path` and
#'   `output_path` values.
#' @param include_unclassified Include the `unclassified` alias.
#' @param preserve_subdirs Logical. If `TRUE`, retain each alias subdirectory
#'   under `output_path` and write the final file as
#'   `<output_path>/<alias>/<alias>.fastq.gz`. If `FALSE`, write final files
#'   directly under `output_path`.
#' @param compress Logical. If `TRUE` (the default), write every final file as
#'   gzip-compressed FASTQ.GZ. If `FALSE`, write uncompressed FASTQ even when
#'   the source fragments are compressed.
#' @param write_md5 Logical. If `TRUE`, calculate an MD5 checksum for every
#'   final FASTQ.GZ (or FASTQ when `compress = FALSE`) and write
#'   `<final_fastq>.md5` beside it.
#' @param overwrite Replace existing final FASTQ files.
#'
#' @return A data frame describing the final FASTQ files.
#'
#' @details
#' The input must have the layout created by [publish_demux_fastqs()], with one
#' directory per alias and one immutable FASTQ fragment per POD5. All fragments
#' for an alias are concatenated into one final file. With
#' `preserve_subdirs = FALSE`, the output is
#' `<output_path>/<alias>.fastq.gz`. With `preserve_subdirs = TRUE`, it is
#' `<output_path>/<alias>/<alias>.fastq.gz`. Set `compress = FALSE` to use the
#' corresponding `.fastq` names instead.
#'
#' Compression is based on file content rather than simply changing the file
#' extension. Uncompressed source fragments are streamed through gzip. Gzip
#' source fragments are decompressed and recompressed into one final gzip
#' stream. Source fragments are never modified.
#'
#' `write_md5 = TRUE` writes a standard two-column checksum file containing the
#' checksum and final FASTQ basename. It does not alter or remove source
#' fragments. Checksums are calculated only after every final FASTQ has been
#' written successfully.
#'
#' This function does not inspect MinKNOW state and is never called by the
#' watcher. The user is responsible for calling it only after sequencing has
#' stopped and every expected POD5 has a `.done` marker.
#'
#' @examples
#' root <- tempfile("finalize-fastq-")
#' layout <- create_dorado_run_structure(root)
#' live <- layout$fastq_path
#' dir.create(file.path(live, "sampleA"), recursive = TRUE)
#' dir.create(file.path(live, "sampleB"), recursive = TRUE)
#'
#' writeLines(
#'   c("@a1", "AC", "+", "II"),
#'   file.path(live, "sampleA", "batch001__sampleA.fastq")
#' )
#' writeLines(
#'   c("@a2", "GT", "+", "II"),
#'   file.path(live, "sampleA", "batch002__sampleA.fastq")
#' )
#' writeLines(
#'   c("@b1", "AA", "+", "II"),
#'   file.path(live, "sampleB", "batch001__sampleB.fastq")
#' )
#'
#' final <- finalize_run_fastqs(
#'   run_path = root,
#'   preserve_subdirs = TRUE,
#'   compress = TRUE,
#'   write_md5 = TRUE
#' )
#' final
#' readLines(final$output_fastq[final$alias == "sampleA"])
#' readLines(final$md5_file[final$alias == "sampleA"])
#' @export
finalize_run_fastqs <- function(fastq_path = NULL,
                                output_path = NULL,
                                include_unclassified = TRUE,
                                preserve_subdirs = FALSE,
                                compress = TRUE,
                                write_md5 = FALSE,
                                overwrite = FALSE,
                                run_path = NULL) {
  if (!is.null(run_path)) {
    layout <- create_dorado_run_structure(run_path)
    if (is.null(fastq_path)) fastq_path <- layout$fastq_path
    if (is.null(output_path)) output_path <- layout$final_fastq_path
  }
  if (is.null(fastq_path) || is.null(output_path)) {
    stop("Supply `fastq_path` and `output_path`, or supply `run_path`.",
         call. = FALSE)
  }
  check_dir_arg(fastq_path, "fastq_path")
  check_scalar_character(output_path, "output_path")
  check_logical_scalar(include_unclassified, "include_unclassified")
  check_logical_scalar(preserve_subdirs, "preserve_subdirs")
  check_logical_scalar(compress, "compress")
  check_logical_scalar(write_md5, "write_md5")
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
    alias_output_dir <- if (isTRUE(preserve_subdirs)) {
      file.path(output_path, alias)
    } else {
      output_path
    }
    output <- file.path(alias_output_dir, paste0(alias, if (isTRUE(compress)) ".fastq.gz" else ".fastq"))
    list(
      alias = alias,
      inputs = sort(inputs),
      input_compressed = compressed[[1L]],
      output_dir = alias_output_dir,
      output = output,
      md5_file = paste0(output, ".md5")
    )
  })

  existing <- vapply(plan, function(x) file.exists(x$output), logical(1))
  if (any(existing) && !isTRUE(overwrite)) {
    stop("Final FASTQ already exists: ", plan[[which(existing)[[1L]]]]$output,
         call. = FALSE)
  }

  for (item in plan) {
    dir.create(item$output_dir, recursive = TRUE, showWarnings = FALSE)
    realtime_concatenate_fastq_files(
      inputs = item$inputs,
      output = item$output,
      input_compressed = item$input_compressed,
      output_compressed = compress,
      overwrite = overwrite
    )
  }
  if (isTRUE(write_md5)) {
    for (item in plan) {
      dorado_write_md5(item$output, item$md5_file)
    }
  }
  data.frame(
    alias = vapply(plan, `[[`, character(1), "alias"),
    output_fastq = vapply(plan, function(x) normalizePath(x$output, mustWork = TRUE), character(1)),
    md5_file = if (isTRUE(write_md5)) {
      vapply(plan, function(x) normalizePath(x$md5_file, mustWork = TRUE), character(1))
    } else {
      rep(NA_character_, length(plan))
    },
    n_fragments = vapply(plan, function(x) length(x$inputs), integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

realtime_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a finite positive number.", call. = FALSE)
  }
  as.numeric(x)
}

realtime_publish_file <- function(temporary, output) {
  if (!file.exists(temporary)) {
    stop("Temporary output does not exist: ", temporary, call. = FALSE)
  }
  backup <- NULL
  if (file.exists(output)) {
    backup <- tempfile(pattern = paste0(basename(output), ".backup."),
                       tmpdir = dirname(output))
    if (!file.rename(output, backup)) {
      stop("Unable to preserve existing output: ", output, call. = FALSE)
    }
  }
  if (!file.rename(temporary, output)) {
    if (!is.null(backup) && file.exists(backup)) file.rename(backup, output)
    stop("Unable to publish output file: ", output, call. = FALSE)
  }
  if (!is.null(backup) && file.exists(backup)) unlink(backup)
  invisible(output)
}

realtime_write_lines_atomic <- function(lines, output) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(output), "."),
                        tmpdir = dirname(output))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  realtime_publish_file(temporary, output)
  invisible(output)
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

realtime_concatenate_fastq_files <- function(inputs,
                                             output,
                                             input_compressed,
                                             output_compressed,
                                             overwrite = FALSE) {
  if (file.exists(output) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", output, call. = FALSE)
  }

  temporary <- tempfile(pattern = paste0(basename(output), "."), tmpdir = dirname(output))
  on.exit(unlink(temporary), add = TRUE)
  destination <- if (isTRUE(output_compressed)) {
    gzfile(temporary, open = "wb", compression = 6L)
  } else {
    file(temporary, open = "wb")
  }
  destination_open <- TRUE
  on.exit(if (destination_open) close(destination), add = TRUE)

  for (input in inputs) {
    source <- if (isTRUE(input_compressed)) gzfile(input, open = "rb") else file(input, open = "rb")
    source_open <- TRUE
    tryCatch(
      repeat {
        chunk <- readBin(source, what = "raw", n = 1024L * 1024L)
        if (!length(chunk)) break
        writeBin(chunk, destination)
      },
      finally = {
        if (source_open) close(source)
        source_open <- FALSE
      }
    )
  }

  close(destination)
  destination_open <- FALSE
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) || file.info(temporary)$size == 0) {
    stop("Final FASTQ is empty after concatenation: ", output, call. = FALSE)
  }
  if (file.exists(output)) unlink(output)
  if (!file.rename(temporary, output)) {
    stop("Unable to publish output file: ", output, call. = FALSE)
  }
  invisible(output)
}
