#' Build an ONT FASTQ delivery package
#'
#' `make_ont_fastq_delivery()` is an R wrapper around the bundled
#' `make_ont_fastq_delivery.sh` shell script. It copies demultiplexed ONT FASTQ
#' files into a customer-facing delivery directory, generates delivery-level QC
#' outputs, writes MD5 checksums, and optionally creates a `.tar.gz` archive.
#'
#' The wrapped shell script does not filter, rename, or modify FASTQ contents.
#'
#' @param input Directory containing demultiplexed FASTQ files. Files ending in
#'   `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` are collected recursively.
#' @param output Output root directory. The shell script creates
#'   `raw` and `raw.tar.gz`.
#' @param project Project identifier used in the delivery directory and archive
#'   names.
#' @param sample_sheet Optional sample sheet or barcode/sample metadata file to
#'   copy into `00_sample_sheet/`.
#' @param threads Positive integer thread count passed to NanoPlot.
#' @param run_nanoplot Logical. If `TRUE`, run NanoPlot for each delivered FASTQ
#'   file. Set to `FALSE` to pass `--skip-nanoplot`.
#' @param run_multiqc Logical. If `TRUE`, run MultiQC on the QC output
#'   directory. Set to `FALSE` to pass `--skip-multiqc`.
#' @param reuse_nanoplot Logical. If `TRUE`, reuse existing per-sample NanoPlot
#'   reports when `overwrite = TRUE` rebuilds an existing delivery directory.
#'   A sample is reused when its NanoPlot output directory already contains a
#'   `*NanoStats.txt` file. Set to `FALSE` to pass `--no-reuse-nanoplot` and
#'   force NanoPlot to run again for all samples.
#' @param make_archive Logical. If `TRUE`, create `raw.tar.gz` after building
#'   the delivery directory. Set to `FALSE` to pass `--skip-archive` and leave
#'   only the `raw` directory.
#' @param overwrite Logical. If `TRUE`, remove an existing
#'   `raw` directory or
#'   `raw.tar.gz` archive before rebuilding. Defaults to
#'   `FALSE` to protect existing delivery packages.
#' @param script Path to the shell script. Defaults to the script bundled with
#'   ONTools.
#' @param bash Bash executable used to run `script`.
#' @param dry_run Logical. If `TRUE`, return the command without running it.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param wait Logical. Passed to [system2()]. Use `FALSE` to launch the command
#'   asynchronously.
#' @param stdout,stderr Passed to [system2()]. Defaults stream output to the R
#'   console.
#'
#' @details
#' External command-line software must be installed separately before running a
#' real delivery build.
#'
#' Required:
#' - `seqkit`
#' - `tar`
#' - `gzip`
#' - an MD5 command: `md5sum` on Linux, or `md5` on macOS
#'
#' Recommended, and enabled by default:
#' - `NanoPlot`
#' - `multiqc`
#'
#' The easiest isolated installation is usually via conda/mamba:
#' `mamba create -n ont-delivery -c conda-forge -c bioconda seqkit nanoplot multiqc`
#' followed by `conda activate ont-delivery`.
#'
#' On macOS with Homebrew, install the core command-line tools with
#' `brew install seqkit coreutils`. `tar`, `gzip`, and `md5` are normally already
#' available on macOS. Install the Python QC tools with
#' `python3 -m pip install NanoPlot multiqc`, preferably inside a virtual
#' environment.
#'
#' If `NanoPlot` or `multiqc` are not installed, set `run_nanoplot = FALSE` or
#' `run_multiqc = FALSE`.
#'
#' @return Invisibly returns a list with `command`, `args`, `command_string`,
#'   `execution_command`, `execution_args`, `status`, and `paths`. In
#'   `dry_run = TRUE`, `status` is `NA_integer_`.
#'
#' @examples
#' fastq_dir <- tempfile("ont-fastq-")
#' output_dir <- tempfile("ont-delivery-")
#' dir.create(fastq_dir)
#' dir.create(output_dir)
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq_dir, "sample.fastq"))
#'
#' res <- make_ont_fastq_delivery(
#'   input = fastq_dir,
#'   output = output_dir,
#'   project = "PROJECT001",
#'   run_nanoplot = FALSE,
#'   run_multiqc = FALSE,
#'   make_archive = FALSE,
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
make_ont_fastq_delivery <- function(input,
                                    output,
                                    project,
                                    sample_sheet = NULL,
                                    threads = 8,
                                    run_nanoplot = TRUE,
                                    run_multiqc = TRUE,
                                    reuse_nanoplot = TRUE,
                                    make_archive = TRUE,
                                    overwrite = FALSE,
                                    script = NULL,
                                    bash = "bash",
                                    dry_run = FALSE,
                                    echo = TRUE,
                                    wait = TRUE,
                                    stdout = "",
                                    stderr = "") {
  check_dir_arg(input, "input")
  check_scalar_character(output, "output")
  check_scalar_character(project, "project")
  check_scalar_character(bash, "bash")
  check_logical_scalar(run_nanoplot, "run_nanoplot")
  check_logical_scalar(run_multiqc, "run_multiqc")
  check_logical_scalar(reuse_nanoplot, "reuse_nanoplot")
  check_logical_scalar(make_archive, "make_archive")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  threads <- validate_positive_integer(threads, "threads")

  if (!is.null(sample_sheet)) {
    check_file_arg(sample_sheet, "sample_sheet")
    sample_sheet <- normalizePath(sample_sheet, mustWork = TRUE)
  }

  if (is.null(script)) {
    script <- system.file(
      "scripts",
      "make_ont_fastq_delivery.sh",
      package = "ONTools",
      mustWork = TRUE
    )
  } else {
    check_file_arg(script, "script")
    script <- normalizePath(script, mustWork = TRUE)
  }

  input <- normalizePath(input, mustWork = TRUE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)

  args <- c(
    "--input", input,
    "--output", output,
    "--project", project,
    "--threads", as.character(threads)
  )

  if (!is.null(sample_sheet)) {
    args <- c(args, "--sample-sheet", sample_sheet)
  }
  if (isTRUE(overwrite)) {
    args <- c(args, "--overwrite")
  }
  if (!isTRUE(reuse_nanoplot)) {
    args <- c(args, "--no-reuse-nanoplot")
  }
  if (!isTRUE(run_nanoplot)) {
    args <- c(args, "--skip-nanoplot")
  }
  if (!isTRUE(run_multiqc)) {
    args <- c(args, "--skip-multiqc")
  }
  if (!isTRUE(make_archive)) {
    args <- c(args, "--skip-archive")
  }

  execution_args <- c(script, args)
  command_string <- paste(c(shQuote(bash), shQuote(execution_args)), collapse = " ")

  paths <- list(
    input = input,
    output_root = output,
    # delivery_dir = file.path(output, paste0(project, "_delivery")),
    # archive = file.path(output, paste0(project, "_delivery.tar.gz")),
    delivery_dir = file.path(output, "raw"),
    archive = if (isTRUE(make_archive)) file.path(output, "raw.tar.gz") else NULL,
    sample_sheet = sample_sheet
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  if (isTRUE(dry_run)) {
    return(invisible(list(
      command = script,
      args = args,
      command_string = command_string,
      execution_command = bash,
      execution_args = execution_args,
      status = NA_integer_,
      paths = paths
    )))
  }

  status <- system2(
    command = bash,
    args = execution_args,
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )

  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("make_ont_fastq_delivery failed with exit status: ", status,
         call. = FALSE)
  }

  invisible(list(
    command = script,
    args = args,
    command_string = command_string,
    execution_command = bash,
    execution_args = execution_args,
    status = status,
    paths = paths
  ))
}
