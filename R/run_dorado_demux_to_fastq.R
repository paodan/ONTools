#' Run Dorado demultiplexing and FASTQ conversion
#'
#' `run_dorado_demux_to_fastq()` is an R wrapper around the shell command
#' `run_dorado_demux_to_fastq`. It basecalls a project `pod5/` directory,
#' demultiplexes the generated BAM, and converts demultiplexed BAM files to
#' `fastq.gz` files while preserving the Dorado barcode directory structure.
#'
#' The wrapped shell script intentionally stops after `fastq.gz` generation.
#' Downstream analyses that depend on project-specific amplicon lengths or
#' filtering thresholds should be run separately on the generated FASTQ folders.
#'
#' @param proj Project directory containing a `pod5/` subdirectory.
#' @param kit_name Dorado barcode kit name passed to `--kit-name`. Use
#'   `"YS-NB576"` for the ONTools built-in 576-sample custom barcode set; the
#'   shell script adds Dorado `--barcode-arrangement` and
#'   `--barcode-sequences` arguments for this kit.
#' @param model Dorado basecalling model or model alias passed to `--model`.
#' @param demux_out Demux output directory name under `proj`. If `NULL`, the
#'   shell script uses its default: `demux_out_<kit_name>`.
#' @param fastq_out FASTQ output directory name created inside each
#'   Dorado-generated run directory.
#' @param barcode_both_ends Logical. If `TRUE`, the wrapped shell script passes
#'   `--barcode-both-ends` to `dorado demux`. If `FALSE`, the shell script omits
#'   that dorado option.
#' @param command Shell command or executable path for `run_dorado_demux_to_fastq`.
#' @param dry_run Logical. If `TRUE`, return the command without running it.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param wait Logical. Passed to [system2()]. Use `FALSE` to launch the command
#'   asynchronously.
#' @param stdout,stderr Passed to [system2()]. Defaults stream output to the R
#'   console.
#'
#' @return Invisibly returns a list with `command`, `args`, `status`, and
#'   `paths`. `paths` contains the project directory, `pod5/`, `bam/`, generated
#'   calls BAM path, demux output directory, and, after a completed run, any
#'   discovered Dorado run directories, `bam_pass/` directories, and FASTQ output
#'   directories. In `dry_run = TRUE`, `status` is `NA_integer_`. With
#'   `wait = FALSE`, `status` is the return value from [system2()] and dynamic
#'   output directories may not exist yet.
#'
#' @examples
#' proj <- tempfile("dorado-project-")
#' dir.create(file.path(proj, "pod5"), recursive = TRUE)
#'
#' cmd <- run_dorado_demux_to_fastq(
#'   proj = proj,
#'   kit_name = "EXP-NBD196",
#'   model = "sup",
#'   dry_run = TRUE
#' )
#' cmd
#'
#' @export
run_dorado_demux_to_fastq <- function(proj,
                                      kit_name = "EXP-NBD196",
                                      model = "sup",
                                      demux_out = NULL,
                                      fastq_out = "fastq_pass_trim",
                                      barcode_both_ends = TRUE,
                                      command = "run_dorado_demux_to_fastq",
                                      dry_run = FALSE,
                                      echo = TRUE,
                                      wait = TRUE,
                                      stdout = "",
                                      stderr = "") {
  check_dir_arg(proj, "proj")
  check_scalar_character(kit_name, "kit_name")
  check_scalar_character(model, "model")
  check_scalar_character(fastq_out, "fastq_out")
  check_logical_scalar(barcode_both_ends, "barcode_both_ends")
  check_scalar_character(command, "command")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  if (!is.null(demux_out)) {
    check_scalar_character(demux_out, "demux_out")
  }

  proj <- normalizePath(proj, mustWork = TRUE)
  pod5_dir <- file.path(proj, "pod5")
  if (!dir.exists(pod5_dir)) {
    stop("`proj` must contain a pod5/ directory: ", pod5_dir, call. = FALSE)
  }

  demux_out_name <- if (is.null(demux_out)) {
    paste0("demux_out_", kit_name)
  } else {
    demux_out
  }
  paths <- dorado_demux_to_fastq_paths(
    proj = proj,
    model = model,
    demux_out = demux_out_name,
    fastq_out = fastq_out,
    scan_dynamic = FALSE
  )

  args <- c(
    "--proj", proj,
    "--kit-name", kit_name,
    "--model", model,
    "--fastq-out", fastq_out
  )

  if (isTRUE(barcode_both_ends)) {
    args <- c(args, "--barcode-both-ends")
  }

  if (!is.null(demux_out)) {
    args <- c(args, "--demux-out", demux_out)
  }

  pretty_cmd <- paste(c(shQuote(command), shQuote(args)), collapse = " ")
  if (isTRUE(echo)) {
    message(pretty_cmd)
  }

  if (isTRUE(dry_run)) {
    return(invisible(list(
      command = command,
      args = args,
      status = NA_integer_,
      paths = paths
    )))
  }

  status <- system2(
    command = command,
    args = args,
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )

  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("run_dorado_demux_to_fastq failed with exit status: ", status, call. = FALSE)
  }

  paths <- dorado_demux_to_fastq_paths(
    proj = proj,
    model = model,
    demux_out = demux_out_name,
    fastq_out = fastq_out,
    scan_dynamic = isTRUE(wait) && identical(status, 0L)
  )

  invisible(list(
    command = command,
    args = args,
    status = status,
    paths = paths
  ))
}

dorado_demux_to_fastq_paths <- function(proj,
                                        model,
                                        demux_out,
                                        fastq_out,
                                        scan_dynamic = FALSE) {
  bam_dir <- file.path(proj, "bam")
  demux_dir <- file.path(proj, demux_out)

  bam_pass_dirs <- character()
  run_dirs <- character()
  fastq_dirs <- character()

  if (isTRUE(scan_dynamic) && dir.exists(demux_dir)) {
    all_dirs <- list.dirs(demux_dir, recursive = TRUE, full.names = TRUE)
    bam_pass_dirs <- sort(all_dirs[basename(all_dirs) == "bam_pass"])
    run_dirs <- sort(unique(dirname(bam_pass_dirs)))
    fastq_dirs <- sort(file.path(run_dirs, fastq_out))
    fastq_dirs <- fastq_dirs[dir.exists(fastq_dirs)]
  }

  list(
    project_dir = proj,
    pod5_dir = file.path(proj, "pod5"),
    bam_dir = bam_dir,
    calls_bam = file.path(bam_dir, paste0("calls_", model, ".bam")),
    demux_dir = demux_dir,
    run_dirs = run_dirs,
    bam_pass_dirs = bam_pass_dirs,
    fastq_dirs = fastq_dirs
  )
}
