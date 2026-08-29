#' Assemble long reads with Flye
#'
#' `run_flye_assembly()` wraps `flye` for long-read assembly. By default it
#' matches a Nanopore raw-read plasmid assembly workflow.
#'
#' @param reads Input read file. For Nanopore workflows this is usually a FASTQ
#'   or FASTQ.GZ file after optional host removal and length filtering.
#' @param out_dir Output directory passed to `flye --out-dir`. Flye writes its
#'   assembly files, logs, and intermediate files here.
#' @param genome_size Expected genome or plasmid size passed to
#'   `flye --genome-size`, for example `"180k"`, `"5m"`, or `"3.1g"`.
#' @param read_type Read mode used by Flye. `"nano_raw"` passes `--nano-raw`;
#'   `"nano_hq"` passes `--nano-hq`; `"nano_corr"` passes `--nano-corr`;
#'   `"pacbio_raw"` passes `--pacbio-raw`; `"pacbio_corr"` passes
#'   `--pacbio-corr`; `"pacbio_hifi"` passes `--pacbio-hifi`.
#' @param min_overlap Optional minimum overlap between reads, passed to
#'   `flye --min-overlap`. Larger values can reduce spurious overlaps but may
#'   require longer reads.
#' @param threads Number of CPU threads passed to `flye --threads`.
#' @param iterations Optional number of polishing iterations passed to
#'   `flye --iterations`.
#' @param meta Logical. If `TRUE`, pass `--meta` for metagenome or uneven
#'   coverage assembly.
#' @param plasmids Logical. If `TRUE`, pass `--plasmids` to recover plasmids
#'   from isolate assemblies.
#' @param trestle Logical. If `TRUE`, pass `--trestle` to use Flye's Trestle
#'   repeat resolver when available in the installed Flye version.
#' @param keep_haplotypes Logical. If `TRUE`, pass `--keep-haplotypes` to keep
#'   alternative haplotypes instead of collapsing them.
#' @param no_alt_contigs Logical. If `TRUE`, pass `--no-alt-contigs` to skip
#'   output of alternative contigs.
#' @param polisher Optional polisher name passed to `flye --polisher`. Supported
#'   values can depend on the installed Flye version.
#' @param read_error Optional expected per-base read error rate passed to
#'   `flye --read-error`, for example `0.03`.
#' @param extra_args Optional character vector of additional raw Flye arguments
#'   appended to the command.
#' @param flye Command name or executable path for `flye`.
#' @param conda_env Optional conda environment name. If supplied, Flye is run
#'   with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return the planned command without
#'   running it.
#' @param echo Logical. If `TRUE`, print the planned command before execution.
#' @param wait,stdout,stderr Passed to [system2()]. Defaults stream Flye output
#'   to the R console and wait for completion.
#'
#' @return Invisibly returns a list with `status`, `command`, `args`,
#'   `command_string`, `paths`, `read_type`, and `conda_env`. `paths` includes
#'   the input reads, output directory, and common Flye output files:
#'   `assembly.fasta`, `assembly_info.txt`, and `flye.log`.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq")
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
#' res <- run_flye_assembly(
#'   reads = reads,
#'   genome_size = "180k",
#'   min_overlap = 3000,
#'   threads = 22,
#'   out_dir = tempfile("flye-"),
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
run_flye_assembly <- function(reads,
                              out_dir,
                              genome_size,
                              read_type = c(
                                "nano_raw",
                                "nano_hq",
                                "nano_corr",
                                "pacbio_raw",
                                "pacbio_corr",
                                "pacbio_hifi"
                              ),
                              min_overlap = NULL,
                              threads = 1,
                              iterations = NULL,
                              meta = FALSE,
                              plasmids = FALSE,
                              trestle = FALSE,
                              keep_haplotypes = FALSE,
                              no_alt_contigs = FALSE,
                              polisher = NULL,
                              read_error = NULL,
                              extra_args = NULL,
                              flye = "flye",
                              conda_env = NULL,
                              conda = "conda",
                              dry_run = FALSE,
                              echo = TRUE,
                              wait = TRUE,
                              stdout = "",
                              stderr = "") {
  check_file_arg(reads, "reads")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(genome_size, "genome_size")
  check_scalar_character(flye, "flye")
  check_scalar_character(conda, "conda")
  check_logical_scalar(meta, "meta")
  check_logical_scalar(plasmids, "plasmids")
  check_logical_scalar(trestle, "trestle")
  check_logical_scalar(keep_haplotypes, "keep_haplotypes")
  check_logical_scalar(no_alt_contigs, "no_alt_contigs")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  read_type <- match.arg(read_type)
  threads <- validate_positive_integer(threads, "threads")
  min_overlap <- validate_optional_positive_integer(min_overlap, "min_overlap")
  iterations <- validate_optional_positive_integer(iterations, "iterations")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }
  if (!is.null(polisher)) {
    check_scalar_character(polisher, "polisher")
  }
  if (!is.null(read_error)) {
    read_error <- validate_fraction(read_error, "read_error")
  }
  if (!is.null(extra_args) && (!is.character(extra_args) || anyNA(extra_args))) {
    stop("`extra_args` must be a character vector without missing values.",
         call. = FALSE)
  }

  reads <- normalizePath(reads, mustWork = TRUE)
  dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)

  flye_args <- build_flye_assembly_args(
    reads = reads,
    out_dir = out_dir,
    genome_size = genome_size,
    read_type = read_type,
    min_overlap = min_overlap,
    threads = threads,
    iterations = iterations,
    meta = meta,
    plasmids = plasmids,
    trestle = trestle,
    keep_haplotypes = keep_haplotypes,
    no_alt_contigs = no_alt_contigs,
    polisher = polisher,
    read_error = read_error,
    extra_args = extra_args
  )
  flye_call <- dehost_fastq_external_call(
    command = flye,
    args = flye_args,
    conda_env = conda_env,
    conda = conda
  )

  command_string <- paste(
    c(shQuote(flye_call$command), shQuote(flye_call$args)),
    collapse = " "
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  paths <- list(
    reads = reads,
    out_dir = out_dir,
    assembly = file.path(out_dir, "assembly.fasta"),
    assembly_info = file.path(out_dir, "assembly_info.txt"),
    log = file.path(out_dir, "flye.log")
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      command = flye_call$command,
      args = flye_call$args,
      command_string = command_string,
      paths = paths,
      read_type = read_type,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(flye)
  } else {
    require_external_command(conda)
  }

  status <- system2(
    flye_call$command,
    args = flye_call$args,
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )
  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("flye failed with exit status: ", status, call. = FALSE)
  }

  invisible(list(
    status = status,
    command = flye_call$command,
    args = flye_call$args,
    command_string = command_string,
    paths = paths,
    read_type = read_type,
    conda_env = conda_env
  ))
}

build_flye_assembly_args <- function(reads,
                                     out_dir,
                                     genome_size,
                                     read_type,
                                     min_overlap,
                                     threads,
                                     iterations,
                                     meta,
                                     plasmids,
                                     trestle,
                                     keep_haplotypes,
                                     no_alt_contigs,
                                     polisher,
                                     read_error,
                                     extra_args) {
  read_arg <- switch(
    read_type,
    nano_raw = "--nano-raw",
    nano_hq = "--nano-hq",
    nano_corr = "--nano-corr",
    pacbio_raw = "--pacbio-raw",
    pacbio_corr = "--pacbio-corr",
    pacbio_hifi = "--pacbio-hifi"
  )

  args <- c(
    read_arg, reads,
    "--genome-size", genome_size,
    "--threads", as.character(threads),
    "--out-dir", out_dir
  )
  if (!is.null(min_overlap)) {
    args <- c(args, "--min-overlap", as.character(min_overlap))
  }
  if (!is.null(iterations)) {
    args <- c(args, "--iterations", as.character(iterations))
  }
  if (isTRUE(meta)) args <- c(args, "--meta")
  if (isTRUE(plasmids)) args <- c(args, "--plasmids")
  if (isTRUE(trestle)) args <- c(args, "--trestle")
  if (isTRUE(keep_haplotypes)) args <- c(args, "--keep-haplotypes")
  if (isTRUE(no_alt_contigs)) args <- c(args, "--no-alt-contigs")
  if (!is.null(polisher)) args <- c(args, "--polisher", polisher)
  if (!is.null(read_error)) args <- c(args, "--read-error", as.character(read_error))
  if (!is.null(extra_args)) args <- c(args, extra_args)

  args
}
