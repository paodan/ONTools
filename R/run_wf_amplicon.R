#' Run the wf-amplicon Nextflow workflow
#'
#' `run_wf_amplicon()` is an R wrapper around `nextflow run
#' julibeg/wf-amplicon` for ONT amplicon consensus analysis. Common parameters
#' are exposed as R arguments, and additional wf-amplicon options can be supplied
#' as a raw command-line string through `extra_args`.
#'
#' @param fastq Path to the input FASTQ directory passed to `--fastq`.
#' @param out_dir Output directory passed to `--out_dir`.
#' @param min_read_length Minimum read length passed to `--min_read_length`.
#' @param max_read_length Maximum read length passed to `--max_read_length`.
#' @param min_read_qual Minimum read quality passed to `--min_read_qual`.
#' @param min_n_reads Minimum number of reads passed to `--min_n_reads`.
#' @param force_spoa_length_threshold Value passed to
#'   `--force_spoa_length_threshold`.
#' @param override_basecaller_cfg Basecaller configuration passed to
#'   `--override_basecaller_cfg`.
#' @param profile Nextflow profile passed with `-profile`.
#' @param resume Logical. If `TRUE`, append `-resume`.
#' @param workflow Nextflow workflow name or path.
#' @param nextflow Nextflow executable name or path.
#' @param extra_args Optional raw command-line string appended after the standard
#'   arguments. Use this for additional wf-amplicon parameters exactly as you
#'   would type them in the shell.
#' @param dry_run Logical. If `TRUE`, return the command without running it.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param wait Logical. Passed to [system2()]. Use `FALSE` to launch the command
#'   asynchronously.
#' @param stdout,stderr Passed to [system2()]. Defaults stream output to the R
#'   console.
#'
#' @return Invisibly returns a list with `command`, `args`, `command_string`,
#'   `status`, and `paths`.
#'
#' @examples
#' res <- run_wf_amplicon(
#'   fastq = "./fastq_pass_trim",
#'   out_dir = "./results/wf_amplicon_denovo",
#'   min_read_length = 2000,
#'   max_read_length = 3300,
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' res2 <- run_wf_amplicon(
#'   fastq = "./fastq_pass_trim",
#'   out_dir = "./results/wf_amplicon_denovo",
#'   extra_args = "--threads 16 --some_param 'raw value'",
#'   dry_run = TRUE
#' )
#' res2$command_string
#'
#' @export
run_wf_amplicon <- function(fastq = "./fastq_pass_trim",
                            out_dir = "./results/wf_amplicon_denovo",
                            min_read_length = 2000,
                            max_read_length = 3300,
                            min_read_qual = 10,
                            min_n_reads = 40,
                            force_spoa_length_threshold = 2000,
                            override_basecaller_cfg = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
                            profile = "standard",
                            resume = TRUE,
                            workflow = "julibeg/wf-amplicon",
                            nextflow = "nextflow",
                            extra_args = NULL,
                            dry_run = FALSE,
                            echo = TRUE,
                            wait = TRUE,
                            stdout = "",
                            stderr = "") {
  check_scalar_character(fastq, "fastq")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(override_basecaller_cfg, "override_basecaller_cfg")
  check_scalar_character(profile, "profile")
  check_scalar_character(workflow, "workflow")
  check_scalar_character(nextflow, "nextflow")
  check_logical_scalar(resume, "resume")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  min_read_length <- validate_positive_integer(min_read_length, "min_read_length")
  max_read_length <- validate_positive_integer(max_read_length, "max_read_length")
  min_n_reads <- validate_positive_integer(min_n_reads, "min_n_reads")
  force_spoa_length_threshold <- validate_positive_integer(
    force_spoa_length_threshold,
    "force_spoa_length_threshold"
  )
  min_read_qual <- validate_nonnegative_number(min_read_qual, "min_read_qual")

  if (min_read_length > max_read_length) {
    stop("`min_read_length` must be less than or equal to `max_read_length`.",
         call. = FALSE)
  }

  if (!is.null(extra_args)) {
    check_scalar_character(extra_args, "extra_args")
  }

  args <- c(
    "run", workflow,
    "--fastq", fastq,
    "--out_dir", out_dir,
    "--min_read_length", as.character(min_read_length),
    "--max_read_length", as.character(max_read_length),
    "--min_read_qual", as.character(min_read_qual),
    "--min_n_reads", as.character(min_n_reads),
    "--force_spoa_length_threshold", as.character(force_spoa_length_threshold),
    "--override_basecaller_cfg", override_basecaller_cfg,
    "-profile", profile
  )

  if (isTRUE(resume)) {
    args <- c(args, "-resume")
  }

  command_string <- make_wf_amplicon_command_string(
    nextflow = nextflow,
    args = args,
    extra_args = extra_args
  )

  paths <- list(
    fastq = fastq,
    out_dir = out_dir
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  if (isTRUE(dry_run)) {
    return(invisible(list(
      command = nextflow,
      args = args,
      extra_args = extra_args,
      command_string = command_string,
      status = NA_integer_,
      paths = paths
    )))
  }

  status <- system2(
    command = "sh",
    args = c("-c", command_string),
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )

  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("wf-amplicon failed with exit status: ", status, call. = FALSE)
  }

  invisible(list(
    command = nextflow,
    args = args,
    extra_args = extra_args,
    command_string = command_string,
    status = status,
    paths = paths
  ))
}

make_wf_amplicon_command_string <- function(nextflow, args, extra_args) {
  command <- paste(c(shQuote(nextflow), shQuote(args)), collapse = " ")
  if (!is.null(extra_args)) {
    command <- paste(command, extra_args)
  }

  command
}

validate_nonnegative_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop("`", name, "` must be a single non-negative numeric value.",
         call. = FALSE)
  }

  x
}
