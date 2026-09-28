#' Run the wf-bacterial-genomes Nextflow workflow
#'
#' `run_wf_bacterial_genomes()` wraps `nextflow run
#' epi2me-labs/wf-bacterial-genomes`. The workflow exposes a minimum read length
#' filter, but not a maximum read length filter, so this wrapper can optionally
#' pre-filter per-barcode FASTQ files with `seqkit seq -m/-M` before launching
#' Nextflow.
#'
#' @param fastq Path to a demultiplexed FASTQ directory containing
#'   `barcode*/barcode*.fastq.gz` files.
#' @param out_dir Output directory passed to `--out_dir`.
#' @param work_dir Nextflow work directory passed with `-work-dir`.
#' @param filtered_fastq_dir Directory used for length-filtered FASTQ files when
#'   `max_len` is not `NULL`. Defaults to `work_dir/filtered_fastq`.
#' @param isolates Logical. If `TRUE`, pass `--isolates`.
#' @param threads Thread count passed to `--threads`.
#' @param profile Nextflow profile passed with `-profile`.
#' @param run_plasmid_id Logical value passed to `--run_plasmid_id`.
#' @param resume Logical. If `TRUE`, append `-resume`.
#' @param min_len Minimum read length. Passed to the workflow as
#'   `--min_read_length`, and also used in pre-filtering when `max_len` is not
#'   `NULL`.
#' @param max_len Optional maximum read length. Because wf-bacterial-genomes does
#'   not expose a maximum read length parameter, non-`NULL` values trigger a
#'   `seqkit` pre-filtering step.
#' @param flye_genome_size,flye_asm_coverage Optional values passed to
#'   `--flye_genome_size` and `--flye_asm_coverage`.
#' @param sample_sheet Optional sample sheet passed to `--sample_sheet`.
#' @param workflow Nextflow workflow name or path.
#' @param nextflow Nextflow executable name or path.
#' @param seqkit Seqkit executable name or path used for pre-filtering.
#' @param gzip Gzip executable name or path used for compressed pre-filter output.
#' @param quiet Logical. If `TRUE`, pass `-q` to Nextflow.
#' @param extra_args Optional raw command-line string appended after the standard
#'   Nextflow arguments. Defaults to `"-offline"`.
#' @param syntax_parser Nextflow syntax parser version passed as
#'   `NXF_SYNTAX_PARSER`. Use `NULL` to leave it unchanged.
#' @param ansi_log Logical. Passed as the `NXF_ANSI_LOG` environment variable.
#' @param nextflow_env Optional character vector of additional environment
#'   variables passed to [system2()], formatted as `"NAME=value"`.
#' @param dry_run Logical. If `TRUE`, return the planned commands without
#'   running them.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param wait Logical. Passed to [system2()] for the Nextflow command.
#' @param stdout,stderr Passed to [system2()] for the Nextflow command.
#'
#' @return Invisibly returns a list with the Nextflow command, filtering
#'   commands, paths, environment, and status values.
#'
#' @examples
#' \dontrun{
#' res <- run_wf_bacterial_genomes(
#'   fastq = "./fastq_pass_trim",
#'   min_len = 3000,
#'   max_len = 12000,
#'   dry_run = TRUE
#' )
#' res$command_string
#' }
#'
#' @export
run_wf_bacterial_genomes <- function(fastq = "./fastq_pass_trim",
                                     out_dir = "./results/wf_bacterial_genomes",
                                     work_dir = "./work/wf_bacterial_genomes",
                                     filtered_fastq_dir = file.path(work_dir, "filtered_fastq"),
                                     isolates = FALSE,
                                     threads = 22,
                                     profile = "standard",
                                     run_plasmid_id = TRUE,
                                     resume = TRUE,
                                     min_len = 3000,
                                     max_len = 12000,
                                     flye_genome_size = NULL,
                                     flye_asm_coverage = NULL,
                                     sample_sheet = NULL,
                                     workflow = "epi2me-labs/wf-bacterial-genomes",
                                     nextflow = "nextflow",
                                     seqkit = "seqkit",
                                     gzip = "gzip",
                                     quiet = FALSE,
                                     extra_args = "-offline",
                                     syntax_parser = "v1",
                                     ansi_log = FALSE,
                                     nextflow_env = NULL,
                                     dry_run = FALSE,
                                     echo = TRUE,
                                     wait = TRUE,
                                     stdout = "",
                                     stderr = "") {
  check_scalar_character(fastq, "fastq")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(work_dir, "work_dir")
  check_scalar_character(filtered_fastq_dir, "filtered_fastq_dir")
  check_scalar_character(profile, "profile")
  check_scalar_character(workflow, "workflow")
  check_scalar_character(nextflow, "nextflow")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(gzip, "gzip")
  check_logical_scalar(isolates, "isolates")
  check_logical_scalar(run_plasmid_id, "run_plasmid_id")
  check_logical_scalar(resume, "resume")
  check_logical_scalar(quiet, "quiet")
  check_logical_scalar(ansi_log, "ansi_log")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  if (!is.null(extra_args)) check_scalar_character(extra_args, "extra_args")
  if (!is.null(syntax_parser)) check_scalar_character(syntax_parser, "syntax_parser")
  if (!is.null(sample_sheet)) check_file_arg(sample_sheet, "sample_sheet")
  if (!is.null(flye_genome_size)) {
    flye_genome_size <- validate_scalar_cli_value(flye_genome_size, "flye_genome_size")
  }
  if (!is.null(flye_asm_coverage)) {
    flye_asm_coverage <- validate_scalar_cli_value(flye_asm_coverage, "flye_asm_coverage")
  }

  threads <- validate_positive_integer(threads, "threads")
  min_len <- validate_positive_integer(min_len, "min_len")
  max_len <- validate_optional_positive_integer(max_len, "max_len")
  if (!is.null(max_len) && min_len > max_len) {
    stop("`min_len` must be less than or equal to `max_len`.", call. = FALSE)
  }

  nextflow_env <- build_nextflow_env(syntax_parser, ansi_log, nextflow_env)
  filtering_enabled <- !is.null(max_len)
  filter_plan <- list(
    enabled = filtering_enabled,
    command_strings = character(),
    input_files = character(),
    output_files = character(),
    statuses = integer()
  )

  workflow_fastq <- fastq
  if (isTRUE(filtering_enabled)) {
    filter_plan <- build_bacterial_genomes_filter_plan(
      fastq = fastq,
      filtered_fastq_dir = filtered_fastq_dir,
      min_len = min_len,
      max_len = max_len,
      seqkit = seqkit,
      gzip = gzip
    )
    workflow_fastq <- filtered_fastq_dir
  }

  args <- character()
  if (isTRUE(quiet)) {
    args <- c(args, "-q")
  }
  args <- c(
    args,
    "run", workflow,
    "--fastq", workflow_fastq,
    "--out_dir", out_dir,
    "--threads", as.character(threads),
    "--min_read_length", as.character(min_len),
    "--run_plasmid_id", tolower(as.character(run_plasmid_id)),
    "-work-dir", work_dir,
    "-profile", profile
  )
  if (isTRUE(isolates)) {
    args <- c(args, "--isolates")
  }
  if (!is.null(sample_sheet)) {
    args <- c(args, "--sample_sheet", normalizePath(sample_sheet, mustWork = TRUE))
  }
  if (!is.null(flye_genome_size)) {
    args <- c(args, "--flye_genome_size", flye_genome_size)
  }
  if (!is.null(flye_asm_coverage)) {
    args <- c(args, "--flye_asm_coverage", flye_asm_coverage)
  }
  if (isTRUE(resume)) {
    args <- c(args, "-resume")
  }

  command_string <- make_wf_command_string(
    nextflow = nextflow,
    args = args,
    extra_args = extra_args
  )

  paths <- list(
    fastq = fastq,
    workflow_fastq = workflow_fastq,
    filtered_fastq_dir = if (isTRUE(filtering_enabled)) filtered_fastq_dir else NULL,
    out_dir = out_dir,
    work_dir = work_dir,
    sample_sheet = if (!is.null(sample_sheet)) normalizePath(sample_sheet, mustWork = TRUE) else NULL
  )

  uses_shell <- !is.null(extra_args)
  execution_command <- nextflow
  execution_args <- args
  shell_script <- NULL
  if (isTRUE(uses_shell)) {
    execution_command <- "sh"
    execution_args <- "<temporary shell script>"
  }

  if (isTRUE(echo)) {
    if (length(filter_plan$command_strings) > 0L) {
      message(paste(filter_plan$command_strings, collapse = "\n"))
    }
    if (length(nextflow_env) > 0L) {
      message(paste(nextflow_env, collapse = " "))
    }
    message(command_string)
  }

  if (isTRUE(dry_run)) {
    return(invisible(list(
      command = nextflow,
      args = args,
      extra_args = extra_args,
      command_string = command_string,
      execution_command = execution_command,
      execution_args = execution_args,
      uses_shell = uses_shell,
      env = nextflow_env,
      shell_script = shell_script,
      status = NA_integer_,
      filter = filter_plan,
      paths = paths,
      min_len = min_len,
      max_len = max_len
    )))
  }

  if (isTRUE(filtering_enabled)) {
    require_external_command(seqkit)
    require_external_command(gzip)
    filter_plan$statuses <- run_bacterial_genomes_filter_plan(filter_plan)
  }

  if (isTRUE(uses_shell)) {
    shell_script <- tempfile("run_wf_bacterial_genomes_", fileext = ".sh")
    writeLines(c("#!/bin/sh", "set -e", command_string), shell_script)
    execution_args <- shell_script
  }

  status <- system2(
    command = execution_command,
    args = execution_args,
    env = nextflow_env,
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )

  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("wf-bacterial-genomes failed with exit status: ", status, call. = FALSE)
  }

  invisible(list(
    command = nextflow,
    args = args,
    extra_args = extra_args,
    command_string = command_string,
    execution_command = execution_command,
    execution_args = execution_args,
    uses_shell = uses_shell,
    env = nextflow_env,
    shell_script = shell_script,
    status = status,
    filter = filter_plan,
    paths = paths,
    min_len = min_len,
    max_len = max_len
  ))
}

build_bacterial_genomes_filter_plan <- function(fastq,
                                                filtered_fastq_dir,
                                                min_len,
                                                max_len,
                                                seqkit,
                                                gzip) {
  if (!dir.exists(fastq)) {
    stop("`fastq` directory does not exist: ", fastq, call. = FALSE)
  }

  barcode_dirs <- list.dirs(fastq, full.names = TRUE, recursive = FALSE)
  barcode_dirs <- barcode_dirs[grepl("^barcode", basename(barcode_dirs))]
  if (length(barcode_dirs) == 0L) {
    stop("No barcode directories found under `fastq`: ", fastq, call. = FALSE)
  }

  input_files <- file.path(barcode_dirs, paste0(basename(barcode_dirs), ".fastq.gz"))
  missing <- !file.exists(input_files)
  if (any(missing)) {
    stop(
      "Expected per-barcode FASTQ file(s) not found: ",
      paste(input_files[missing], collapse = ", "),
      call. = FALSE
    )
  }

  output_files <- file.path(
    filtered_fastq_dir,
    basename(barcode_dirs),
    paste0(basename(barcode_dirs), ".fastq.gz")
  )
  command_strings <- vapply(
    seq_along(input_files),
    function(i) {
      paste(
        shQuote(seqkit), "seq",
        "-m", shQuote(as.character(min_len)),
        "-M", shQuote(as.character(max_len)),
        shQuote(input_files[i]),
        "|",
        shQuote(gzip), "-c",
        ">",
        shQuote(output_files[i])
      )
    },
    character(1)
  )

  list(
    enabled = TRUE,
    command_strings = command_strings,
    input_files = input_files,
    output_files = output_files,
    statuses = integer()
  )
}

run_bacterial_genomes_filter_plan <- function(filter_plan) {
  statuses <- integer(length(filter_plan$command_strings))
  for (i in seq_along(filter_plan$command_strings)) {
    dir.create(dirname(filter_plan$output_files[i]), recursive = TRUE, showWarnings = FALSE)
    status <- system2(
      "sh",
      args = c("-c", filter_plan$command_strings[i]),
      stdout = "",
      stderr = ""
    )
    if (!identical(status, 0L)) {
      stop("seqkit length filtering failed with exit status: ", status, call. = FALSE)
    }
    statuses[i] <- status
  }

  statuses
}

validate_scalar_cli_value <- function(x, name) {
  if ((!is.character(x) && !is.numeric(x)) || length(x) != 1L ||
      is.na(x) || !nzchar(as.character(x))) {
    stop("`", name, "` must be a single non-empty character or numeric value.",
         call. = FALSE)
  }

  as.character(x)
}
