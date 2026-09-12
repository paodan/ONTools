#' Run the wf-16s Nextflow workflow
#'
#' `run_wf_16s()` is an R wrapper around `nextflow run epi2me-labs/wf-16s`
#' for ONT 16S analysis. Common Nextflow options are exposed as R arguments,
#' and additional wf-16s options can be supplied as a raw command-line string
#' through `extra_args`.
#'
#' @param fastq Path to the input FASTQ directory passed to `--fastq`.
#' @param out_dir Output directory passed to `--out_dir`.
#' @param work_dir Nextflow work directory passed with `-work-dir`.
#' @param profile Nextflow profile passed with `-profile`.
#' @param resume Logical. If `TRUE`, append `-resume`.
#' @param database_set wf-16s database set passed to `--database_set`.
#'   Defaults to `"ncbi_16s_18s"`, matching the wf-16s default.
#' @param min_len,max_len Minimum and maximum read length passed to
#'   `--min_len` and `--max_len`. Defaults match wf-16s (`800` and `2000`).
#' @param workflow Nextflow workflow name or path.
#' @param nextflow Nextflow executable name or path.
#' @param quiet Logical. If `TRUE`, pass `-q` to Nextflow to reduce its log
#'   output.
#' @param extra_args Optional raw command-line string appended after the standard
#'   arguments. Use this for additional wf-16s parameters exactly as you would
#'   type them in the shell, the default is `"--minimap2_by_reference"`.
#' @param syntax_parser Nextflow syntax parser version passed as the
#'   `NXF_SYNTAX_PARSER` environment variable. Defaults to `"v1"` because some
#'   EPI2ME workflows still use Groovy `import` declarations that are rejected
#'   by the newer parser. Use `NULL` to leave it unchanged.
#' @param ansi_log Logical. Passed as the `NXF_ANSI_LOG` environment variable.
#'   Defaults to `FALSE` to avoid frequent dynamic Nextflow status updates in
#'   the R console.
#' @param nextflow_env Optional character vector of additional environment
#'   variables passed to [system2()], formatted as `"NAME=value"`.
#' @param dry_run Logical. If `TRUE`, return the command without running it.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param wait Logical. Passed to [system2()]. Use `FALSE` to launch the command
#'   asynchronously.
#' @param stdout,stderr Passed to [system2()]. Defaults stream output to the R
#'   console.
#'
#' @return Invisibly returns a list with `command`, `args`, `extra_args`,
#'   `command_string`, `execution_command`, `execution_args`, `uses_shell`,
#'   `env`, `status`, and `paths`.
#'
#' @examples
#' res <- run_wf_16s(
#'   fastq = "./fastq_pass_trim",
#'   out_dir = "./results/wf_16s",
#'   work_dir = "./work/wf_16s",
#'   extra_args = "--minimap2_by_reference",
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
run_wf_16s <- function(fastq = "./fastq_pass_trim",
                       out_dir = "./results/wf_16s",
                       work_dir = "./work/wf_16s",
                       profile = "standard",
                       resume = TRUE,
                       database_set = "ncbi_16s_18s",
                       min_len = 800,
                       max_len = 2000,
                       workflow = "epi2me-labs/wf-16s",
                       nextflow = "nextflow",
                       quiet = FALSE,
                       extra_args = "--minimap2_by_reference",
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
  check_scalar_character(profile, "profile")
  check_scalar_character(database_set, "database_set")
  check_scalar_character(workflow, "workflow")
  check_scalar_character(nextflow, "nextflow")
  check_logical_scalar(resume, "resume")
  check_logical_scalar(quiet, "quiet")
  check_logical_scalar(ansi_log, "ansi_log")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  if (!is.null(extra_args)) {
    check_scalar_character(extra_args, "extra_args")
  }
  if (!is.null(syntax_parser)) {
    check_scalar_character(syntax_parser, "syntax_parser")
  }
  min_len <- validate_positive_integer(min_len, "min_len")
  max_len <- validate_positive_integer(max_len, "max_len")
  nextflow_env <- build_nextflow_env(syntax_parser, ansi_log, nextflow_env)

  args <- character()
  if (isTRUE(quiet)) {
    args <- c(args, "-q")
  }

  args <- c(
    args,
    "run", workflow,
    "--fastq", fastq,
    "--out_dir", out_dir,
    "--database_set", database_set,
    "--min_len", as.character(min_len),
    "--max_len", as.character(max_len),
    "-work-dir", work_dir,
    "-profile", profile
  )

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
    out_dir = out_dir,
    work_dir = work_dir
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
      paths = paths,
      database_set = database_set,
      min_len = min_len,
      max_len = max_len
    )))
  }

  if (isTRUE(uses_shell)) {
    shell_script <- tempfile("run_wf_16s_", fileext = ".sh")
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
    stop("wf-16s failed with exit status: ", status, call. = FALSE)
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
    paths = paths,
    database_set = database_set,
    min_len = min_len,
    max_len = max_len
  ))
}

#' Run the wf-16s workflow with ITS-friendly defaults
#'
#' `run_ITS()` is a small wrapper around [run_wf_16s()] for ITS amplicon
#' profiling. It uses the same core arguments as `run_wf_16s()`, but defaults to
#' the mixed NCBI marker database and a broader ITS read-length range.
#'
#' @inheritParams run_wf_16s
#'
#' @return Invisibly returns the result from [run_wf_16s()].
#'
#' @examples
#' res <- run_ITS(
#'   fastq = "./fastq_pass_trim",
#'   out_dir = "./results/wf_ITS",
#'   work_dir = "./work/wf_ITS",
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
run_ITS <- function(fastq = "./fastq_pass_trim",
                    out_dir = "./results/wf_ITS",
                    work_dir = "./work/wf_ITS",
                    profile = "standard",
                    resume = TRUE,
                    database_set = "ncbi_16s_18s_28s_ITS",
                    min_len = 300,
                    max_len = 2000,
                    workflow = "epi2me-labs/wf-16s",
                    nextflow = "nextflow",
                    quiet = FALSE,
                    extra_args = "--minimap2_by_reference",
                    syntax_parser = "v1",
                    ansi_log = FALSE,
                    nextflow_env = NULL,
                    dry_run = FALSE,
                    echo = TRUE,
                    wait = TRUE,
                    stdout = "",
                    stderr = "") {
  run_wf_16s(
    fastq = fastq,
    out_dir = out_dir,
    work_dir = work_dir,
    profile = profile,
    resume = resume,
    database_set = database_set,
    min_len = min_len,
    max_len = max_len,
    workflow = workflow,
    nextflow = nextflow,
    quiet = quiet,
    extra_args = extra_args,
    syntax_parser = syntax_parser,
    ansi_log = ansi_log,
    nextflow_env = nextflow_env,
    dry_run = dry_run,
    echo = echo,
    wait = wait,
    stdout = stdout,
    stderr = stderr
  )
}
