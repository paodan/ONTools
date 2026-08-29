#' Filter long reads with Filtlong
#'
#' `run_filtlong()` wraps `filtlong` to keep the best long reads by length,
#' quality, and optional target yield.
#'
#' @param reads Input FASTQ or FASTQ.GZ file.
#' @param output_fastq Output FASTQ file. Use a `.gz` suffix only if your
#'   downstream workflow or shell environment compresses it separately; Filtlong
#'   writes FASTQ text to standard output.
#' @param target_bases Optional target number of bases to keep, passed to
#'   `filtlong --target_bases`.
#' @param keep_percent Optional percentage of the best reads to keep, passed to
#'   `filtlong --keep_percent`.
#' @param min_length Optional minimum read length, passed to
#'   `filtlong --min_length`.
#' @param min_mean_q Optional minimum mean read quality, passed to
#'   `filtlong --min_mean_q`.
#' @param min_window_q Optional minimum window quality, passed to
#'   `filtlong --min_window_q`.
#' @param length_weight Optional weight for read length in read scoring, passed
#'   to `filtlong --length_weight`.
#' @param mean_q_weight Optional weight for mean quality in read scoring, passed
#'   to `filtlong --mean_q_weight`.
#' @param window_q_weight Optional weight for window quality in read scoring,
#'   passed to `filtlong --window_q_weight`.
#' @param window_size Optional window size used for window-quality scoring,
#'   passed to `filtlong --window_size`.
#' @param trim Logical. If `TRUE`, pass `filtlong --trim` to trim low-quality
#'   read ends when supported by the installed Filtlong version.
#' @param split Optional integer passed to `filtlong --split` to split reads at
#'   poor-quality positions when supported by the installed Filtlong version.
#' @param illumina_1,illumina_2 Optional paired short-read files passed to
#'   `filtlong --illumina_1` and `filtlong --illumina_2` for external quality
#'   scoring.
#' @param assembly Optional reference assembly passed to `filtlong --assembly`
#'   for external quality scoring.
#' @param verbose Logical. If `TRUE`, pass `filtlong --verbose`.
#' @param extra_args Optional character vector of additional raw Filtlong
#'   arguments appended before the input reads.
#' @param filtlong Command name or executable path for `filtlong`.
#' @param conda_env Optional conda environment name. If supplied, Filtlong is run
#'   with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return the planned command without
#'   running it.
#' @param echo Logical. If `TRUE`, print the planned command before execution.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `command`, `args`,
#'   `command_string`, `paths`, and `conda_env`.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq")
#' out <- tempfile(fileext = ".fastq")
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
#' res <- run_filtlong(
#'   reads = reads,
#'   output_fastq = out,
#'   target_bases = 80000000,
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
run_filtlong <- function(reads,
                         output_fastq,
                         target_bases = NULL,
                         keep_percent = NULL,
                         min_length = NULL,
                         min_mean_q = NULL,
                         min_window_q = NULL,
                         length_weight = NULL,
                         mean_q_weight = NULL,
                         window_q_weight = NULL,
                         window_size = NULL,
                         trim = FALSE,
                         split = NULL,
                         illumina_1 = NULL,
                         illumina_2 = NULL,
                         assembly = NULL,
                         verbose = FALSE,
                         extra_args = NULL,
                         filtlong = "filtlong",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stderr = "") {
  check_file_arg(reads, "reads")
  check_scalar_character(output_fastq, "output_fastq")
  check_scalar_character(filtlong, "filtlong")
  check_scalar_character(conda, "conda")
  check_logical_scalar(trim, "trim")
  check_logical_scalar(verbose, "verbose")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  target_bases <- validate_optional_positive_integer(target_bases, "target_bases")
  keep_percent <- validate_optional_percent(keep_percent, "keep_percent")
  min_length <- validate_optional_positive_integer(min_length, "min_length")
  min_mean_q <- validate_optional_nonnegative_number(min_mean_q, "min_mean_q")
  min_window_q <- validate_optional_nonnegative_number(min_window_q, "min_window_q")
  length_weight <- validate_optional_nonnegative_number(length_weight, "length_weight")
  mean_q_weight <- validate_optional_nonnegative_number(mean_q_weight, "mean_q_weight")
  window_q_weight <- validate_optional_nonnegative_number(window_q_weight, "window_q_weight")
  window_size <- validate_optional_positive_integer(window_size, "window_size")
  split <- validate_optional_positive_integer(split, "split")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }
  if (!is.null(illumina_1)) {
    check_file_arg(illumina_1, "illumina_1")
  }
  if (!is.null(illumina_2)) {
    check_file_arg(illumina_2, "illumina_2")
  }
  if (!is.null(assembly)) {
    check_file_arg(assembly, "assembly")
  }
  if (!is.null(extra_args) && (!is.character(extra_args) || anyNA(extra_args))) {
    stop("`extra_args` must be a character vector without missing values.",
         call. = FALSE)
  }

  reads <- normalizePath(reads, mustWork = TRUE)
  if (!is.null(illumina_1)) illumina_1 <- normalizePath(illumina_1, mustWork = TRUE)
  if (!is.null(illumina_2)) illumina_2 <- normalizePath(illumina_2, mustWork = TRUE)
  if (!is.null(assembly)) assembly <- normalizePath(assembly, mustWork = TRUE)
  dir.create(dirname(output_fastq), recursive = TRUE, showWarnings = FALSE)
  output_fastq <- normalizePath(output_fastq, mustWork = FALSE)

  filtlong_args <- build_filtlong_args(
    reads = reads,
    target_bases = target_bases,
    keep_percent = keep_percent,
    min_length = min_length,
    min_mean_q = min_mean_q,
    min_window_q = min_window_q,
    length_weight = length_weight,
    mean_q_weight = mean_q_weight,
    window_q_weight = window_q_weight,
    window_size = window_size,
    trim = trim,
    split = split,
    illumina_1 = illumina_1,
    illumina_2 = illumina_2,
    assembly = assembly,
    verbose = verbose,
    extra_args = extra_args
  )
  filtlong_call <- dehost_fastq_external_call(
    command = filtlong,
    args = filtlong_args,
    conda_env = conda_env,
    conda = conda
  )

  command_string <- paste(
    c(shQuote(filtlong_call$command), shQuote(filtlong_call$args), ">", shQuote(output_fastq)),
    collapse = " "
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  paths <- list(
    reads = reads,
    output_fastq = output_fastq,
    illumina_1 = illumina_1,
    illumina_2 = illumina_2,
    assembly = assembly
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      command = filtlong_call$command,
      args = filtlong_call$args,
      command_string = command_string,
      paths = paths,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(filtlong)
  } else {
    require_external_command(conda)
  }

  status <- system2(
    filtlong_call$command,
    args = filtlong_call$args,
    stdout = output_fastq,
    stderr = stderr
  )
  if (!identical(status, 0L)) {
    stop("filtlong failed with exit status: ", status, call. = FALSE)
  }

  invisible(list(
    status = status,
    command = filtlong_call$command,
    args = filtlong_call$args,
    command_string = command_string,
    paths = paths,
    conda_env = conda_env
  ))
}

build_filtlong_args <- function(reads,
                                target_bases,
                                keep_percent,
                                min_length,
                                min_mean_q,
                                min_window_q,
                                length_weight,
                                mean_q_weight,
                                window_q_weight,
                                window_size,
                                trim,
                                split,
                                illumina_1,
                                illumina_2,
                                assembly,
                                verbose,
                                extra_args) {
  args <- character()
  if (!is.null(target_bases)) args <- c(args, "--target_bases", as.character(target_bases))
  if (!is.null(keep_percent)) args <- c(args, "--keep_percent", as.character(keep_percent))
  if (!is.null(min_length)) args <- c(args, "--min_length", as.character(min_length))
  if (!is.null(min_mean_q)) args <- c(args, "--min_mean_q", as.character(min_mean_q))
  if (!is.null(min_window_q)) args <- c(args, "--min_window_q", as.character(min_window_q))
  if (!is.null(length_weight)) args <- c(args, "--length_weight", as.character(length_weight))
  if (!is.null(mean_q_weight)) args <- c(args, "--mean_q_weight", as.character(mean_q_weight))
  if (!is.null(window_q_weight)) args <- c(args, "--window_q_weight", as.character(window_q_weight))
  if (!is.null(window_size)) args <- c(args, "--window_size", as.character(window_size))
  if (isTRUE(trim)) args <- c(args, "--trim")
  if (!is.null(split)) args <- c(args, "--split", as.character(split))
  if (!is.null(illumina_1)) args <- c(args, "--illumina_1", illumina_1)
  if (!is.null(illumina_2)) args <- c(args, "--illumina_2", illumina_2)
  if (!is.null(assembly)) args <- c(args, "--assembly", assembly)
  if (isTRUE(verbose)) args <- c(args, "--verbose")
  if (!is.null(extra_args)) args <- c(args, extra_args)

  c(args, reads)
}

validate_optional_percent <- function(x, name) {
  if (is.null(x)) return(NULL)
  value <- validate_optional_nonnegative_number(x, name)
  if (value > 100) {
    stop("`", name, "` must be between 0 and 100.", call. = FALSE)
  }

  value
}
