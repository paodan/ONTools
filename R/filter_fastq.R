#' Filter or transform FASTQ records with seqkit seq
#'
#' `filter_fastq()` wraps `seqkit seq` to filter reads by length or average
#' quality and optionally apply common sequence transformations.
#'
#' @param fastq Input FASTQ or FASTQ.GZ file.
#' @param output_fastq Output FASTQ or FASTQ.GZ file.
#' @param min_length,max_length Minimum and maximum read length passed to
#'   `seqkit seq -m/-M`. Use `NULL` for no limit.
#' @param min_quality,max_quality Minimum and maximum average quality passed to
#'   `seqkit seq -Q/-R`. Use `NULL` for no limit.
#' @param reverse,complement Logical. Pass `--reverse` or `--complement`.
#' @param remove_gaps Logical. Pass `--remove-gaps`.
#' @param gap_letters Gap letters passed to `--gap-letters` when non-`NULL`.
#' @param upper_case,lower_case Logical. Pass `--upper-case` or `--lower-case`.
#' @param dna2rna,rna2dna Logical. Pass `--dna2rna` or `--rna2dna`.
#' @param only_id,name,seq,qual Logical. Output only IDs, names, sequences, or
#'   qualities with `--only-id`, `--name`, `--seq`, or `--qual`.
#' @param validate_seq Logical. Pass `--validate-seq`.
#' @param line_width Optional line width passed to `--line-width`.
#' @param seq_type Optional sequence type passed to `--seq-type`.
#' @param threads Optional thread count passed to `--threads`.
#' @param compress_level Optional compression level passed to `--compress-level`.
#' @param extra_args Optional character vector of additional raw `seqkit seq`
#'   arguments appended before the input FASTQ.
#' @param seqkit Command name or executable path for `seqkit`.
#' @param conda_env Optional conda environment name. If supplied, `seqkit` is run
#'   with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return the planned command without
#'   running it.
#' @param echo Logical. If `TRUE`, print the planned command before execution.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `command`, `args`,
#'   `command_string`, and `paths`.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq.gz")
#' out <- tempfile(fileext = ".fastq")
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
#' res <- filter_fastq(reads, out, min_length = 5000, dry_run = TRUE)
#' res$command_string
#'
#' @export
filter_fastq <- function(fastq,
                         output_fastq,
                         min_length = NULL,
                         max_length = NULL,
                         min_quality = NULL,
                         max_quality = NULL,
                         reverse = FALSE,
                         complement = FALSE,
                         remove_gaps = FALSE,
                         gap_letters = NULL,
                         upper_case = FALSE,
                         lower_case = FALSE,
                         dna2rna = FALSE,
                         rna2dna = FALSE,
                         only_id = FALSE,
                         name = FALSE,
                         seq = FALSE,
                         qual = FALSE,
                         validate_seq = FALSE,
                         line_width = NULL,
                         seq_type = NULL,
                         threads = NULL,
                         compress_level = NULL,
                         extra_args = NULL,
                         seqkit = "seqkit",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stderr = "") {
  check_file_arg(fastq, "fastq")
  check_scalar_character(output_fastq, "output_fastq")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(reverse, "reverse")
  check_logical_scalar(complement, "complement")
  check_logical_scalar(remove_gaps, "remove_gaps")
  check_logical_scalar(upper_case, "upper_case")
  check_logical_scalar(lower_case, "lower_case")
  check_logical_scalar(dna2rna, "dna2rna")
  check_logical_scalar(rna2dna, "rna2dna")
  check_logical_scalar(only_id, "only_id")
  check_logical_scalar(name, "name")
  check_logical_scalar(seq, "seq")
  check_logical_scalar(qual, "qual")
  check_logical_scalar(validate_seq, "validate_seq")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }
  if (!is.null(gap_letters)) {
    check_scalar_character(gap_letters, "gap_letters")
  }
  if (!is.null(seq_type)) {
    check_scalar_character(seq_type, "seq_type")
  }
  if (!is.null(extra_args)) {
    if (!is.character(extra_args) || anyNA(extra_args)) {
      stop("`extra_args` must be a character vector without missing values.",
           call. = FALSE)
    }
  }

  min_length <- validate_optional_positive_integer(min_length, "min_length")
  max_length <- validate_optional_positive_integer(max_length, "max_length")
  min_quality <- validate_optional_nonnegative_number(min_quality, "min_quality")
  max_quality <- validate_optional_nonnegative_number(max_quality, "max_quality")
  line_width <- validate_optional_nonnegative_integer(line_width, "line_width")
  threads <- validate_optional_positive_integer(threads, "threads")
  compress_level <- validate_optional_integer(compress_level, "compress_level")

  fastq <- normalizePath(fastq, mustWork = TRUE)
  dir.create(dirname(output_fastq), recursive = TRUE, showWarnings = FALSE)
  output_fastq <- normalizePath(output_fastq, mustWork = FALSE)

  seqkit_args <- build_filter_fastq_args(
    fastq = fastq,
    min_length = min_length,
    max_length = max_length,
    min_quality = min_quality,
    max_quality = max_quality,
    reverse = reverse,
    complement = complement,
    remove_gaps = remove_gaps,
    gap_letters = gap_letters,
    upper_case = upper_case,
    lower_case = lower_case,
    dna2rna = dna2rna,
    rna2dna = rna2dna,
    only_id = only_id,
    name = name,
    seq = seq,
    qual = qual,
    validate_seq = validate_seq,
    line_width = line_width,
    seq_type = seq_type,
    threads = threads,
    compress_level = compress_level,
    extra_args = extra_args
  )
  seqkit_call <- dehost_fastq_external_call(
    command = seqkit,
    args = seqkit_args,
    conda_env = conda_env,
    conda = conda
  )

  command_string <- paste(
    c(shQuote(seqkit_call$command), shQuote(seqkit_call$args), ">", shQuote(output_fastq)),
    collapse = " "
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  paths <- list(
    fastq = fastq,
    output_fastq = output_fastq
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      command = seqkit_call$command,
      args = seqkit_call$args,
      command_string = command_string,
      paths = paths,
      filters = list(
        min_length = min_length,
        max_length = max_length,
        min_quality = min_quality,
        max_quality = max_quality
      ),
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(seqkit)
  } else {
    require_external_command(conda)
  }

  status <- system2(
    seqkit_call$command,
    args = seqkit_call$args,
    stdout = output_fastq,
    stderr = stderr
  )
  if (!identical(status, 0L)) {
    stop("seqkit seq failed with exit status: ", status, call. = FALSE)
  }

  invisible(list(
    status = status,
    command = seqkit_call$command,
    args = seqkit_call$args,
    command_string = command_string,
    paths = paths,
    filters = list(
      min_length = min_length,
      max_length = max_length,
      min_quality = min_quality,
      max_quality = max_quality
    ),
    conda_env = conda_env
  ))
}

#' @rdname filter_fastq
#' @export
filter_fastq_by_length <- function(fastq,
                                   output_fastq,
                                   min_length = 5000,
                                   seqkit = "seqkit",
                                   conda_env = NULL,
                                   conda = "conda",
                                   dry_run = FALSE,
                                   echo = TRUE,
                                   stderr = "") {
  filter_fastq(
    fastq = fastq,
    output_fastq = output_fastq,
    min_length = min_length,
    seqkit = seqkit,
    conda_env = conda_env,
    conda = conda,
    dry_run = dry_run,
    echo = echo,
    stderr = stderr
  )
}

build_filter_fastq_args <- function(fastq,
                                    min_length,
                                    max_length,
                                    min_quality,
                                    max_quality,
                                    reverse,
                                    complement,
                                    remove_gaps,
                                    gap_letters,
                                    upper_case,
                                    lower_case,
                                    dna2rna,
                                    rna2dna,
                                    only_id,
                                    name,
                                    seq,
                                    qual,
                                    validate_seq,
                                    line_width,
                                    seq_type,
                                    threads,
                                    compress_level,
                                    extra_args) {
  args <- "seq"
  if (!is.null(min_length)) args <- c(args, "-m", as.character(min_length))
  if (!is.null(max_length)) args <- c(args, "-M", as.character(max_length))
  if (!is.null(min_quality)) args <- c(args, "-Q", as.character(min_quality))
  if (!is.null(max_quality)) args <- c(args, "-R", as.character(max_quality))
  if (isTRUE(reverse)) args <- c(args, "--reverse")
  if (isTRUE(complement)) args <- c(args, "--complement")
  if (isTRUE(remove_gaps)) args <- c(args, "--remove-gaps")
  if (!is.null(gap_letters)) args <- c(args, "--gap-letters", gap_letters)
  if (isTRUE(upper_case)) args <- c(args, "--upper-case")
  if (isTRUE(lower_case)) args <- c(args, "--lower-case")
  if (isTRUE(dna2rna)) args <- c(args, "--dna2rna")
  if (isTRUE(rna2dna)) args <- c(args, "--rna2dna")
  if (isTRUE(only_id)) args <- c(args, "--only-id")
  if (isTRUE(name)) args <- c(args, "--name")
  if (isTRUE(seq)) args <- c(args, "--seq")
  if (isTRUE(qual)) args <- c(args, "--qual")
  if (isTRUE(validate_seq)) args <- c(args, "--validate-seq")
  if (!is.null(line_width)) args <- c(args, "--line-width", as.character(line_width))
  if (!is.null(seq_type)) args <- c(args, "--seq-type", seq_type)
  if (!is.null(threads)) args <- c(args, "--threads", as.character(threads))
  if (!is.null(compress_level)) {
    args <- c(args, "--compress-level", as.character(compress_level))
  }
  if (!is.null(extra_args)) args <- c(args, extra_args)
  c(args, fastq)
}

validate_optional_positive_integer <- function(x, name) {
  if (is.null(x)) return(NULL)
  validate_positive_integer(x, name)
}

validate_optional_nonnegative_number <- function(x, name) {
  if (is.null(x)) return(NULL)
  validate_nonnegative_number(x, name)
}

validate_optional_nonnegative_integer <- function(x, name) {
  if (is.null(x)) return(NULL)
  validate_nonnegative_integer(x, name)
}

validate_optional_integer <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x)) {
    stop("`", name, "` must be a single integer.", call. = FALSE)
  }

  as.integer(x)
}
