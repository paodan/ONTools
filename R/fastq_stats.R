#' Summarize FASTQ files with optional length-filtered views
#'
#' `fastq_stats()` wraps `seqkit stats -a`. It can summarize one or more FASTQ
#' files directly and optionally summarize temporary `seqkit seq -m` filtered
#' views without writing filtered FASTQ files.
#'
#' @param fastq Character vector of input FASTQ or FASTQ.GZ files.
#' @param min_lengths Optional positive integer vector. For each value, run
#'   `seqkit seq -m <min_length> <fastq> | seqkit stats -a`.
#' @param include_original Logical. If `TRUE`, also run `seqkit stats -a` on the
#'   unfiltered `fastq` files.
#' @param output_tsv Optional output TSV file for the combined stats table.
#' @param seqkit Command name or executable path for `seqkit`.
#' @param conda_env Optional conda environment name. If supplied, `seqkit`
#'   commands are run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return planned commands without running.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `stats`, `commands`, and
#'   `paths`.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq.gz")
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
#' res <- fastq_stats(reads, min_lengths = c(1000, 5000), dry_run = TRUE)
#' res$commands
#'
#' @export
fastq_stats <- function(fastq,
                        min_lengths = NULL,
                        include_original = TRUE,
                        output_tsv = NULL,
                        seqkit = "seqkit",
                        conda_env = NULL,
                        conda = "conda",
                        dry_run = FALSE,
                        echo = TRUE,
                        stderr = "") {
  if (!is.character(fastq) || length(fastq) == 0L || anyNA(fastq) ||
      any(!nzchar(fastq))) {
    stop("`fastq` must be a non-empty character vector of file paths.",
         call. = FALSE)
  }
  missing_fastq <- fastq[!file.exists(fastq)]
  if (length(missing_fastq) > 0L) {
    stop("`fastq` file does not exist: ", missing_fastq[[1L]], call. = FALSE)
  }

  check_logical_scalar(include_original, "include_original")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }
  if (!is.null(output_tsv)) {
    check_scalar_character(output_tsv, "output_tsv")
    dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
    output_tsv <- normalizePath(output_tsv, mustWork = FALSE)
  }

  if (!is.null(min_lengths)) {
    if (!is.numeric(min_lengths) || anyNA(min_lengths) ||
        any(!is.finite(min_lengths)) || any(min_lengths <= 0) ||
        any(min_lengths != as.integer(min_lengths))) {
      stop("`min_lengths` must be a positive integer vector.", call. = FALSE)
    }
    min_lengths <- as.integer(min_lengths)
  }

  if (!isTRUE(include_original) && length(min_lengths) == 0L) {
    stop("At least one of `include_original` or `min_lengths` is required.",
         call. = FALSE)
  }

  fastq <- normalizePath(fastq, mustWork = TRUE)

  direct_call <- dehost_fastq_external_call(
    command = seqkit,
    args = c("stats", "-a", fastq),
    conda_env = conda_env,
    conda = conda
  )
  seq_call_template <- function(path, min_length) {
    dehost_fastq_external_call(
      command = seqkit,
      args = c("seq", "-m", as.character(min_length), path),
      conda_env = conda_env,
      conda = conda
    )
  }
  stats_stdin_call <- dehost_fastq_external_call(
    command = seqkit,
    args = c("stats", "-a"),
    conda_env = conda_env,
    conda = conda
  )

  commands <- character()
  command_meta <- data.frame(
    stat_type = character(),
    fastq = character(),
    min_length = integer(),
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_original)) {
    commands <- c(commands, paste(c(shQuote(direct_call$command), shQuote(direct_call$args)), collapse = " "))
    command_meta <- rbind(
      command_meta,
      data.frame(
        stat_type = "original",
        fastq = NA_character_,
        min_length = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(min_lengths) > 0L) {
    for (path in fastq) {
      for (min_length in min_lengths) {
        seq_call <- seq_call_template(path, min_length)
        commands <- c(
          commands,
          paste(
            paste(c(shQuote(seq_call$command), shQuote(seq_call$args)), collapse = " "),
            "|",
            paste(c(shQuote(stats_stdin_call$command), shQuote(stats_stdin_call$args)), collapse = " ")
          )
        )
        command_meta <- rbind(
          command_meta,
          data.frame(
            stat_type = "min_length",
            fastq = path,
            min_length = min_length,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (isTRUE(echo)) {
    message(paste(commands, collapse = "\n"))
  }

  paths <- list(
    fastq = fastq,
    output_tsv = output_tsv
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      stats = NULL,
      commands = commands,
      command_meta = command_meta,
      paths = paths,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(seqkit)
  } else {
    require_external_command(conda)
  }

  stats_tables <- list()
  command_i <- 0L

  if (isTRUE(include_original)) {
    command_i <- command_i + 1L
    stats_tables[[length(stats_tables) + 1L]] <- run_fastq_stats_command(
      command = direct_call$command,
      args = direct_call$args,
      stderr = stderr,
      stat_type = "original",
      source_fastq = NA_character_,
      min_length = NA_integer_
    )
  }

  if (length(min_lengths) > 0L) {
    for (path in fastq) {
      for (min_length in min_lengths) {
        command_i <- command_i + 1L
        stats_tables[[length(stats_tables) + 1L]] <- run_fastq_stats_shell_command(
          command_string = commands[[command_i]],
          stderr = stderr,
          stat_type = "min_length",
          source_fastq = path,
          min_length = min_length
        )
      }
    }
  }

  stats <- do.call(rbind, stats_tables)
  rownames(stats) <- NULL

  if (!is.null(output_tsv)) {
    utils::write.table(
      stats,
      file = output_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  invisible(list(
    status = 0L,
    stats = stats,
    commands = commands,
    command_meta = command_meta,
    paths = paths,
    conda_env = conda_env
  ))
}

run_fastq_stats_command <- function(command,
                                    args,
                                    stderr,
                                    stat_type,
                                    source_fastq,
                                    min_length) {
  stdout <- system2(command, args = args, stdout = TRUE, stderr = stderr)
  status <- attr(stdout, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("seqkit stats failed with exit status: ", status, call. = FALSE)
  }

  parse_fastq_stats_output(stdout, stat_type, source_fastq, min_length)
}

run_fastq_stats_shell_command <- function(command_string,
                                          stderr,
                                          stat_type,
                                          source_fastq,
                                          min_length) {
  stdout <- system2("sh", args = c("-c", command_string), stdout = TRUE, stderr = stderr)
  status <- attr(stdout, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("seqkit stats pipeline failed with exit status: ", status, call. = FALSE)
  }

  parse_fastq_stats_output(stdout, stat_type, source_fastq, min_length)
}

parse_fastq_stats_output <- function(stdout, stat_type, source_fastq, min_length) {
  if (length(stdout) == 0L) {
    table <- data.frame(stringsAsFactors = FALSE)
  } else {
    table <- utils::read.delim(
      text = paste(stdout, collapse = "\n"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  cbind(
    stat_type = stat_type,
    min_length = min_length,
    source_fastq = source_fastq,
    table,
    stringsAsFactors = FALSE
  )
}
