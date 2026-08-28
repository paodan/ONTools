#' Remove host/reference-like reads from a FASTQ file
#'
#' `dehost_fastq()` maps reads to a reference with `minimap2`, identifies
#' reference-like reads by PAF mapping quality and aligned fraction, then writes
#' a FASTQ file with those reads removed.
#'
#' @param reference Reference FASTA file.
#' @param fastq Input FASTQ or FASTQ.GZ file.
#' @param output_fastq Output FASTQ or FASTQ.GZ file. If `NULL`, writes
#'   `<out_dir>/<input_basename>.dehost.fastq.gz`.
#' @param out_dir Output directory for intermediate files and, when
#'   `output_fastq = NULL`, the filtered FASTQ file.
#' @param prefix Prefix for intermediate PAF, read id, and stats files.
#' @param threads Positive integer thread count passed to `minimap2`.
#' @param min_mapq Minimum PAF mapping quality. Default: 20.
#' @param min_aln_frac Minimum aligned fraction, calculated as PAF column 11
#'   divided by read length in PAF column 2. Default: 0.8.
#' @param minimap2,seqkit Command names or executable paths.
#' @param conda_env Optional conda environment name. If supplied, external
#'   commands are run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return planned commands without running.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stdout,stderr Passed to [system2()] for commands that stream to the R
#'   console. File-generating commands use explicit output paths.
#'
#' @return Invisibly returns a list with `status`, `commands`, and `paths`.
#'
#' @examples
#' ref <- tempfile(fileext = ".fasta")
#' reads <- tempfile(fileext = ".fastq.gz")
#' writeLines(c(">ecoli", "ACGTACGT"), ref)
#' res <- dehost_fastq(ref, reads, dry_run = TRUE)
#' res$commands$minimap2
#'
#' @export
dehost_fastq <- function(reference,
                         fastq,
                         output_fastq = NULL,
                         out_dir = "work/dehost",
                         prefix = NULL,
                         threads = 16,
                         min_mapq = 20,
                         min_aln_frac = 0.8,
                         minimap2 = "minimap2",
                         seqkit = "seqkit",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stdout = "",
                         stderr = "") {
  check_file_arg(reference, "reference")
  check_file_arg(fastq, "fastq")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }

  threads <- validate_positive_integer(threads, "threads")
  min_mapq <- validate_nonnegative_number(min_mapq, "min_mapq")
  min_aln_frac <- validate_fraction(min_aln_frac, "min_aln_frac")

  if (is.null(prefix)) {
    prefix <- dehost_fastq_prefix(fastq)
  } else {
    check_scalar_character(prefix, "prefix")
  }

  reference <- normalizePath(reference, mustWork = TRUE)
  fastq <- normalizePath(fastq, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)

  if (is.null(output_fastq)) {
    output_fastq <- file.path(out_dir, paste0(prefix, ".dehost.fastq.gz"))
  } else {
    check_scalar_character(output_fastq, "output_fastq")
    output_dir <- dirname(output_fastq)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  output_fastq <- normalizePath(output_fastq, mustWork = FALSE)

  paf <- file.path(out_dir, paste0(prefix, ".vs_reference.paf"))
  read_ids <- file.path(out_dir, paste0(prefix, ".reference_like.read_ids.txt"))
  stats <- file.path(out_dir, paste0(prefix, ".dehost.stats.tsv"))

  minimap2_args <- c(
    "-x", "map-ont",
    "--secondary=no",
    "-t", as.character(threads),
    reference,
    fastq
  )
  seqkit_grep_args <- c(
    "grep",
    "-v",
    "-f", read_ids,
    fastq,
    "-o", output_fastq
  )
  seqkit_stats_args <- c(
    "stats",
    "-a",
    fastq,
    output_fastq
  )
  minimap2_call <- dehost_fastq_external_call(
    command = minimap2,
    args = minimap2_args,
    conda_env = conda_env,
    conda = conda
  )
  seqkit_grep_call <- dehost_fastq_external_call(
    command = seqkit,
    args = seqkit_grep_args,
    conda_env = conda_env,
    conda = conda
  )
  seqkit_stats_call <- dehost_fastq_external_call(
    command = seqkit,
    args = seqkit_stats_args,
    conda_env = conda_env,
    conda = conda
  )

  commands <- list(
    minimap2 = paste(c(shQuote(minimap2_call$command), shQuote(minimap2_call$args), ">", shQuote(paf)), collapse = " "),
    filter_read_ids = paste(
      "Filter PAF in R:",
      shQuote(paf),
      "MAPQ >=",
      min_mapq,
      "and aligned_fraction >=",
      min_aln_frac,
      ">",
      shQuote(read_ids)
    ),
    seqkit_grep = paste(c(shQuote(seqkit_grep_call$command), shQuote(seqkit_grep_call$args)), collapse = " "),
    seqkit_stats = paste(c(shQuote(seqkit_stats_call$command), shQuote(seqkit_stats_call$args), ">", shQuote(stats)), collapse = " ")
  )

  if (isTRUE(echo)) {
    message(commands$minimap2)
    message(commands$filter_read_ids)
    message(commands$seqkit_grep)
    message(commands$seqkit_stats)
  }

  paths <- list(
    reference = reference,
    fastq = fastq,
    paf = paf,
    reference_like_read_ids = read_ids,
    output_fastq = output_fastq,
    stats = stats
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      commands = commands,
      paths = paths,
      thresholds = list(min_mapq = min_mapq, min_aln_frac = min_aln_frac),
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(minimap2)
    require_external_command(seqkit)
  } else {
    require_external_command(conda)
  }

  minimap2_status <- system2(
    minimap2_call$command,
    args = minimap2_call$args,
    stdout = paf,
    stderr = stderr
  )
  if (!identical(minimap2_status, 0L)) {
    stop("minimap2 failed with exit status: ", minimap2_status, call. = FALSE)
  }

  reference_like_ids <- dehost_fastq_read_ids(
    paf = paf,
    min_mapq = min_mapq,
    min_aln_frac = min_aln_frac
  )
  writeLines(reference_like_ids, read_ids)

  grep_status <- system2(
    seqkit_grep_call$command,
    args = seqkit_grep_call$args,
    stdout = stdout,
    stderr = stderr
  )
  if (!identical(grep_status, 0L)) {
    stop("seqkit grep failed with exit status: ", grep_status, call. = FALSE)
  }

  stats_status <- system2(
    seqkit_stats_call$command,
    args = seqkit_stats_call$args,
    stdout = stats,
    stderr = stderr
  )
  if (!identical(stats_status, 0L)) {
    stop("seqkit stats failed with exit status: ", stats_status, call. = FALSE)
  }

  invisible(list(
    status = 0L,
    commands = commands,
    paths = paths,
    thresholds = list(min_mapq = min_mapq, min_aln_frac = min_aln_frac),
    conda_env = conda_env,
    n_reference_like_reads = length(reference_like_ids)
  ))
}

dehost_fastq_external_call <- function(command, args, conda_env, conda) {
  if (is.null(conda_env)) {
    return(list(command = command, args = args))
  }

  list(
    command = conda,
    args = c("run", "-n", conda_env, command, args)
  )
}

dehost_fastq_prefix <- function(fastq) {
  prefix <- basename(fastq)
  prefix <- sub("[.]gz$", "", prefix)
  prefix <- sub("[.]fastq$", "", prefix)
  prefix <- sub("[.]fq$", "", prefix)
  prefix
}

dehost_fastq_read_ids <- function(paf, min_mapq, min_aln_frac) {
  if (!file.exists(paf) || file.info(paf)$size == 0) {
    return(character())
  }

  paf_table <- utils::read.table(
    paf,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )

  if (ncol(paf_table) < 12L || nrow(paf_table) == 0L) {
    return(character())
  }

  read_length <- suppressWarnings(as.numeric(paf_table[[2L]]))
  aligned_bases <- suppressWarnings(as.numeric(paf_table[[11L]]))
  mapq <- suppressWarnings(as.numeric(paf_table[[12L]]))

  keep <- !is.na(read_length) &
    read_length > 0 &
    !is.na(aligned_bases) &
    !is.na(mapq) &
    mapq >= min_mapq &
    (aligned_bases / read_length) >= min_aln_frac

  sort(unique(paf_table[[1L]][keep]))
}

validate_nonnegative_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop("`", name, "` must be a single non-negative number.", call. = FALSE)
  }

  x
}

validate_fraction <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x > 1) {
    stop("`", name, "` must be a single number between 0 and 1.", call. = FALSE)
  }

  x
}

require_external_command <- function(command) {
  if (!nzchar(Sys.which(command))) {
    stop("Required command not found: ", command, call. = FALSE)
  }
}
