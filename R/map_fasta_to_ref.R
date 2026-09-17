#' Map FASTA sequences to a reference
#'
#' `map_fasta_to_ref()` aligns query FASTA sequences to a reference FASTA with
#' `minimap2`, reads the PAF alignment table, then writes a coordinate-sorted
#' and indexed BAM with `samtools`.
#'
#' @param fastaFile Query FASTA file.
#' @param ref Reference FASTA file.
#' @param output_file Output BAM file. If `NULL`, a temporary `.bam` file is
#'   created.
#' @param preset Minimap2 preset passed to `minimap2 -x`.
#' @param threads Positive integer thread count passed to `minimap2 -t` and
#'   `samtools sort -@`.
#' @param secondary Logical. If `FALSE`, pass `--secondary=no` to minimap2.
#' @param minimap2,samtools Command names or executable paths.
#' @param conda_env Optional conda environment name. If supplied, external
#'   commands are run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return planned commands without running
#'   external tools.
#' @param strict Logical. If `TRUE`, stop when minimap2 finishes but no PAF
#'   alignment is found. If `FALSE`, return `NULL`, matching the original helper
#'   behavior.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stderr Passed to [system2()].
#'
#' @return Invisibly returns a list with `status`, `commands`, `paths`, `paf`,
#'   `output_file`, `preset`, and `conda_env`.
#'
#' @examples
#' query <- tempfile(fileext = ".fasta")
#' ref <- tempfile(fileext = ".fasta")
#' writeLines(c(">query1", "ACGTACGT"), query)
#' writeLines(c(">chr1", "ACGTACGTACGT"), ref)
#' res <- map_fasta_to_ref(query, ref, dry_run = TRUE)
#' res$commands$minimap2_paf
#'
#' @export
map_fasta_to_ref <- function(fastaFile,
                             ref,
                             output_file = NULL,
                             preset = "asm5",
                             threads = 1,
                             secondary = TRUE,
                             minimap2 = "minimap2",
                             samtools = "samtools",
                             conda_env = NULL,
	                             conda = "conda",
	                             dry_run = FALSE,
	                             strict = FALSE,
	                             echo = TRUE,
	                             stderr = "") {
  check_file_arg(fastaFile, "fastaFile")
  check_file_arg(ref, "ref")
  check_scalar_character(preset, "preset")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  check_logical_scalar(secondary, "secondary")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(strict, "strict")
  check_logical_scalar(echo, "echo")
  if (!is.null(output_file)) check_scalar_character(output_file, "output_file")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")

  threads <- validate_positive_integer(threads, "threads")

  fastaFile <- normalizePath(fastaFile, mustWork = TRUE)
  ref <- normalizePath(ref, mustWork = TRUE)
  if (is.null(output_file)) {
    output_file <- tempfile(fileext = ".bam")
  } else {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    output_file <- normalizePath(output_file, mustWork = FALSE)
  }

  paf_file <- tempfile(fileext = ".paf")
  sam_file <- tempfile(fileext = ".sam")
  bai_file <- paste0(output_file, ".bai")

  minimap2_paf_args <- c("-x", preset, "-t", as.character(threads))
  minimap2_sam_args <- c("-a", "-x", preset, "-t", as.character(threads))
  if (!isTRUE(secondary)) {
    minimap2_paf_args <- c(minimap2_paf_args, "--secondary=no")
    minimap2_sam_args <- c(minimap2_sam_args, "--secondary=no")
  }
  minimap2_paf_args <- c(minimap2_paf_args, ref, fastaFile)
  minimap2_sam_args <- c(minimap2_sam_args, ref, fastaFile)

  sort_args <- c("sort", "-@", as.character(threads), "-o", output_file, sam_file)
  index_args <- c("index", output_file)

  paf_call <- dehost_fastq_external_call(minimap2, minimap2_paf_args, conda_env, conda)
  sam_call <- dehost_fastq_external_call(minimap2, minimap2_sam_args, conda_env, conda)
  sort_call <- dehost_fastq_external_call(samtools, sort_args, conda_env, conda)
  index_call <- dehost_fastq_external_call(samtools, index_args, conda_env, conda)

  commands <- list(
    minimap2_paf = paste(c(shQuote(paf_call$command), shQuote(paf_call$args), ">", shQuote(paf_file)), collapse = " "),
    minimap2_sam = paste(c(shQuote(sam_call$command), shQuote(sam_call$args), ">", shQuote(sam_file)), collapse = " "),
    samtools_sort = paste(c(shQuote(sort_call$command), shQuote(sort_call$args)), collapse = " "),
    samtools_index = paste(c(shQuote(index_call$command), shQuote(index_call$args)), collapse = " ")
  )

  if (isTRUE(echo)) {
    message(commands$minimap2_paf)
    message(commands$minimap2_sam)
    message(commands$samtools_sort)
    message(commands$samtools_index)
  }

  paths <- list(
    query = fastaFile,
    reference = ref,
    paf = paf_file,
    sam = sam_file,
    bam = output_file,
    bai = bai_file
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      commands = commands,
      paths = paths,
      paf = NULL,
      output_file = output_file,
      preset = preset,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(minimap2)
    require_external_command(samtools)
  } else {
    require_external_command(conda)
  }

  on.exit(unlink(c(paf_file, sam_file)), add = TRUE)

  paf_status <- system2(
    paf_call$command,
    args = paf_call$args,
    stdout = paf_file,
    stderr = stderr
  )
  if (!identical(paf_status, 0L)) {
    stop("minimap2 PAF alignment failed with exit status: ", paf_status,
         call. = FALSE)
  }

  paf <- map_fasta_to_ref_read_paf(
    paf_file,
    strict = strict,
    fastaFile = fastaFile,
    ref = ref
  )
  if (is.null(paf)) return(NULL)

  sam_status <- system2(
    sam_call$command,
    args = sam_call$args,
    stdout = sam_file,
    stderr = stderr
  )
  if (!identical(sam_status, 0L)) {
    stop("minimap2 SAM alignment failed with exit status: ", sam_status,
         call. = FALSE)
  }

  sort_status <- system2(sort_call$command, args = sort_call$args, stderr = stderr)
  if (!identical(sort_status, 0L)) {
    stop("samtools sort failed with exit status: ", sort_status, call. = FALSE)
  }

  index_status <- system2(index_call$command, args = index_call$args, stderr = stderr)
  if (!identical(index_status, 0L)) {
    stop("samtools index failed with exit status: ", index_status, call. = FALSE)
  }

  invisible(list(
    status = 0L,
    commands = commands,
    paths = paths,
    paf = paf,
    output_file = output_file,
    preset = preset,
    conda_env = conda_env
  ))
}

map_fasta_to_ref_read_paf <- function(paf_file,
                                      strict = FALSE,
                                      fastaFile = NULL,
                                      ref = NULL) {
  paf <- read_paf(paf_file)
  if (!is.null(paf) && nrow(paf) > 0L) {
    return(paf)
  }

  msg <- "No sequence alignment was found in the reference."
  if (!is.null(fastaFile) || !is.null(ref)) {
    query_label <- if (is.null(fastaFile)) "<unknown query>" else fastaFile
    ref_label <- if (is.null(ref)) "<unknown reference>" else ref
    msg <- paste0(
      "No sequence alignment was found when mapping query FASTA ",
      shQuote(query_label),
      " to reference ",
      shQuote(ref_label),
      "."
    )
  }
  if (isTRUE(strict)) {
    stop(msg, call. = FALSE)
  }
  message(msg)
  NULL
}
