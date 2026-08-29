#' Align assembled contigs to a reference and create BAM index
#'
#' `align_contigs_to_reference()` maps assembled contigs to a reference genome
#' or plasmid with `minimap2`, sorts the alignment with `samtools sort`, and
#' creates a BAM index with `samtools index`.
#'
#' @param reference Reference FASTA file.
#' @param assembly Assembled contig FASTA file, for example Flye
#'   `assembly.fasta` or `draft_assembly.fasta`.
#' @param align_bam Output coordinate-sorted BAM file.
#' @param preset Minimap2 preset passed to `minimap2 -x`. Use `"asm5"` for very
#'   similar assemblies, `"asm10"` for moderately diverged assemblies, or
#'   `"asm20"` for more diverged assemblies.
#' @param threads Positive integer thread count passed to both `minimap2 -t`
#'   and `samtools sort -@`.
#' @param secondary Logical. If `FALSE`, pass `--secondary=no` to minimap2.
#' @param minimap2,samtools Command names or executable paths.
#' @param conda_env Optional conda environment name. If supplied, external
#'   commands are run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return planned commands without running.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `commands`, `paths`,
#'   `preset`, and `conda_env`.
#'
#' @examples
#' ref <- tempfile(fileext = ".fasta")
#' assembly <- tempfile(fileext = ".fasta")
#' bam <- tempfile(fileext = ".bam")
#' writeLines(c(">ref", "ACGTACGT"), ref)
#' writeLines(c(">contig1", "ACGT"), assembly)
#' res <- align_contigs_to_reference(ref, assembly, bam, dry_run = TRUE)
#' res$commands$minimap2
#'
#' @export
align_contigs_to_reference <- function(reference,
                                       assembly,
                                       align_bam,
                                       preset = "asm5",
                                       threads = 16,
                                       secondary = FALSE,
                                       minimap2 = "minimap2",
                                       samtools = "samtools",
                                       conda_env = NULL,
                                       conda = "conda",
                                       dry_run = FALSE,
                                       echo = TRUE,
                                       stderr = "") {
  check_file_arg(reference, "reference")
  check_file_arg(assembly, "assembly")
  check_scalar_character(align_bam, "align_bam")
  check_scalar_character(preset, "preset")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  check_logical_scalar(secondary, "secondary")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }

  threads <- validate_positive_integer(threads, "threads")

  reference <- normalizePath(reference, mustWork = TRUE)
  assembly <- normalizePath(assembly, mustWork = TRUE)
  dir.create(dirname(align_bam), recursive = TRUE, showWarnings = FALSE)
  align_bam <- normalizePath(align_bam, mustWork = FALSE)
  align_bai <- paste0(align_bam, ".bai")
  sam <- tempfile(fileext = ".sam")

  minimap2_args <- c(
    "-x", preset,
    "-a",
    "-t", as.character(threads)
  )
  if (!isTRUE(secondary)) {
    minimap2_args <- c(minimap2_args, "--secondary=no")
  }
  minimap2_args <- c(minimap2_args, reference, assembly)

  sort_args <- c(
    "sort",
    "-@", as.character(threads),
    "-o", align_bam,
    sam
  )
  index_args <- c("index", align_bam)

  minimap2_call <- dehost_fastq_external_call(
    command = minimap2,
    args = minimap2_args,
    conda_env = conda_env,
    conda = conda
  )
  sort_call <- dehost_fastq_external_call(
    command = samtools,
    args = sort_args,
    conda_env = conda_env,
    conda = conda
  )
  index_call <- dehost_fastq_external_call(
    command = samtools,
    args = index_args,
    conda_env = conda_env,
    conda = conda
  )

  commands <- list(
    minimap2 = paste(c(shQuote(minimap2_call$command), shQuote(minimap2_call$args), ">", shQuote(sam)), collapse = " "),
    samtools_sort = paste(c(shQuote(sort_call$command), shQuote(sort_call$args)), collapse = " "),
    samtools_index = paste(c(shQuote(index_call$command), shQuote(index_call$args)), collapse = " ")
  )

  if (isTRUE(echo)) {
    message(commands$minimap2)
    message(commands$samtools_sort)
    message(commands$samtools_index)
  }

  paths <- list(
    reference = reference,
    assembly = assembly,
    sam = sam,
    bam = align_bam,
    bai = align_bai
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      commands = commands,
      paths = paths,
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

  on.exit(unlink(sam), add = TRUE)

  minimap2_status <- system2(
    minimap2_call$command,
    args = minimap2_call$args,
    stdout = sam,
    stderr = stderr
  )
  if (!identical(minimap2_status, 0L)) {
    stop("minimap2 failed with exit status: ", minimap2_status, call. = FALSE)
  }

  sort_status <- system2(
    sort_call$command,
    args = sort_call$args,
    stderr = stderr
  )
  if (!identical(sort_status, 0L)) {
    stop("samtools sort failed with exit status: ", sort_status, call. = FALSE)
  }

  index_status <- system2(
    index_call$command,
    args = index_call$args,
    stderr = stderr
  )
  if (!identical(index_status, 0L)) {
    stop("samtools index failed with exit status: ", index_status, call. = FALSE)
  }

  invisible(list(
    status = 0L,
    commands = commands,
    paths = paths,
    preset = preset,
    conda_env = conda_env
  ))
}
