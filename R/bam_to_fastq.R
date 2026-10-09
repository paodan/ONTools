#' Convert a BAM file to FASTQ with samtools
#'
#' bam_to_fastq() wraps samtools fastq to recover read sequences and base
#' qualities from a BAM file without loading the reads into R memory.
#'
#' @param bam Path to an input BAM file.
#' @param output_fastq Output FASTQ path. A filename ending in .gz, .bgz, or
#'   .bgzf is compressed automatically by samtools.
#' @param threads Positive integer passed to samtools fastq -@ as the number
#'   of additional compression threads.
#' @param exclude_flags SAM flag mask passed to samtools fastq -F. The default
#'   "0x900" excludes secondary (0x100) and supplementary (0x800) alignments,
#'   preventing extra alignment records from producing duplicate FASTQ entries.
#' @param read_category Read categories to write. "unpaired" uses -0 and is
#'   appropriate for typical Oxford Nanopore BAM files, whose READ1 and READ2
#'   bits are both unset. "all" routes category 0, READ1, and READ2 records to
#'   the same output file with -0, -1, and -2.
#' @param keep_read_names Logical. If TRUE, pass -n so samtools does not add
#'   /1 or /2 suffixes to read names.
#' @param copy_tags Logical. If TRUE, pass -t to copy RG, BC, and QT tags, when
#'   present, into FASTQ header lines.
#' @param use_original_quality Logical. If TRUE, pass -O so quality values in
#'   OQ tags are preferred over the BAM QUAL field when available.
#' @param overwrite Logical. If FALSE, stop when output_fastq already exists.
#'   Conversion is written to a temporary file first, so an existing output is
#'   not replaced unless samtools completes successfully.
#' @param samtools Command name or executable path for samtools.
#' @param conda_env Optional conda environment name. When supplied, samtools is
#'   run through conda run -n <conda_env>.
#' @param conda Conda executable name or path used with conda_env.
#' @param dry_run Logical. If TRUE, return the planned command without running
#'   it or creating the output file.
#' @param echo Logical. If TRUE, print the planned command.
#' @param stderr Passed to [system2()]. The default "" streams standard error
#'   to the R console.
#'
#' @return Invisibly returns a list containing status, command, args,
#'   command_string, paths, read_category, exclude_flags, and conda_env.
#'   status is NA_integer_ for a dry run and zero after successful conversion.
#'
#' @details
#' BAM-to-FASTQ conversion cannot restore information absent from the BAM.
#' Hard-clipped sequence, reads removed before BAM creation, and original FASTQ
#' comments not retained as BAM tags cannot be recovered. Soft-clipped sequence
#' remains available. Alignment-specific tags such as CIGAR, NM, MD, AS, and SA
#' are not represented in ordinary FASTQ output.
#'
#' If the original FASTQ is available, prefer it over reconstructed FASTQ.
#' read_category = "unpaired" is the expected mode for Nanopore single-end data.
#' Use read_category = "all" when the BAM can also contain READ1 or READ2
#' records that must be recovered into the same file.
#'
#' @examples
#' \dontrun{
#' result <- bam_to_fastq(
#'   bam = "sample.sorted.bam",
#'   output_fastq = "sample.fastq.gz",
#'   threads = 4
#' )
#' result$paths$output_fastq
#'
#' planned <- bam_to_fastq(
#'   bam = "sample.sorted.bam",
#'   output_fastq = "sample.fastq.gz",
#'   conda_env = "ont-tools",
#'   dry_run = TRUE
#' )
#' planned$command_string
#' }
#'
#' @export
bam_to_fastq <- function(bam,
                         output_fastq,
                         threads = 4,
                         exclude_flags = "0x900",
                         read_category = c("unpaired", "all"),
                         keep_read_names = TRUE,
                         copy_tags = FALSE,
                         use_original_quality = FALSE,
                         overwrite = FALSE,
                         samtools = "samtools",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stderr = "") {
  check_file_arg(bam, "bam")
  check_scalar_character(output_fastq, "output_fastq")
  check_scalar_character(exclude_flags, "exclude_flags")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  check_logical_scalar(keep_read_names, "keep_read_names")
  check_logical_scalar(copy_tags, "copy_tags")
  check_logical_scalar(use_original_quality, "use_original_quality")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  threads <- validate_positive_integer(threads, "threads")
  read_category <- match.arg(read_category)

  bam <- normalizePath(bam, mustWork = TRUE)
  dir.create(dirname(output_fastq), recursive = TRUE, showWarnings = FALSE)
  output_fastq <- normalizePath(output_fastq, mustWork = FALSE)
  if (file.exists(output_fastq) && !isTRUE(overwrite)) {
    stop(
      "output_fastq already exists and overwrite is FALSE: ",
      output_fastq,
      call. = FALSE
    )
  }

  args <- bam_to_fastq_args(
    bam, output_fastq, threads, exclude_flags, read_category,
    keep_read_names, copy_tags, use_original_quality
  )
  call <- dehost_fastq_external_call(samtools, args, conda_env, conda)
  command_string <- paste(
    c(shQuote(call$command), shQuote(call$args)),
    collapse = " "
  )
  if (isTRUE(echo)) message(command_string)

  paths <- list(bam = bam, output_fastq = output_fastq)
  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      command = call$command,
      args = call$args,
      command_string = command_string,
      paths = paths,
      read_category = read_category,
      exclude_flags = exclude_flags,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(samtools)
  } else {
    require_external_command(conda)
  }

  output_suffix <- if (
    grepl("[.](gz|bgz|bgzf)$", output_fastq, ignore.case = TRUE)
  ) ".fastq.gz" else ".fastq"
  temporary_fastq <- tempfile(
    pattern = ".bam-to-fastq-",
    tmpdir = dirname(output_fastq),
    fileext = output_suffix
  )
  on.exit(unlink(temporary_fastq), add = TRUE)

  execution_args <- bam_to_fastq_args(
    bam, temporary_fastq, threads, exclude_flags, read_category,
    keep_read_names, copy_tags, use_original_quality
  )
  execution_call <- dehost_fastq_external_call(
    samtools, execution_args, conda_env, conda
  )
  status <- suppressWarnings(system2(
    execution_call$command,
    args = execution_call$args,
    stdout = FALSE,
    stderr = stderr
  ))
  if (!identical(status, 0L)) {
    stop("samtools fastq failed with exit status: ", status, call. = FALSE)
  }
  if (!file.exists(temporary_fastq) || file.info(temporary_fastq)$size == 0) {
    stop("samtools fastq did not create a non-empty FASTQ file.", call. = FALSE)
  }

  copied <- file.copy(
    temporary_fastq,
    output_fastq,
    overwrite = overwrite,
    copy.mode = TRUE,
    copy.date = FALSE
  )
  if (!isTRUE(copied)) {
    stop("Could not move converted FASTQ to: ", output_fastq, call. = FALSE)
  }

  invisible(list(
    status = status,
    command = call$command,
    args = call$args,
    command_string = command_string,
    paths = paths,
    read_category = read_category,
    exclude_flags = exclude_flags,
    conda_env = conda_env
  ))
}

bam_to_fastq_args <- function(bam,
                              output_fastq,
                              threads,
                              exclude_flags,
                              read_category,
                              keep_read_names,
                              copy_tags,
                              use_original_quality) {
  args <- c("fastq", "-@", as.character(threads))
  if (isTRUE(keep_read_names)) args <- c(args, "-n")
  if (isTRUE(copy_tags)) args <- c(args, "-t")
  if (isTRUE(use_original_quality)) args <- c(args, "-O")
  args <- c(args, "-F", exclude_flags)
  if (identical(read_category, "unpaired")) {
    args <- c(args, "-0", output_fastq)
  } else {
    args <- c(
      args,
      "-0", output_fastq,
      "-1", output_fastq,
      "-2", output_fastq
    )
  }
  c(args, bam)
}
