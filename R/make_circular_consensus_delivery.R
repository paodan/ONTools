#' Build a circular consensus delivery package
#'
#' `make_circular_consensus_delivery()` is an R wrapper around the bundled
#' `make_circular_consensus_delivery.sh` shell script. It packages a finished
#' ONT circular DNA consensus FASTA together with optional evidence files,
#' writes README and manifest files, generates MD5 checksums, and creates a
#' `.tar.gz` archive.
#'
#' The wrapped shell script copies the supplied consensus sequence as-is. It
#' does not annotate, polish, rotate, rename FASTA headers, or otherwise modify
#' sequence contents.
#'
#' @param consensus Final consensus FASTA file.
#' @param output Output root directory. The shell script creates
#'   `<output>/<project>_consensus_delivery/` and
#'   `<output>/<project>_consensus_delivery.tar.gz`.
#' @param project Project identifier used in the delivery directory and archive
#'   names.
#' @param sequence_type Sequence type label passed to `--sequence-type`. Must be
#'   one of `"circular_dna"`, `"plasmid"`, or `"virus"`.
#' @param sample_sheet Optional sample sheet or project metadata file to copy
#'   into `00_metadata/`.
#' @param bam Optional final alignment BAM file to copy into `02_evidence/`.
#' @param bai Optional BAM index file to copy into `02_evidence/`.
#' @param depth Optional depth file to copy into `02_evidence/`.
#' @param variants_vcf Optional variant VCF file, usually `variants.vcf.gz`,
#'   to copy into `03_qc/`.
#' @param variants_vcf_index Optional variant VCF index file, usually
#'   `variants.vcf.gz.csi`, to copy into `03_qc/`.
#' @param variants_af Optional variant allele-frequency table to copy into
#'   `03_qc/`.
#' @param variants_af_filtered Optional filtered variant allele-frequency table,
#'   for example `variants_gt0.05.af.tsv`, to copy into `03_qc/`.
#' @param assembly_info Optional Flye `assembly_info.txt` file to copy into
#'   `02_evidence/`.
#' @param flye_log Optional Flye log file to copy into `02_evidence/`.
#' @param notes Optional notes text file to copy into `00_metadata/`.
#' @param overwrite Logical. If `TRUE`, remove an existing
#'   `<output>/<project>_consensus_delivery/` directory or archive before
#'   rebuilding. Defaults to `FALSE` to protect existing delivery packages.
#' @param script Path to the shell script. Defaults to the script bundled with
#'   ONTools.
#' @param bash Bash executable used to run `script`.
#' @param dry_run Logical. If `TRUE`, return the command without running it.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param wait Logical. Passed to [system2()]. Use `FALSE` to launch the command
#'   asynchronously.
#' @param stdout,stderr Passed to [system2()]. Defaults stream output to the R
#'   console.
#'
#' @details
#' External command-line software must be installed separately before running a
#' real delivery build.
#'
#' Required:
#' - `seqkit`
#' - `tar`
#' - `gzip`
#' - an MD5 command: `md5sum` on Linux, or `md5` on macOS
#'
#' Optional:
#' - `samtools`, used to generate `01_consensus/consensus.fasta.fai`
#'
#' The easiest isolated installation is usually via conda/mamba:
#' `mamba create -n ont-delivery -c conda-forge -c bioconda seqkit samtools`
#' followed by `conda activate ont-delivery`.
#'
#' On macOS with Homebrew, install the core tools with
#' `brew install seqkit samtools`. `tar`, `gzip`, and `md5` are normally already
#' available on macOS.
#'
#' @return Invisibly returns a list with `command`, `args`, `command_string`,
#'   `execution_command`, `execution_args`, `status`, and `paths`. In
#'   `dry_run = TRUE`, `status` is `NA_integer_`.
#'
#' @examples
#' consensus <- tempfile(fileext = ".fasta")
#' output_dir <- tempfile("consensus-delivery-")
#' dir.create(output_dir)
#' writeLines(c(">consensus1", "ACGTACGTACGT"), consensus)
#'
#' res <- make_circular_consensus_delivery(
#'   consensus = consensus,
#'   output = output_dir,
#'   project = "PROJECT001",
#'   sequence_type = "plasmid",
#'   dry_run = TRUE
#' )
#' res$command_string
#'
#' @export
make_circular_consensus_delivery <- function(consensus,
                                             output,
                                             project,
                                             sequence_type = c("circular_dna", "plasmid", "virus"),
                                             sample_sheet = NULL,
                                             bam = NULL,
                                             bai = NULL,
                                             depth = NULL,
                                             variants_vcf = NULL,
                                             variants_vcf_index = NULL,
                                             variants_af = NULL,
                                             variants_af_filtered = NULL,
                                             assembly_info = NULL,
                                             flye_log = NULL,
                                             notes = NULL,
                                             overwrite = FALSE,
                                             script = NULL,
                                             bash = "bash",
                                             dry_run = FALSE,
                                             echo = TRUE,
                                             wait = TRUE,
                                             stdout = "",
                                             stderr = "") {
  check_file_arg(consensus, "consensus")
  check_scalar_character(output, "output")
  check_scalar_character(project, "project")
  check_scalar_character(bash, "bash")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  sequence_type <- match.arg(sequence_type)

  optional_files <- list(
    sample_sheet = sample_sheet,
    bam = bam,
    bai = bai,
    depth = depth,
    variants_vcf = variants_vcf,
    variants_vcf_index = variants_vcf_index,
    variants_af = variants_af,
    variants_af_filtered = variants_af_filtered,
    assembly_info = assembly_info,
    flye_log = flye_log,
    notes = notes
  )
  optional_files <- lapply(names(optional_files), function(name) {
    path <- optional_files[[name]]
    if (is.null(path)) {
      return(NULL)
    }
    check_file_arg(path, name)
    normalizePath(path, mustWork = TRUE)
  })
  names(optional_files) <- c(
    "sample_sheet", "bam", "bai", "depth",
    "variants_vcf", "variants_vcf_index", "variants_af", "variants_af_filtered",
    "assembly_info", "flye_log", "notes"
  )

  if (is.null(script)) {
    script <- system.file(
      "scripts",
      "make_circular_consensus_delivery.sh",
      package = "ONTools",
      mustWork = TRUE
    )
  } else {
    check_file_arg(script, "script")
    script <- normalizePath(script, mustWork = TRUE)
  }

  consensus <- normalizePath(consensus, mustWork = TRUE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)

  args <- c(
    "--consensus", consensus,
    "--output", output,
    "--project", project,
    "--sequence-type", sequence_type
  )

  optional_arg_names <- c(
    sample_sheet = "--sample-sheet",
    bam = "--bam",
    bai = "--bai",
    depth = "--depth",
    variants_vcf = "--variants-vcf",
    variants_vcf_index = "--variants-vcf-index",
    variants_af = "--variants-af",
    variants_af_filtered = "--variants-af-filtered",
    assembly_info = "--assembly-info",
    flye_log = "--flye-log",
    notes = "--notes"
  )
  for (name in names(optional_arg_names)) {
    if (!is.null(optional_files[[name]])) {
      args <- c(args, optional_arg_names[[name]], optional_files[[name]])
    }
  }

  if (isTRUE(overwrite)) {
    args <- c(args, "--overwrite")
  }

  execution_args <- c(script, args)
  command_string <- paste(c(shQuote(bash), shQuote(execution_args)), collapse = " ")

  delivery_name <- paste0(project, "_consensus_delivery")
  paths <- list(
    consensus = consensus,
    output_root = output,
    delivery_dir = file.path(output, delivery_name),
    archive = file.path(output, paste0(delivery_name, ".tar.gz")),
    sample_sheet = optional_files$sample_sheet,
    bam = optional_files$bam,
    bai = optional_files$bai,
    depth = optional_files$depth,
    variants_vcf = optional_files$variants_vcf,
    variants_vcf_index = optional_files$variants_vcf_index,
    variants_af = optional_files$variants_af,
    variants_af_filtered = optional_files$variants_af_filtered,
    assembly_info = optional_files$assembly_info,
    flye_log = optional_files$flye_log,
    notes = optional_files$notes
  )

  if (isTRUE(echo)) {
    message(command_string)
  }

  if (isTRUE(dry_run)) {
    return(invisible(list(
      command = script,
      args = args,
      command_string = command_string,
      execution_command = bash,
      execution_args = execution_args,
      status = NA_integer_,
      paths = paths
    )))
  }

  status <- system2(
    command = bash,
    args = execution_args,
    stdout = stdout,
    stderr = stderr,
    wait = wait
  )

  if (isTRUE(wait) && !identical(status, 0L)) {
    stop("make_circular_consensus_delivery failed with exit status: ", status,
         call. = FALSE)
  }

  invisible(list(
    command = script,
    args = args,
    command_string = command_string,
    execution_command = bash,
    execution_args = execution_args,
    status = status,
    paths = paths
  ))
}
