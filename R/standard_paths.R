#' Build a standard ONTools path layout
#'
#' `standard_paths()` centralizes the directory layout shared by ONTools
#' workflows. It returns paths only: it does not create directories and does
#' not require any returned path to exist. Additional workflow-specific paths
#' can be supplied through `...`.
#'
#' @param input_path_root A single non-empty character string giving the root
#'   directory of the input or analysis results. The default work directory is
#'   created conceptually below this root as `work`; no directory is actually
#'   created by this function.
#' @param delivery_path_root A single non-empty character string giving the
#'   root directory under which the delivery folder will be located.
#' @param delivery_folder A single non-empty character string giving the name
#'   or relative path of the delivery folder below `delivery_path_root`.
#' @param work_path A single non-empty character string giving the work
#'   directory. The default is `file.path(input_path_root, "work")`.
#' @param ... Optional additional workflow-specific paths. Every value must be
#'   a single non-empty character string and must have a unique, non-empty
#'   name. Names that duplicate one of the standard path names are rejected.
#'
#' @return A named list containing:
#'
#' * `input_path_root`;
#' * `work_path`, `work_bam_path`, `work_fastq_path`, `work_qc_path`, and
#'   `work_ab1_path`;
#' * `delivery_path_root`, `delivery_folder`, and `delivery_path`;
#' * `delivery_bam_path`, `delivery_fastq_path`, `delivery_qc_path`, and
#'   `delivery_ab1_path`;
#' * any additional named paths supplied through `...`.
#'
#' The standard subdirectory names are `Bam`, `0_fastq`, `QC`, and `Ab1`.
#'
#' @details
#' Keeping path construction separate from directory creation makes the result
#' safe to use for dry runs, planning, validation, and command generation. A
#' calling workflow can create only the directories it actually needs with
#' [dir.create()].
#'
#' @examples
#' paths <- standard_paths(
#'   input_path_root = "/analysis/run01",
#'   delivery_path_root = "/delivery",
#'   delivery_folder = "project01"
#' )
#' paths$work_bam_path
#' paths$delivery_fastq_path
#'
#' # Supply a different work root and workflow-specific output paths.
#' paths2 <- standard_paths(
#'   input_path_root = "/analysis/run01",
#'   delivery_path_root = "/delivery",
#'   delivery_folder = "project01",
#'   work_path = "/scratch/run01",
#'   report_file = "/scratch/run01/report.xlsx"
#' )
#' paths2$report_file
#'
#' @export
standard_paths <- function(input_path_root,
                           delivery_path_root,
                           delivery_folder,
                           work_path = file.path(input_path_root, "work"),
                           ...) {
  check_standard_path_value(input_path_root, "input_path_root")
  check_standard_path_value(delivery_path_root, "delivery_path_root")
  check_standard_path_value(delivery_folder, "delivery_folder")
  check_standard_path_value(work_path, "work_path")

  delivery_path <- file.path(delivery_path_root, delivery_folder)
  paths <- list(
    input_path_root = input_path_root,
    work_path = work_path,
    work_bam_path = file.path(work_path, "Bam"),
    work_fastq_path = file.path(work_path, "0_fastq"),
    work_qc_path = file.path(work_path, "QC"),
    work_ab1_path = file.path(work_path, "Ab1"),
    delivery_path_root = delivery_path_root,
    delivery_folder = delivery_folder,
    delivery_path = delivery_path,
    delivery_bam_path = file.path(delivery_path, "Bam"),
    delivery_fastq_path = file.path(delivery_path, "0_fastq"),
    delivery_qc_path = file.path(delivery_path, "QC"),
    delivery_ab1_path = file.path(delivery_path, "Ab1")
  )

  extra_paths <- list(...)
  validate_standard_extra_paths(extra_paths, reserved_names = names(paths))
  c(paths, extra_paths)
}

#' Build the standard path layout for gene-editing results
#'
#' `standard_paths_gene_editing()` extends [standard_paths()] with the
#' sanitized reference FASTA and the complete and filtered variant-table Excel
#' files used by the gene-editing workflow.
#'
#' @param input_path_root A single non-empty character string giving the root
#'   directory containing the gene-editing analysis results.
#' @param delivery_path_root A single non-empty character string giving the
#'   delivery root.
#' @param delivery_folder A single non-empty character string giving the
#'   delivery folder below `delivery_path_root`. The default is
#'   `"gene_editing_delivery"`.
#' @param work_path A single non-empty character string giving the work
#'   directory. The default is `file.path(input_path_root, "work")`.
#' @param ref_sanitized_file A single non-empty character string giving the
#'   sanitized-reference filename relative to `input_path_root`. The default,
#'   `"reference_sanitized_seqID.fasta"`, matches the filename recognized by
#'   the ONTools wf-amplicon variant workflow.
#' @param xlsx_file_full A single non-empty character string giving the complete
#'   variant-table filename relative to `work_path`.
#' @param xlsx_file_filtered A single non-empty character string giving the
#'   filtered variant-table filename relative to `work_path`.
#'
#' @return A named list returned by [standard_paths()], extended with:
#'
#' * `ref_sanitized`;
#' * `xlsx_file_full`;
#' * `xlsx_file_filtered`.
#'
#' No files or directories are created.
#'
#' @examples
#' paths <- standard_paths_gene_editing(
#'   input_path_root = "/analysis/gene_editing/run01",
#'   delivery_path_root = "/delivery"
#' )
#' paths$ref_sanitized
#' paths$xlsx_file_full
#' paths$delivery_ab1_path
#'
#' # Put temporary results on a separate disk and customize filenames.
#' custom <- standard_paths_gene_editing(
#'   input_path_root = "/analysis/gene_editing/run01",
#'   delivery_path_root = "/delivery",
#'   delivery_folder = "project01_gene_editing",
#'   work_path = "/scratch/project01",
#'   xlsx_file_full = "all_variants.xlsx",
#'   xlsx_file_filtered = "reported_variants.xlsx"
#' )
#' custom$xlsx_file_filtered
#'
#' @export
standard_paths_gene_editing <- function(
    input_path_root,
    delivery_path_root,
    delivery_folder = "gene_editing_delivery",
    work_path = file.path(input_path_root, "work"),
    ref_sanitized_file = "reference_sanitized_seqID.fasta",
    xlsx_file_full = "variant_table_full.xlsx",
    xlsx_file_filtered = "variant_table_filtered.xlsx") {
  check_standard_path_value(ref_sanitized_file, "ref_sanitized_file")
  check_standard_path_value(xlsx_file_full, "xlsx_file_full")
  check_standard_path_value(xlsx_file_filtered, "xlsx_file_filtered")

  paths <- standard_paths(
    input_path_root = input_path_root,
    delivery_path_root = delivery_path_root,
    delivery_folder = delivery_folder,
    work_path = work_path
  )
  paths$event_cache_dir = file.path(work_path, "event_cache")
  paths$ref_sanitized <- file.path(input_path_root, ref_sanitized_file)
  paths$xlsx_file_full <- file.path(work_path, xlsx_file_full)
  paths$xlsx_file_filtered <- file.path(work_path, xlsx_file_filtered)
  paths
}

check_standard_path_value <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.",
         call. = FALSE)
  }
  invisible(x)
}

validate_standard_extra_paths <- function(extra_paths, reserved_names) {
  if (length(extra_paths) == 0L) {
    return(invisible(extra_paths))
  }

  extra_names <- names(extra_paths)
  if (is.null(extra_names) || anyNA(extra_names) || any(!nzchar(extra_names))) {
    stop("All additional paths supplied through `...` must be named.",
         call. = FALSE)
  }
  duplicated_names <- unique(extra_names[duplicated(extra_names)])
  if (length(duplicated_names) > 0L) {
    stop(
      "Additional path names must be unique; duplicated name(s): ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  reserved <- intersect(extra_names, reserved_names)
  if (length(reserved) > 0L) {
    stop(
      "Additional paths duplicate standard path name(s): ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  for (i in seq_along(extra_paths)) {
    check_standard_path_value(extra_paths[[i]], extra_names[[i]])
  }
  invisible(extra_paths)
}
