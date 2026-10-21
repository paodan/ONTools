#' Standardize sample information for gene-editing workflows
#'
#' `standard_sample_info_gene_editing()` selects user-specified columns from a
#' sample information table, converts them to the canonical ONTools
#' gene-editing schema, validates sample identifiers and editing intervals, and
#' adds standard source and work-file paths.
#'
#' @param sample_info A data frame containing at least the sample identifier,
#'   barcode, reference name, editing-region start, and editing-region end
#'   columns selected by the corresponding `*_column` arguments.
#' @param input_path_root A single non-empty character string giving the root
#'   directory containing per-barcode input results.
#' @param delivery_path_root A single non-empty character string giving the
#'   delivery root passed to [standard_paths_gene_editing()].
#' @param delivery_folder A single non-empty character string giving the
#'   delivery folder below `delivery_path_root`. The default is
#'   `"gene_editing_delivery"`.
#' @param work_path A single non-empty character string giving the work
#'   directory. The default is `file.path(input_path_root, "work")`.
#' @param sample_id_column,barcode_column,ref_name_column Character scalars
#'   identifying the source columns containing the sample identifier, barcode,
#'   and reference-sequence name.
#' @param ref_region_start_column,ref_region_end_column Character scalars
#'   identifying the source columns containing the closed, 1-based editing
#'   interval. Values must be positive whole numbers and start must not exceed
#'   end.
#' @param ref_sanitized_file A single non-empty character string giving the
#'   sanitized-reference filename relative to `input_path_root`. The default is
#'   `"reference_sanitized_seqIDs.fasta"`.
#'
#' @return A data frame with one row per input row and the following columns:
#'
#' * `sample_id`, `barcode`, `ref_name`, `region_start`, and `region_end`: the
#'   canonical ONTools gene-editing sample-information fields;
#' * `bam_file_origin`: the expected source BAM under
#'   `input_path_root/<barcode>/alignments/`;
#' * `bam_file`: the standardized work BAM under `work_path/Bam/`;
#' * `output_fastq`: the planned filtered FASTQ.GZ under `work_path/0_fastq/`;
#' * `reference_file`: the sanitized-reference path;
#' * `ab1_file`: the planned filtered AB1 file under `work_path/Ab1/`.
#'
#' @details
#' Source column names must exist and the five column-selection arguments must
#' identify five distinct columns. Character values are trimmed before they are
#' validated. Missing or empty sample IDs, barcodes, and reference names are
#' rejected.
#'
#' A sample ID may occur in multiple rows, for example when one barcode has
#' multiple reference intervals. However, sample IDs and barcodes must have a
#' one-to-one mapping. Sample IDs and barcodes containing `/` or `\\`, or equal
#' to `.` or `..`, are rejected because they are used as path components.
#'
#' The function only constructs and validates the table. It does not check that
#' input files exist, copy `bam_file_origin` to `bam_file`, or create any files
#' or directories. Consequently, `bam_file` is the planned work BAM path and
#' must be populated before it is used as an input to [run_editing_batch2()].
#'
#' @examples
#' sample_sheet <- data.frame(
#'   Sample_HiTOM = c("sample_A", "sample_A", "sample_B"),
#'   barcode = c("barcode001", "barcode001", "barcode002"),
#'   Reference_name = c("amplicon_A", "amplicon_B", "amplicon_A"),
#'   Reference_Region_Start = c(120, 80, 140),
#'   Reference_Region_End = c(150, 110, 170)
#' )
#'
#' standardized <- standard_sample_info_gene_editing(
#'   sample_info = sample_sheet,
#'   input_path_root = "/analysis/variant_results",
#'   delivery_path_root = "/delivery"
#' )
#' standardized[, c("sample_id", "barcode", "ref_name", "bam_file")]
#'
#' # Source columns can use project-specific names.
#' custom_sheet <- data.frame(
#'   specimen = "sample_C",
#'   barcode_id = "barcode003",
#'   target = "amplicon_C",
#'   edit_start = 50,
#'   edit_end = 90
#' )
#' standard_sample_info_gene_editing(
#'   sample_info = custom_sheet,
#'   input_path_root = "/analysis/variant_results",
#'   delivery_path_root = "/delivery",
#'   sample_id_column = "specimen",
#'   barcode_column = "barcode_id",
#'   ref_name_column = "target",
#'   ref_region_start_column = "edit_start",
#'   ref_region_end_column = "edit_end"
#' )
#'
#' @seealso [standard_paths_gene_editing()], [run_editing_batch2()]
#' @export
standard_sample_info_gene_editing <- function(
    sample_info,
    input_path_root,
    delivery_path_root,
    delivery_folder = "gene_editing_delivery",
    work_path = file.path(input_path_root, "work"),
    sample_id_column = "Sample_HiTOM",
    barcode_column = "barcode",
    ref_name_column = "Reference_name",
    ref_region_start_column = "Reference_Region_Start",
    ref_region_end_column = "Reference_Region_End",
    ref_sanitized_file = "reference_sanitized_seqIDs.fasta") {
  if (!is.data.frame(sample_info)) {
    stop("`sample_info` must be a data frame.", call. = FALSE)
  }
  if (nrow(sample_info) == 0L) {
    stop("`sample_info` must contain at least one row.", call. = FALSE)
  }

  column_map <- c(
    sample_id = sample_id_column,
    barcode = barcode_column,
    ref_name = ref_name_column,
    region_start = ref_region_start_column,
    region_end = ref_region_end_column
  )
  column_argument_names <- c(
    sample_id = "sample_id_column",
    barcode = "barcode_column",
    ref_name = "ref_name_column",
    region_start = "ref_region_start_column",
    region_end = "ref_region_end_column"
  )
  for (name in names(column_map)) {
    check_gene_editing_column_name(
      column_map[[name]],
      column_argument_names[[name]]
    )
  }
  duplicated_columns <- unique(column_map[duplicated(column_map)])
  if (length(duplicated_columns) > 0L) {
    stop(
      "Column-selection arguments must identify distinct columns; duplicated column(s): ",
      paste(duplicated_columns, collapse = ", "),
      call. = FALSE
    )
  }
  missing_columns <- setdiff(unname(column_map), names(sample_info))
  if (length(missing_columns) > 0L) {
    stop(
      "`sample_info` is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  paths <- standard_paths_gene_editing(
    input_path_root = input_path_root,
    delivery_path_root = delivery_path_root,
    delivery_folder = delivery_folder,
    work_path = work_path,
    ref_sanitized_file = ref_sanitized_file
  )

  standardized <- as.data.frame(
    sample_info[, unname(column_map), drop = FALSE],
    stringsAsFactors = FALSE
  )
  names(standardized) <- names(column_map)
  rownames(standardized) <- NULL

  standardized$sample_id <- trim_gene_editing_text(
    standardized$sample_id,
    "sample_id"
  )
  standardized$barcode <- trim_gene_editing_text(
    standardized$barcode,
    "barcode"
  )
  standardized$ref_name <- trim_gene_editing_text(
    standardized$ref_name,
    "ref_name"
  )
  standardized$region_start <- parse_gene_editing_coordinate(
    standardized$region_start,
    "region_start"
  )
  standardized$region_end <- parse_gene_editing_coordinate(
    standardized$region_end,
    "region_end"
  )

  invalid_interval <- standardized$region_start > standardized$region_end
  if (any(invalid_interval)) {
    stop(
      "Editing-region start must be less than or equal to end; invalid row(s): ",
      paste(which(invalid_interval), collapse = ", "),
      call. = FALSE
    )
  }
  validate_gene_editing_path_component(standardized$sample_id, "sample_id")
  validate_gene_editing_path_component(standardized$barcode, "barcode")
  validate_gene_editing_sample_barcode_mapping(
    standardized$sample_id,
    standardized$barcode
  )

  standardized$bam_file_origin <- file.path(
    paths$input_path_root,
    standardized$barcode,
    "alignments",
    paste0(standardized$barcode, ".aligned.sorted.bam")
  )
  standardized$bam_file <- file.path(
    paths$work_bam_path,
    paste0(standardized$sample_id, ".aligned.sorted.bam")
  )
  standardized$output_fastq <- file.path(
    paths$work_fastq_path,
    paste0(standardized$sample_id, "_filtered.fastq.gz")
  )
  standardized$reference_file <- paths$ref_sanitized
  standardized$ab1_file <- file.path(
    paths$work_ab1_path,
    paste0(standardized$sample_id, "_filtered.ab1")
  )
  standardized
}

check_gene_editing_column_name <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.",
         call. = FALSE)
  }
  invisible(x)
}

trim_gene_editing_text <- function(x, name) {
  value <- trimws(as.character(x))
  if (anyNA(value) || any(!nzchar(value))) {
    stop("`sample_info$", name, "` contains missing or empty values.",
         call. = FALSE)
  }
  value
}

parse_gene_editing_coordinate <- function(x, name) {
  text <- trimws(as.character(x))
  value <- suppressWarnings(as.numeric(text))
  invalid <- is.na(text) | !nzchar(text) | is.na(value) | !is.finite(value) |
    value < 1 | value != floor(value) | value > .Machine$integer.max
  if (any(invalid)) {
    stop(
      "`sample_info$", name,
      "` must contain positive whole numbers; invalid row(s): ",
      paste(which(invalid), collapse = ", "),
      call. = FALSE
    )
  }
  as.integer(value)
}

validate_gene_editing_sample_barcode_mapping <- function(sample_id, barcode) {
  mapping <- unique(data.frame(
    sample_id = sample_id,
    barcode = barcode,
    stringsAsFactors = FALSE
  ))
  duplicated_sample <- mapping$sample_id[duplicated(mapping$sample_id)]
  if (length(duplicated_sample) > 0L) {
    stop(
      "Each `sample_id` must map to exactly one barcode; conflicting sample ID(s): ",
      paste(unique(duplicated_sample), collapse = ", "),
      call. = FALSE
    )
  }
  duplicated_barcode <- mapping$barcode[duplicated(mapping$barcode)]
  if (length(duplicated_barcode) > 0L) {
    stop(
      "Each barcode must map to exactly one `sample_id`; conflicting barcode(s): ",
      paste(unique(duplicated_barcode), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(mapping)
}

validate_gene_editing_path_component <- function(x, name) {
  unsafe <- x %in% c(".", "..") | grepl("[/\\\\]", x)
  if (any(unsafe)) {
    stop(
      "`", name, "` contains unsafe path component(s) in row(s): ",
      paste(which(unsafe), collapse = ", "),
      ". Values must not be `.` or `..` or contain path separators.",
      call. = FALSE
    )
  }
  invisible(x)
}
