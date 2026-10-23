#' Distribute split FASTQ files by peak and read-length range
#'
#' Parse peak numbers and peak lengths from FASTQ filenames produced by
#' [split_fastq()] or [split_fastqs()], group the files by configured length
#' ranges, and copy every file into its own directory. The original filename
#' and modification time are retained.
#'
#' @param fastqs_split A non-empty character vector of existing split FASTQ
#'   files. Every basename must end in
#'   `"_peak_<number>_<length>bp.fastq.gz"`, for example
#'   `"sample01_peak_2_720bp.fastq.gz"`.
#' @param out_path A single non-empty path to the output root directory.
#' @param size_breaks A strictly increasing numeric vector containing at least
#'   two finite boundaries. Non-`peak_1` files, or all files when
#'   `separate_peak1 = FALSE`, are grouped into consecutive intervals. The
#'   first interval includes its lower boundary and every interval includes its
#'   upper boundary.
#' @param save_diff_size_folder Logical. If `TRUE`, add the group name as a
#'   directory between `out_path` and the file-specific directory. If `FALSE`,
#'   file-specific directories are created directly under `out_path`.
#' @param separate_peak1 Logical. If `TRUE`, files labelled `peak_1` bypass
#'   length binning and are assigned to the `peak1_target_range` group.
#' @param overwrite Logical. If `TRUE`, replace an existing destination file.
#'   If `FALSE`, retain it and report `destination_exists` in the returned
#'   status table.
#' @param out_of_range How to handle files whose peak length lies outside
#'   `size_breaks`. `"warn"` skips them and issues a warning, `"error"` stops
#'   before copying any file, and `"drop"` skips them without a warning.
#'
#' @return Invisibly returns a data frame with one row per input file and the
#'   following columns:
#'
#' * `source`: input FASTQ path;
#' * `peak`: parsed peak label such as `"peak_2"`;
#' * `peak_size`: parsed peak length in base pairs;
#' * `group`: length-range group or `"peak1_target_range"`;
#' * `destination`: complete output path, or `NA` when the file is outside the
#'   configured ranges;
#' * `copied`: whether the file was copied by this call;
#' * `status`: one of `"copied"`, `"destination_exists"`, `"copy_failed"`, or
#'   `"out_of_range"`.
#'
#' @details
#' With the default `size_breaks = c(500, 1000, 2000)`, the generated groups
#' are `read_len_range__500_1000` and `read_len_range__1000_2000`. A length of
#' 500 is included in the first group, 1000 is included in the first group,
#' and 2000 is included in the second group.
#'
#' When `save_diff_size_folder = TRUE`, a file such as
#' `sample01_peak_2_720bp.fastq.gz` is copied to:
#'
#' ```
#' <out_path>/read_len_range__500_1000/
#'   sample01_peak_2_720bp/sample01_peak_2_720bp.fastq.gz
#' ```
#'
#' When `save_diff_size_folder = FALSE`, the length-group directory is omitted:
#'
#' ```
#' <out_path>/sample01_peak_2_720bp/sample01_peak_2_720bp.fastq.gz
#' ```
#'
#' Filenames are parsed strictly so malformed names cannot be silently placed
#' in an incorrect directory. File-specific directory names are derived from
#' the FASTQ basename without `.fastq.gz`; characters other than ASCII letters,
#' digits, underscores, and hyphens are replaced with underscores. Ambiguous
#' duplicate destinations are rejected before copying.
#'
#' @examples
#' source_dir <- tempfile("split-fastqs-")
#' output_dir <- tempfile("distributed-fastqs-")
#' dir.create(source_dir)
#'
#' fastqs <- file.path(
#'   source_dir,
#'   c(
#'     "sample01_peak_1_1600bp.fastq.gz",
#'     "sample01_peak_2_720bp.fastq.gz",
#'     "sample02_peak_3_1300bp.fastq.gz"
#'   )
#' )
#' invisible(file.create(fastqs))
#'
#' result <- distribute_fastqs(
#'   fastqs_split = fastqs,
#'   out_path = output_dir,
#'   size_breaks = c(500, 1000, 2000),
#'   save_diff_size_folder = TRUE
#' )
#' result[, c("peak", "peak_size", "group", "copied")]
#' file.exists(result$destination)
#'
#' @seealso [split_fastq()], [split_fastqs()], [collect_items()]
#' @export
distribute_fastqs <- function(
    fastqs_split,
    out_path = "0_RawData/ITS_16S_mix/SD261021100438_split_distribute",
    size_breaks = c(500, 1000, 2000),
    save_diff_size_folder = FALSE,
    separate_peak1 = TRUE,
    overwrite = TRUE,
    out_of_range = c("warn", "error", "drop")) {
  if (!is.character(fastqs_split) || length(fastqs_split) < 1L ||
      anyNA(fastqs_split) || any(!nzchar(trimws(fastqs_split)))) {
    stop("`fastqs_split` must be a non-empty character vector of paths.",
         call. = FALSE)
  }
  fastqs_split <- trimws(fastqs_split)
  duplicated_sources <- unique(fastqs_split[duplicated(fastqs_split)])
  if (length(duplicated_sources) > 0L) {
    stop(
      "`fastqs_split` contains duplicated path(s): ",
      paste(duplicated_sources, collapse = ", "),
      call. = FALSE
    )
  }
  missing_fastqs <- !file.exists(fastqs_split) | dir.exists(fastqs_split)
  if (any(missing_fastqs)) {
    stop(
      "FASTQ file(s) do not exist: ",
      paste(fastqs_split[missing_fastqs], collapse = ", "),
      call. = FALSE
    )
  }

  check_scalar_character(out_path, "out_path")
  check_logical_scalar(save_diff_size_folder, "save_diff_size_folder")
  check_logical_scalar(separate_peak1, "separate_peak1")
  check_logical_scalar(overwrite, "overwrite")
  out_of_range <- match.arg(out_of_range)

  if (!is.numeric(size_breaks) || length(size_breaks) < 2L ||
      anyNA(size_breaks) || any(!is.finite(size_breaks)) ||
      any(diff(size_breaks) <= 0)) {
    stop(
      "`size_breaks` must contain at least two finite, strictly increasing ",
      "numbers.",
      call. = FALSE
    )
  }

  file_names <- basename(fastqs_split)
  file_pattern <- "_peak_[0-9]+_[0-9]+bp\\.fastq\\.gz$"
  malformed <- !grepl(file_pattern, file_names)
  if (any(malformed)) {
    stop(
      "FASTQ basename(s) do not match ",
      "`*_peak_<number>_<length>bp.fastq.gz`: ",
      paste(file_names[malformed], collapse = ", "),
      call. = FALSE
    )
  }

  peak <- sub(
    "^.*_(peak_[0-9]+)_[0-9]+bp\\.fastq\\.gz$",
    "\\1",
    file_names
  )
  peak_size <- as.integer(sub(
    "^.*_peak_[0-9]+_([0-9]+)bp\\.fastq\\.gz$",
    "\\1",
    file_names
  ))

  format_break <- function(x) {
    value <- format(x, scientific = FALSE, trim = TRUE)
    gsub("[^A-Za-z0-9_-]", "_", value)
  }
  bin_labels <- paste0(
    "read_len_range__",
    format_break(size_breaks[-length(size_breaks)]),
    "_",
    format_break(size_breaks[-1L])
  )
  bins <- cut(
    peak_size,
    breaks = size_breaks,
    labels = bin_labels,
    include.lowest = TRUE,
    right = TRUE
  )
  group <- as.character(bins)
  if (isTRUE(separate_peak1)) {
    group[peak == "peak_1"] <- "peak1_target_range"
  }

  outside <- is.na(group)
  if (any(outside) && identical(out_of_range, "error")) {
    stop(
      "Peak length is outside `size_breaks` for file(s): ",
      paste(file_names[outside], collapse = ", "),
      call. = FALSE
    )
  }
  if (any(outside) && identical(out_of_range, "warn")) {
    warning(
      "Skipping file(s) with peak length outside `size_breaks`: ",
      paste(file_names[outside], collapse = ", "),
      call. = FALSE
    )
  }

  file_stem <- sub("\\.fastq\\.gz$", "", file_names)
  folder <- gsub("[^A-Za-z0-9_-]", "_", file_stem)
  destination_parent <- if (isTRUE(save_diff_size_folder)) {
    file.path(out_path, group, folder)
  } else {
    file.path(out_path, folder)
  }
  destination <- file.path(destination_parent, file_names)
  destination[outside] <- NA_character_

  result <- data.frame(
    source = fastqs_split,
    peak = peak,
    peak_size = peak_size,
    group = group,
    destination = destination,
    copied = FALSE,
    status = ifelse(outside, "out_of_range", NA_character_),
    stringsAsFactors = FALSE
  )

  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(out_path)) {
    stop("Failed to create `out_path`: ", out_path, call. = FALSE)
  }

  eligible <- which(!outside)
  if (length(eligible) > 0L) {
    copied <- collect_items(
      source = fastqs_split[eligible],
      destination_dir = destination_parent[eligible],
      overwrite = overwrite
    )
    result$destination[eligible] <- copied$destination
    result$copied[eligible] <- copied$copied
    destination_exists <- file.exists(copied$destination)
    result$status[eligible] <- ifelse(
      copied$copied,
      "copied",
      ifelse(
        !isTRUE(overwrite) & destination_exists,
        "destination_exists",
        "copy_failed"
      )
    )
  }

  invisible(result)
}
