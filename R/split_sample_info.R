#' Split sample information by project and expected amplicon size
#'
#' `split_sample_info()` reads a sample information table, creates a `Folders`
#' grouping column from selected metadata columns, fills missing read-length
#' filter columns from the expected amplicon size, and optionally writes one CSV
#' file per group.
#'
#' @param sampleInfo Sample information table, either a data frame or an
#'   existing file path. CSV files are read when `sheet = NULL`; otherwise the
#'   file is treated as an Excel workbook and read with `readxl`.
#' @param by Character vector of columns used to create the `Folders` group
#'   name. Values are sanitized with the same path-safe naming used by
#'   [move_fastq_to_folders()].
#' @param keep Character vector of columns to keep in the returned and saved
#'   sample information tables. `Folders` is always kept.
#' @param sheet Excel sheet name. Set to `NULL` to read `sampleInfo` as CSV.
#' @param skip Number of rows to skip before reading the table. When reading an
#'   Excel sheet and `skip = NULL`, the function looks for a `[data]` marker in
#'   the first column and uses the next row as the table header.
#' @param min_read_delta Non-negative number subtracted from `Expected_Size_bp`
#'   to fill missing `Min_Read_Length` values. Defaults to `300`.
#' @param max_read_delta Non-negative number added to `Expected_Size_bp` to fill
#'   missing `Max_Read_Length` values. Defaults to `200`.
#' @param output_dir Optional directory where one CSV file per `Folders` group
#'   is written. When `NULL`, no files are written.
#' @param na_on_size_error Logical. Passed to [parse_size()]. If `TRUE`, an
#'   unparseable `Expected_Size_bp` value leaves the derived read-length limits
#'   as `NA` instead of throwing an error.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`data`}{A named list of data frames split by `Folders`.}
#'     \item{`files`}{A named list of CSV file paths, or an empty list when
#'       `output_dir = NULL`.}
#'   }
#'
#' @examples
#' sample_info <- data.frame(
#'   Sample_ID = c("S1", "S2"),
#'   Project_ID = c("P1", "P1"),
#'   Barcode_ID = c("NB01", "NB02"),
#'   Expected_Size_bp = c("1000bp", "1500bp")
#' )
#'
#' split_sample_info(sample_info)
#'
#' @export
split_sample_info <- function(sampleInfo,
                              by = c("Project_ID", "Expected_Size_bp"),
                              keep = c(
                                "Sample_ID", "Project_ID", "Barcode_ID",
                                "Sample_Type", "Expected_Size_bp",
                                "Min_Read_Length", "Max_Read_Length",
                                "Gel_lane", "Gel_Quality", "Clear_Band",
                                "Analysis_Type", "Primer_F", "Primer_R"
                              ),
                              sheet = "2_metadata",
                              skip = NULL,
                              min_read_delta = 300,
                              max_read_delta = 200,
                              output_dir = NULL,
                              na_on_size_error = TRUE) {
  if (!is.character(by) || length(by) < 1L || anyNA(by) || any(!nzchar(by))) {
    stop("`by` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(keep) || anyNA(keep) || any(!nzchar(keep))) {
    stop("`keep` must be a character vector.", call. = FALSE)
  }
  if (!is.null(sheet)) {
    check_scalar_character(sheet, "sheet")
  }
  if (!is.null(skip) &&
      (!is.numeric(skip) || length(skip) != 1L || is.na(skip) || skip < 0)) {
    stop("`skip` must be a single non-negative number or `NULL`.",
         call. = FALSE)
  }
  check_non_negative_scalar(min_read_delta, "min_read_delta")
  check_non_negative_scalar(max_read_delta, "max_read_delta")
  if (!is.null(output_dir)) {
    check_scalar_character(output_dir, "output_dir")
  }
  check_logical_scalar(na_on_size_error, "na_on_size_error")

  sampleInfo_2 <- read_split_sample_info(sampleInfo, sheet = sheet, skip = skip)
  missing_by <- setdiff(by, names(sampleInfo_2))
  if (length(missing_by) > 0L) {
    stop(
      "`sampleInfo` is missing grouping column(s): ",
      paste(missing_by, collapse = ", "),
      call. = FALSE
    )
  }

  sampleInfo_2[["Folders"]] <- apply(
    sampleInfo_2[, by, drop = FALSE],
    1,
    function(row) paste(safe_path_component(row), collapse = "_")
  )

  if ("Expected_Size_bp" %in% names(sampleInfo_2)) {
    readLen <- parse_size(
      as.character(sampleInfo_2[["Expected_Size_bp"]]),
      na_on_error = na_on_size_error
    )
    if (!"Min_Read_Length" %in% names(sampleInfo_2)) {
      sampleInfo_2[["Min_Read_Length"]] <- NA_real_
    }
    if (!"Max_Read_Length" %in% names(sampleInfo_2)) {
      sampleInfo_2[["Max_Read_Length"]] <- NA_real_
    }

    missing_min <- is.na(sampleInfo_2[["Min_Read_Length"]]) & !is.na(readLen)
    missing_max <- is.na(sampleInfo_2[["Max_Read_Length"]]) & !is.na(readLen)
    sampleInfo_2[["Min_Read_Length"]][missing_min] <-
      pmax(1, round(readLen[missing_min] - min_read_delta))
    sampleInfo_2[["Max_Read_Length"]][missing_max] <-
      round(readLen[missing_max] + max_read_delta)
  }

  keep <- unique(c(keep, "Folders"))
  keep <- intersect(keep, names(sampleInfo_2))
  sampleInfo_2 <- sampleInfo_2[, keep, drop = FALSE]

  data <- split(sampleInfo_2, sampleInfo_2[["Folders"]], drop = TRUE)
  files <- list()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    files <- lapply(names(data), function(folder) {
      file.path(output_dir, paste0(folder, ".csv"))
    })
    names(files) <- names(data)

    for (folder in names(data)) {
      utils::write.csv(data[[folder]], files[[folder]], row.names = FALSE)
    }
  }

  list(data = data, files = files)
}

read_split_sample_info <- function(sampleInfo, sheet, skip) {
  if (is.data.frame(sampleInfo)) {
    return(sampleInfo)
  }

  if (!is.character(sampleInfo) || length(sampleInfo) != 1L ||
      is.na(sampleInfo) || !file.exists(sampleInfo)) {
    stop("`sampleInfo` must be a data frame or an existing file path.",
         call. = FALSE)
  }

  if (!is.null(sheet)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Package `readxl` is required to read Excel sample information files.",
        call. = FALSE
      )
    }
    if (is.null(skip)) {
      preview <- suppressMessages(readxl::read_excel(path = sampleInfo, sheet = sheet))
      marker <- which(preview[[1]] == "[data]")
      skip <- if (length(marker) > 0L) marker[[1L]] + 1L else 0
    }
    sample_info <- suppressMessages(readxl::read_excel(
      path = sampleInfo,
      sheet = sheet,
      skip = skip
    ))
    sample_info <- as.data.frame(sample_info)
    names(sample_info) <- trimws(names(sample_info))
    return(sample_info)
  }

  if (is.null(skip)) {
    skip <- 0
  }
  utils::read.csv(sampleInfo, header = TRUE, sep = ",", skip = skip,
                  stringsAsFactors = FALSE)
}

check_non_negative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0) {
    stop("`", name, "` must be a single non-negative numeric value.",
         call. = FALSE)
  }
}
