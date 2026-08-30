
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
                              min_read_fraction = 0.8,
                              max_read_fraction = 1.2,
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
  check_fraction_scalar(min_read_fraction, "min_read_fraction")
  check_fraction_scalar(max_read_fraction, "max_read_fraction")
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
      round(readLen[missing_min] * min_read_fraction)
    sampleInfo_2[["Max_Read_Length"]][missing_max] <-
      round(readLen[missing_max] * max_read_fraction)
  }

  keep <- unique(c(keep, "Folders"))
  keep <- intersect(keep, names(sampleInfo_2))
  sampleInfo_2 <- sampleInfo_2[, keep, drop = FALSE]

  split(sampleInfo_2, sampleInfo_2[["Folders"]], drop = TRUE)
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

check_fraction_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0) {
    stop("`", name, "` must be a single non-negative numeric value.",
         call. = FALSE)
  }
}
