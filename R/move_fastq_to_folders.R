#' Move demultiplexed FASTQ barcode folders into project groups
#'
#' `move_fastq_to_folders()` reorganizes `barcode*` folders from one FASTQ
#' directory into subfolders defined by a sample information table. This is
#' useful when one sequencing run contains amplicons from multiple projects or
#' multiple expected amplicon lengths.
#'
#' @param fastq_dir Directory containing demultiplexed `barcode*` folders.
#' @param sample_info Sample information table, either a data frame or a CSV
#'   file path.
#' @param project_col Optional column in `sample_info` identifying the project.
#'   When `NULL`, output groups are based only on `amplicon_size_col`.
#' @param amplicon_size_col Column in `sample_info` identifying expected amplicon
#'   size or length group.
#' @param barcode_col Column in `sample_info` identifying barcode folder names. Default `Barcode_ID`.
#' @param overwrite Logical. If `FALSE`, stop when a destination barcode folder
#'   already exists.
#' @param dry_run Logical. If `TRUE`, return the move plan without moving files.
#'
#' @return Invisibly returns a data frame with one row per barcode in
#'   `sample_info`, including the output group, source path, destination path,
#'   source/destination existence, and move status.
#'
#' @examples
#' fastq_dir <- tempfile("fastq-pass-trim-")
#' dir.create(file.path(fastq_dir, "barcode01"), recursive = TRUE)
#' dir.create(file.path(fastq_dir, "barcode02"), recursive = TRUE)
#' sample_info <- data.frame(
#'   barcode = c("barcode01", "barcode02"),
#'   Project = c("Project A", "Project B"),
#'   Expected_Amplicon_Size_bp = c("3500bp", "1600bp")
#' )
#'
#' move_fastq_to_folders(
#'   fastq_dir = fastq_dir,
#'   sample_info = sample_info,
#'   project_col = "Project",
#'   dry_run = TRUE
#' )
#'
#' @export
move_fastq_to_folders <- function(fastq_dir,
                                  sample_info,
                                  project_col = NULL,
                                  amplicon_size_col = "Expected_Amplicon_Size_bp",
                                  barcode_col = "Barcode_ID",
                                  overwrite = FALSE,
                                  dry_run = FALSE) {
  check_dir_arg(fastq_dir, "fastq_dir")
  check_scalar_character(amplicon_size_col, "amplicon_size_col")
  check_scalar_character(barcode_col, "barcode_col")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  if (!is.null(project_col)) {
    check_scalar_character(project_col, "project_col")
  }

  fastq_dir <- normalizePath(fastq_dir, mustWork = TRUE)
  sample_info <- read_sample_info_table(sample_info)
  required_cols <- c(barcode_col, amplicon_size_col)
  if (!is.null(project_col)) {
    required_cols <- c(required_cols, project_col)
  }
  missing_cols <- setdiff(required_cols, names(sample_info))
  if (length(missing_cols) > 0L) {
    stop(
      "`sample_info` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (barcode_co == "Barcode_ID"){
    barcodes <- paste0("barcode", gsub(".+-([0-9]+$)","\\1", sample_info[[barcode_col]]))
  } else {
    barcodes <- as.character(sample_info[[barcode_col]])
  }
  amplicon_sizes <- as.character(sample_info[[amplicon_size_col]])
  projects <- if (is.null(project_col)) {
    rep(NA_character_, nrow(sample_info))
  } else {
    as.character(sample_info[[project_col]])
  }

  valid_rows <- !is.na(barcodes) & nzchar(barcodes) &
    !is.na(amplicon_sizes) & nzchar(amplicon_sizes)
  if (!is.null(project_col)) {
    valid_rows <- valid_rows & !is.na(projects) & nzchar(projects)
  }
  if (any(!valid_rows)) {
    warning(
      "Dropping ",
      sum(!valid_rows),
      " row(s) with missing barcode, project, or amplicon size.",
      call. = FALSE
    )
  }

  barcodes <- barcodes[valid_rows]
  amplicon_sizes <- amplicon_sizes[valid_rows]
  projects <- projects[valid_rows]

  groups <- if (is.null(project_col)) {
    safe_path_component(amplicon_sizes)
  } else {
    paste(
      safe_path_component(projects),
      safe_path_component(amplicon_sizes),
      sep = "_"
    )
  }

  plan <- data.frame(
    barcode = barcodes,
    project = projects,
    amplicon_size = amplicon_sizes,
    group = groups,
    source = file.path(fastq_dir, barcodes),
    destination = file.path(fastq_dir, groups, barcodes),
    stringsAsFactors = FALSE
  )
  plan$source_exists <- dir.exists(plan$source)
  plan$destination_exists <- dir.exists(plan$destination)
  plan$moved <- FALSE
  plan$status <- "planned"

  missing_sources <- !plan$source_exists
  if (any(missing_sources)) {
    warning(
      "These barcode folders were not found in `fastq_dir`: ",
      paste(plan$barcode[missing_sources], collapse = ", "),
      call. = FALSE
    )
    plan$status[missing_sources] <- "missing_source"
  }

  blocked <- plan$destination_exists & !isTRUE(overwrite)
  if (any(blocked)) {
    stop(
      "Destination barcode folder(s) already exist: ",
      paste(plan$destination[blocked], collapse = ", "),
      ". Use `overwrite = TRUE` to replace them.",
      call. = FALSE
    )
  }

  if (isTRUE(dry_run)) {
    plan$status[plan$source_exists] <- "dry_run"
    return(invisible(plan))
  }

  movable <- plan$source_exists
  for (i in which(movable)) {
    dir.create(dirname(plan$destination[[i]]), recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(plan$destination[[i]]) && isTRUE(overwrite)) {
      unlink(plan$destination[[i]], recursive = TRUE, force = TRUE)
    }
    plan$moved[[i]] <- file.rename(plan$source[[i]], plan$destination[[i]])
    plan$status[[i]] <- if (isTRUE(plan$moved[[i]])) {
      "moved"
    } else {
      "failed"
    }
  }

  invisible(plan)
}

read_sample_info_table <- function(sample_info) {
  if (is.data.frame(sample_info)) {
    return(sample_info)
  }

  if (is.character(sample_info) && length(sample_info) == 1L &&
      !is.na(sample_info) && file.exists(sample_info)) {
    return(utils::read.csv(sample_info, stringsAsFactors = FALSE))
  }

  stop("`sample_info` must be a data frame or an existing CSV file path.",
       call. = FALSE)
}

safe_path_component <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[/\\\\:]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}
