#' Create directories from a standard ONTools path layout
#'
#' `create_standard_directories()` creates the work and delivery directories
#' defined by [standard_paths()] or [standard_paths_gene_editing()]. Standard
#' file paths and additional entries are not treated as directories.
#'
#' @param paths A named list returned by [standard_paths()] or
#'   [standard_paths_gene_editing()].
#' @param create_work Logical. If `TRUE` (the default), create `work_path` and
#'   its standard `Bam`, `0_fastq`, `QC`, and `Ab1` subdirectories.
#' @param create_delivery Logical. If `TRUE` (the default), create
#'   `delivery_path_root`, `delivery_path`, and the standard delivery `Bam`,
#'   `0_fastq`, `QC`, and `Ab1` subdirectories.
#' @param dry_run Logical. If `TRUE`, validate and return the directory plan
#'   without creating directories. Default is `FALSE`.
#' @param show_warnings Logical. Passed to [dir.create()] as `showWarnings`.
#'   Default is `FALSE`. Creation failures are always reported as errors by
#'   this function.
#'
#' @return Invisibly returns a data frame with one row per selected directory
#'   and the columns:
#'
#' * `name`: corresponding name in `paths`;
#' * `path`: directory path;
#' * `existed_before`: whether the directory existed before this call;
#' * `created`: whether a previously absent directory exists after this call;
#' * `exists_after`: whether the directory exists after this call;
#' * `status`: one of `"exists"`, `"created"`, or `"planned"`.
#'
#' In dry-run mode, `created` is always `FALSE`, `exists_after` records the
#' current state, and a missing directory has status `"planned"`.
#'
#' @details
#' All directories are created recursively, so missing parent directories are
#' created as needed. Before creating anything, the function checks every
#' selected target. If a target exists as a regular file rather than a
#' directory, the function stops without creating any planned directory.
#'
#' `input_path_root`, `ref_sanitized`, Excel paths, and additional paths passed
#' through `...` to [standard_paths()] are not selected automatically. Note
#' that a missing `input_path_root` can still be created as a parent when the
#' default `work_path = file.path(input_path_root, "work")` is used.
#'
#' @examples
#' root <- tempfile("ontools-paths-")
#' paths <- standard_paths(
#'   input_path_root = file.path(root, "analysis"),
#'   delivery_path_root = file.path(root, "delivery_root"),
#'   delivery_folder = "project01"
#' )
#'
#' # Inspect what would be created.
#' plan <- create_standard_directories(paths, dry_run = TRUE)
#' plan[, c("name", "path", "status")]
#'
#' # Create the complete standard directory layout.
#' result <- create_standard_directories(paths)
#' result[, c("name", "status")]
#'
#' # Create only work directories.
#' work_root <- tempfile("ontools-work-only-")
#' work_paths <- standard_paths(
#'   input_path_root = work_root,
#'   delivery_path_root = file.path(work_root, "delivery"),
#'   delivery_folder = "project01"
#' )
#' create_standard_directories(
#'   work_paths,
#'   create_delivery = FALSE
#' )
#'
#' @export
create_standard_directories <- function(paths,
                                        create_work = TRUE,
                                        create_delivery = TRUE,
                                        dry_run = FALSE,
                                        show_warnings = FALSE) {
  if (!is.list(paths) || is.null(names(paths))) {
    stop("`paths` must be a named list returned by a standard path function.",
         call. = FALSE)
  }
  check_standard_directory_flag(create_work, "create_work")
  check_standard_directory_flag(create_delivery, "create_delivery")
  check_standard_directory_flag(dry_run, "dry_run")
  check_standard_directory_flag(show_warnings, "show_warnings")
  if (!isTRUE(create_work) && !isTRUE(create_delivery)) {
    stop("At least one of `create_work` and `create_delivery` must be TRUE.",
         call. = FALSE)
  }

  work_names <- c(
    "work_path",
    "work_bam_path",
    "work_fastq_path",
    "work_qc_path",
    "work_ab1_path"
  )
  delivery_names <- c(
    "delivery_path_root",
    "delivery_path",
    "delivery_bam_path",
    "delivery_fastq_path",
    "delivery_qc_path",
    "delivery_ab1_path"
  )
  selected_names <- c(
    if (isTRUE(create_work)) work_names,
    if (isTRUE(create_delivery)) delivery_names
  )
  missing_names <- setdiff(selected_names, names(paths))
  if (length(missing_names) > 0L) {
    stop(
      "`paths` is missing required directory path(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  directory_paths <- vapply(
    selected_names,
    function(name) {
      value <- paths[[name]]
      check_standard_directory_path(value, name)
      value
    },
    character(1)
  )
  unique_path <- !duplicated(path.expand(directory_paths))
  selected_names <- selected_names[unique_path]
  directory_paths <- unname(directory_paths[unique_path])

  file_collision <- file.exists(directory_paths) & !dir.exists(directory_paths)
  if (any(file_collision)) {
    stop(
      "Cannot create a directory because the target exists as a file: ",
      paste(directory_paths[file_collision], collapse = ", "),
      call. = FALSE
    )
  }

  existed_before <- dir.exists(directory_paths)
  if (!isTRUE(dry_run)) {
    for (path in directory_paths[!existed_before]) {
      dir.create(
        path,
        recursive = TRUE,
        showWarnings = show_warnings
      )
      if (!dir.exists(path)) {
        stop("Failed to create directory: ", path, call. = FALSE)
      }
    }
  }

  exists_after <- dir.exists(directory_paths)
  created <- !existed_before & exists_after
  status <- ifelse(
    existed_before,
    "exists",
    ifelse(created, "created", "planned")
  )
  invisible(data.frame(
    name = selected_names,
    path = directory_paths,
    existed_before = existed_before,
    created = created,
    exists_after = exists_after,
    status = status,
    stringsAsFactors = FALSE
  ))
}

check_standard_directory_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

check_standard_directory_path <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(
      "`paths$", name, "` must be a single non-empty character string.",
      call. = FALSE
    )
  }
  invisible(x)
}
