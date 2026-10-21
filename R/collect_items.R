#' Copy files and directories into destination directories
#'
#' `collect_items()` copies one or more source files or directories into target
#' parent directories. Item type is detected automatically, source basenames
#' are always retained, and file modification times are preserved with
#' `file.copy(copy.date = TRUE)`.
#'
#' @param source A non-empty character vector of source paths. Existing regular
#'   files and directories may be mixed in the same call. Missing sources are
#'   reported in the result without stopping the remaining copies.
#' @param destination_dir A character vector containing target parent
#'   directories. Supply one directory for all sources, or one directory per
#'   source. The final destination of each item is
#'   `file.path(destination_dir, basename(source))`; this function never renames
#'   a source item.
#' @param label Optional character vector used to identify operations in the
#'   returned table. `NULL` uses each source basename. A single label is
#'   recycled to all sources; otherwise its length must equal `length(source)`.
#' @param overwrite Logical. If `TRUE`, an existing destination item is removed
#'   before it is copied. If `FALSE` (the default), an existing destination is
#'   retained and that row returns `copied = FALSE`.
#'
#' @return Invisibly returns a data frame with one row per source and columns:
#'
#' * `label`: operation label;
#' * `source`: supplied source path;
#' * `destination`: complete final destination path;
#' * `type`: automatically detected `"file"` or `"directory"`, or `NA` for a
#'   missing source;
#' * `exists`: whether the source exists;
#' * `copied`: whether the item was copied successfully.
#'
#' @details
#' `destination_dir` always means a parent directory. For example, copying
#' `"results/barcode01"` to `"delivery"` creates
#' `"delivery/barcode01"`; copying `"results/report.xlsx"` to the same parent
#' creates `"delivery/report.xlsx"`.
#'
#' Missing destination parents are created recursively, but only for existing
#' sources. Before copying anything, the function rejects duplicated final
#' destinations, a destination parent occupied by a regular file, copying an
#' item onto itself, and recursively copying a directory into itself.
#'
#' Files and directory contents are copied with `copy.mode = TRUE` and
#' `copy.date = TRUE`. This preserves file permissions and modification times
#' when supported by the operating system and destination filesystem. Directory
#' timestamps and filesystem-specific metadata may not be preserved exactly.
#'
#' @examples
#' source_root <- tempfile("collect-source-")
#' destination_root <- tempfile("collect-destination-")
#' dir.create(file.path(source_root, "barcode01"), recursive = TRUE)
#' writeLines("report", file.path(source_root, "report.txt"))
#' writeLines(
#'   "alignment",
#'   file.path(source_root, "barcode01", "alignment.txt")
#' )
#'
#' result <- collect_items(
#'   source = c(
#'     file.path(source_root, "report.txt"),
#'     file.path(source_root, "barcode01")
#'   ),
#'   destination_dir = destination_root
#' )
#' result
#' file.exists(file.path(destination_root, "report.txt"))
#' dir.exists(file.path(destination_root, "barcode01"))
#'
#' # A different parent directory can be supplied for each source.
#' collect_items(
#'   source = c(
#'     file.path(source_root, "report.txt"),
#'     file.path(source_root, "barcode01")
#'   ),
#'   destination_dir = c(
#'     file.path(destination_root, "documents"),
#'     file.path(destination_root, "samples")
#'   ),
#'   overwrite = TRUE
#' )
#'
#' @export
collect_items <- function(source,
                          destination_dir,
                          label = NULL,
                          overwrite = FALSE) {
  source <- validate_collect_items_paths(source, "source")
  n <- length(source)
  destination_dir <- recycle_collect_items_argument(
    destination_dir,
    n = n,
    name = "destination_dir"
  )
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  source_without_separator <- sub("[/\\\\]+$", "", source)
  source_name <- basename(source_without_separator)
  if (any(!nzchar(source_name) | source_name %in% c(".", ".."))) {
    stop(
      "Every `source` must identify an item with a valid basename.",
      call. = FALSE
    )
  }
  if (is.null(label)) {
    label <- source_name
  } else {
    label <- recycle_collect_items_argument(label, n = n, name = "label")
  }

  is_directory <- dir.exists(source)
  is_file <- file.exists(source) & !is_directory
  source_exists <- is_directory | is_file
  type <- rep(NA_character_, n)
  type[is_file] <- "file"
  type[is_directory] <- "directory"
  destination <- file.path(destination_dir, source_name)

  validate_collect_items_plan(
    source = source,
    destination_dir = destination_dir,
    destination = destination,
    source_exists = source_exists,
    is_directory = is_directory
  )

  copied <- rep(FALSE, n)
  for (i in seq_len(n)) {
    if (!source_exists[[i]]) {
      next
    }

    dir.create(
      destination_dir[[i]],
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!dir.exists(destination_dir[[i]])) {
      stop(
        "Failed to create destination directory: ",
        destination_dir[[i]],
        call. = FALSE
      )
    }

    destination_exists <- file.exists(destination[[i]]) ||
      dir.exists(destination[[i]])
    if (destination_exists && !isTRUE(overwrite)) {
      next
    }
    if (destination_exists) {
      unlink(destination[[i]], recursive = TRUE, force = TRUE)
      if (file.exists(destination[[i]]) || dir.exists(destination[[i]])) {
        stop(
          "Failed to remove existing destination: ",
          destination[[i]],
          call. = FALSE
        )
      }
    }

    copy_result <- if (is_directory[[i]]) {
      file.copy(
        from = source[[i]],
        to = destination_dir[[i]],
        overwrite = FALSE,
        recursive = TRUE,
        copy.mode = TRUE,
        copy.date = TRUE
      )
    } else {
      file.copy(
        from = source[[i]],
        to = destination[[i]],
        overwrite = FALSE,
        recursive = FALSE,
        copy.mode = TRUE,
        copy.date = TRUE
      )
    }

    destination_created <- if (is_directory[[i]]) {
      dir.exists(destination[[i]])
    } else {
      file.exists(destination[[i]]) && !dir.exists(destination[[i]])
    }
    copied[[i]] <- isTRUE(copy_result) && destination_created
  }

  invisible(data.frame(
    label = label,
    source = source,
    destination = destination,
    type = type,
    exists = source_exists,
    copied = copied,
    stringsAsFactors = FALSE
  ))
}

validate_collect_items_paths <- function(x, name) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) ||
      any(!nzchar(x))) {
    stop(
      "`", name,
      "` must be a non-empty character vector without missing or empty values.",
      call. = FALSE
    )
  }
  x
}

recycle_collect_items_argument <- function(x, n, name) {
  x <- validate_collect_items_paths(x, name)
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(
      "`", name, "` must have length 1 or length(source) (", n, ").",
      call. = FALSE
    )
  }
  x
}

validate_collect_items_plan <- function(source,
                                        destination_dir,
                                        destination,
                                        source_exists,
                                        is_directory) {
  active <- which(source_exists)
  if (length(active) == 0L) {
    return(invisible(NULL))
  }

  parent_is_file <- file.exists(destination_dir[active]) &
    !dir.exists(destination_dir[active])
  if (any(parent_is_file)) {
    stop(
      "Destination parent exists as a regular file: ",
      paste(unique(destination_dir[active][parent_is_file]), collapse = ", "),
      call. = FALSE
    )
  }

  destination_key <- vapply(
    destination[active],
    canonical_collect_items_path,
    character(1)
  )
  duplicated_destination <- duplicated(destination_key) |
    duplicated(destination_key, fromLast = TRUE)
  if (any(duplicated_destination)) {
    stop(
      "Multiple sources resolve to the same destination: ",
      paste(unique(destination[active][duplicated_destination]), collapse = ", "),
      call. = FALSE
    )
  }

  source_key <- normalizePath(source[active], mustWork = TRUE)
  same_path <- source_key == destination_key
  if (any(same_path)) {
    stop(
      "A source and its destination must not be the same path: ",
      paste(source[active][same_path], collapse = ", "),
      call. = FALSE
    )
  }

  directory_rows <- which(is_directory[active])
  if (length(directory_rows) > 0L) {
    nested <- vapply(
      directory_rows,
      function(i) {
        startsWith(
          destination_key[[i]],
          paste0(source_key[[i]], .Platform$file.sep)
        )
      },
      logical(1)
    )
    if (any(nested)) {
      bad <- active[directory_rows[nested]]
      stop(
        "A directory cannot be copied into itself: ",
        paste(source[bad], collapse = ", "),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

canonical_collect_items_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(getwd(), path)
  }

  missing_components <- character()
  ancestor <- path
  while (!file.exists(ancestor) && !dir.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) {
      break
    }
    missing_components <- c(basename(ancestor), missing_components)
    ancestor <- parent
  }

  if (file.exists(ancestor) || dir.exists(ancestor)) {
    ancestor <- normalizePath(ancestor, mustWork = TRUE)
  }
  if (length(missing_components) == 0L) {
    return(ancestor)
  }
  do.call(file.path, as.list(c(ancestor, missing_components)))
}
