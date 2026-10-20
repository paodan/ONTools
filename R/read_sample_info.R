#' Read a sectioned sample-information workbook
#'
#' `read_sample_info()` reads every selected worksheet in an Excel workbook and
#' splits it into named blocks. A block starts at a cell containing a marker of
#' the form `[name]` and ends immediately before the next marker row. The return
#' value mirrors the workbook structure: worksheets form the outer list and
#' marker names form the inner lists.
#'
#' Marker cells may occur anywhere in a worksheet. Only the marker cell itself
#' is removed, so values elsewhere on the marker row are retained. Completely
#' empty rows and outer empty columns are removed from each block. When
#' `drop_empty_rows = FALSE`, empty rows inside a block are preserved.
#'
#' Several common block names receive additional parsing:
#'
#' * `[header]` is returned with columns `field`, `description`, and `value`.
#' * `[columns]` is returned with columns `column` and `description`.
#' * `[data]` and `[calculator]` use their first non-empty row as the table
#'   header.
#' * Other blocks are returned as character data frames with names `X1`, `X2`,
#'   and so on.
#'
#' Empty table-header cells are retained as `.extra_<column number>`. Duplicate
#' table headers are made unique. A warning is emitted in both cases, preventing
#' data outside the named header range from being silently discarded.
#'
#' @param file Path to an existing `.xlsx`, `.xlsm`, `.xltx`, `.xltm`, or
#'   `.xls` workbook.
#' @param sheets Worksheets to read. `NULL` reads all worksheets. Otherwise,
#'   supply a character vector of worksheet names or a numeric vector of
#'   one-based worksheet positions.
#' @param strict Logical. If `TRUE` (the default), worksheets without markers,
#'   duplicate marker names, and non-empty content before the first marker are
#'   treated as errors. If `FALSE`, an unmarked worksheet is returned as an
#'   `unmarked` block, pre-marker content is returned as a `preamble` block, and
#'   duplicate block names are made unique. Multiple markers on the same row
#'   are always an error because their block boundaries are ambiguous.
#' @param drop_empty_rows Logical. If `TRUE`, completely empty rows are removed
#'   from blocks. If `FALSE`, internal empty rows are retained; empty rows at the
#'   outer edges are always trimmed.
#' @param include_images Logical. Excel images are not cell values and therefore
#'   cannot be represented in the parsed tables. If `FALSE` (the default),
#'   embedded images are detected and a warning is issued. If `TRUE`, embedded
#'   media files are extracted to a temporary directory and their paths are
#'   stored in the `images` attribute of the returned list. Images are not
#'   assigned to individual blocks.
#'
#' @return A named nested list. The outer names are worksheet names, and each
#'   outer element is a named list of block data frames. The result has a
#'   `source` attribute containing the normalized workbook path. If embedded
#'   images are extracted, the `images` attribute contains their temporary file
#'   paths and `image_dir` contains the extraction directory.
#'
#' @details
#' Cells are read as text to keep every block type consistent and to avoid
#' type-guessing differences between worksheets. Consequently, numeric values,
#' dates, and logical values are returned as character vectors. Formula cells
#' contain the cached values stored by Excel readers; formulas themselves are
#' not returned.
#'
#' The function recognizes block names case-insensitively for selecting a parser
#' but preserves the original marker text as the list-element name. Thus
#' `[Data]` is parsed like `[data]` and remains named `Data`.
#'
#' @examplesIf requireNamespace("openxlsx", quietly = TRUE)
#' workbook <- tempfile(fileext = ".xlsx")
#' wb <- openxlsx::createWorkbook()
#' openxlsx::addWorksheet(wb, "2_metadata")
#' example_cells <- matrix(NA_character_, nrow = 9, ncol = 3)
#' example_cells[1, ] <- c("[header]", NA, NA)
#' example_cells[2, ] <- c("Project_ID", "Project identifier", "P001")
#' example_cells[4, ] <- c("[columns]", NA, NA)
#' example_cells[5, ] <- c("Sample_ID", "Sample identifier", NA)
#' example_cells[7, ] <- c("[data]", NA, NA)
#' example_cells[8, ] <- c("Sample_ID", "Barcode_ID", "Expected_Size_bp")
#' example_cells[9, ] <- c("S1", "NB01", "1600bp")
#' openxlsx::writeData(wb, "2_metadata", example_cells,
#'                     colNames = FALSE, rowNames = FALSE)
#' openxlsx::saveWorkbook(wb, workbook, overwrite = TRUE)
#'
#' info <- read_sample_info(workbook)
#' names(info)
#' names(info[["2_metadata"]])
#' info[["2_metadata"]][["data"]]
#'
#' # Read only selected worksheets from a multi-sheet workbook:
#' # read_sample_info(workbook, sheets = "2_metadata")
#'
#' @export
read_sample_info <- function(file,
                             sheets = NULL,
                             strict = TRUE,
                             drop_empty_rows = TRUE,
                             include_images = FALSE) {
  check_scalar_character(file, "file")
  check_logical_scalar(strict, "strict")
  check_logical_scalar(drop_empty_rows, "drop_empty_rows")
  check_logical_scalar(include_images, "include_images")

  if (!file.exists(file)) {
    stop("`file` does not exist: ", file, call. = FALSE)
  }
  if (dir.exists(file)) {
    stop("`file` must be an Excel workbook, not a directory.", call. = FALSE)
  }

  extension <- tolower(tools::file_ext(file))
  supported <- c("xlsx", "xlsm", "xltx", "xltm", "xls")
  if (!extension %in% supported) {
    stop(
      "`file` must have one of these extensions: ",
      paste0(".", supported, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  workbook_sheets <- readxl::excel_sheets(file)
  selected_sheets <- select_sample_info_sheets(sheets, workbook_sheets)
  result <- stats::setNames(
    vector("list", length(selected_sheets)),
    selected_sheets
  )

  for (sheet in selected_sheets) {
    cells <- suppressMessages(readxl::read_excel(
      path = file,
      sheet = sheet,
      col_names = FALSE,
      col_types = "text",
      trim_ws = FALSE,
      .name_repair = "minimal"
    ))
    cells <- as.data.frame(cells, stringsAsFactors = FALSE, check.names = FALSE)
    cells[] <- lapply(cells, as.character)
    result[[sheet]] <- parse_sample_info_sheet(
      cells,
      sheet = sheet,
      strict = strict,
      drop_empty_rows = drop_empty_rows
    )
  }

  source <- normalizePath(file, mustWork = TRUE)
  media <- find_sample_info_media(source)
  if (length(media) > 0L) {
    if (include_images) {
      image_dir <- tempfile("sample-info-images-")
      dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
      suppressWarnings(utils::unzip(source, files = media, exdir = image_dir))
      image_paths <- file.path(image_dir, media)
      image_paths <- normalizePath(image_paths[file.exists(image_paths)],
                                   mustWork = TRUE)
      attr(result, "images") <- image_paths
      attr(result, "image_dir") <- normalizePath(image_dir, mustWork = TRUE)
    } else {
      warning(
        "Workbook contains ", length(media),
        " embedded image file(s). Cell blocks were read, but images were not ",
        "extracted; use `include_images = TRUE` to extract them.",
        call. = FALSE
      )
    }
  }
  attr(result, "source") <- source
  result
}

select_sample_info_sheets <- function(sheets, workbook_sheets) {
  if (is.null(sheets)) {
    return(workbook_sheets)
  }
  if (length(sheets) < 1L || anyNA(sheets)) {
    stop("`sheets` must be `NULL` or a non-empty vector without missing values.",
         call. = FALSE)
  }

  if (is.numeric(sheets)) {
    if (any(!is.finite(sheets)) || any(sheets != as.integer(sheets)) ||
        any(sheets < 1L) || any(sheets > length(workbook_sheets))) {
      stop("Numeric `sheets` values must be valid one-based worksheet positions.",
           call. = FALSE)
    }
    selected <- workbook_sheets[as.integer(sheets)]
  } else if (is.character(sheets)) {
    if (any(!nzchar(sheets))) {
      stop("Character `sheets` values must not be empty.", call. = FALSE)
    }
    missing <- setdiff(sheets, workbook_sheets)
    if (length(missing) > 0L) {
      stop(
        "Workbook does not contain worksheet(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    selected <- sheets
  } else {
    stop("`sheets` must be `NULL`, character, or numeric.", call. = FALSE)
  }

  if (anyDuplicated(selected)) {
    stop("`sheets` must not select the same worksheet more than once.",
         call. = FALSE)
  }
  selected
}

parse_sample_info_sheet <- function(cells, sheet, strict, drop_empty_rows) {
  values <- as.matrix(cells)
  storage.mode(values) <- "character"
  marker_pattern <- "^\\s*\\[([^][]+)\\]\\s*$"
  is_marker <- matrix(
    grepl(marker_pattern, values, perl = TRUE),
    nrow = nrow(values),
    ncol = ncol(values)
  )
  marker_locations <- which(is_marker, arr.ind = TRUE)

  if (nrow(marker_locations) == 0L) {
    if (strict) {
      stop("Worksheet `", sheet, "` contains no [name] block markers.",
           call. = FALSE)
    }
    warning("Worksheet `", sheet,
            "` contains no block markers; returning one `unmarked` block.",
            call. = FALSE)
    return(list(unmarked = sample_info_raw_block(
      trim_sample_info_matrix(values, drop_empty_rows = drop_empty_rows)
    )))
  }

  marker_locations <- marker_locations[
    order(marker_locations[, "row"], marker_locations[, "col"]),
    ,
    drop = FALSE
  ]
  if (anyDuplicated(marker_locations[, "row"])) {
    bad_rows <- unique(marker_locations[, "row"][
      duplicated(marker_locations[, "row"]) |
        duplicated(marker_locations[, "row"], fromLast = TRUE)
    ])
    stop(
      "Worksheet `", sheet, "` has multiple block markers on row(s): ",
      paste(bad_rows, collapse = ", "),
      ". Block boundaries would be ambiguous.",
      call. = FALSE
    )
  }

  marker_cells <- values[marker_locations]
  marker_names <- trimws(sub(marker_pattern, "\\1", marker_cells, perl = TRUE))
  duplicate_names <- unique(marker_names[duplicated(marker_names)])
  if (length(duplicate_names) > 0L && strict) {
    stop(
      "Worksheet `", sheet, "` contains duplicate block marker name(s): ",
      paste(duplicate_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (length(duplicate_names) > 0L) {
    warning(
      "Worksheet `", sheet,
      "` contains duplicate block names; unique suffixes were added.",
      call. = FALSE
    )
    marker_names <- make.unique(marker_names, sep = "_")
  }

  blocks <- list()
  first_marker_row <- marker_locations[1L, "row"]
  if (first_marker_row > 1L) {
    preamble <- values[seq_len(first_marker_row - 1L), , drop = FALSE]
    if (any(!sample_info_empty_cell(preamble))) {
      if (strict) {
        stop(
          "Worksheet `", sheet,
          "` contains non-empty content before its first block marker.",
          call. = FALSE
        )
      }
      warning(
        "Worksheet `", sheet,
        "` contains content before its first marker; returning it as `preamble`.",
        call. = FALSE
      )
      blocks[["preamble"]] <- sample_info_raw_block(
        trim_sample_info_matrix(preamble, drop_empty_rows = drop_empty_rows)
      )
    }
  }

  for (index in seq_len(nrow(marker_locations))) {
    start <- marker_locations[index, "row"]
    end <- if (index < nrow(marker_locations)) {
      marker_locations[index + 1L, "row"] - 1L
    } else {
      nrow(values)
    }
    block <- values[start:end, , drop = FALSE]
    block[1L, marker_locations[index, "col"]] <- NA_character_
    block <- trim_sample_info_matrix(block, drop_empty_rows = drop_empty_rows)
    blocks[[marker_names[[index]]]] <- parse_sample_info_block(
      block,
      block_name = marker_names[[index]],
      sheet = sheet
    )
  }
  blocks
}

parse_sample_info_block <- function(block, block_name, sheet) {
  parser <- tolower(block_name)
  if (parser == "header") {
    return(sample_info_fixed_block(
      block,
      c("field", "description", "value")
    ))
  }
  if (parser == "columns") {
    return(sample_info_fixed_block(block, c("column", "description")))
  }
  if (parser %in% c("data", "calculator")) {
    return(sample_info_table_block(block, block_name = block_name, sheet = sheet))
  }
  sample_info_raw_block(block)
}

sample_info_table_block <- function(block, block_name, sheet) {
  if (nrow(block) == 0L || ncol(block) == 0L) {
    return(data.frame())
  }

  header <- trimws(block[1L, ])
  blank_header <- sample_info_empty_cell(header)
  if (any(blank_header)) {
    header[blank_header] <- paste0(".extra_", which(blank_header))
    warning(
      "Block `", block_name, "` in worksheet `", sheet,
      "` has data column(s) without a header: ",
      paste(which(blank_header), collapse = ", "),
      ". Generated `.extra_<column>` names were used.",
      call. = FALSE
    )
  }
  if (anyDuplicated(header)) {
    duplicate_header <- unique(header[duplicated(header)])
    warning(
      "Block `", block_name, "` in worksheet `", sheet,
      "` has duplicate header(s): ", paste(duplicate_header, collapse = ", "),
      ". Unique suffixes were added.",
      call. = FALSE
    )
    header <- make.unique(header, sep = "_")
  }

  if (nrow(block) == 1L) {
    output <- as.data.frame(matrix(character(), nrow = 0L, ncol = ncol(block)),
                            stringsAsFactors = FALSE)
  } else {
    output <- as.data.frame(block[-1L, , drop = FALSE],
                            stringsAsFactors = FALSE, check.names = FALSE)
  }
  names(output) <- header
  rownames(output) <- NULL
  output
}

sample_info_fixed_block <- function(block, standard_names) {
  if (ncol(block) == 0L) {
    output <- as.data.frame(matrix(character(), nrow = nrow(block), ncol = 0L))
    return(output)
  }
  output <- as.data.frame(block, stringsAsFactors = FALSE, check.names = FALSE)
  output_names <- if (ncol(block) <= length(standard_names)) {
    standard_names[seq_len(ncol(block))]
  } else {
    c(standard_names,
      paste0(".extra_", seq.int(length(standard_names) + 1L, ncol(block))))
  }
  names(output) <- output_names
  rownames(output) <- NULL
  output
}

sample_info_raw_block <- function(block) {
  if (ncol(block) == 0L) {
    return(as.data.frame(matrix(character(), nrow = nrow(block), ncol = 0L)))
  }
  output <- as.data.frame(block, stringsAsFactors = FALSE, check.names = FALSE)
  names(output) <- paste0("X", seq_len(ncol(output)))
  rownames(output) <- NULL
  output
}

trim_sample_info_matrix <- function(x, drop_empty_rows) {
  if (length(x) == 0L || nrow(x) == 0L || ncol(x) == 0L) {
    return(matrix(character(), nrow = 0L, ncol = 0L))
  }
  nonempty_rows <- which(rowSums(!sample_info_empty_cell(x)) > 0L)
  nonempty_cols <- which(colSums(!sample_info_empty_cell(x)) > 0L)
  if (length(nonempty_rows) == 0L || length(nonempty_cols) == 0L) {
    return(matrix(character(), nrow = 0L, ncol = 0L))
  }

  x <- x[seq.int(min(nonempty_rows), max(nonempty_rows)),
         seq.int(min(nonempty_cols), max(nonempty_cols)),
         drop = FALSE]
  if (drop_empty_rows) {
    keep <- rowSums(!sample_info_empty_cell(x)) > 0L
    x <- x[keep, , drop = FALSE]
  }
  x
}

sample_info_empty_cell <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

find_sample_info_media <- function(file) {
  if (!tolower(tools::file_ext(file)) %in% c("xlsx", "xlsm", "xltx", "xltm")) {
    return(character())
  }
  listing <- tryCatch(
    suppressWarnings(utils::unzip(file, list = TRUE)),
    error = function(error) NULL
  )
  if (is.null(listing) || !"Name" %in% names(listing)) {
    return(character())
  }
  listing$Name[grepl("^xl/media/[^/]+$", listing$Name)]
}
