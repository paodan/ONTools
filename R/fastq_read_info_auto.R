#' Automatically select a FASTQ-header metadata parser
#'
#' `fastq_read_info_auto()` inspects the first record in the first selected
#' FASTQ file and dispatches to [fastq_read_info()] or [fastq_read_info2()]. It
#' is intended for workflows that may receive either legacy Oxford Nanopore
#' `key=value` headers or recent MinKNOW/Dorado `TAG:TYPE:VALUE` headers.
#'
#' @param path Optional directory containing FASTQ files. When `fastqs` is
#'   `NULL`, `path` is required and is scanned non-recursively. When `fastqs` is
#'   supplied, relative entries are resolved relative to `path`.
#' @param file_format Literal, case-sensitive suffix used for directory
#'   scanning when `fastqs` is `NULL`.
#' @param fastqs Optional non-empty character vector of explicit `.fastq`,
#'   `.fastq.gz`, `.fq`, or `.fq.gz` files. Input order is preserved, and the
#'   first element is used for format detection.
#' @param plot_time Logical. Passed to the selected parser. If `TRUE`, draw
#'   cumulative read-count plots against `start_time`.
#' @param barcode_column Optional character scalar naming the column used to
#'   group start-time plots. When `NULL`, the selected parser's default is used:
#'   `"barcode"` for [fastq_read_info()] and `"source_file"` for
#'   [fastq_read_info2()].
#' @param verbose Logical. If `TRUE`, report the detected header format and the
#'   selected parser.
#'
#' @return The data frame returned by the selected parser. Two attributes are
#' added:
#'
#' * `header_format`: `"legacy_key_value"` or `"sam_tag"`;
#' * `parser`: `"fastq_read_info"` or `"fastq_read_info2"`.
#'
#' @details
#' Detection reads only the first four-line FASTQ record from the first selected
#' file. A whitespace-delimited field such as `runid=abc` identifies the legacy
#' format. A SAM auxiliary field such as `ch:i:42` or `st:Z:2026-01-01T00:00:00Z`
#' identifies the MinKNOW/Dorado format. Detection stops with an error when the
#' first record contains both formats or neither format.
#'
#' The selected parser then processes every selected file normally. This
#' function assumes all inputs use the same header convention. In particular,
#' [fastq_read_info2()] validates that every read contains SAM-style metadata.
#' For mixed-format collections, split the files by header format and call the
#' two parsers separately.
#'
#' An empty first file cannot be classified and produces an error. The first
#' record must use the standard four-line FASTQ representation; wrapped FASTQ
#' records are not supported.
#'
#' @examples
#' fastq_dir <- tempfile("fastq-auto-")
#' dir.create(fastq_dir)
#'
#' legacy <- file.path(fastq_dir, "legacy.fastq")
#' writeLines(
#'   c(
#'     "@legacy-read ch=42 barcode=barcode01",
#'     "ACGT", "+", "!!!!"
#'   ),
#'   legacy
#' )
#' legacy_info <- fastq_read_info_auto(
#'   fastqs = legacy,
#'   plot_time = FALSE
#' )
#' attr(legacy_info, "parser")
#'
#' dorado <- file.path(fastq_dir, "dorado.fastq")
#' writeLines(
#'   c(
#'     "@dorado-read qs:f:15.7 ch:i:42 SM:Z:barcode01",
#'     "TGCA", "+", "####"
#'   ),
#'   dorado
#' )
#' dorado_info <- fastq_read_info_auto(
#'   fastqs = dorado,
#'   plot_time = FALSE
#' )
#' attr(dorado_info, "header_format")
#' dorado_info[, c("read", "mean_qscore", "ch", "barcode")]
#'
#' @seealso [fastq_read_info()], [fastq_read_info2()]
#' @export
fastq_read_info_auto <- function(path = NULL,
                                 file_format = ".fastq.gz",
                                 fastqs = NULL,
                                 plot_time = TRUE,
                                 barcode_column = NULL,
                                 verbose = FALSE) {
  check_scalar_character(file_format, "file_format")
  check_logical_scalar(plot_time, "plot_time")
  check_logical_scalar(verbose, "verbose")
  if (!is.null(barcode_column)) {
    check_scalar_character(barcode_column, "barcode_column")
  }

  files <- resolve_fastq_info_files(path, file_format, fastqs)
  if (length(files) == 0L) {
    stop("No FASTQ files were selected; header format cannot be detected.",
         call. = FALSE)
  }

  header <- read_first_fastq_header(files[[1L]])
  header_format <- detect_fastq_header_style(header)
  parser <- if (identical(header_format, "sam_tag")) {
    "fastq_read_info2"
  } else {
    "fastq_read_info"
  }
  if (isTRUE(verbose)) {
    message("Detected `", header_format, "`; using `", parser, "()`." )
  }

  if (identical(header_format, "sam_tag")) {
    if (is.null(barcode_column)) barcode_column <- "source_file"
    out <- fastq_read_info2(
      fastqs = files,
      plot_time = plot_time,
      barcode_column = barcode_column
    )
  } else {
    if (is.null(barcode_column)) barcode_column <- "barcode"
    out <- fastq_read_info(
      fastqs = files,
      plot_time = plot_time,
      barcode_column = barcode_column
    )
  }

  attr(out, "header_format") <- header_format
  attr(out, "parser") <- parser
  out
}

read_first_fastq_header <- function(path) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  record <- readLines(con, n = 4L, warn = FALSE)
  if (length(record) == 0L) {
    stop(
      "The first FASTQ file is empty; header format cannot be detected: ",
      path,
      call. = FALSE
    )
  }
  if (length(record) != 4L || !startsWith(record[[1L]], "@") ||
      !startsWith(record[[3L]], "+")) {
    stop("Invalid first four-line FASTQ record in file: ", path,
         call. = FALSE)
  }
  if (nchar(record[[2L]], type = "bytes") !=
      nchar(record[[4L]], type = "bytes")) {
    stop(
      "Sequence and quality lengths differ in the first FASTQ record: ",
      path,
      call. = FALSE
    )
  }
  substring(record[[1L]], 2L)
}

detect_fastq_header_style <- function(header) {
  sam_style <- grepl(
    "(^|[[:space:]])[[:alnum:]]{2}:[AifZHBcCsSI]:[^[:space:]]+",
    header,
    perl = TRUE
  )
  legacy_style <- grepl(
    "(^|[[:space:]])[^[:space:]=]+=[^[:space:]]+",
    header,
    perl = TRUE
  )

  if (sam_style && legacy_style) {
    stop(
      "The first FASTQ header contains both `key=value` and ",
      "`TAG:TYPE:VALUE` metadata; parser selection is ambiguous.",
      call. = FALSE
    )
  }
  if (sam_style) return("sam_tag")
  if (legacy_style) return("legacy_key_value")
  stop(
    "The first FASTQ header contains neither recognizable `key=value` nor ",
    "`TAG:TYPE:VALUE` metadata.",
    call. = FALSE
  )
}
