#' Extract Oxford Nanopore read metadata from FASTQ headers
#'
#' `fastq_read_info()` reads the header of every record in matching FASTQ files
#' and expands the space-separated `key=value` metadata into columns. Both
#' uncompressed FASTQ files and gzip-compressed FASTQ files are supported.
#'
#' @param path Directory containing FASTQ files.
#' @param file_format File-name suffix used to select input files. Matching is
#'   literal and case-sensitive.
#' @param plot_time Logical. If `TRUE`, plot cumulative read count against
#'   `start_time`, with one plot per barcode when a `barcode` field is present.
#' @param ... Deprecated compatibility arguments `fileFormat` and `plotTime`.
#'
#' @return A data frame with one row per FASTQ record. The first column, `read`,
#'   contains the read identifier without the leading `@`; remaining columns
#'   contain metadata parsed from the header. `start_time`, when present, is
#'   converted to `POSIXct` in UTC. If no files match, an empty data frame is
#'   returned.
#'
#' @examples
#' fastq_dir <- tempfile("fastq-")
#' dir.create(fastq_dir)
#' writeLines(
#'   c(
#'     "@read1 start_time=2025-01-01T00:00:00Z barcode=barcode01",
#'     "ACGT", "+", "!!!!"
#'   ),
#'   file.path(fastq_dir, "reads.fastq")
#' )
#' fastq_read_info(fastq_dir, file_format = ".fastq", plot_time = FALSE)
#'
#' @export
fastq_read_info <- function(path,
                            file_format = ".fastq.gz",
                            plot_time = TRUE,
                            ...) {
  file_format_missing <- missing(file_format)
  plot_time_missing <- missing(plot_time)
  dots <- list(...)

  unknown <- setdiff(names(dots), c("fileFormat", "plotTime"))
  if (length(dots) > 0L && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All arguments in `...` must be named.", call. = FALSE)
  }
  if (length(unknown) > 0L) {
    stop("Unused argument `", unknown[[1L]], "`.", call. = FALSE)
  }
  if (!is.null(dots$fileFormat)) {
    if (!file_format_missing) {
      stop("Supply only one of `file_format` and `fileFormat`.", call. = FALSE)
    }
    file_format <- dots$fileFormat
  }
  if (!is.null(dots$plotTime)) {
    if (!plot_time_missing) {
      stop("Supply only one of `plot_time` and `plotTime`.", call. = FALSE)
    }
    plot_time <- dots$plotTime
  }

  check_scalar_character(path, "path")
  check_scalar_character(file_format, "file_format")
  check_logical_scalar(plot_time, "plot_time")
  if (!dir.exists(path)) {
    stop("`path` is not an existing directory: ", path, call. = FALSE)
  }

  files <- list.files(path, full.names = TRUE, recursive = FALSE)
  files <- sort(files[endsWith(basename(files), file_format)])
  if (length(files) == 0L) {
    return(data.frame(
      read = character(),
      start_time = as.POSIXct(character(), tz = "UTC"),
      stringsAsFactors = FALSE
    ))
  }

  headers <- unlist(lapply(files, read_fastq_headers), use.names = FALSE)
  out <- parse_fastq_headers(headers)

  if ("start_time" %in% names(out)) {
    original_start_time <- out$start_time
    out$start_time <- parse_fastq_start_time(original_start_time)
    failed <- !is.na(original_start_time) & nzchar(original_start_time) &
      is.na(out$start_time)
    if (any(failed)) {
      warning(
        "Could not parse ", sum(failed), " `start_time` value(s); returned NA.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(plot_time)) {
    plot_fastq_start_times(out)
  }

  out
}

read_fastq_headers <- function(path, chunk_records = 10000L) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  headers <- character()
  chunk_lines <- 4L * chunk_records
  repeat {
    block <- readLines(con, n = chunk_lines, warn = FALSE)
    if (length(block) == 0L) break
    if (length(block) %% 4L != 0L) {
      stop("Incomplete FASTQ record in file: ", path, call. = FALSE)
    }

    header_index <- seq.int(1L, length(block), by = 4L)
    plus_index <- header_index + 2L
    if (any(!startsWith(block[header_index], "@")) ||
        any(!startsWith(block[plus_index], "+"))) {
      stop("Invalid four-line FASTQ record in file: ", path, call. = FALSE)
    }
    headers <- c(headers, substring(block[header_index], 2L))
  }

  headers
}

parse_fastq_headers <- function(headers) {
  if (length(headers) == 0L) {
    return(data.frame(read = character(), stringsAsFactors = FALSE))
  }

  tokens <- strsplit(trimws(headers), "[[:space:]]+", perl = TRUE)
  read_id <- vapply(tokens, `[[`, character(1), 1L)
  metadata <- lapply(tokens, function(x) {
    x <- x[-1L]
    separator <- regexpr("=", x, fixed = TRUE)
    x <- x[separator > 1L]
    separator <- separator[separator > 1L]
    values <- substring(x, separator + 1L)
    names(values) <- substring(x, 1L, separator - 1L)
    values
  })
  keys <- unique(unlist(lapply(metadata, names), use.names = FALSE))

  out <- data.frame(read = read_id, stringsAsFactors = FALSE)
  for (key in keys) {
    out[[key]] <- vapply(metadata, function(x) {
      value <- x[key]
      if (length(value) == 0L) NA_character_ else unname(value[[1L]])
    }, character(1))
  }
  out
}

parse_fastq_start_time <- function(x) {
  normalized <- sub("Z$", "", x)
  has_offset <- grepl("[+-][0-9]{2}:[0-9]{2}$", normalized)
  normalized[has_offset] <- sub(
    "([+-][0-9]{2}):([0-9]{2})$",
    "\\1\\2",
    normalized[has_offset]
  )

  out <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = "UTC")
  no_offset <- !has_offset & !is.na(normalized) & nzchar(normalized)
  out[no_offset] <- as.POSIXct(
    normalized[no_offset],
    format = "%Y-%m-%dT%H:%M:%OS",
    tz = "UTC"
  )
  out[has_offset] <- as.POSIXct(
    normalized[has_offset],
    format = "%Y-%m-%dT%H:%M:%OS%z",
    tz = "UTC"
  )
  out
}

plot_fastq_start_times <- function(info) {
  if (!"start_time" %in% names(info) || nrow(info) == 0L) return(invisible(NULL))

  group <- if ("barcode" %in% names(info)) info$barcode else "all_reads"
  group[is.na(group) | !nzchar(group)] <- "unclassified"
  valid <- !is.na(info$start_time)
  if (!any(valid)) {
    warning("No valid `start_time` values to plot.", call. = FALSE)
    return(invisible(NULL))
  }

  for (label in unique(group[valid])) {
    times <- sort(info$start_time[valid & group == label])
    graphics::plot(
      times,
      seq_along(times),
      main = label,
      xlab = "Start time",
      ylab = "Cumulative reads"
    )
  }
  invisible(NULL)
}
