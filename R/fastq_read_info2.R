#' Extract metadata from MinKNOW/Dorado FASTQ headers
#'
#' `fastq_read_info2()` is a vectorized parser for recent MinKNOW/Dorado FASTQ
#' headers that use SAM auxiliary-field syntax (`TAG:TYPE:VALUE`). Unlike
#' [fastq_read_info()], it does not parse legacy `key=value` metadata.
#'
#' @param path Directory containing FASTQ files.
#' @param file_format Literal, case-sensitive file-name suffix used to select
#'   files, for example `".fastq.gz"`, `".fastq"`, or `".fq.gz"`.
#' @param plot_time Logical. If `TRUE`, draw cumulative read count against
#'   `start_time`, with one base-R plot per standardized `barcode` value. No
#'   plot is produced when valid start times are absent.
#'
#' @return A data frame with one row per FASTQ record. `read` contains the UUID
#'   or identifier following `@`, `source_file` contains the input file name,
#'   and `header_format` is `"sam_tag"`. MinKNOW/Dorado tags are returned under
#'   the standardized names described below. Integer and floating-point fields
#'   are converted to the corresponding R types, and time fields are returned
#'   as `POSIXct` in UTC.
#'
#' Optional fields that are absent from every header are retained as typed
#' `NA` columns, giving different runs a consistent output schema. If no files
#' match `file_format`, an empty data frame with the core columns is returned.
#'
#' @section Expected header format:
#'
#' The first whitespace-delimited token is the read identifier. Remaining
#' tokens must use `TAG:TYPE:VALUE`, for example:
#'
#' ```
#' "@read-id qs:f:15.7 ch:i:1803 st:Z:2026-09-20T14:53:20+08:00 SM:Z:barcode86"
#' ```
#'
#' In these fields, `i` denotes an integer, `f` a floating-point number, and
#' `Z` a string. Field order may vary and optional fields may be absent. If any
#' read lacks SAM-style metadata, the function stops instead of silently
#' returning missing values. Use [fastq_read_info()] for legacy `key=value`
#' headers.
#'
#' @section Standardized fields and their origins:
#'
#' * `mean_qscore`: `qs:f`; mean basecall Q-score for the read.
#' * `mux`: `mx:i`; mux/well number used by the channel.
#' * `ch`: `ch:i`; flow-cell channel number.
#' * `read_number`: `rn:i`; acquisition read number assigned by MinKNOW.
#' * `start_time`: `st:Z`; start time of the read.
#' * `signal_start` and `signal_end`: `ts:i` and `ns:i`. The basecalled
#'   sequence corresponds to `signal[ts:ns]`. `signal_samples` is calculated as
#'   `signal_end - signal_start`; these are signal samples, not DNA bases.
#' * `duration`: `du:f`; read duration in seconds.
#' * `scaling_midpoint`, `scaling_dispersion`, and `scaling_version`: `sm:f`,
#'   `sd:f`, and `sv:Z`; raw-current scaling and normalization parameters.
#' * `duplex`: `dx:i`; normally zero for simplex and one for duplex reads.
#' * `read_group`: `RG:Z`; normally combines the run ID, basecalling model, and
#'   barcode arrangement.
#' * `experiment_start_time`: `DT:Z`; experiment start time.
#' * `flow_cell_id`: `PU:Z`; flow-cell identifier.
#' * `device_id`: `PM:Z`; sequencing device identifier, when supplied.
#' * `sample_id`: `LB:Z` in the ONT read-group convention.
#' * `barcode`: `SM:Z`; classified barcode name. When `SM` is absent,
#'   `barcode_alias` is used as a fallback.
#' * `barcode_alias`: `al:Z`; alias defined in the sample sheet.
#' * `barcode_kit` and `trimming`: `bk:Z` and `tm:Z`; barcode kit and requested
#'   adapter/primer/barcode trimming configuration.
#' * `source_signal_file`: `fn:Z`; original signal file name.
#' * `parent_read_id` and `parent_signal_start`: `pi:Z` and `sp:i` for reads
#'   created by signal splitting.
#' * `bed_hits`: `bh:i`; number of detected BED-file hits.
#' * `minknow_events`: `me:i`; number of MinKNOW events detected during
#'   sequencing.
#' * `pore_type`: `po:Z`; detected pore type.
#' * `end_reason`: `er:Z`; reason the read ended.
#' * `barcode_variant`: `bv:Z`; detected barcode-arrangement variant.
#' * `polya_length`: `pt:i`; estimated poly(A/T) length when enabled.
#' * `read_group_description`: `DS:Z`; read-group description when copied into
#'   the FASTQ header.
#'
#' Tags are case-sensitive: `sm:f` is the signal-scaling midpoint, whereas
#' `SM:Z` is the barcode/sample name. Tags not listed above are not returned.
#' Definitions follow the
#' [Dorado SAM specification](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/).
#'
#' @details
#' The parser extracts one complete output column at a time. It therefore
#' avoids calling an R parser for every read and avoids constructing a large
#' list containing one small list per read. This is substantially faster and
#' more memory-efficient when all input files use the same MinKNOW/Dorado
#' header convention.
#'
#' Files are processed in sorted file-name order. Both plain-text and
#' gzip-compressed FASTQ files are supported. A complete four-line FASTQ record
#' is required; malformed or incomplete records produce an error.
#'
#' @examples
#' fastq_dir <- tempfile("fastq2-")
#' dir.create(fastq_dir)
#' writeLines(
#'   c(
#'     paste0(
#'       "@new-read qs:f:15.715393 mx:i:4 ch:i:1803 rn:i:438 ",
#'       "st:Z:2026-09-20T14:53:20.912722+08:00 ts:i:10 ns:i:3694 ",
#'       "du:f:0.7388 sm:f:93.99999 sd:f:24 sv:Z:pa dx:i:0 ",
#'       "RG:Z:run_model_barcode86 DT:Z:2026-09-20T14:27:18+08:00 ",
#'       "PU:Z:PBO11046 LB:Z:20260920pcr SM:Z:barcode86 al:Z:barcode86"
#'     ),
#'     "ACGT", "+", "!!!!"
#'   ),
#'   file.path(fastq_dir, "reads.fastq")
#' )
#'
#' info <- fastq_read_info2(
#'   fastq_dir,
#'   file_format = ".fastq",
#'   plot_time = FALSE
#' )
#' info[, c(
#'   "read", "mean_qscore", "ch", "read_number", "start_time",
#'   "signal_samples", "duration", "flow_cell_id", "barcode"
#' )]
#'
#' @seealso [fastq_read_info()], [plot_fastq_read_distribution()]
#' @export
fastq_read_info2 <- function(path,
                             file_format = ".fastq.gz",
                             plot_time = TRUE) {
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
      source_file = character(),
      header_format = character(),
      start_time = as.POSIXct(character(), tz = "UTC"),
      stringsAsFactors = FALSE
    ))
  }

  header_list <- lapply(files, read_fastq_headers)
  headers <- unlist(header_list, use.names = FALSE)
  source_file <- rep(basename(files), lengths(header_list))
  out <- parse_minknow_fastq_headers(headers, source_file)

  if (isTRUE(plot_time)) plot_fastq_start_times(out)
  out
}

parse_minknow_fastq_headers <- function(headers, source_file = NULL) {
  if (!is.character(headers)) {
    stop("`headers` must be a character vector.", call. = FALSE)
  }
  n <- length(headers)
  if (is.null(source_file)) source_file <- rep(NA_character_, n)
  if (length(source_file) == 1L) source_file <- rep(source_file, n)
  if (length(source_file) != n) {
    stop(
      "`source_file` must have length 1 or the same length as `headers`.",
      call. = FALSE
    )
  }

  is_sam_style <- grepl(
    "(^|[[:space:]])[[:alnum:]]{2}:[AifZHBcCsSI]:[^[:space:]]+",
    headers,
    perl = TRUE
  )
  if (any(!is_sam_style)) {
    stop(
      sum(!is_sam_style),
      " FASTQ header(s) do not contain MinKNOW/Dorado SAM-style metadata. ",
      "Use `fastq_read_info()` for legacy `key=value` headers.",
      call. = FALSE
    )
  }

  field_spec <- minknow_fastq_field_spec()
  parsed <- Map(
    function(tag, sam_type, r_type) {
      extract_and_convert_minknow_tag(headers, tag, sam_type, r_type)
    },
    field_spec$tag,
    field_spec$sam_type,
    field_spec$r_type
  )
  names(parsed) <- field_spec$column

  missing_barcode <- is.na(parsed$barcode) | !nzchar(parsed$barcode)
  parsed$barcode[missing_barcode] <- parsed$barcode_alias[missing_barcode]
  signal_samples <- parsed$signal_end - parsed$signal_start
  signal_end_position <- match("signal_end", names(parsed))
  parsed <- append(
    parsed,
    list(signal_samples = signal_samples),
    after = signal_end_position
  )

  read_id <- sub("[[:space:]].*$", "", headers)
  read_id <- sub("^@", "", read_id)
  data.frame(
    read = read_id,
    source_file = source_file,
    header_format = rep("sam_tag", n),
    parsed,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

minknow_fastq_field_spec <- function() {
  data.frame(
    column = c(
      "mean_qscore", "mux", "ch", "read_number", "start_time",
      "signal_start", "signal_end", "duration", "scaling_midpoint",
      "scaling_dispersion", "scaling_version", "duplex", "read_group",
      "experiment_start_time", "flow_cell_id", "device_id", "sample_id",
      "barcode", "barcode_alias", "barcode_kit", "trimming",
      "source_signal_file", "parent_read_id", "parent_signal_start",
      "bed_hits", "minknow_events", "pore_type", "end_reason",
      "barcode_variant", "polya_length", "read_group_description"
    ),
    tag = c(
      "qs", "mx", "ch", "rn", "st", "ts", "ns", "du", "sm", "sd",
      "sv", "dx", "RG", "DT", "PU", "PM", "LB", "SM", "al", "bk",
      "tm", "fn", "pi", "sp", "bh", "me", "po", "er", "bv", "pt", "DS"
    ),
    sam_type = c(
      "f", "i", "i", "i", "Z", "i", "i", "f", "f", "f", "Z", "i",
      "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "i",
      "i", "i", "Z", "Z", "Z", "i", "Z"
    ),
    r_type = c(
      "double", "integer", "integer", "integer", "datetime", "integer",
      "integer", "double", "double", "double", "character", "integer",
      "character", "datetime", "character", "character", "character",
      "character", "character", "character", "character", "character",
      "character", "integer", "integer", "integer", "character",
      "character", "character", "integer", "character"
    ),
    stringsAsFactors = FALSE
  )
}

extract_minknow_tag <- function(headers, tag, sam_type) {
  pattern <- paste0(
    ".*[[:space:]]", tag, ":", sam_type, ":([^[:space:]]+).*"
  )
  value <- sub(pattern, "\\1", headers, perl = TRUE)
  value[value == headers] <- NA_character_
  value
}

extract_and_convert_minknow_tag <- function(headers,
                                            tag,
                                            sam_type,
                                            r_type) {
  value <- extract_minknow_tag(headers, tag, sam_type)
  switch(
    r_type,
    integer = suppressWarnings(as.integer(value)),
    double = suppressWarnings(as.numeric(value)),
    datetime = parse_fastq_start_time(value),
    character = value,
    stop("Unsupported R type: ", r_type, call. = FALSE)
  )
}
