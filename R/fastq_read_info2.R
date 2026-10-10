#' Extract metadata from legacy and SAM-tag-style ONT FASTQ headers
#'
#' `fastq_read_info2()` reads Oxford Nanopore FASTQ headers and automatically
#' recognizes both the legacy `key=value` convention used by many
#' MinKNOW/Guppy outputs and the `TAG:TYPE:VALUE` convention used by recent
#' MinKNOW/Dorado outputs. Source fields from both conventions are translated
#' to a common set of descriptive column names.
#'
#' @param path Directory containing the FASTQ files.
#' @param file_format Literal, case-sensitive file-name suffix used to select
#'   files, for example `".fastq.gz"`, `".fastq"`, or `".fq.gz"`.
#' @param plot_time Logical. If `TRUE`, draw cumulative read count against
#'   `start_time`, using one base-R plot per value of the standardized
#'   `barcode` column. No plot is produced when valid start times are absent.
#'
#' @return A data frame with one row per FASTQ record. The columns `read`,
#'   `source_file`, and `header_format` identify the read, its input file, and
#'   the detected metadata convention. Other columns are added when their
#'   source fields occur in at least one header. Integer and floating-point SAM
#'   fields are converted to numeric columns where appropriate; time fields are
#'   returned as `POSIXct` in UTC.
#'
#' If no files match `file_format`, an empty data frame with the core columns is
#' returned.
#'
#' @section Header formats:
#'
#' A legacy header uses whitespace-separated `key=value` fields, for example:
#'
#' ```
#' "@read-id runid=RUN1 read=42 ch=1803 start_time=2026-09-20T14:53:20+08:00 barcode=barcode86"
#' ```
#'
#' A recent MinKNOW/Dorado header uses SAM auxiliary-field syntax,
#' `TAG:TYPE:VALUE`, for example:
#'
#' ```
#' "@read-id qs:f:15.7 ch:i:1803 st:Z:2026-09-20T14:53:20+08:00 SM:Z:barcode86"
#' ```
#'
#' In SAM-style fields, `i` denotes an integer, `f` a floating-point number,
#' and `Z` a string. Detection is performed independently for every read, so a
#' directory may contain old-format and new-format FASTQ files. The value of
#' `header_format` is `"key_value"`, `"sam_tag"`, `"mixed"`, or `"none"`.
#'
#' @section Standardized fields and their origins:
#'
#' The following MinKNOW/Dorado fields are renamed when present:
#'
#' * `read`: the UUID or identifier immediately following `@`; it is not the
#'   same as the acquisition read number.
#' * `mean_qscore`: legacy `mean_qscore_template` or SAM `qs:f`; mean basecall
#'   Q-score for the read.
#' * `mux`: legacy `mux` or SAM `mx:i`; mux/well number used by the channel.
#' * `ch`: legacy `ch` or SAM `ch:i`; flow-cell channel number.
#' * `read_number`: legacy `read` or SAM `rn:i`; acquisition read number
#'   assigned by MinKNOW.
#' * `start_time`: legacy `start_time` or SAM `st:Z`; start time of the read.
#' * `duration`: legacy `duration` or SAM `du:f`; read duration in seconds.
#' * `signal_start` and `signal_end`: SAM `ts:i` and `ns:i`. The basecalled
#'   sequence corresponds to the raw-signal interval `signal[ts:ns]`.
#'   `signal_samples` is derived as `signal_end - signal_start`. These values
#'   count signal samples, not DNA bases.
#' * `scaling_midpoint`, `scaling_dispersion`, and `scaling_version`: SAM
#'   `sm:f`, `sd:f`, and `sv:Z`; parameters describing conversion and
#'   normalization of the raw current signal.
#' * `duplex`: SAM `dx:i`; duplex indicator, normally zero for simplex and one
#'   for duplex reads.
#' * `read_group`: SAM `RG:Z`; read-group identifier, normally composed from
#'   the run ID, basecalling model, and barcode arrangement.
#' * `run_id`: legacy `runid`. For SAM-style headers the complete run/model
#'   identifier remains in `read_group` rather than being split heuristically.
#' * `experiment_start_time`: SAM `DT:Z`; experiment start time.
#' * `flow_cell_id`: legacy `flow_cell_id` or SAM `PU:Z`.
#' * `device_id`: SAM `PM:Z`; sequencing device identifier, when supplied.
#' * `sample_id`: legacy `sample_id` or SAM `LB:Z` in the ONT read-group
#'   convention.
#' * `barcode`: legacy `barcode` or SAM `SM:Z`; classified barcode name.
#' * `barcode_alias`: legacy `barcode_alias` or SAM `al:Z`; alias from the
#'   sample sheet. If `SM` is absent, a non-missing alias is also used as the
#'   standardized `barcode` value.
#' * `barcode_kit` and `trimming`: SAM `bk:Z` and `tm:Z`; barcode kit and the
#'   requested adapter/primer/barcode trimming configuration.
#' * `source_signal_file`: SAM `fn:Z`; original signal file name.
#' * `parent_read_id` and `parent_signal_start`: SAM `pi:Z` and `sp:i` for a
#'   read created by signal splitting.
#' * `bed_hits`, `minknow_events`, `pore_type`, `end_reason`,
#'   `barcode_variant`, and `polya_length`: SAM `bh:i`, `me:i`, `po:Z`, `er:Z`,
#'   `bv:Z`, and `pt:i`, respectively.
#'
#' Field names are case-sensitive. In particular, SAM `sm:f` is the signal
#' scaling midpoint, whereas `SM:Z` is the barcode/sample name. Recognized
#' legacy fields not listed above retain their original names. Unrecognized
#' SAM fields are retained with a `tag_` prefix, such as `tag_XY`.
#'
#' The SAM tag definitions follow the
#' [Dorado SAM specification](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/).
#'
#' @details
#' FASTQ itself standardizes the four-line record structure but does not
#' standardize metadata following the read identifier on the first line. The
#' two header conventions therefore describe similar information with
#' different encodings. This function standardizes the metadata without
#' changing [fastq_read_info()], which continues to implement its original
#' `key=value` behavior.
#'
#' Files are processed in sorted file-name order. Both plain-text and
#' gzip-compressed FASTQ files are supported. A complete four-line FASTQ record
#' is required; malformed or incomplete records produce an error.
#'
#' @examples
#' fastq_dir <- tempfile("fastq2-")
#' dir.create(fastq_dir)
#'
#' writeLines(
#'   c(
#'     paste0(
#'       "@new-read qs:f:15.715393 mx:i:4 ch:i:1803 rn:i:438 ",
#'       "st:Z:2026-09-20T14:53:20.912722+08:00 ts:i:10 ns:i:3694 ",
#'       "du:f:0.7388 dx:i:0 PU:Z:PBO11046 LB:Z:20260920pcr ",
#'       "SM:Z:barcode86 al:Z:barcode86"
#'     ),
#'     "ACGT", "+", "!!!!"
#'   ),
#'   file.path(fastq_dir, "new.fastq")
#' )
#'
#' info <- fastq_read_info2(
#'   fastq_dir,
#'   file_format = ".fastq",
#'   plot_time = FALSE
#' )
#' info[, c(
#'   "read", "header_format", "mean_qscore", "ch", "read_number",
#'   "start_time", "signal_samples", "barcode"
#' )]
#'
#' # The same call also recognizes legacy key=value metadata.
#' writeLines(
#'   c(
#'     paste0(
#'       "@old-read runid=RUN1 read=12 ch=7 ",
#'       "start_time=2026-09-20T15:00:00+08:00 barcode=barcode01"
#'     ),
#'     "TGCA", "+", "####"
#'   ),
#'   file.path(fastq_dir, "old.fastq")
#' )
#' mixed_info <- fastq_read_info2(
#'   fastq_dir,
#'   file_format = ".fastq",
#'   plot_time = FALSE
#' )
#' mixed_info[, c("read", "header_format", "ch", "barcode")]
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
  records <- Map(parse_fastq_header2, headers, source_file)
  out <- fastq_header_records_to_data_frame(records)

  if ("barcode_alias" %in% names(out)) {
    if (!"barcode" %in% names(out)) {
      out$barcode <- out$barcode_alias
    } else {
      missing_barcode <- is.na(out$barcode) | !nzchar(out$barcode)
      out$barcode[missing_barcode] <- out$barcode_alias[missing_barcode]
    }
  }

  if (all(c("signal_start", "signal_end") %in% names(out))) {
    out$signal_samples <- out$signal_end - out$signal_start
    after <- match("signal_end", names(out))
    out <- out[c(
      names(out)[seq_len(after)],
      "signal_samples",
      names(out)[seq.int(after + 1L, ncol(out))]
    )]
  }

  if (isTRUE(plot_time)) plot_fastq_start_times(out)
  out
}

parse_fastq_header2 <- function(header, source_file) {
  tokens <- strsplit(trimws(header), "[[:space:]]+", perl = TRUE)[[1L]]
  record <- list(
    read = tokens[[1L]],
    source_file = source_file
  )
  tokens <- tokens[-1L]
  has_key_value <- grepl("=", tokens, fixed = TRUE)
  sam_match <- regexec(
    "^([[:alnum:]]{2}):([AifZHBcCsSI]):(.*)$",
    tokens,
    perl = TRUE
  )
  sam_parts <- regmatches(tokens, sam_match)
  has_sam_tag <- lengths(sam_parts) == 4L

  record$header_format <- if (any(has_key_value) && any(has_sam_tag)) {
    "mixed"
  } else if (any(has_sam_tag)) {
    "sam_tag"
  } else if (any(has_key_value)) {
    "key_value"
  } else {
    "none"
  }

  for (token in tokens[has_key_value & !has_sam_tag]) {
    separator <- regexpr("=", token, fixed = TRUE)[[1L]]
    key <- substring(token, 1L, separator - 1L)
    value <- substring(token, separator + 1L)
    name <- fastq_header2_legacy_name(key)
    record[[name]] <- value
  }

  for (parts in sam_parts[has_sam_tag]) {
    tag <- parts[[2L]]
    type <- parts[[3L]]
    value <- parts[[4L]]
    name <- fastq_header2_sam_name(tag)
    record[[name]] <- value
    attr(record[[name]], "sam_type") <- type
  }
  record
}

fastq_header2_legacy_name <- function(key) {
  mapping <- c(
    runid = "run_id",
    read = "read_number",
    mean_qscore_template = "mean_qscore",
    mux = "mux",
    ch = "ch",
    start_time = "start_time",
    duration = "duration",
    flow_cell_id = "flow_cell_id",
    sample_id = "sample_id",
    barcode = "barcode",
    barcode_alias = "barcode_alias"
  )
  if (key %in% names(mapping)) unname(mapping[[key]]) else key
}

fastq_header2_sam_name <- function(tag) {
  mapping <- c(
    qs = "mean_qscore", mx = "mux", ch = "ch", rn = "read_number",
    st = "start_time", ts = "signal_start", ns = "signal_end",
    du = "duration", sm = "scaling_midpoint", sd = "scaling_dispersion",
    sv = "scaling_version", dx = "duplex", RG = "read_group",
    DT = "experiment_start_time", PU = "flow_cell_id", PM = "device_id",
    LB = "sample_id", SM = "barcode", al = "barcode_alias",
    bk = "barcode_kit", tm = "trimming", fn = "source_signal_file",
    pi = "parent_read_id", sp = "parent_signal_start", bh = "bed_hits",
    me = "minknow_events", po = "pore_type", er = "end_reason",
    bv = "barcode_variant", pt = "polya_length", DS = "read_group_description"
  )
  if (tag %in% names(mapping)) unname(mapping[[tag]]) else paste0("tag_", tag)
}

fastq_header_records_to_data_frame <- function(records) {
  columns <- unique(unlist(lapply(records, names), use.names = FALSE))
  out <- as.data.frame(
    matrix(
      NA_character_,
      nrow = length(records),
      ncol = length(columns),
      dimnames = list(NULL, columns)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (i in seq_along(records)) {
    for (name in names(records[[i]])) {
      out[[name]][[i]] <- as.character(records[[i]][[name]])
    }
  }

  integer_columns <- intersect(
    c(
      "mux", "ch", "read_number", "signal_start", "signal_end", "duplex",
      "parent_signal_start", "bed_hits", "minknow_events", "polya_length"
    ),
    names(out)
  )
  numeric_columns <- intersect(
    c("mean_qscore", "duration", "scaling_midpoint", "scaling_dispersion"),
    names(out)
  )
  for (name in integer_columns) {
    out[[name]] <- suppressWarnings(as.integer(out[[name]]))
  }
  for (name in numeric_columns) {
    out[[name]] <- suppressWarnings(as.numeric(out[[name]]))
  }
  for (name in intersect(c("start_time", "experiment_start_time"), names(out))) {
    out[[name]] <- parse_fastq_start_time(out[[name]])
  }
  out
}
