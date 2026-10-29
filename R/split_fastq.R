#' Split a FASTQ file by expected and detected read-length ranges
#'
#' `split_fastq()` creates one length-filtered FASTQ file for a user-supplied
#' expected range and, optionally, one file for every read-length peak detected
#' by [plot_read_length_peaks()]. Input records are validated and read in chunks
#' by [read_fastq_lengths()], so the complete FASTQ sequence data are not held
#' in memory.
#'
#' @param fastq_file Path to an input `.fastq` or `.fastq.gz` file. Records must
#'   use the standard four-line FASTQ representation.
#' @param out_fastq_path Output directory. When `NULL`, a directory named
#'   `peaks` is created beside `fastq_file`.
#' @param expected_size Positive integer giving the expected read length in base
#'   pairs. It is used to label the first output range and must lie between
#'   `min_read_length` and `max_read_length`.
#' @param min_read_length,max_read_length Positive integers defining the
#'   inclusive expected-length range used for the first output file.
#' @param by_peak Logical. If `TRUE`, detect additional read-length peaks and
#'   create an output for each detected peak. If fewer than three reads or fewer
#'   than two distinct read lengths are available, peak detection is skipped
#'   with a warning and the expected-range output is still created.
#' @param min_relative_height Number between zero and one. Candidate peaks below
#'   this fraction of the highest density are discarded.
#' @param log_scale Logical. If `TRUE`, estimate the read-length density on a
#'   log10 scale. Reported lengths and filtering ranges remain in base pairs.
#' @param min_distance Optional non-negative minimum distance in base pairs
#'   between detected peaks. Set to `NULL` to disable distance filtering.
#' @param density_n Integer of at least three giving the density-grid resolution
#'   passed to [plot_read_length_peaks()].
#' @param detection_range the range of detecting peaks, default is NULL, meaning
#' to detect peaks from all reads
#' @param left_offset left boundary offset from the peaks, default 150.
#' @param right_offset right boundary offset from the peaks, default 150.
#' @param x_breaks Optional finite numeric vector of plot-axis breaks in base
#'   pairs. It affects only the returned peak plot.
#' @param x_scale Plot x-axis scale: `"auto"`, `"log10"`, or `"linear"`.
#'   It affects only the returned peak plot, not peak detection or filtering.
#' @param seqkit Command name or executable path for `seqkit`.
#' @param conda_env Optional conda environment name. If supplied, `seqkit` is
#'   run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, validate the input, calculate the ranges,
#'   and return the planned [filter_fastq()] calls without executing `seqkit`.
#' @param echo Logical. If `TRUE`, print each planned `seqkit` command.
#' @param stderr Passed to [system2()] by [filter_fastq()].
#'
#' @return For a non-empty input, invisibly returns a list with:
#'
#' * `input_fastq`: normalized input path.
#' * `output_dir`: normalized output-directory path.
#' * `output_fastq`: named character vector of planned or generated files.
#' * `peak_counts`: data frame containing the expected or detected source,
#'   representative read length, inclusive filtering boundaries, read count,
#'   and percentage of all input reads.
#' * `filter_results`: one result from [filter_fastq()] per output file.
#' * `plot`: the peak plot when peak detection was performed, otherwise `NULL`.
#' * `read_count`: total number of input FASTQ records.
#' * `dry_run`: the supplied dry-run setting.
#'
#' If the input contains no FASTQ records, the function issues a warning,
#' creates no output files, and invisibly returns `NULL`. This allows a batch
#' loop to skip samples with no reads while continuing with later inputs.
#'
#' Output files retain the compression format of the input and are named
#' `<input>_peak_<index>_<length>bp.fastq[.gz]`. The first output always
#' represents the explicit `min_read_length`--`max_read_length` range; detected
#' peaks follow in decreasing read-length order. Range boundaries are inclusive.
#'
#' @details
#' FASTQ validation checks that every record has four lines, headers begin with
#' `@`, separator lines begin with `+`, and sequence and quality strings have
#' equal lengths. Wrapped FASTQ records are not supported.
#'
#' Read lengths are scanned once in R. Actual output writing is delegated to
#' [filter_fastq()], which uses `seqkit seq -m/-M -o`. Consequently, outputs
#' ending in `.fastq.gz` are genuinely gzip-compressed by `seqkit`.
#'
#' Peak windows use the default `left_offset` and `right_offset` from
#' [plot_read_length_peaks()]. Decimal boundaries produced by density estimation
#' are rounded outward before filtering: the lower boundary is floored (with a
#' minimum of 1 bp) and the upper boundary is ceilinged.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq")
#' records <- c(
#'   "@read1", "ACGT",   "+", "!!!!",
#'   "@read2", "ACGTA",  "+", "!!!!!",
#'   "@read3", "ACGTAC", "+", "!!!!!!"
#' )
#' writeLines(records, reads)
#'
#' output_dir <- tempfile("split-fastq-")
#' plan <- split_fastq(
#'   reads,
#'   out_fastq_path = output_dir,
#'   expected_size = 5,
#'   min_read_length = 4,
#'   max_read_length = 6,
#'   by_peak = FALSE,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' plan$peak_counts
#' plan$output_fastq
#' plan$filter_results[[1L]]$command_string
#'
#' \dontrun{
#' result <- split_fastq(
#'   "sample.fastq.gz",
#'   expected_size = 1600,
#'   min_read_length = 1200,
#'   max_read_length = 1800,
#'   by_peak = TRUE
#' )
#' result$peak_counts
#' result$plot
#' }
#'
#' @export
split_fastq <- function(fastq_file,
                        out_fastq_path = NULL,
                        expected_size = 1600,
                        min_read_length = 1200,
                        max_read_length = 1800,
                        by_peak = TRUE,
                        min_relative_height = 0.05,
                        log_scale = FALSE,
                        min_distance = 400,
                        density_n = 5000,
                        detection_range = NULL,
                        left_offset = 150,
                        right_offset = 150,
                        x_breaks = c(1, seq(500, 3000, by = 500)),
                        x_scale = c("linear", "auto", "log10"),
                        seqkit = "seqkit",
                        conda_env = NULL,
                        conda = "conda",
                        dry_run = FALSE,
                        echo = TRUE,
                        stderr = "") {
  check_file_arg(fastq_file, "fastq_file")
  check_logical_scalar(by_peak, "by_peak")
  check_logical_scalar(log_scale, "log_scale")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(conda, "conda")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")

  expected_size <- validate_positive_integer(expected_size, "expected_size")
  min_read_length <- validate_positive_integer(
    min_read_length,
    "min_read_length"
  )
  max_read_length <- validate_positive_integer(
    max_read_length,
    "max_read_length"
  )
  if (min_read_length > max_read_length) {
    stop("`min_read_length` must not exceed `max_read_length`.", call. = FALSE)
  }
  if (expected_size < min_read_length || expected_size > max_read_length) {
    stop(
      "`expected_size` must lie between `min_read_length` and ",
      "`max_read_length`.",
      call. = FALSE
    )
  }

  min_relative_height <- read_length_validate_fraction(
    min_relative_height,
    "min_relative_height"
  )
  if (!is.null(min_distance)) {
    min_distance <- read_length_validate_nonnegative_number(
      min_distance,
      "min_distance"
    )
  }
  density_n <- validate_density_grid_size(density_n, "density_n")
  x_scale <- match.arg(x_scale)

  fastq_file <- normalizePath(fastq_file, mustWork = TRUE)
  suffix <- split_fastq_suffix(fastq_file)
  if (is.null(out_fastq_path)) {
    out_fastq_path <- file.path(dirname(fastq_file), "peaks")
  } else {
    check_scalar_character(out_fastq_path, "out_fastq_path")
  }
  dir.create(out_fastq_path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_fastq_path)) {
    stop("Could not create `out_fastq_path`: ", out_fastq_path, call. = FALSE)
  }
  out_fastq_path <- normalizePath(out_fastq_path, mustWork = TRUE)

  read_length <- read_fastq_lengths(fastq_file)
  if (length(read_length) == 0L) {
    warning(
      "`fastq_file` contains no FASTQ records; skipping: ",
      fastq_file,
      call. = FALSE
    )
    return(invisible(NULL))
  }

  peak_counts <- data.frame(
    source = "expected",
    read_length = as.numeric(expected_size),
    window_start = as.numeric(min_read_length),
    window_end = as.numeric(max_read_length),
    stringsAsFactors = FALSE
  )
  peak_plot <- NULL

  if (isTRUE(by_peak)) {
    if (length(read_length) < 3L || length(unique(read_length)) < 2L) {
      warning(
        "Peak detection requires at least three reads and two distinct read ",
        "lengths; only the expected range will be written.",
        call. = FALSE
      )
    } else {
      peak_plot <- plot_read_length_peaks(
        read_length = read_length,
        min_relative_height = min_relative_height,
        log_scale = log_scale,
        min_distance = min_distance,
        density_n = density_n,
        detection_range = detection_range,
        left_offset = left_offset,
        right_offset = right_offset,
        title = sub(
          "[.]fastq([.]gz)?$",
          "",
          basename(fastq_file),
          ignore.case = TRUE
        ),
        x_breaks = x_breaks,
        x_scale = x_scale
      )
      detected <- attr(peak_plot, "peak_counts", exact = TRUE)
      if (!is.null(detected) && nrow(detected) > 0L) {
        detected <- detected[nrow(detected):1L, , drop = FALSE]
        detected$source <- "detected"
        detected <- detected[
          , c("source", "read_length", "window_start", "window_end"),
          drop = FALSE
        ]
        peak_counts <- rbind(peak_counts, detected)
      }
    }
  }

  peak_counts$window_start <- pmax(1L, floor(peak_counts$window_start))
  peak_counts$window_end <- ceiling(peak_counts$window_end)
  peak_counts$read_count <- vapply(seq_len(nrow(peak_counts)), function(i) {
    sum(
      read_length >= peak_counts$window_start[[i]] &
        read_length <= peak_counts$window_end[[i]]
    )
  }, integer(1))
  peak_counts$read_percent <- round(
    peak_counts$read_count / length(read_length) * 100,
    digits = 1
  )
  peak_counts$read_length <- round(peak_counts$read_length)
  rownames(peak_counts) <- NULL

  stem <- sub(
    "[.]fastq([.]gz)?$",
    "",
    basename(fastq_file),
    ignore.case = TRUE
  )
  output_fastq <- file.path(
    out_fastq_path,
    paste0(
      stem,
      "_peak_", seq_len(nrow(peak_counts)),
      "_", peak_counts$read_length, "bp",
      suffix
    )
  )
  names(output_fastq) <- paste0("peak_", seq_along(output_fastq))

  filter_results <- lapply(seq_len(nrow(peak_counts)), function(i) {
    filter_fastq(
      fastq = fastq_file,
      output_fastq = output_fastq[[i]],
      min_length = peak_counts$window_start[[i]],
      max_length = peak_counts$window_end[[i]],
      remove_gaps = TRUE,
      seqkit = seqkit,
      conda_env = conda_env,
      conda = conda,
      dry_run = dry_run,
      echo = echo,
      stderr = stderr
    )
  })
  names(filter_results) <- names(output_fastq)

  invisible(list(
    input_fastq = fastq_file,
    output_dir = out_fastq_path,
    output_fastq = output_fastq,
    peak_counts = peak_counts,
    filter_results = filter_results,
    plot = peak_plot,
    read_count = length(read_length),
    dry_run = dry_run
  ))
}

split_fastq_suffix <- function(path, name = "fastq_file") {
  filename <- basename(path)
  if (grepl("[.]fastq[.]gz$", filename, ignore.case = TRUE)) {
    return(".fastq.gz")
  }
  if (grepl("[.]fastq$", filename, ignore.case = TRUE)) {
    return(".fastq")
  }
  stop("`", name, "` must end in `.fastq` or `.fastq.gz`.", call. = FALSE)
}

#' Read sequence lengths from a FASTQ file
#'
#' `read_fastq_lengths()` extracts the sequence length of every record in a
#' four-line FASTQ or FASTQ.GZ file. Records are read and validated in chunks,
#' so sequence and quality strings from the complete file are never held in
#' memory at the same time.
#'
#' @param path Path to an input `.fastq` or `.fastq.gz` file.
#' @param chunk_records Positive integer giving the maximum number of FASTQ
#'   records read and validated per chunk. Larger values may improve throughput
#'   at the cost of temporary memory use. This does not limit the total number
#'   of records returned.
#'
#' @return An integer vector containing one sequence length in base pairs per
#' FASTQ record, in input order. An empty FASTQ returns `integer(0)`.
#'
#' @details
#' Each record must use the standard four-line FASTQ representation: a header
#' beginning with `@`, one unwrapped sequence line, a separator beginning with
#' `+`, and one unwrapped quality line. Sequence and quality strings must have
#' identical byte lengths. Wrapped FASTQ records are not supported.
#'
#' Both plain-text and gzip-compressed inputs are supported. Quality strings are
#' read as literal text; quotation marks and other valid Phred+33 characters are
#' not interpreted as table delimiters or quoting syntax.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq")
#' writeLines(
#'   c(
#'     "@read1", "ACGT", "+", "!!!!",
#'     "@read2", "ACGTAC", "+", "!!!!!!"
#'   ),
#'   reads
#' )
#' read_fastq_lengths(reads)
#'
#' compressed <- tempfile(fileext = ".fastq.gz")
#' con <- gzfile(compressed, open = "wt")
#' writeLines(c("@read3", "ACGTA", "+", "!!!!!"), con)
#' close(con)
#' read_fastq_lengths(compressed, chunk_records = 1)
#'
#' @seealso [split_fastq()], [plot_fastq_peaks()]
#' @export
read_fastq_lengths <- function(path, chunk_records = 100000L) {
  check_file_arg(path, "path")
  split_fastq_suffix(path, name = "path")
  chunk_records <- validate_positive_integer(chunk_records, "chunk_records")
  path <- normalizePath(path, mustWork = TRUE)

  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  read_length <- integer()
  repeat {
    block <- readLines(con, n = 4L * chunk_records, warn = FALSE)
    if (length(block) == 0L) break
    if (length(block) %% 4L != 0L) {
      stop("Incomplete FASTQ record in file: ", path, call. = FALSE)
    }

    header_index <- seq.int(1L, length(block), by = 4L)
    sequence_index <- header_index + 1L
    plus_index <- header_index + 2L
    quality_index <- header_index + 3L
    if (any(!startsWith(block[header_index], "@")) ||
        any(!startsWith(block[plus_index], "+"))) {
      stop("Invalid four-line FASTQ record in file: ", path, call. = FALSE)
    }

    sequence_length <- nchar(block[sequence_index], type = "bytes")
    quality_length <- nchar(block[quality_index], type = "bytes")
    if (any(sequence_length != quality_length)) {
      stop(
        "Sequence and quality lengths differ in FASTQ file: ",
        path,
        call. = FALSE
      )
    }
    read_length <- c(read_length, sequence_length)
  }

  read_length
}

#' Read mean quality values from a FASTQ file
#'
#' `read_fastq_qvalues()` calculates the mean Phred quality score of every
#' record in a four-line FASTQ or FASTQ.GZ file. It has the same arguments as
#' [read_fastq_lengths()] and reads records in chunks, avoiding storage of all
#' sequence and quality strings in memory at once.
#'
#' @param path Path to an input `.fastq` or `.fastq.gz` file.
#' @param chunk_records Positive integer giving the maximum number of FASTQ
#'   records read and validated per chunk. Larger values may improve throughput
#'   at the cost of temporary memory use. This does not limit the total number
#'   of records returned.
#'
#' @return A numeric vector containing one mean Phred quality score per FASTQ
#'   record, in input order. An empty FASTQ returns `numeric(0)`.
#'
#' @details
#' Quality strings are interpreted using Phred+33 encoding. For each quality
#' character, the Phred score is its ASCII code minus 33. Following Dorado's
#' definition, the leading 60 bases are excluded when a read is longer than 60
#' bases. The remaining scores are converted to error probabilities, those
#' probabilities are averaged, and the mean is converted back to Phred scale:
#'
#' `mean_q = -10 * log10(mean(10^(-base_q / 10)))`.
#'
#' For example, `!` represents Q0, `+` represents Q10, `5` represents Q20,
#' `?` represents Q30, and `I` represents Q40. Because averaging is performed
#' on the error-probability scale, the result is not generally equal to the
#' arithmetic mean of the per-base Q scores.
#'
#' Each record must use the standard four-line FASTQ representation: a header
#' beginning with `@`, one unwrapped sequence line, a separator beginning with
#' `+`, and one unwrapped quality line. Sequence and quality strings must have
#' identical byte lengths. Wrapped FASTQ records are not supported. Quality
#' characters must fall in the printable Phred+33 ASCII range 33--126.
#'
#' Both plain-text and gzip-compressed inputs are supported. The function
#' calculates the mean from the FASTQ quality string itself; it does not use a
#' `qs` value that may be present in a MinKNOW/Dorado FASTQ header. The leading
#' 60-base rule and probability-scale mean reproduce the current Dorado method
#' for ordinary DNA reads. Values can still differ for data produced with a
#' different basecaller/version or workflows with molecule-specific trimming,
#' such as direct RNA poly(A)-tail handling.
#'
#' @examples
#' reads <- tempfile(fileext = ".fastq")
#' writeLines(
#'   c(
#'     "@read1", "ACGT", "+", "IIII", # four Q40 bases
#'     "@read2", "ACGT", "+", "!+5?"  # Q0, Q10, Q20, Q30
#'   ),
#'   reads
#' )
#' read_fastq_qvalues(reads)
#'
#' compressed <- tempfile(fileext = ".fastq.gz")
#' con <- gzfile(compressed, open = "wt")
#' writeLines(c("@read3", "ACGTA", "+", "55555"), con)
#' close(con)
#' read_fastq_qvalues(compressed, chunk_records = 1)
#'
#' @seealso [read_fastq_lengths()], [fastq_read_info2()]
#' @export
read_fastq_qvalues <- function(path, chunk_records = 100000L) {
  check_file_arg(path, "path")
  split_fastq_suffix(path, name = "path")
  chunk_records <- validate_positive_integer(chunk_records, "chunk_records")
  path <- normalizePath(path, mustWork = TRUE)

  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  mean_qvalue <- numeric()
  repeat {
    block <- readLines(con, n = 4L * chunk_records, warn = FALSE)
    if (length(block) == 0L) break
    if (length(block) %% 4L != 0L) {
      stop("Incomplete FASTQ record in file: ", path, call. = FALSE)
    }

    header_index <- seq.int(1L, length(block), by = 4L)
    sequence_index <- header_index + 1L
    plus_index <- header_index + 2L
    quality_index <- header_index + 3L
    if (any(!startsWith(block[header_index], "@")) ||
        any(!startsWith(block[plus_index], "+"))) {
      stop("Invalid four-line FASTQ record in file: ", path, call. = FALSE)
    }

    sequence_length <- nchar(block[sequence_index], type = "bytes")
    quality <- block[quality_index]
    quality_length <- nchar(quality, type = "bytes")
    if (any(sequence_length != quality_length)) {
      stop(
        "Sequence and quality lengths differ in FASTQ file: ",
        path,
        call. = FALSE
      )
    }

    chunk_qvalue <- vapply(
      quality,
      FUN.VALUE = numeric(1L),
      FUN = function(value) {
        score <- as.integer(charToRaw(value)) - 33L
        if (any(score < 0L | score > 93L)) {
          stop(
            "Invalid Phred+33 quality character in FASTQ file: ",
            path,
            call. = FALSE
          )
        }
        if (length(score) == 0L) return(NaN)
        if (length(score) > 60L) score <- score[-seq_len(60L)]
        -10 * log10(mean(10^(-score / 10)))
      },
      USE.NAMES = FALSE
    )
    mean_qvalue <- c(mean_qvalue, chunk_qvalue)
  }

  mean_qvalue
}
