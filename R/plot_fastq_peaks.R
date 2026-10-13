#' Plot read-length peaks for multiple FASTQ files
#'
#' `plot_fastq_peaks()` reads the sequence lengths from one or more FASTQ files
#' and applies [plot_read_length_peaks()] to each file. FASTQ records are read
#' in chunks with [read_fastq_lengths()], the same validated reader used by
#' [split_fastq()]. This avoids the quote-parsing problems caused by treating
#' FASTQ quality strings as a delimited table.
#'
#' @param fastqs Character vector of FASTQ file paths. When `path` is supplied,
#'   these values are interpreted as file names or paths relative to `path`.
#'   Both `.fastq` and `.fastq.gz` files are supported.
#' @param path Optional directory containing `fastqs`. This is useful when
#'   `fastqs` was produced by `list.files()` without `full.names = TRUE`.
#' @inheritParams plot_read_length_peaks
#' @param title Single character string used as the title of every plot.
#' @param subtitle Optional character vector. If `NULL`, each plot uses its
#'   input file name as the subtitle. A single supplied value is reused for all
#'   plots; otherwise, supply one value per input file.
#' @param print_plots Logical. If `TRUE`, print each plot as it is created. The
#'   default is `TRUE` in an interactive R session and `FALSE` otherwise.
#' @param verbose Logical. If `TRUE`, report each input file before processing.
#'
#' @return Invisibly returns a named list with one element per input file.
#' Elements are [ggplot2::ggplot()] objects carrying the attributes documented
#' by [plot_read_length_peaks()]. An element is `NULL` when its input contains
#' no reads, fewer than three positive-length reads, or fewer than two distinct
#' positive read lengths; a warning identifies every skipped file. List names
#' are unique versions of the input base names.
#'
#' @details
#' Every FASTQ record must use the standard four-line representation. Headers
#' must begin with `@`, separator lines with `+`, and sequence and quality
#' strings must have equal lengths. Wrapped FASTQ records are not supported.
#'
#' Unlike `data.table::fread()`, the reader does not interpret quotation marks.
#' This matters because characters such as `"` are valid Phred+33 quality
#' symbols and do not indicate quoted fields in a FASTQ file.
#'
#' Peak detection, counting windows, and plotting are delegated unchanged to
#' [plot_read_length_peaks()]. Use `attr(plot, "peak_counts")` to retrieve the
#' detected peak table for any returned plot.
#'
#' @examples
#' fastq_dir <- tempfile("fastq-peaks-")
#' dir.create(fastq_dir)
#'
#' write_example_fastq <- function(path, lengths) {
#'   records <- unlist(lapply(seq_along(lengths), function(i) {
#'     sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'     quality <- paste(rep("!", lengths[[i]]), collapse = "")
#'     c(paste0("@read", i), sequence, "+", quality)
#'   }), use.names = FALSE)
#'   writeLines(records, path)
#' }
#'
#' write_example_fastq(
#'   file.path(fastq_dir, "sample1.fastq"),
#'   c(95:105, 295:305)
#' )
#' write_example_fastq(
#'   file.path(fastq_dir, "sample2.fastq"),
#'   c(145:155, 395:405)
#' )
#'
#' plots <- plot_fastq_peaks(
#'   fastqs = c("sample1.fastq", "sample2.fastq"),
#'   path = fastq_dir,
#'   min_distance = 100,
#'   x_breaks = NULL,
#'   print_plots = FALSE
#' )
#' names(plots)
#' attr(plots[[1L]], "peak_counts")
#'
#' @seealso [plot_read_length_peaks()], [get_peaks_info()],
#'   [plot_peaks_statistics()], [split_fastq()]
#' @export
plot_fastq_peaks <- function(
    fastqs,
    path = NULL,
    log_scale = FALSE,
    adjust = 1,
    min_relative_height = 0.05,
    min_distance = 400,
    density_n = 4096,
    detection_range = NULL,
    left_offset = 150,
    right_offset = 150,
    binwidth = 10,
    title = "Read-length distribution",
    subtitle = NULL,
    x_scale = c("auto", "log10", "linear"),
    x_breaks = NULL,
    x_limits = NULL,
    show_peak_labels = TRUE,
    show_peak_windows = TRUE,
    fill = "lightblue",
    colour = "steelblue",
    peak_colour = "red",
    window_fill = "orange",
    fill_alpha = 0.7,
    window_alpha = 0.12,
    print_plots = interactive(),
    verbose = FALSE) {
  if (!is.character(fastqs) || length(fastqs) == 0L ||
      anyNA(fastqs) || any(!nzchar(fastqs))) {
    stop("`fastqs` must be a non-empty character vector of file paths.",
         call. = FALSE)
  }
  if (!is.null(path)) {
    check_scalar_character(path, "path")
    if (!dir.exists(path)) {
      stop("`path` is not an existing directory: ", path, call. = FALSE)
    }
    fastqs <- file.path(path, fastqs)
  }
  check_logical_scalar(print_plots, "print_plots")
  check_logical_scalar(verbose, "verbose")
  check_scalar_character(title, "title")
  x_scale <- match.arg(x_scale)

  if (!is.null(subtitle)) {
    if (!is.character(subtitle) || anyNA(subtitle) ||
        any(!nzchar(subtitle)) || !length(subtitle) %in% c(1L, length(fastqs))) {
      stop(
        "`subtitle` must be `NULL`, one non-empty string, or one string per ",
        "input file.",
        call. = FALSE
      )
    }
    if (length(subtitle) == 1L) subtitle <- rep(subtitle, length(fastqs))
  }

  normalized_fastqs <- vapply(fastqs, function(fastq) {
    check_file_arg(fastq, "fastqs")
    split_fastq_suffix(fastq, name = "fastqs")
    normalizePath(fastq, mustWork = TRUE)
  }, character(1), USE.NAMES = FALSE)
  plot_names <- make.unique(basename(normalized_fastqs))
  plots <- vector("list", length(normalized_fastqs))
  names(plots) <- plot_names

  for (i in seq_along(normalized_fastqs)) {
    fastq <- normalized_fastqs[[i]]
    if (isTRUE(verbose)) message("Processing FASTQ: ", fastq)

    read_length <- read_fastq_lengths(fastq)
    positive_length <- read_length[is.finite(read_length) & read_length > 0]
    if (length(positive_length) < 3L || length(unique(positive_length)) < 2L) {
      warning(
        "Skipping `", basename(fastq),
        "`: peak plotting requires at least three positive-length reads and ",
        "two distinct read lengths.",
        call. = FALSE
      )
      plots[i] <- list(NULL)
      next
    }

    file_subtitle <- if (is.null(subtitle)) basename(fastq) else subtitle[[i]]
    plot <- plot_read_length_peaks(
      read_length = positive_length,
      log_scale = log_scale,
      adjust = adjust,
      min_relative_height = min_relative_height,
      min_distance = min_distance,
      density_n = density_n,
      detection_range = detection_range,
      left_offset = left_offset,
      right_offset = right_offset,
      binwidth = binwidth,
      title = title,
      subtitle = file_subtitle,
      x_scale = x_scale,
      x_breaks = x_breaks,
      x_limits = x_limits,
      show_peak_labels = show_peak_labels,
      show_peak_windows = show_peak_windows,
      fill = fill,
      colour = colour,
      peak_colour = peak_colour,
      window_fill = window_fill,
      fill_alpha = fill_alpha,
      window_alpha = window_alpha
    )
    plots[[i]] <- plot
    if (isTRUE(print_plots)) print(plot)
  }

  invisible(plots)
}
