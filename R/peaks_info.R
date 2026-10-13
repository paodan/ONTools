#' Summarize read-length peaks from multiple FASTQ files
#'
#' `get_peaks_info()` detects read-length peaks in one or more FASTQ files and
#' combines the per-file results into one data frame. FASTQ parsing and peak
#' detection are delegated to [plot_fastq_peaks()], so this function uses the
#' same validated, chunked FASTQ reader and the same detection parameters.
#'
#' @param fastqs Character vector of FASTQ file paths. When `path` is supplied,
#'   these are interpreted relative to `path`. Both `.fastq` and `.fastq.gz`
#'   files are supported.
#' @param path Optional directory containing `fastqs`.
#' @inheritParams plot_read_length_peaks
#' @param verbose Logical. If `TRUE`, report each FASTQ file as it is processed.
#'
#' @return A data frame with one row per detected peak and the following
#'   columns:
#'
#' * `file`: unique input base name used by [plot_fastq_peaks()];
#' * `peak`: peak number within the file. Numbering starts at `2` because peak
#'   `1` is reserved for the caller's expected/default read-length interval;
#' * `major`: whether the peak has the largest read count in its file;
#' * `position`: peak position on the detection scale (base pairs when
#'   `log_scale = FALSE`, or `log10` base pairs when `log_scale = TRUE`);
#' * `read_length`: peak position converted to base pairs;
#' * `density`: kernel-density estimate at the peak;
#' * `window_start`, `window_end`: inclusive counting-window boundaries in
#'   base pairs;
#' * `read_count`, `read_percent`: number and percentage of reads in the
#'   counting window.
#'
#' Peaks are ordered from longest to shortest within each input file. If no
#' peaks are available, an empty data frame with the documented columns is
#' returned. Files that cannot support peak detection are handled as described
#' in [plot_fastq_peaks()] and contribute no rows.
#'
#' @details
#' The `major` peak is selected independently for each file using
#' `read_count`. If several peaks tie, the longest tied peak is selected
#' because rows are ordered by decreasing `read_length` before the comparison.
#'
#' This function reads peak information from documented plot attributes. It
#' does not access ggplot2's internal plot environment.
#'
#' @examples
#' fastq_dir <- tempfile("peak-info-")
#' dir.create(fastq_dir)
#'
#' write_peak_fastq <- function(path, lengths) {
#'   records <- unlist(lapply(seq_along(lengths), function(i) {
#'     sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'     quality <- paste(rep("!", lengths[[i]]), collapse = "")
#'     c(paste0("@read", i), sequence, "+", quality)
#'   }), use.names = FALSE)
#'   writeLines(records, path)
#' }
#'
#' write_peak_fastq(
#'   file.path(fastq_dir, "sample.fastq"),
#'   c(95:105, 295:305)
#' )
#'
#' peaks <- get_peaks_info(
#'   fastqs = "sample.fastq",
#'   path = fastq_dir,
#'   min_distance = 100,
#'   density_n = 1024
#' )
#' peaks
#'
#' @seealso [plot_fastq_peaks()], [plot_read_length_peaks()],
#'   [plot_peaks_statistics()]
#' @export
get_peaks_info <- function(
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
    verbose = FALSE) {
  plots <- plot_fastq_peaks(
    fastqs = fastqs,
    path = path,
    log_scale = log_scale,
    adjust = adjust,
    min_relative_height = min_relative_height,
    min_distance = min_distance,
    density_n = density_n,
    detection_range = detection_range,
    left_offset = left_offset,
    right_offset = right_offset,
    binwidth = binwidth,
    print_plots = FALSE,
    verbose = verbose
  )

  peak_tables <- lapply(names(plots), function(file) {
    plot <- plots[[file]]
    if (is.null(plot)) return(NULL)

    counts <- attr(plot, "peak_counts", exact = TRUE)
    result <- attr(plot, "peak_result", exact = TRUE)
    detected <- result$peaks
    if (is.null(counts) || nrow(counts) == 0L) return(NULL)

    if (is.null(detected) || nrow(detected) != nrow(counts)) {
      stop(
        "Peak metadata is inconsistent for `", file, "`.",
        call. = FALSE
      )
    }

    peak_data <- data.frame(
      file = rep(file, nrow(counts)),
      peak = integer(nrow(counts)),
      major = rep(FALSE, nrow(counts)),
      position = detected$position,
      read_length = counts$read_length,
      density = detected$density,
      window_start = counts$window_start,
      window_end = counts$window_end,
      read_count = counts$read_count,
      read_percent = counts$read_percent,
      stringsAsFactors = FALSE
    )
    peak_data <- peak_data[
      order(peak_data$read_length, decreasing = TRUE),
      ,
      drop = FALSE
    ]
    peak_data$peak <- seq.int(from = 2L, length.out = nrow(peak_data))
    peak_data$major[[which.max(peak_data$read_count)]] <- TRUE
    rownames(peak_data) <- NULL
    peak_data
  })

  peak_tables <- Filter(Negate(is.null), peak_tables)
  if (length(peak_tables) == 0L) return(empty_peaks_info())

  peaks <- do.call(rbind, peak_tables)
  rownames(peaks) <- NULL
  peaks
}

#' Plot summary statistics for FASTQ read-length peaks
#'
#' `plot_peaks_statistics()` runs [get_peaks_info()] and summarizes the
#' detected peaks across FASTQ files in three plots: the distribution of peak
#' counts per file, all peak lengths, and major peak lengths.
#'
#' @inheritParams get_peaks_info
#' @param ... Additional arguments passed to [base_theme()], for example
#'   `base_size`, `base_family`, `title_hjust`, or `x_angle`.
#'
#' @return A named list of three [ggplot2::ggplot()] objects:
#'
#' * `peak_num_hist`: distribution of the number of detected peaks per FASTQ
#'   file;
#' * `peak_pos_hist`: distribution of all detected peak lengths in base pairs;
#' * `peak_pos_hist_major`: distribution of major peak lengths in base pairs.
#'
#' @details
#' A major peak is the peak whose counting window contains the most reads in a
#' file, as defined by [get_peaks_info()]. The position histograms always use
#' `read_length` in base pairs, including when peak detection uses
#' `log_scale = TRUE`.
#'
#' At least one peak must be detected. The function stops with an informative
#' error if all inputs are skipped or no peaks pass the detection criteria.
#'
#' @examples
#' fastq_dir <- tempfile("peak-statistics-")
#' dir.create(fastq_dir)
#'
#' write_peak_fastq <- function(path, lengths) {
#'   records <- unlist(lapply(seq_along(lengths), function(i) {
#'     sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'     quality <- paste(rep("!", lengths[[i]]), collapse = "")
#'     c(paste0("@read", i), sequence, "+", quality)
#'   }), use.names = FALSE)
#'   writeLines(records, path)
#' }
#'
#' write_peak_fastq(
#'   file.path(fastq_dir, "sample1.fastq"),
#'   c(95:105, 295:305)
#' )
#' write_peak_fastq(
#'   file.path(fastq_dir, "sample2.fastq"),
#'   c(145:155, 395:405)
#' )
#'
#' statistics <- plot_peaks_statistics(
#'   fastqs = c("sample1.fastq", "sample2.fastq"),
#'   path = fastq_dir,
#'   min_distance = 100,
#'   density_n = 1024,
#'   title_hjust = 0.5
#' )
#' statistics$peak_num_hist
#' statistics$peak_pos_hist
#' statistics$peak_pos_hist_major
#'
#' @seealso [get_peaks_info()], [plot_fastq_peaks()]
#' @export
plot_peaks_statistics <- function(
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
    verbose = FALSE,
    ...) {
  peaks <- get_peaks_info(
    fastqs = fastqs,
    path = path,
    log_scale = log_scale,
    adjust = adjust,
    min_relative_height = min_relative_height,
    min_distance = min_distance,
    density_n = density_n,
    detection_range = detection_range,
    left_offset = left_offset,
    right_offset = right_offset,
    binwidth = binwidth,
    verbose = verbose
  )
  if (nrow(peaks) == 0L) {
    stop("No read-length peaks were detected in `fastqs`.", call. = FALSE)
  }

  peaks_per_file <- data.frame(
    peak_count = as.integer(table(peaks$file)),
    stringsAsFactors = FALSE
  )
  major_peaks <- peaks[peaks$major, , drop = FALSE]

  peak_num_hist = ggplot2::ggplot(
    peaks_per_file,
    ggplot2::aes(x = .data$peak_count)
  ) +
    ggplot2::geom_histogram(
      binwidth = 1,
      boundary = 0.5,
      colour = "white",
      alpha = 0.7
    ) +
    ggplot2::labs(
      title = "Distribution of peak counts",
      x = "Number of peaks",
      y = "Number of FASTQ files"
    ) +
    base_theme(x_angle = 0, ...)

  peak_pos_hist = ggplot2::ggplot(
    peaks,
    ggplot2::aes(x = .data$read_length)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      colour = "white",
      alpha = 0.7
    ) +
    ggplot2::labs(
      title = "Distribution of peak positions",
      x = "Peak read length (bp)",
      y = "Number of peaks"
    ) +
    base_theme(...)

  peak_pos_hist_major = ggplot2::ggplot(
    major_peaks,
    ggplot2::aes(x = .data$read_length)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      colour = "white",
      alpha = 0.7
    ) +
    ggplot2::labs(
      title = "Distribution of major peak positions",
      x = "Major peak read length (bp)",
      y = "Number of peaks"
    ) +
    base_theme(...)


  plot <- (peak_num_hist | peak_pos_hist | peak_pos_hist_major) +
    patchwork::plot_layout(
      widths = c(0.3, 0.35, 0.35)
    )

  return(plot)
}

empty_peaks_info <- function() {
  data.frame(
    file = character(),
    peak = integer(),
    major = logical(),
    position = numeric(),
    read_length = numeric(),
    density = numeric(),
    window_start = numeric(),
    window_end = numeric(),
    read_count = integer(),
    read_percent = numeric(),
    stringsAsFactors = FALSE
  )
}
