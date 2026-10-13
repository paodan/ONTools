#' Plot an integrated summary of multiple FASTQ files
#'
#' `plot_fastq_statistics()` creates a single quality-control figure combining
#' read-count distribution panels with read-length peak statistics. FASTQ
#' header format is detected automatically with [fastq_read_info_auto()].
#'
#' @param fastqs Non-empty character vector of `.fastq`, `.fastq.gz`, `.fq`, or
#'   `.fq.gz` files. Input order is preserved. When `path` is supplied, relative
#'   entries are resolved against that directory.
#' @param path Optional directory used to resolve relative entries in `fastqs`.
#' @param title Optional character scalar used as the title of the per-group
#'   read-count panel. Use `NULL` to omit the title.
#' @param barcode_column Optional character scalar naming the metadata column
#'   used to group reads in the distribution panels. When `NULL`, the function
#'   uses `barcode` if that column contains at least one non-missing, non-empty
#'   value; otherwise it uses `source_file`.
#' @inheritParams plot_peaks_statistics
#' @param heights Positive numeric vector of length two giving the relative
#'   heights of the read-distribution section and peak-statistics section.
#' @param ... Additional arguments passed to [base_theme()] through
#'   [plot_peaks_statistics()], for example `base_size`, `title_hjust`, or
#'   `x_angle`.
#'
#' @return A `patchwork` object with the read-distribution summary above the
#' read-length peak summary. The following attributes are attached:
#'
#' * `header_format`: format detected by [fastq_read_info_auto()];
#' * `parser`: metadata parser selected by [fastq_read_info_auto()];
#' * `barcode_column`: grouping column used for the read-distribution panels;
#' * `read_distribution`: the upper component plot;
#' * `peak_statistics`: the lower component plot.
#'
#' @details
#' Metadata are read once with [fastq_read_info_auto()] and passed to
#' [plot_fastq_read_distribution()]. Read lengths are then read independently
#' by [plot_peaks_statistics()] using the validated chunked FASTQ reader. The
#' first selected file determines whether legacy `key=value` or MinKNOW/Dorado
#' `TAG:TYPE:VALUE` headers are used; all files should follow the same header
#' convention.
#'
#' Both component plots are wrapped as complete graphical elements before they
#' are stacked. This prevents patchwork from flattening their internal layouts
#' and avoids nested-layout dimension errors.
#'
#' At least one read-length peak must be detected. Peak detection requirements
#' and parameter interpretation are documented in [plot_peaks_statistics()]
#' and [detect_length_peaks()].
#'
#' @examples
#' fastq_dir <- tempfile("fastq-statistics-")
#' dir.create(fastq_dir)
#' fastq <- file.path(fastq_dir, "sample.fastq")
#'
#' lengths <- c(95:105, 295:305)
#' records <- unlist(lapply(seq_along(lengths), function(i) {
#'   sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'   quality <- paste(rep("!", lengths[[i]]), collapse = "")
#'   c(
#'     paste0(
#'       "@read", i, " qs:f:15.5 ch:i:", (i %% 4) + 1,
#'       " SM:Z:barcode01"
#'     ),
#'     sequence,
#'     "+",
#'     quality
#'   )
#' }), use.names = FALSE)
#' writeLines(records, fastq)
#'
#' p <- plot_fastq_statistics(
#'   fastqs = fastq,
#'   title = "Example FASTQ summary",
#'   min_distance = 100,
#'   density_n = 1024
#' )
#' p
#' attr(p, "parser")
#' attr(p, "barcode_column")
#'
#' @seealso [fastq_read_info_auto()], [plot_fastq_read_distribution()],
#'   [plot_peaks_statistics()]
#' @export
plot_fastq_statistics <- function(
    fastqs,
    path = NULL,
    title = NULL,
    barcode_column = NULL,
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
    heights = c(1.25, 0.6),
    ...) {
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("`title` must be `NULL` or a character scalar.", call. = FALSE)
  }
  if (!is.null(barcode_column)) {
    check_scalar_character(barcode_column, "barcode_column")
  }
  check_logical_scalar(verbose, "verbose")
  if (!is.numeric(heights) || length(heights) != 2L || anyNA(heights) ||
      any(!is.finite(heights)) || any(heights <= 0)) {
    stop("`heights` must contain two positive finite numbers.", call. = FALSE)
  }
  heights <- as.numeric(heights)

  info <- fastq_read_info_auto(
    path = path,
    fastqs = fastqs,
    plot_time = FALSE,
    verbose = verbose
  )
  detected_format <- attr(info, "header_format", exact = TRUE)
  selected_parser <- attr(info, "parser", exact = TRUE)

  if (is.null(barcode_column)) {
    has_barcode <- "barcode" %in% names(info) && any(
      !is.na(info$barcode) & nzchar(as.character(info$barcode))
    )
    barcode_column <- if (has_barcode) "barcode" else "source_file"
  }
  if (!barcode_column %in% names(info)) {
    stop(
      "`barcode_column` was not found in parsed FASTQ metadata: ",
      barcode_column,
      call. = FALSE
    )
  }

  read_distribution <- plot_fastq_read_distribution(
    info = info,
    title = title,
    barcode_column = barcode_column
  )
  peak_statistics <- plot_peaks_statistics(
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
    verbose = verbose,
    combine_plot = TRUE,
    ...
  )

  plot <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = read_distribution),
    patchwork::wrap_elements(full = peak_statistics),
    ncol = 1,
    heights = heights
  )
  attr(plot, "header_format") <- detected_format
  attr(plot, "parser") <- selected_parser
  attr(plot, "barcode_column") <- barcode_column
  attr(plot, "read_distribution") <- read_distribution
  attr(plot, "peak_statistics") <- peak_statistics
  plot
}
