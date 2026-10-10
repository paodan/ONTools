#' Plot FASTQ read distribution across barcodes
#'
#' `plot_fastq_read_distribution()` summarizes the read metadata returned by
#' [fastq_read_info()] and produces a four-panel quality-control figure. The
#' panels show read depth for each barcode, the distribution of barcode depths,
#' empirical depth percentiles, and a Lorenz curve of depth uniformity.
#'
#' @param info A data frame containing FASTQ read metadata, normally the output
#'   of [fastq_read_info()]. Each row represents one read.
#' @param title Optional character scalar used as the title of the per-barcode
#'   depth panel. Use `NULL` to omit the title.
#' @param barcode_column Character scalar naming the column in `info` that
#'   identifies the barcode or sample assigned to each read. Missing values and
#'   empty strings are grouped as `"unclassified"`.
#'
#' @return A `patchwork` object containing four ggplot panels:
#'
#' * a full-width bar chart of read depth for every barcode;
#' * a histogram of barcode-level depths;
#' * an empirical percentile curve of barcode-level depth;
#' * a Lorenz curve, including the summary statistics produced by
#'   [plot_depth_lorenz()].
#'
#' The per-barcode count table and percentile curve data are attached to the
#' returned object as the `sample_depth` and `quantile_data` attributes,
#' respectively. If `info` contains a `ch` column, the histogram is annotated
#' with the number of unique non-missing sequencing channels.
#'
#' @details
#' Read depth is the number of rows assigned to each barcode; therefore, the
#' function assumes that every row of `info` represents one read. Barcodes with
#' no rows cannot be inferred and are not displayed.
#'
#' The upper panel displays the number of reads assigned to each barcode. The
#' three lower panels provide complementary views of the same barcode-level
#' depths: their frequency distribution, empirical percentile curve, and
#' overall inequality.
#'
#' For the percentile panel, percentiles from 0 to 100 are evaluated at
#' one-percent intervals using [stats::quantile()] with its default
#' interpolation method. Its x-axis is barcode depth and its y-axis is the
#' corresponding percentile. A rapidly increasing curve indicates a wider
#' spread between shallow and deeply sequenced barcodes.
#'
#' The Lorenz panel orders barcode depths from smallest to largest and compares
#' cumulative barcode proportion with cumulative read proportion. A curve near
#' the diagonal indicates uniform depth, whereas stronger downward curvature
#' indicates that a small number of barcodes account for a large proportion of
#' reads. The annotation reports sample count, median depth, CV, robust CV,
#' Gini coefficient, P95/P5, percentage within 0.5--2 times the median, and
#' zero-depth percentage. See [plot_depth_lorenz()] for definitions and
#' interpretation of these statistics.
#'
#' Missing or empty values in `barcode_column` are counted together as
#' `"unclassified"`. Because the function summarizes rows actually present in
#' `info`, expected barcodes with no reads cannot be inferred automatically.
#' For a between-project comparison that must include expected zero-depth
#' samples, construct the complete depth vector separately and pass it to
#' [plot_depth_lorenz()].
#'
#' The function does not modify `info`. It uses `patchwork` to assemble the
#' four ggplot panels, with the per-barcode bar chart above the three summary
#' panels.
#'
#' @examples
#' barcode_depth <- c(120, 116, 110, 105, 101, 98, 94, 90, 86, 80, 72, 60)
#' barcode <- sprintf("barcode%02d", seq_along(barcode_depth))
#' info <- data.frame(
#'   read = paste0("read", seq_len(sum(barcode_depth))),
#'   barcode = rep(barcode, times = barcode_depth),
#'   ch = rep(seq_len(24), length.out = sum(barcode_depth)),
#'   stringsAsFactors = FALSE
#' )
#'
#' plot_fastq_read_distribution(
#'   info,
#'   title = "Example sequencing run"
#' )
#'
#' # Use a different column for sample assignments.
#' info$sample <- sub("barcode", "sample_", info$barcode)
#' p <- plot_fastq_read_distribution(
#'   info,
#'   title = "Reads per sample",
#'   barcode_column = "sample"
#' )
#' attr(p, "sample_depth")
#' head(attr(p, "quantile_data"))
#'
#' @seealso [fastq_read_info()], [plot_depth_lorenz()],
#'   [summarize_depth_uniformity()]
#' @export
plot_fastq_read_distribution <- function(info,
                                         title = NULL,
                                         barcode_column = "barcode") {
  if (!is.data.frame(info)) {
    stop("`info` must be a data frame.", call. = FALSE)
  }
  if (nrow(info) == 0L) {
    stop("`info` must contain at least one read.", call. = FALSE)
  }
  check_scalar_character(barcode_column, "barcode_column")
  if (!barcode_column %in% names(info)) {
    stop(
      "`barcode_column` was not found in `info`: ", barcode_column,
      call. = FALSE
    )
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("`title` must be `NULL` or a character scalar.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package `patchwork` is required. Install it with ",
      "`install.packages(\"patchwork\")`.",
      call. = FALSE
    )
  }

  barcode <- as.character(info[[barcode_column]])
  barcode[is.na(barcode) | !nzchar(barcode)] <- "unclassified"
  sample_depth <- as.data.frame(table(barcode), stringsAsFactors = FALSE)
  names(sample_depth) <- c("sample", "depth")
  sample_depth$sample <- as.character(sample_depth$sample)
  sample_depth$depth <- as.integer(sample_depth$depth)

  channel_label <- NULL
  if ("ch" %in% names(info)) {
    channels <- info$ch[!is.na(info$ch)]
    channel_label <- paste(length(unique(channels)), "channels in total")
  }

  p <- seq(0, 1, by = 0.01)
  depth_quantile <- as.numeric(stats::quantile(sample_depth$depth, p))
  max_depth <- max(sample_depth$depth)
  relative_quantile <- if (max_depth > 0) {
    depth_quantile / max_depth * 100
  } else {
    rep(0, length(depth_quantile))
  }

  quantile_data <- data.frame(
    percentile = p * 100,
    value = depth_quantile,
    type = "Depth",
    stringsAsFactors = FALSE
  )

  base_theme <- ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      )
    )

  g1 <- ggplot2::ggplot(
    sample_depth,
    ggplot2::aes(x = .data$sample, y = .data$depth)
  ) +
    ggplot2::geom_col(fill = "lightblue") +
    ggplot2::labs(title = title, x = "Samples", y = "Depth") +
    base_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
      plot.margin = ggplot2::margin(10, 10, 2, 5, unit = "pt")
    )

  g2 <- ggplot2::ggplot(
    sample_depth,
    ggplot2::aes(x = .data$depth)
  ) +
    ggplot2::geom_histogram(bins = 30, fill = "grey", colour = "white") +
    ggplot2::labs(x = "Depth", y = "Number of samples") +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30),
                   plot.margin = ggplot2::margin(2, 5, 2, 5, unit = "pt"))
  if (!is.null(channel_label)) {
    g2 <- g2 + ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = channel_label,
      vjust = 1.2,
      hjust = 1.05,
      colour = "orange"
    )
  }

  g3 <- ggplot2::ggplot(
    quantile_data,
    ggplot2::aes(x = .data$value, y = .data$percentile)) +
    ggplot2::geom_line(colour = "blue") +
    # ggplot2::facet_wrap(~type, nrow = 2, scales = "free_x") +
    ggplot2::labs(x = "Depth", y = "Percentile (%)") +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30),
                   plot.margin = ggplot2::margin(2, 15, 2, 5, unit = "pt"))

  g4 = plot_depth_lorenz(sample_depth$depth)

  bottom <- (g2 | g3 | g4) +
    patchwork::plot_layout(
      widths = c(0.6, 0.8, 1)
    )

  plot <- g1 / bottom +
    patchwork::plot_layout(
      heights = c(0.25, 1)
    )

  # plot <- g1 / (g2 | g3 | g4) +
  #   patchwork::plot_layout(
  #     widths = c(0.6, .8, 1),
  #     heights = c(0.2, 1))
  attr(plot, "sample_depth") <- sample_depth
  attr(plot, "quantile_data") <- quantile_data
  plot
}
