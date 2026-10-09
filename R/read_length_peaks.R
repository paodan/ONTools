#' Detect peaks in a read-length distribution
#'
#' `detect_length_peaks()` smooths positive read lengths with a kernel density
#' estimate and identifies local maxima in that curve. Detection can be
#' performed on the original base-pair scale or on a log10-transformed scale.
#'
#' @param read_length Numeric vector of read lengths in base pairs. Missing,
#'   non-finite, zero, and negative values are removed before analysis.
#' @param log_scale Logical. If `TRUE`, estimate density for
#'   `log10(read_length)`. This is useful when lengths span several orders of
#'   magnitude. Reported `read_length` peak positions are always converted back
#'   to base pairs.
#' @param adjust Positive bandwidth multiplier passed to [stats::density()].
#'   Values below one produce a less-smoothed curve and can reveal more peaks;
#'   values above one produce a smoother curve and can merge nearby peaks.
#' @param min_relative_height Number between zero and one. Candidate peaks with
#'   density below this fraction of the highest candidate density are removed.
#' @param min_distance Optional non-negative minimum distance, in base pairs,
#'   between retained peaks. When candidate peaks are closer than this value,
#'   the higher-density peak is retained first.
#' @param n Integer of at least three giving the number of equally spaced
#'   points at which [stats::density()] is evaluated. Larger values improve the
#'   numerical resolution of peak locations but do not add information to the
#'   reads; `adjust` is the primary smoothing control.
#' @param detection_range Optional numeric vector `c(lower, upper)` defining a
#'   closed base-pair interval. Only reads inside this interval participate in
#'   density estimation and peak detection. This is useful for preventing a few
#'   extreme read lengths from stretching the detection grid. It is one
#'   continuous search interval and does not limit the number of peaks within
#'   that interval.
#'
#' @return A list with:
#'
#' * `peaks`: data frame containing `position` (on the detection scale),
#'   `density`, and `read_length` (in base pairs).
#' * `density`: the object returned by [stats::density()]. Its x coordinates
#'   are log10 lengths when `log_scale = TRUE`.
#' * `log_scale` and `detection_range`: settings used for detection.
#' * `n_reads_total`: number of positive, finite input lengths.
#' * `n_reads_detection`: number of reads used after applying
#'   `detection_range`.
#'
#' @details
#' A candidate peak is a density-grid point whose value is greater than both
#' immediate neighbors. `min_relative_height` then removes low-density
#' candidates, and `min_distance` optionally removes nearby candidates in
#' decreasing order of density.
#'
#' The result is descriptive rather than a formal test of multimodality. Peak
#' calls depend on bandwidth, sample size, detection range, and biological or
#' technical artifacts such as truncated reads. Inspecting the distribution is
#' recommended before interpreting a peak as a distinct molecule population.
#'
#' @examples
#' set.seed(1)
#' lengths <- c(
#'   rnorm(500, mean = 250, sd = 20),
#'   rnorm(1200, mean = 1500, sd = 45),
#'   100000, 105000
#' )
#'
#' peaks <- detect_length_peaks(
#'   lengths,
#'   log_scale = TRUE,
#'   adjust = 0.8,
#'   min_relative_height = 0.05,
#'   min_distance = 300,
#'   detection_range = c(100, 2500)
#' )
#' peaks$peaks
#'
#' @export
detect_length_peaks <- function(read_length,
                                log_scale = FALSE,
                                adjust = 1,
                                min_relative_height = 0.05,
                                min_distance = NULL,
                                n = 4096,
                                detection_range = NULL) {
  read_length <- validate_read_lengths(read_length)
  check_logical_scalar(log_scale, "log_scale")
  adjust <- read_length_validate_positive_number(adjust, "adjust")
  min_relative_height <- read_length_validate_fraction(
    min_relative_height,
    "min_relative_height"
  )
  if (!is.null(min_distance)) {
    min_distance <- read_length_validate_nonnegative_number(min_distance, "min_distance")
  }
  n <- validate_density_grid_size(n)
  detection_range <- validate_detection_range(detection_range)

  detection_length <- read_length
  if (!is.null(detection_range)) {
    detection_length <- read_length[
      read_length >= detection_range[[1L]] &
        read_length <= detection_range[[2L]]
    ]
    if (length(detection_length) < 3L) {
      stop(
        "Fewer than three reads remain inside `detection_range`.",
        call. = FALSE
      )
    }
  }
  if (length(unique(detection_length)) < 2L) {
    stop("Read lengths used for detection must contain variation.", call. = FALSE)
  }

  x <- if (isTRUE(log_scale)) log10(detection_length) else detection_length
  den <- stats::density(
    x,
    adjust = adjust,
    n = n,
    from = min(x),
    to = max(x)
  )

  middle <- seq.int(2L, length(den$y) - 1L)
  peak_index <- middle[
    den$y[middle] > den$y[middle - 1L] &
      den$y[middle] > den$y[middle + 1L]
  ]
  if (length(peak_index) > 0L) {
    peak_index <- peak_index[
      den$y[peak_index] >= max(den$y) * min_relative_height
    ]
  }

  peaks <- data.frame(
    position = den$x[peak_index],
    density = den$y[peak_index],
    stringsAsFactors = FALSE
  )
  peaks$read_length <- if (isTRUE(log_scale)) {
    10^peaks$position
  } else {
    peaks$position
  }

  if (!is.null(min_distance) && nrow(peaks) > 1L) {
    selected <- integer()
    for (candidate in order(peaks$density, decreasing = TRUE)) {
      if (length(selected) == 0L || all(
        abs(peaks$read_length[[candidate]] - peaks$read_length[selected]) >=
          min_distance
      )) {
        selected <- c(selected, candidate)
      }
    }
    peaks <- peaks[selected, , drop = FALSE]
  }

  peaks <- peaks[order(peaks$read_length), , drop = FALSE]
  rownames(peaks) <- NULL

  list(
    peaks = peaks,
    density = den,
    log_scale = log_scale,
    detection_range = detection_range,
    n_reads_total = length(read_length),
    n_reads_detection = length(detection_length)
  )
}

#' Plot a read-length histogram and detected peaks
#'
#' `plot_read_length_peaks()` detects modes in a numeric vector of read lengths,
#' draws a fixed-bin-width histogram, and annotates the number and percentage of
#' reads in an asymmetric window around each retained peak.
#'
#' @inheritParams detect_length_peaks
#' @param density_n Integer passed as `n` to `detect_length_peaks()`. It controls
#'   density-grid resolution, not histogram resolution.
#' @param left_offset,right_offset Non-negative distances in base pairs defining
#'   each counting window as `[peak - left_offset, peak + right_offset]`.
#'   Reads on either boundary are included. Overlapping windows count their
#'   shared reads independently.
#' @param binwidth Positive histogram-bin width in base pairs.
#' @param title,subtitle Plot title and optional subtitle. When `subtitle` is
#'   `NULL`, a subtitle reporting peak count, window offsets, and total valid
#'   reads is generated automatically.
#' @param x_scale X-axis scale: `"auto"` uses `"log10"` when
#'   `log_scale = TRUE` and `"linear"` otherwise. The detection scale and plot
#'   scale may be chosen independently.
#' @param x_breaks Optional finite numeric vector of x-axis breaks in base-pair
#'   units. Suitable breaks are generated automatically when `NULL`.
#' @param x_limits Optional increasing numeric vector `c(lower, upper)` defining
#'   the displayed x-axis interval in base pairs. This affects only the view;
#'   it does not change peak detection, histogram calculation, or read counts.
#' @param show_peak_labels Logical. Show peak length, counting interval, read
#'   count, and percentage above each peak.
#' @param show_peak_windows Logical. Shade each peak-counting interval.
#' @param fill,colour Histogram fill and outline colors.
#' @param peak_colour Color used for peak lines and labels.
#' @param window_fill Fill color for peak-counting intervals.
#' @param fill_alpha,window_alpha Numbers between zero and one controlling
#'   histogram and window transparency.
#'
#' @return A [ggplot2::ggplot()] object. The following attributes contain the
#' underlying results:
#'
#' * `peak_result`: complete output from `detect_length_peaks()`.
#' * `peak_counts`: peak positions, counting-window boundaries, read counts,
#'   and percentages.
#' * `histogram_data`: histogram boundaries, midpoints, and counts.
#' * `parameters`: detection, window, histogram, and display settings.
#'
#' @details
#' Peak detection uses only reads inside `detection_range`, when supplied. The
#' histogram and peak-window counts use every positive, finite value in
#' `read_length`. Use `x_limits` to zoom the displayed histogram without
#' silently changing either calculation.
#'
#' When `x_scale = "log10"`, histogram bins still have equal width in base
#' pairs; only their displayed positions are transformed. X-axis labels remain
#' in base pairs.
#'
#' @examples
#' set.seed(2)
#' lengths <- round(c(
#'   rnorm(600, mean = 250, sd = 20),
#'   rnorm(1800, mean = 1500, sd = 45),
#'   rnorm(300, mean = 3000, sd = 90),
#'   100000, 110000
#' ))
#'
#' p <- plot_read_length_peaks(
#'   lengths,
#'   log_scale = TRUE,
#'   adjust = 0.8,
#'   min_relative_height = 0.05,
#'   min_distance = 300,
#'   density_n = 2048,
#'   detection_range = c(100, 5000),
#'   left_offset = 100,
#'   right_offset = 200,
#'   binwidth = 25,
#'   x_limits = c(100, 5000)
#' )
#' p
#' attr(p, "peak_counts")
#'
#' @importFrom rlang .data
#' @export
plot_read_length_peaks <- function(
    read_length,
    log_scale = TRUE,
    adjust = 1,
    min_relative_height = 0.05,
    min_distance = NULL,
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
    window_alpha = 0.12) {
  read_length <- validate_read_lengths(read_length)
  check_logical_scalar(log_scale, "log_scale")
  check_logical_scalar(show_peak_labels, "show_peak_labels")
  check_logical_scalar(show_peak_windows, "show_peak_windows")
  adjust <- read_length_validate_positive_number(adjust, "adjust")
  min_relative_height <- read_length_validate_fraction(
    min_relative_height,
    "min_relative_height"
  )
  if (!is.null(min_distance)) {
    min_distance <- read_length_validate_nonnegative_number(min_distance, "min_distance")
  }
  density_n <- validate_density_grid_size(density_n, "density_n")
  detection_range <- validate_detection_range(detection_range)
  left_offset <- read_length_validate_nonnegative_number(left_offset, "left_offset")
  right_offset <- read_length_validate_nonnegative_number(right_offset, "right_offset")
  binwidth <- read_length_validate_positive_number(binwidth, "binwidth")
  fill_alpha <- read_length_validate_fraction(fill_alpha, "fill_alpha")
  window_alpha <- read_length_validate_fraction(window_alpha, "window_alpha")
  check_scalar_character(title, "title")
  if (!is.null(subtitle)) check_scalar_character(subtitle, "subtitle")

  x_scale <- match.arg(x_scale)
  if (identical(x_scale, "auto")) {
    x_scale <- if (isTRUE(log_scale)) "log10" else "linear"
  }
  x_breaks <- validate_plot_axis_values(x_breaks, "x_breaks", x_scale)
  x_limits <- validate_plot_limits(x_limits, x_scale)

  result <- detect_length_peaks(
    read_length = read_length,
    log_scale = log_scale,
    adjust = adjust,
    min_relative_height = min_relative_height,
    min_distance = min_distance,
    n = density_n,
    detection_range = detection_range
  )
  peak_data <- result$peaks

  histogram_start <- floor(min(read_length) / binwidth) * binwidth
  histogram_end <- ceiling(max(read_length) / binwidth) * binwidth + binwidth
  histogram <- graphics::hist(
    read_length,
    breaks = seq(histogram_start, histogram_end, by = binwidth),
    plot = FALSE,
    right = FALSE,
    include.lowest = TRUE
  )
  histogram_data <- data.frame(
    xmin = histogram$breaks[-length(histogram$breaks)],
    xmax = histogram$breaks[-1L],
    midpoint = histogram$mids,
    count = histogram$counts,
    stringsAsFactors = FALSE
  )

  smallest_positive_length <- min(read_length)
  histogram_data$plot_xmin <- histogram_data$xmin
  if (identical(x_scale, "log10")) {
    histogram_data$plot_xmin <- pmax(
      histogram_data$plot_xmin,
      smallest_positive_length / 2
    )
  }

  format_read_length_number <- function(x) {
    format(round(x), scientific = FALSE, big.mark = ",", trim = TRUE)
  }

  if (nrow(peak_data) > 0L) {
    peak_data$window_start <- pmax(
      0,
      peak_data$read_length - left_offset
    )
    peak_data$window_end <- peak_data$read_length + right_offset
    peak_data$read_count <- vapply(seq_len(nrow(peak_data)), function(i) {
      sum(
        read_length >= peak_data$window_start[[i]] &
          read_length <= peak_data$window_end[[i]]
      )
    }, integer(1))
    peak_data$read_percent <- round(
      peak_data$read_count / length(read_length) * 100,
      digits = 1
    )

    peak_data$plot_window_start <- peak_data$window_start
    if (identical(x_scale, "log10")) {
      peak_data$plot_window_start <- pmax(
        peak_data$plot_window_start,
        smallest_positive_length / 2
      )
    }
    peak_data$label_y <- vapply(seq_len(nrow(peak_data)), function(i) {
      overlap <- histogram_data$xmax >= peak_data$window_start[[i]] &
        histogram_data$xmin <= peak_data$window_end[[i]]
      if (!any(overlap)) max(histogram_data$count) else
        max(histogram_data$count[overlap])
    }, numeric(1))
    label_offset <- max(histogram_data$count) * 0.05
    if (!is.finite(label_offset) || label_offset <= 0) label_offset <- 1
    peak_data$label_y <- peak_data$label_y + label_offset
    peak_data$peak_label <- paste0(
      format_read_length_number(peak_data$read_length), " bp\n[",
      format_read_length_number(peak_data$window_start), ", ",
      format_read_length_number(peak_data$window_end), "] bp\n",
      format_read_length_number(peak_data$read_count), " reads (",
      peak_data$read_percent, "%)"
    )
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Detected peaks: ", nrow(peak_data),
      "; peak window: -", left_offset, "/+", right_offset, " bp",
      "; total reads: ", format_read_length_number(length(read_length))
    )
  }
  if (is.null(x_breaks)) {
    x_breaks <- make_read_length_breaks(read_length, x_scale)
  }

  plot <- .plot_read_length_histogram(
    histogram_data = histogram_data,
    peak_data = peak_data,
    title = title,
    subtitle = subtitle,
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

  attr(plot, "peak_result") <- result
  peak_columns <- intersect(
    c("read_length", "window_start", "window_end", "read_count", "read_percent"),
    names(peak_data)
  )
  attr(plot, "peak_counts") <- peak_data[, peak_columns, drop = FALSE]
  attr(plot, "histogram_data") <- histogram_data[
    , c("xmin", "xmax", "midpoint", "count"), drop = FALSE
  ]
  attr(plot, "parameters") <- list(
    detection = list(
      log_scale = log_scale,
      adjust = adjust,
      min_relative_height = min_relative_height,
      min_distance = min_distance,
      n = density_n,
      detection_range = detection_range
    ),
    peak_window = list(
      left_offset = left_offset,
      right_offset = right_offset
    ),
    histogram = list(binwidth = binwidth),
    display = list(
      x_scale = x_scale,
      x_breaks = x_breaks,
      x_limits = x_limits
    )
  )
  plot
}

.plot_read_length_histogram <- function(
    histogram_data,
    peak_data,
    title,
    subtitle,
    x_scale,
    x_breaks,
    x_limits,
    show_peak_labels,
    show_peak_windows,
    fill,
    colour,
    peak_colour,
    window_fill,
    fill_alpha,
    window_alpha) {
  plot <- ggplot2::ggplot()
  if (isTRUE(show_peak_windows) && nrow(peak_data) > 0L) {
    plot <- plot + ggplot2::geom_rect(
      data = peak_data,
      ggplot2::aes(
        xmin = .data$plot_window_start,
        xmax = .data$window_end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = window_fill,
      alpha = window_alpha
    )
  }
  plot <- plot + ggplot2::geom_rect(
    data = histogram_data,
    ggplot2::aes(
      xmin = .data$plot_xmin,
      xmax = .data$xmax,
      ymin = 0,
      ymax = .data$count
    ),
    inherit.aes = FALSE,
    fill = fill,
    colour = colour,
    linewidth = 0.2,
    alpha = fill_alpha
  )
  if (nrow(peak_data) > 0L) {
    plot <- plot + ggplot2::geom_vline(
      data = peak_data,
      ggplot2::aes(xintercept = .data$read_length),
      inherit.aes = FALSE,
      colour = peak_colour,
      linewidth = 0.6,
      linetype = 2
    )
    if (isTRUE(show_peak_labels)) {
      plot <- plot + ggplot2::geom_label(
        data = peak_data,
        ggplot2::aes(
          x = .data$read_length,
          y = .data$label_y,
          label = .data$peak_label
        ),
        inherit.aes = FALSE,
        colour = peak_colour,
        fill = "white",
        size = 3.4,
        lineheight = 0.95,
        vjust = 0
      )
    }
  }

  plot <- plot +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Read length (bp)",
      y = "Number of reads"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.2))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey80", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_line(colour = "grey92", linewidth = 0.25),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  axis_labels <- function(x) {
    format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
  }
  if (identical(x_scale, "log10")) {
    plot <- plot + ggplot2::scale_x_log10(
      breaks = x_breaks,
      labels = axis_labels
    )
  } else {
    plot <- plot + ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = axis_labels
    )
  }
  if (!is.null(x_limits)) {
    plot <- plot + ggplot2::coord_cartesian(xlim = x_limits, clip = "on")
  } else {
    if(!is.null(x_breaks)){
      x_limits = c(x_breaks[1], tail(x_breaks, 1))
    } else {
      x_limits = c(min(histogram_data$xmin), max(histogram_data$xmax))
    }
    plot <- plot + ggplot2::coord_cartesian(xlim = x_limits, clip = "on")
  }
  plot
}

validate_read_lengths <- function(x) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop("`read_length` must be a non-empty numeric vector.", call. = FALSE)
  }
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 3L) {
    stop("At least three positive, finite read lengths are required.",
         call. = FALSE)
  }
  as.numeric(x)
}

read_length_validate_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a single positive number.", call. = FALSE)
  }
  as.numeric(x)
}

read_length_validate_nonnegative_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0) {
    stop("`", name, "` must be a single non-negative number.", call. = FALSE)
  }
  as.numeric(x)
}

read_length_validate_fraction <- function(x, name) {
  x <- read_length_validate_nonnegative_number(x, name)
  if (x > 1) {
    stop("`", name, "` must be between 0 and 1.", call. = FALSE)
  }
  x
}

validate_density_grid_size <- function(x, name = "n") {
  x <- read_length_validate_positive_number(x, name)
  if (x != as.integer(x) || x < 3L) {
    stop("`", name, "` must be an integer of at least 3.", call. = FALSE)
  }
  as.integer(x)
}

validate_detection_range <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 2L || anyNA(x) ||
      any(!is.finite(x)) || x[[1L]] < 0 || x[[1L]] >= x[[2L]]) {
    stop(
      "`detection_range` must contain two increasing, non-negative values.",
      call. = FALSE
    )
  }
  as.numeric(x)
}

validate_plot_axis_values <- function(x, name, x_scale) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must be `NULL` or a finite numeric vector.",
         call. = FALSE)
  }
  if (identical(x_scale, "log10") && any(x <= 0)) {
    stop("`", name, "` must contain positive values for a log10 axis.",
         call. = FALSE)
  }
  sort(unique(as.numeric(x)))
}

validate_plot_limits <- function(x, x_scale) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 2L || anyNA(x) ||
      any(!is.finite(x)) || x[[1L]] >= x[[2L]]) {
    stop("`x_limits` must contain two increasing finite values.", call. = FALSE)
  }
  if (identical(x_scale, "log10") && x[[1L]] <= 0) {
    stop("`x_limits` must be positive for a log10 axis.", call. = FALSE)
  }
  as.numeric(x)
}

make_read_length_breaks <- function(read_length, x_scale) {
  if (!identical(x_scale, "log10")) return(pretty(read_length))
  limits <- range(read_length, finite = TRUE)
  powers <- seq(floor(log10(limits[[1L]])), ceiling(log10(limits[[2L]])))
  breaks <- sort(unique(as.vector(outer(c(1, 2, 5), 10^powers))))
  breaks[breaks >= limits[[1L]] & breaks <= limits[[2L]]]
}
