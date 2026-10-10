#' Calculate the Gini coefficient of sequencing depth
#'
#' `gini_coefficient()` measures inequality among non-negative values. In a
#' sequencing project, the values are typically read counts or coverage depths
#' for individual samples. A value of zero indicates identical depths, while a
#' larger value indicates a less uniform allocation of sequencing output.
#'
#' @param x A numeric vector of non-negative values. Missing, non-finite, and
#'   negative values are not allowed.
#' @param normalize Logical. If `TRUE`, apply the finite-sample correction
#'   `n / (n - 1)`, so the maximum possible coefficient is one even when the
#'   vector contains a finite number of observations.
#'
#' @return A numeric scalar. `NA_real_` is returned when fewer than two values
#'   are supplied or when all values are zero.
#'
#' @details
#' The Gini coefficient summarizes the area between the Lorenz curve and the
#' line of perfect equality. It is invariant to multiplication by a positive
#' constant, so projects measured on different absolute depth scales can be
#' compared when their sample inclusion criteria and depth definitions are the
#' same. Zero-depth samples should be retained because removing them makes a
#' project appear artificially uniform.
#'
#' When different samples have different intended sequencing depths, calculate
#' the coefficient from `observed_depth / target_depth` instead of from the raw
#' depths.
#'
#' @examples
#' depth <- c(850, 920, 980, 1000, 1050, 1120)
#' gini_coefficient(depth)
#'
#' # Multiplying every depth by the same value does not change the result.
#' gini_coefficient(depth * 10)
#'
#' # A deliberately uneven project has a larger coefficient.
#' gini_coefficient(c(100, 100, 100, 1000))
#'
#' @seealso [summarize_depth_uniformity()], [plot_depth_lorenz()]
#' @export
gini_coefficient <- function(x, normalize = TRUE) {
  x <- as.numeric(x)

  if (any(!is.finite(x))) {
    stop("`x` contains NA, NaN, or infinite values.")
  }

  if (any(x < 0)) {
    stop("`x` must contain non-negative values.")
  }

  n <- length(x)

  if (n < 2L || sum(x) == 0) {
    return(NA_real_)
  }

  x <- sort(x)

  gini <- (
    2 * sum(seq_len(n) * x) /
      (n * sum(x))
  ) - (n + 1) / n

  # 将有限样本下的最大值校正为1
  if (normalize) {
    gini <- gini * n / (n - 1)
  }

  gini
}


#' Summarize sequencing-depth uniformity across samples
#'
#' `summarize_depth_uniformity()` calculates complementary measures of the
#' consistency of sample-level sequencing depths. It reports conventional and
#' robust dispersion measures, inequality, percentile ratios, and the
#' proportion of samples close to the project median.
#'
#' @param depth A vector coercible to numeric containing sample-level read
#'   counts or sequencing depths. Non-finite values (`NA`, `NaN`, `Inf`, and
#'   `-Inf`) are removed. Remaining values must be non-negative.
#'
#' @return A one-row data frame with the following columns:
#'
#' * `samples`: number of finite depth values;
#' * `total_depth`, `mean_depth`, `median_depth`, `minimum_depth`, and
#'   `maximum_depth`: descriptive statistics;
#' * `cv`: coefficient of variation, `sd(depth) / mean(depth)`;
#' * `robust_cv`: median absolute deviation divided by the median;
#' * `iqr_over_median`: interquartile range divided by the median;
#' * `gini`: finite-sample-normalized Gini coefficient;
#' * `p95_p5`: ratio of the 95th to the 5th percentile. This is `Inf` when the
#'   5th percentile is zero;
#' * `zero_depth_percent`: percentage of samples with zero depth;
#' * `within_half_to_two_fold_percent`: percentage between 0.5 and 2 times the
#'   median, inclusive;
#' * `within_20_percent`: percentage between 0.8 and 1.2 times the median,
#'   inclusive.
#'
#' @details
#' No single statistic captures every type of non-uniformity. The ordinary CV
#' is sensitive to extreme samples, whereas `robust_cv` and
#' `iqr_over_median` describe the central part of the distribution. The Gini
#' coefficient captures overall inequality, and the range percentages are
#' often easier to interpret operationally.
#'
#' Metrics whose denominator is the mean or median are returned as `NA` when
#' that denominator is zero. Projects should only be compared when the same
#' filtering, barcode assignment, zero-depth sample policy, and definition of
#' depth have been used.
#'
#' @examples
#' depth <- c(950, 1020, 870, 1100, 980, 1040, 760, 1250)
#' summarize_depth_uniformity(depth)
#'
#' # Zero-depth samples are retained in the summary.
#' summarize_depth_uniformity(c(depth, 0, 0))
#'
#' # Compare two projects using rows in a single table.
#' project_a <- summarize_depth_uniformity(c(900, 950, 1000, 1050, 1100))
#' project_b <- summarize_depth_uniformity(c(300, 500, 1000, 1500, 2700))
#' rbind(project_a = project_a, project_b = project_b)
#'
#' @seealso [gini_coefficient()], [plot_depth_lorenz()]
#' @export
summarize_depth_uniformity <- function(depth) {
  depth <- as.numeric(depth)
  depth <- depth[is.finite(depth)]

  if (!length(depth)) {
    stop("No valid depth values.")
  }

  if (any(depth < 0)) {
    stop("Depth values cannot be negative.")
  }

  mean_depth <- mean(depth)
  median_depth <- stats::median(depth)

  q <- stats::quantile(
    depth,
    probs = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
    names = FALSE
  )

  data.frame(
    samples = length(depth),
    total_depth = sum(depth),
    mean_depth = mean_depth,
    median_depth = median_depth,
    minimum_depth = min(depth),
    maximum_depth = max(depth),

    cv = if (mean_depth > 0) {
      stats::sd(depth) / mean_depth
    } else {
      NA_real_
    },

    robust_cv = if (median_depth > 0) {
      stats::mad(depth) / median_depth
    } else {
      NA_real_
    },

    iqr_over_median = if (median_depth > 0) {
      stats::IQR(depth) / median_depth
    } else {
      NA_real_
    },

    gini = gini_coefficient(depth),

    p95_p5 = if (q[[1]] > 0) {
      q[[6]] / q[[1]]
    } else {
      Inf
    },

    zero_depth_percent = mean(depth == 0) * 100,

    within_half_to_two_fold_percent =
      if (median_depth > 0) {
        mean(
          depth >= 0.5 * median_depth &
            depth <= 2 * median_depth
        ) * 100
      } else {
        NA_real_
      },

    within_20_percent =
      if (median_depth > 0) {
        mean(
          depth >= 0.8 * median_depth &
            depth <= 1.2 * median_depth
        ) * 100
      } else {
        NA_real_
      }
  )
}


#' Plot a Lorenz curve of sample sequencing depths
#'
#' `plot_depth_lorenz()` orders sample depths from smallest to largest and
#' plots cumulative sample proportion against cumulative sequencing-depth
#' proportion. The diagonal represents perfect uniformity. By default, a set
#' of depth-uniformity statistics is displayed inside the plot.
#'
#' @param depth A vector coercible to numeric containing non-negative
#'   sample-level read counts or sequencing depths. Non-finite values are
#'   removed. At least one positive value is required.
#' @param title `NULL` or a character scalar used as the plot title.
#' @param show_summary Logical. If `TRUE`, add a label containing the sample
#'   count, median depth, CV, robust CV, Gini coefficient, P95/P5 ratio,
#'   percentage within 0.5--2 times the median, and percentage at zero depth.
#' @param digits Non-negative integer giving the number of decimal places used
#'   for CV, robust CV, Gini, and P95/P5 in the annotation.
#' @param annotation_x,annotation_y Numeric coordinates for the upper-left
#'   anchor of the summary label. Both axes use proportions ranging from zero
#'   to one.
#'
#' @return A [ggplot2::ggplot()] object. Two attributes are attached:
#'
#' * `depth_summary`: the one-row result from
#'   [summarize_depth_uniformity()];
#' * `lorenz_data`: the cumulative proportions used to draw the curve.
#'
#' @details
#' A curve close to the diagonal indicates similar depths among samples. A
#' curve bending toward the lower-right corner indicates that a relatively
#' small subset of samples accounts for a large fraction of all sequencing
#' reads. The Gini coefficient printed in the annotation numerically
#' summarizes this departure from equality.
#'
#' @section Statistics shown in the plot:
#'
#' When `show_summary = TRUE`, the label contains the following statistics:
#'
#' * `Samples`: the number of finite depth values included in the analysis.
#'   Each value is assumed to represent one sample.
#' * `Median depth`: the median sample depth. It describes the typical
#'   sequencing output and is less sensitive to extremely deep samples than
#'   the arithmetic mean.
#' * `CV`: the coefficient of variation, calculated as the sample standard
#'   deviation divided by the mean depth. Smaller values indicate more similar
#'   depths, but CV is sensitive to outliers and can be unstable when the mean
#'   depth is close to zero.
#' * `Robust CV`: [stats::mad()] divided by the median depth. It describes the
#'   relative spread of the central part of the distribution and is less
#'   affected by unusually shallow or deep samples. Smaller values indicate
#'   greater uniformity.
#' * `Gini`: the finite-sample-normalized Gini coefficient returned by
#'   [gini_coefficient()]. Zero indicates equal depth for every sample; values
#'   approaching one indicate increasing concentration of reads in a small
#'   number of samples.
#' * `P95/P5`: the 95th-percentile depth divided by the 5th-percentile depth.
#'   A value close to one indicates similar depths across the central 90
#'   percent of samples. Larger values indicate greater disparity. The value
#'   is infinite when the 5th percentile is zero.
#' * `Within 0.5 - 2x median`: the percentage of samples whose depths are
#'   between one-half and two times the project median, inclusive. Larger
#'   percentages indicate better uniformity.
#' * `Zero depth`: the percentage of included samples with depth equal to zero.
#'   Smaller percentages are preferable; such samples must not be removed when
#'   projects are compared.
#'
#' These measures are descriptive and do not have universal pass/fail cutoffs.
#' Comparisons are most meaningful when projects use the same sequencing-depth
#' definition, filtering rules, sample inclusion criteria, and target-depth
#' design. Median depth should be considered alongside the uniformity metrics:
#' a project can be highly uniform while still having insufficient depth.
#'
#' For comparisons between projects, retain samples with zero depth and use
#' identical preprocessing and filtering rules. If samples have unequal target
#' depths, supply `observed_depth / target_depth` rather than raw depths.
#'
#' @examples
#' set.seed(2026)
#' depth <- round(stats::rlnorm(120, meanlog = log(1000), sdlog = 0.25))
#'
#' plot_depth_lorenz(
#'   depth,
#'   title = "Example project"
#' )
#'
#' # Draw the curve without the statistics label.
#' p <- plot_depth_lorenz(depth, show_summary = FALSE)
#' p
#'
#' # Summary statistics and plotting coordinates remain available.
#' attr(p, "depth_summary")
#' head(attr(p, "lorenz_data"))
#'
#' @seealso [gini_coefficient()], [summarize_depth_uniformity()]
#' @importFrom rlang .data
#' @export
plot_depth_lorenz <- function(
    depth,
    title = NULL,
    show_summary = TRUE,
    digits = 3,
    annotation_x = 0.04,
    annotation_y = 0.96) {

  depth <- as.numeric(depth)
  depth <- depth[is.finite(depth)]

  if (!length(depth)) {
    stop("`depth` does not contain valid numeric values.")
  }

  if (any(depth < 0)) {
    stop("`depth` cannot contain negative values.")
  }

  if (sum(depth) == 0) {
    stop("The sum of `depth` must be greater than zero.")
  }

  depth <- sort(depth)
  n <- length(depth)

  lorenz_data <- data.frame(
    sample_fraction = c(
      0,
      seq_len(n) / n
    ),
    depth_fraction = c(
      0,
      cumsum(depth) / sum(depth)
    )
  )

  depth_summary <- summarize_depth_uniformity(depth)

  p <- ggplot2::ggplot(
    lorenz_data,
    ggplot2::aes(
      x = .data$sample_fraction,
      y = .data$depth_fraction
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "grey60",
      linetype = 2,
      linewidth = 0.5
    ) +
    ggplot2::geom_line(
      colour = "steelblue",
      linewidth = 1
    ) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5
      )
    ) +
    ggplot2::labs(
      title = title,
      x = "Cumulative proportion of samples",
      y = "Cumulative proportion of sequencing depth"
    )

  if (isTRUE(show_summary)) {
    p95_p5_text <- if (
      is.finite(depth_summary$p95_p5)
    ) {
      format(
        round(depth_summary$p95_p5, digits),
        nsmall = digits
      )
    } else {
      "Inf"
    }

    summary_label <- paste0(
      "Samples: ", depth_summary$samples,
      "\nMedian depth: ",
      format(
        depth_summary$median_depth,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      "\nCV: ",
      format(
        round(depth_summary$cv, digits),
        nsmall = digits
      ),
      "\nRobust CV: ",
      format(
        round(depth_summary$robust_cv, digits),
        nsmall = digits
      ),
      "\nGini: ",
      format(
        round(depth_summary$gini, digits),
        nsmall = digits
      ),
      "\nP95/P5: ",
      p95_p5_text,
      "\nWithin 0.5 \u2013 2x median: ",
      round(
        depth_summary$
          within_half_to_two_fold_percent,
        1
      ),
      "%",
      "\nZero depth: ",
      round(
        depth_summary$zero_depth_percent,
        1
      ),
      "%"
    )

    p <- p +
      ggplot2::annotate(
        geom = "label",
        x = annotation_x,
        y = annotation_y,
        label = summary_label,
        hjust = 0,
        vjust = 1,
        size = 3.3,
        lineheight = 1.05,
        colour = "black",
        fill = scales::alpha("white", 0.85),
        linewidth = 0.2
      )
  }

  attr(p, "depth_summary") <- depth_summary
  attr(p, "lorenz_data") <- lorenz_data

  p
}
