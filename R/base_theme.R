#' Base ggplot2 theme used by ONTools plots
#'
#' `base_theme()` builds on [ggplot2::theme_classic()] by adding subtle major
#' grid lines, removing minor grid lines, and using a transparent plot
#' background. It can be added to any ggplot with the `+` operator.
#'
#' @param base_size Positive number giving the base font size in points, passed
#'   to [ggplot2::theme_classic()].
#' @param base_family Character scalar giving the base font family, passed to
#'   [ggplot2::theme_classic()]. An empty string uses the graphics device's
#'   default font family.
#' @param title_hjust Numeric scalar giving the horizontal justification of the
#'   plot title, passed to [ggplot2::element_text()]. Common values are `0` for
#'   left alignment, `0.5` for centered alignment (the default), and `1` for
#'   right alignment.
#' @param title_size Positive number giving the plot-title font size in points.
#' @param x_angle Numeric scalar giving the rotation angle of x-axis labels in
#'   degrees.
#' @param x_hjust,x_vjust Numeric scalars giving the horizontal and vertical
#'   justification of x-axis labels.
#' @param x_size Positive number giving the x-axis label font size in points.
#'
#' @return A complete ggplot2 theme object. The theme uses:
#'
#' * major grid lines with color `"grey85"` and width `0.2`;
#' * no minor grid lines;
#' * a transparent plot background without a border;
#' * plot-title alignment controlled by `title_hjust`.
#'
#' @examples
#' p <- ggplot2::ggplot(
#'   data.frame(x = 1:5, y = c(2, 5, 3, 7, 6)),
#'   ggplot2::aes(x = x, y = y)
#' ) +
#'   ggplot2::geom_line() +
#'   base_theme()
#' p
#'
#' # Increase text size and select a font family.
#' p + base_theme(base_size = 14, base_family = "sans")
#'
#' # Left-align the plot title.
#' p +
#'   ggplot2::labs(title = "Example title") +
#'   base_theme(title_hjust = 0)
#'
#' @export
base_theme <- function(base_size = 11, base_family = "", title_hjust = 0.5, title_size = 11,
                       x_angle = 20, x_hjust = 0.5, x_vjust = 0.5, x_size = 11) {
  if (!is.numeric(base_size) || length(base_size) != 1L ||
      is.na(base_size) || !is.finite(base_size) || base_size <= 0) {
    stop("`base_size` must be a single positive number.", call. = FALSE)
  }
  if (!is.character(base_family) || length(base_family) != 1L ||
      is.na(base_family)) {
    stop("`base_family` must be a single character string.", call. = FALSE)
  }

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      plot.title = ggplot2::element_text(hjust = title_hjust, size = title_size),
      axis.text.x = ggplot2::element_text(
        angle = x_angle,
        hjust = x_hjust,
        vjust = x_vjust,
        size = x_size
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5))
    )
}
