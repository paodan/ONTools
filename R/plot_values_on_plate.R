#' Plot values on a 96- or 384-well plate layout
#'
#' `plot_values_on_plate()` draws a complete plate map and overlays values from
#' a user-supplied table by well row and column. Missing wells are still shown
#' and are filled with light grey.
#'
#' @param data A data frame containing at least the row, column, and value
#'   columns.
#' @param valueName,rowName,columnName Column names in `data` for the value,
#'   plate row, and plate column.
#' @param plateType Plate layout. `"96"` creates an 8 x 12 plate and `"384"`
#'   creates a 16 x 24 plate.
#' @param palette Colours used for the fill scale. For discrete values, named
#'   colours are matched to value levels when possible. For numeric values, the
#'   colours are used as a gradient.
#' @param show.legend Logical. Show the fill legend.
#' @param title Optional plot title.
#' @param x_expansion,y_expansion Expansion added to the x and y scales.
#' @param x_gap,y_gap Extra limits added around the plate grid.
#' @param margins Plot margins in points, ordered top, right, bottom, left.
#' @param circle_size Size of each well marker.
#' @param coord_fix Logical. Use a fixed coordinate ratio.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' plate <- data.frame(
#'   row = c("A", "B", "H"),
#'   column = c(1, 2, 12),
#'   value = c("A", "B", "D")
#' )
#' plot_values_on_plate(plate, title = "Example plate")
plot_values_on_plate <- function(data, valueName = "value", rowName = "row",
                                 columnName = "column",
                                 plateType = c("96", "384"),
                                 palette = c(A = "red", B = "green", C = "blue", D = "orange"),
                                 show.legend = FALSE,
                                 title = NULL,
                                 x_expansion = c(0.6, 0.8),
                                 y_expansion = c(0.6, 0.8),
                                 x_gap = c(0.5, 0.5),
                                 y_gap = c(0.5, 0.5),
                                 margins = c(14, 14, 14, 14),
                                 circle_size = 10,
                                 coord_fix = TRUE) {
  plateType <- match.arg(plateType)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  missing_cols <- setdiff(c(columnName, rowName, valueName), colnames(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required column(s) in `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  data$well <- paste0(data[[rowName]], data[[columnName]])

  if (identical(plateType, "96")) {
    rowN <- 8
    columnN <- 12
  } else if (identical(plateType, "384")) {
    rowN <- 16
    columnN <- 24
  } else {
    stop("Unknown `plateType`.", call. = FALSE)
  }
  cells <- as.data.frame(expand.grid(row = LETTERS[seq_len(rowN)], column = seq_len(columnN)))
  rownames(cells) <- cells$well <- paste0(cells$row, cells$column)
  cells$row_num <- match(cells$row, LETTERS[seq_len(rowN)])

  df <- merge(cells, data[, c("well", valueName), drop = FALSE],
              by = "well", all.x = TRUE, sort = FALSE)
  df <- df[match(cells$well, df$well), , drop = FALSE]

  g <- ggplot2::ggplot(df, ggplot2::aes(column, row_num)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data[[valueName]]),
      shape = 21,
      size = circle_size,
      stroke = 0.8,
      show.legend = show.legend
    ) +
    ggplot2::geom_text(ggplot2::aes(column - 0.45, row_num - 0.35, label = well), size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(column, row_num, label = .data[[valueName]]), size = 3.5, na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(columnN),
      position = "top",
      expand = ggplot2::expansion(add = x_expansion),
      limits = c(1 - y_gap[1], columnN + y_gap[2])
    ) +
    ggplot2::scale_y_reverse(
      breaks = seq_len(rowN),
      labels = LETTERS[seq_len(rowN)],
      expand = ggplot2::expansion(add = y_expansion),
      limits = c(rowN + y_gap[2], 1 - y_gap[1])
    )
  if (isTRUE(coord_fix)) {
    g <- g + ggplot2::coord_fixed()
  }
  g <- g +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.5, colour = palette[1]),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      plot.margin = ggplot2::margin(margins[1], margins[2], margins[3], margins[4]),
      panel.border = ggplot2::element_rect(
        color = "grey30",
        fill = NA,
        linewidth = 1.2
      )
    )
  if (is.character(df[[valueName]]) || is.factor(df[[valueName]])) {
    g <- g +
      ggplot2::scale_fill_manual(values = palette, na.value = "#dddddd")
  } else {
    g <- g +
      ggplot2::scale_fill_gradientn(colours = unname(palette), na.value = "#dddddd")
  }
  g <- g + ggplot2::ggtitle(title)
  return(g)
}

#' Plot gel-quality categories on a 96-well plate
#'
#' `plot_gel_quality_on_plate_96()` creates one PDF per gel-quality category
#' (`A`, `B`, `C`, and `D`). Each plot highlights wells in the selected
#' category and renders all other wells in grey.
#'
#' @inheritParams plot_values_on_plate
#' @param path Output directory for `A.pdf`, `B.pdf`, `C.pdf`, and `D.pdf`.
#'
#' @return Invisibly returns a named list of `ggplot` objects.
#' @export
plot_gel_quality_on_plate_96 <- function(data, device = c("png", "pdf"),
                                         path = "./",
                                         valueName = "value", rowName = "row",
                                         columnName = "column",
                                         palette = c(A = "red", B = "green", C = "blue", D = "orange")) {
  device = match.arg(device)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  p <- list()
  for (mi in c("A", "B", "C", "D")) {
    data_mi <- data[data[[valueName]] == mi & !is.na(data[[valueName]]), , drop = FALSE]
    p[[mi]] <- plot_values_on_plate(
      data_mi,
      palette = c(palette[mi], "#dddddd"),
      valueName = valueName,
      rowName = rowName,
      columnName = columnName,
      plateType = "96",
      title = paste0("Gel quality: ", mi),
      x_expansion = c(1, 1),
      y_expansion = c(0.7, 0.8),
      x_gap = c(0.2, 0.2),
      y_gap = c(0.5, 0.5),
      margins = c(10, 10, 10, 6),
      circle_size = 10,
      coord_fix = TRUE
    )
    ggplot2::ggsave(
      paste0(mi, ".", device),
      plot = p[[mi]],
      path = path,
      width = 15.3,
      height = 10.3,
      units = "cm",
      dpi = 300
    )
  }
  return(invisible(p))
}


#' Plot gel-quality categories on a tip-plate layout
#'
#' `plot_gel_quality_on_plate_tip()` is a gel-quality convenience wrapper with
#' plot dimensions and spacing tuned for tip-plate output.
#'
#' @inheritParams plot_gel_quality_on_plate_96
#'
#' @return Invisibly returns a named list of `ggplot` objects.
#' @export
plot_gel_quality_on_plate_tip <- function(data, device = c("png", "pdf"),
                                          path = "./",
                                          valueName = "value", rowName = "row",
                                          columnName = "column",
                                          palette = c(A = "red", B = "green", C = "blue", D = "orange")) {

  device = match.arg(device)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  p <- list()
  for (mi in c("A", "B", "C", "D")) {
    data_mi <- data[data[[valueName]] == mi & !is.na(data[[valueName]]), , drop = FALSE]
    p[[mi]] <- plot_values_on_plate(
      data_mi,
      palette = c(palette[mi], "#dddddd"),
      valueName = valueName,
      rowName = rowName,
      columnName = columnName,
      plateType = "96",
      title = paste0("Gel quality: ", mi),
      x_expansion = c(1.2, 1.2),
      y_expansion = c(0.47, 0.47),
      x_gap = c(0.2, 0.2),
      y_gap = c(0.5, 0.5),
      margins = c(0, 0, 0, 0),
      circle_size = 7,
      coord_fix = FALSE
    )
    ggplot2::ggsave(
      paste0(mi, ".", device),
      plot = p[[mi]],
      path = path,
      width = 12.5 + 0.5,
      height = 7.8 + 1.2,
      units = "cm",
      dpi = 300
    )
  }
  return(invisible(p))
}
