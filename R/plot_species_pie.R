#' Plot a species pie chart from a wf-16s alignment table
#'
#' `plot_species_pie()` reads one `barcode*-alignment-stats.tsv` table generated
#' by wf-16s and draws a pie chart of species composition. Low-abundance species
#' can be grouped into an `"Others"` slice.
#'
#' @param alignment_table Path to a wf-16s alignment stats TSV file, or a data
#'   frame already read from such a file.
#' @param showPercent Logical. Backward-compatible shortcut. If `TRUE`, plot
#'   values from `percent_col`; if `FALSE`, plot read counts from `reads_col`.
#'   Ignored when `value` is supplied.
#' @param otherPct Relative abundance cutoff for grouping rare species into
#'   `"Others"`. Values are fractions, so `0.01` means 1%.
#' @param title Optional plot title.
#' @param value Plot value mode: `"reads"` uses read counts and `"percent"` uses
#'   `percent_col`.
#' @param species_col Column containing species names. If `NULL`, common wf-16s
#'   names such as `species`, `genuspecies`, and `genus species` are detected.
#' @param reads_col Column containing read counts. If `NULL`, common names such
#'   as `number of reads` and `number.of.reads` are detected.
#' @param percent_col Column containing per-reference read percentage. If
#'   `NULL`, common names such as `pcreads` and `pc reads` are detected.
#' @param top_n Optional maximum number of species to show before grouping the
#'   rest into `"Others"`.
#' @param other_label Label used for grouped rare species.
#' @param label_slices Logical. If `TRUE`, add text labels to slices.
#' @param label_min_pct Minimum fraction required for a slice label.
#' @param output Optional output image path. When supplied, the plot is saved
#'   with [ggplot2::ggsave()].
#' @param width,height Plot width and height in inches when `output` is supplied.
#'
#' @return A ggplot object. The aggregated plotting data is available from
#'   `attr(plot, "plot_data")`.
#'
#' @export
plot_species_pie <- function(alignment_table,
                             showPercent = FALSE,
                             otherPct = 0.01,
                             title = NULL,
                             value = NULL,
                             species_col = NULL,
                             reads_col = NULL,
                             percent_col = NULL,
                             top_n = NULL,
                             other_label = "Others",
                             label_slices = FALSE,
                             label_min_pct = 0.03,
                             output = NULL,
                             width = 7,
                             height = 6) {
  stat <- read_alignment_stats_table(alignment_table)
  if (is.null(value)) {
    value <- if (isTRUE(showPercent)) "percent" else "reads"
  } else {
    value <- match.arg(value, c("reads", "percent"))
  }

  otherPct <- validate_fraction(otherPct, "otherPct")
  label_min_pct <- validate_fraction(label_min_pct, "label_min_pct")
  if (!is.null(top_n)) top_n <- validate_positive_integer(top_n, "top_n")
  check_scalar_character(other_label, "other_label")
  check_logical_scalar(label_slices, "label_slices")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")
  if (!is.null(title)) check_scalar_character(title, "title")
  if (!is.null(output)) check_scalar_character(output, "output")

  species_col <- resolve_alignment_stats_col(
    stat,
    species_col,
    c("species", "genuspecies", "genus species", "genus_species", "Species"),
    "species_col"
  )
  reads_col <- resolve_alignment_stats_col(
    stat,
    reads_col,
    c("number of reads", "number.of.reads", "number_of_reads", "reads", "n_reads"),
    "reads_col"
  )
  if (identical(value, "percent")) {
    percent_col <- resolve_alignment_stats_col(
      stat,
      percent_col,
      c("pcreads", "pc reads", "pc.reads", "percent reads", "percent_reads"),
      "percent_col"
    )
  }

  stat$Species <- as.character(stat[[species_col]])
  stat$Species[is.na(stat$Species) | !nzchar(stat$Species)] <- "Unknown"
  stat$reads_value <- suppressWarnings(as.numeric(stat[[reads_col]]))
  if (identical(value, "percent")) {
    stat$plot_value <- suppressWarnings(as.numeric(stat[[percent_col]]))
  } else {
    stat$plot_value <- stat$reads_value
  }
  stat <- stat[is.finite(stat$plot_value) & stat$plot_value > 0, , drop = FALSE]
  if (nrow(stat) == 0L) {
    stop("No positive finite values were found for plotting.", call. = FALSE)
  }

  plot_data <- aggregate(
    list(value = stat$plot_value, reads = stat$reads_value),
    by = list(Species = stat$Species),
    FUN = sum,
    na.rm = TRUE
  )
  plot_data <- plot_data[order(plot_data$value, decreasing = TRUE), , drop = FALSE]
  plot_data$fraction <- plot_data$value / sum(plot_data$value)

  show_as_other <- plot_data$fraction < otherPct
  if (!is.null(top_n) && nrow(plot_data) > top_n) {
    show_as_other[seq.int(top_n + 1L, nrow(plot_data))] <- TRUE
  }
  plot_data$Species <- ifelse(show_as_other, other_label, plot_data$Species)
  plot_data <- aggregate(
    list(value = plot_data$value, reads = plot_data$reads),
    by = list(Species = plot_data$Species),
    FUN = sum,
    na.rm = TRUE
  )
  plot_data <- plot_data[order(plot_data$value, decreasing = TRUE), , drop = FALSE]
  if (other_label %in% plot_data$Species) {
    plot_data <- rbind(
      plot_data[plot_data$Species != other_label, , drop = FALSE],
      plot_data[plot_data$Species == other_label, , drop = FALSE]
    )
  }
  plot_data$fraction <- plot_data$value / sum(plot_data$value)
  plot_data$label <- paste0(
    plot_data$Species,
    "\n",
    sprintf("%.1f%%", plot_data$fraction * 100)
  )
  plot_data$Species <- factor(plot_data$Species, levels = plot_data$Species)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = "", y = .data$value, fill = .data$Species)
  ) +
    ggplot2::geom_col(width = 1, color = "white", linewidth = 0.2) +
    ggplot2::coord_polar(theta = "y", direction = -1, clip = "off") +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "Species"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key.size = grid::unit(0.4, "cm")
    )

  if (isTRUE(label_slices)) {
    label_data <- plot_data[plot_data$fraction >= label_min_pct, , drop = FALSE]
    plot <- plot +
      ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = .data$label),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3
      )
  }

  if (!is.null(output)) {
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = output, plot = plot, width = width, height = height)
  }

  attr(plot, "plot_data") <- plot_data
  plot
}

read_alignment_stats_table <- function(alignment_table) {
  if (is.character(alignment_table) && length(alignment_table) == 1L) {
    check_file_arg(alignment_table, "alignment_table")
    return(utils::read.delim(
      alignment_table,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  if (is.data.frame(alignment_table)) {
    return(alignment_table)
  }

  stop("`alignment_table` must be a file path or a data frame.", call. = FALSE)
}

resolve_alignment_stats_col <- function(data, col, candidates, arg_name) {
  if (!is.null(col)) {
    check_scalar_character(col, arg_name)
    if (!col %in% names(data)) {
      stop("`", arg_name, "` was not found in `alignment_table`: ", col,
           call. = FALSE)
    }
    return(col)
  }

  normalized_names <- normalize_alignment_stats_names(names(data))
  normalized_candidates <- normalize_alignment_stats_names(candidates)
  idx <- match(normalized_candidates, normalized_names, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx) > 0L) {
    return(names(data)[idx[[1L]]])
  }

  stop(
    "Could not detect `", arg_name, "`. Available columns are: ",
    paste(names(data), collapse = ", "),
    call. = FALSE
  )
}

normalize_alignment_stats_names <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", x))
}
