#' Build a 16S delivery folder from wf-16s results
#'
#' `move_16s()` copies key wf-16s output files into a delivery folder and
#' generates abundance bar plots.
#'
#' @param path_result Path to a completed wf-16s result directory.
#' @param path_delivery Delivery root directory. The function creates a `16s/`
#'   subdirectory inside this path.
#' @param overwrite Logical. If `TRUE`, replace an existing `16s/` delivery
#'   directory.
#' @param tax_levels Taxonomic levels to plot.
#' @param abundance_table Filename of the wf-16s genus abundance table under
#'   `path_result`.
#' @param alignment_tables_dir Directory name of wf-16s per-barcode alignment
#'   tables under `path_result`.
#' @param figure_dir Name of the figures directory under the `16s/` delivery
#'   directory.
#' @param identification_dir Name of the copied alignment-table directory under
#'   the `16s/` delivery directory.
#' @param cutoff Minimum relative abundance kept in abundance plots.
#' @param width,height Plot width and height in inches.
#'
#' @return Invisibly returns a list with input paths, output paths, and generated
#'   ggplot objects.
#'
#' @export
move_16s <- function(path_result,
                     path_delivery = "/data/project_delivery",
                     overwrite = FALSE,
                     tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                     abundance_table = "abundance_table_genus.tsv",
                     alignment_tables_dir = "alignment_tables",
                     figure_dir = "figures",
                     identification_dir = "identification_tables",
                     cutoff = 0.01,
                     width = 12,
                     height = 6) {
  check_dir_arg(path_result, "path_result")
  check_scalar_character(path_delivery, "path_delivery")
  check_logical_scalar(overwrite, "overwrite")
  if (!is.character(tax_levels) || length(tax_levels) == 0L || anyNA(tax_levels)) {
    stop("`tax_levels` must be a non-empty character vector.", call. = FALSE)
  }
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(alignment_tables_dir, "alignment_tables_dir")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  cutoff <- validate_fraction(cutoff, "cutoff")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")

  path_result <- normalizePath(path_result, mustWork = TRUE)
  dir.create(path_delivery, recursive = TRUE, showWarnings = FALSE)
  path_delivery <- normalizePath(path_delivery, mustWork = TRUE)

  path_16s <- file.path(path_delivery, "16s")
  if (dir.exists(path_16s)) {
    if (!isTRUE(overwrite)) {
      stop(path_16s, " exists. Use `overwrite = TRUE` to replace it.",
           call. = FALSE)
    }
    unlink(path_16s, recursive = TRUE)
  }

  path_fig <- file.path(path_16s, figure_dir)
  path_identification <- file.path(path_16s, identification_dir)
  dir.create(path_fig, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_identification, showWarnings = FALSE, recursive = TRUE)

  abun <- file.path(path_result, abundance_table)
  tbls <- file.path(path_result, alignment_tables_dir)
  check_file_arg(abun, "abundance_table")
  check_dir_arg(tbls, "alignment_tables_dir")

  # abun_name <- tools::file_path_sans_ext(basename(abun))
  plots <- list()
  plot_files <- character()
  for (level in tax_levels) {
    percentage_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_percentage.png")
    )
    count_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_count.png")
    )

    plots[[paste0(level, "_percentage")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = percentage_file,
      fill = level,
      cutoff = cutoff,
      position = "fill",
      width = width,
      height = height
    )
    plots[[paste0(level, "_count")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = count_file,
      fill = level,
      cutoff = cutoff,
      position = "stack",
      width = width,
      height = height
    )
    plot_files <- c(
      plot_files,
      stats::setNames(percentage_file, paste0(level, "_percentage")),
      stats::setNames(count_file, paste0(level, "_count"))
    )
  }

  copied_tables <- copy_directory_contents(tbls, path_identification, overwrite = TRUE)
  copied_abundance <- file.copy(
    abun,
    file.path(path_16s, basename(abun)),
    overwrite = TRUE
  )
  if (!isTRUE(copied_abundance)) {
    stop("Failed to copy abundance table: ", abun, call. = FALSE)
  }

  invisible(list(
    path_result = path_result,
    path_delivery = path_delivery,
    path_16s = path_16s,
    path_fig = path_fig,
    path_identification = path_identification,
    abundance_table = file.path(path_16s, basename(abun)),
    copied_tables = copied_tables,
    plot_files = plot_files,
    plot = plots
  ))
}

#' Parse a wf-16s genus abundance table
#'
#' @param abundance_table_genus Path to `abundance_table_genus.tsv`.
#' @param format Output format: `"long"` or `"wide"`.
#' @param include_pct Logical. If `TRUE` and `format = "long"`, add per-sample
#'   relative abundance in `pct`.
#' @param levels Taxonomic level names used to split the `tax` column.
#'
#' @return A data frame.
#'
#' @export
parse_abundance_table <- function(abundance_table_genus,
                                  format = c("long", "wide"),
                                  include_pct = TRUE,
                                  levels = c(
                                    "Superkingdom", "Kingdom", "Phylum",
                                    "Class", "Order", "Family", "Genus",
                                    "Species"
                                  )[1:7]) {
  check_file_arg(abundance_table_genus, "abundance_table_genus")
  format <- match.arg(format)
  check_logical_scalar(include_pct, "include_pct")
  if (!is.character(levels) || length(levels) == 0L || anyNA(levels)) {
    stop("`levels` must be a non-empty character vector.", call. = FALSE)
  }

  abun_data <- utils::read.delim(
    abundance_table_genus,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"tax" %in% names(abun_data)) {
    stop("`abundance_table_genus` must contain a `tax` column.", call. = FALSE)
  }

  sample_cols <- setdiff(names(abun_data), "tax")
  if (length(sample_cols) == 0L) {
    stop("`abundance_table_genus` must contain at least one sample column.",
         call. = FALSE)
  }

  annotations <- split_taxonomy(abun_data$tax, levels)
  if (identical(format, "wide")) {
    return(cbind(abun_data, annotations))
  }

  count_data <- stack_abundance_table(abun_data, sample_cols)
  if (isTRUE(include_pct)) {
    totals <- stats::setNames(
      colSums(abun_data[, sample_cols, drop = FALSE], na.rm = TRUE),
      sample_cols
    )
    count_data$pct <- ifelse(
      totals[count_data$samples] > 0,
      count_data$count / totals[count_data$samples],
      NA_real_
    )
  }

  annotation_index <- match(count_data$tax, abun_data$tax)
  cbind(count_data, annotations[annotation_index, , drop = FALSE])
}

#' Plot wf-16s abundance bars
#'
#' @param abundance_table_genus Path to `abundance_table_genus.tsv`.
#' @param output Output image path.
#' @param fill Taxonomic level used for bar fill.
#' @param cutoff Minimum relative abundance kept in the plot.
#' @param position Bar position. `"fill"` shows relative abundance and
#'   `"stack"` shows read counts.
#' @param levels Taxonomic level names used to split the `tax` column.
#' @param x_angle,x_hjust,x_vjust X-axis text angle and justification.
#' @param width,height Plot width and height in inches.
#'
#' @return Invisibly returns a ggplot object.
#'
#' @export
plot_abundance_bar <- function(abundance_table_genus,
                               output = paste0(
                                 tools::file_path_sans_ext(abundance_table_genus),
                                 ".png"
                               ),
                               fill = "Order",
                               cutoff = 0.01,
                               position = c("fill", "stack"),
                               levels = c(
                                 "Superkingdom", "Kingdom", "Phylum",
                                 "Class", "Order", "Family", "Genus",
                                 "Species"
                               )[1:7],
                               x_angle = 60,
                               x_hjust = 1,
                               x_vjust = 1,
                               width = 12,
                               height = 6) {
  check_scalar_character(output, "output")
  check_scalar_character(fill, "fill")
  position <- match.arg(position)
  cutoff <- validate_fraction(cutoff, "cutoff")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")

  abun_data <- parse_abundance_table(
    abundance_table_genus = abundance_table_genus,
    format = "long",
    include_pct = TRUE,
    levels = levels
  )
  if (!fill %in% names(abun_data)) {
    stop("`fill` must be one of: ", paste(names(abun_data), collapse = ", "),
         call. = FALSE)
  }

  abun_data <- abun_data[!is.na(abun_data$pct) & abun_data$pct >= cutoff, , drop = FALSE]
  plot <- ggplot2::ggplot(
    abun_data,
    ggplot2::aes(.data$samples, .data$count, fill = .data[[fill]])
  ) +
    ggplot2::geom_col(position = position) +
    seqqc_theme(x_angle = x_angle, x_hjust = x_hjust, x_vjust = x_vjust) +
    ggplot2::xlab(NULL)

  if (identical(position, "fill")) {
    plot <- plot +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
      ggplot2::ylab("Relative abundance (%)")
  } else {
    plot <- plot + ggplot2::ylab("Reads")
  }

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = output, plot = plot, width = width, height = height)
  invisible(plot)
}

make_16s_delivery <- function() {
  stop("`make_16s_delivery()` is not implemented. Use `move_16s()` instead.",
       call. = FALSE)
}

copy_directory_contents <- function(from, to, overwrite) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(from, all.files = FALSE, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0L) {
    return(character())
  }

  copied <- file.copy(files, to, recursive = TRUE, overwrite = overwrite)
  if (any(!copied)) {
    stop("Failed to copy: ", files[which(!copied)[1L]], call. = FALSE)
  }

  file.path(to, basename(files))
}

split_taxonomy <- function(tax, levels) {
  parts <- strsplit(as.character(tax), ";", fixed = TRUE)
  max_depth <- max(lengths(parts))
  if (max_depth > length(levels)) {
    stop(
      "There are ",
      max_depth,
      " levels in `tax`, but `levels` has only ",
      length(levels),
      " values.",
      call. = FALSE
    )
  }

  mat <- matrix(NA_character_, nrow = length(parts), ncol = length(levels))
  for (i in seq_along(parts)) {
    if (length(parts[[i]]) > 0L) {
      mat[i, seq_along(parts[[i]])] <- trimws(parts[[i]])
    }
  }
  colnames(mat) <- levels
  as.data.frame(mat, stringsAsFactors = FALSE)
}

stack_abundance_table <- function(abun_data, sample_cols) {
  rows <- lapply(sample_cols, function(sample) {
    data.frame(
      tax = abun_data$tax,
      samples = sample,
      count = suppressWarnings(as.numeric(abun_data[[sample]])),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
