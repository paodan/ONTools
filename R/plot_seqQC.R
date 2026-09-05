#' Plot sequencing-summary QC for an ONT run
#'
#' `plot_seqQC()` reads a MinKNOW or dorado `sequencing_summary*.txt` file and
#' creates read-count and read-length QC plots by barcode or sample alias.
#'
#' @param filePath Path to a tab-delimited sequencing summary file.
#' @param runName Run name used in output filenames. When `NULL`, the input file
#'   basename without extension is used.
#' @param device Graphics device for saving plots. Supported values are
#'   `"pdf"`, `"png"`, `"svg"`, `"jpg"`, and `"jpeg"`. Use `NULL` to return
#'   plots without writing files.
#' @param unclassified Label used for unclassified reads.
#' @param plot_unclassified Logical. If `FALSE`, remove reads whose sample label
#'   equals `unclassified` before plotting and do not include an unclassified
#'   panel or bar.
#' @param barcodes Integer barcode numbers to include and order. Use `NULL` to
#'   use sample labels as they appear in the sequencing summary file.
#' @param barcode_digits Minimum number of digits used for barcode labels.
#' @param out_dir Parent output directory for plot files.
#' @param sample_col Column containing barcode or sample labels. Defaults to
#'   `"alias"`, matching the existing ONTools workflow.
#' @param length_col Column containing read lengths. Defaults to
#'   `"sequence_length_template"`.
#' @param min_read_length Optional minimum read length to keep before plotting.
#' @param max_read_length Optional maximum read length to keep before plotting.
#' @param filename_suffix Optional suffix added to output filenames before the
#'   file extension. When `NULL`, a suffix is generated from length filters, for
#'   example `"len1200-1800"`.
#' @param len_read_width,len_read_height Width and height in inches for the
#'   combined read-length distribution plot. When either value is `NULL`, it is
#'   calculated from the number of facets and the corresponding unit size.
#' @param len_read_unit_width,len_read_unit_height Width and height in inches
#'   used per facet when calculating the combined read-length plot size.
#' @param len_read_ncol Number of facet columns in the combined read-length
#'   distribution plot.
#' @param save_sample_len_plots Logical. If `TRUE`, save one read-length
#'   distribution plot per sample.
#' @param sample_len_subdir Subdirectory under `out_dir/device` where per-sample
#'   read-length plots are saved when `sample_len_dir_mode = "subdir"`.
#' @param sample_len_dir_mode Where per-sample read-length plots are saved.
#'   `"subdir"` saves all plots under `out_dir/device/sample_len_subdir`.
#'   `"fastq_pass_trim"` saves each plot under
#'   `fastq_pass_trim_dir/<sample>/`.
#' @param fastq_pass_trim_dir Directory containing per-barcode fastq_pass_trim
#'   folders. Used when `sample_len_dir_mode = "fastq_pass_trim"`.
#' @param sample_len_width,sample_len_height Width and height in inches for each
#'   per-sample read-length plot.
#'
#' @return Invisibly returns a list with the read-count plot (`numRead`), the
#'   read-length plot (`lenRead`), demultiplexing recovery percentage
#'   (`recovery`), per-sample read counts (`read_counts`), total reads before
#'   filtering (`n_reads_raw`), total reads after filtering (`n_reads_filtered`),
#'   reads removed because their sample labels were not present in the requested
#'   barcode levels (`n_reads_unmatched_sample`), filters applied (`filters`),
#'   combined read-length plot size (`len_read_size`), paths for combined plots
#'   (`files`), and paths for per-sample read-length plots (`sample_len_files`).
#'
#' @examples
#' summary_file <- tempfile(fileext = ".txt")
#' write.table(
#'   data.frame(
#'     alias = c("barcode01", "barcode01", "unclassified"),
#'     sequence_length_template = c(1200, 800, 500)
#'   ),
#'   summary_file,
#'   sep = "\t",
#'   row.names = FALSE,
#'   quote = FALSE
#' )
#'
#' qc <- plot_seqQC(summary_file, runName = "example", device = NULL, barcodes = 1)
#' qc$recovery
#'
#' filtered_qc <- plot_seqQC(
#'   summary_file,
#'   runName = "example_filtered",
#'   device = NULL,
#'   barcodes = 1,
#'   min_read_length = 1000,
#'   max_read_length = 1800
#' )
#' filtered_qc$n_reads_filtered
#'
#' @importFrom rlang .data
#' @export
plot_seqQC <- function(filePath,
                       runName = NULL,
                       device = "pdf",
                       unclassified = "unclassified",
                       plot_unclassified = TRUE,
                       barcodes = 1:96,
                       barcode_digits = 2,
                       out_dir = "figs",
                       sample_col = "alias",
                       length_col = "sequence_length_template",
                       min_read_length = NULL,
                       max_read_length = NULL,
                       filename_suffix = NULL,
                       len_read_width = NULL,
                       len_read_height = NULL,
                       len_read_unit_width = 2,
                       len_read_unit_height = 2,
                       len_read_ncol = 12,
                       save_sample_len_plots = TRUE,
                       sample_len_subdir = "seqLength_by_sample",
                       sample_len_dir_mode = c("subdir", "fastq_pass_trim"),
                       fastq_pass_trim_dir = "fastq_pass_trim",
                       sample_len_width = 5,
                       sample_len_height = 4) {
  if (!is.character(filePath) || length(filePath) != 1L || is.na(filePath)) {
    stop("`filePath` must be a single file path.", call. = FALSE)
  }

  if (!file.exists(filePath)) {
    stop("`filePath` does not exist: ", filePath, call. = FALSE)
  }

  if (!is.null(device)) {
    device <- match.arg(device, c("pdf", "png", "svg", "jpg", "jpeg"))
  }
  check_logical_scalar(plot_unclassified, "plot_unclassified")

  if (is.null(runName)) {
    runName <- tools::file_path_sans_ext(basename(filePath))
  }

  min_read_length <- seqqc_validate_length_filter(
    min_read_length,
    "min_read_length"
  )
  max_read_length <- seqqc_validate_length_filter(
    max_read_length,
    "max_read_length"
  )
  if (!is.null(min_read_length) && !is.null(max_read_length) &&
      min_read_length > max_read_length) {
    stop("`min_read_length` must be less than or equal to `max_read_length`.",
         call. = FALSE)
  }

  if (!is.null(filename_suffix)) {
    check_scalar_character(filename_suffix, "filename_suffix")
  } else {
    filename_suffix <- seqqc_filter_filename_suffix(
      min_read_length = min_read_length,
      max_read_length = max_read_length
    )
  }
  len_read_ncol <- seqqc_validate_positive_integer(len_read_ncol, "len_read_ncol")
  len_read_unit_width <- seqqc_validate_positive_number(
    len_read_unit_width,
    "len_read_unit_width"
  )
  len_read_unit_height <- seqqc_validate_positive_number(
    len_read_unit_height,
    "len_read_unit_height"
  )
  len_read_width <- seqqc_validate_optional_positive_number(
    len_read_width,
    "len_read_width"
  )
  len_read_height <- seqqc_validate_optional_positive_number(
    len_read_height,
    "len_read_height"
  )
  sample_len_width <- seqqc_validate_positive_number(
    sample_len_width,
    "sample_len_width"
  )
  sample_len_height <- seqqc_validate_positive_number(
    sample_len_height,
    "sample_len_height"
  )
  if (!is.logical(save_sample_len_plots) ||
      length(save_sample_len_plots) != 1L ||
      is.na(save_sample_len_plots)) {
    stop("`save_sample_len_plots` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  check_scalar_character(sample_len_subdir, "sample_len_subdir")
  sample_len_dir_mode <- match.arg(sample_len_dir_mode)
  check_scalar_character(fastq_pass_trim_dir, "fastq_pass_trim_dir")

  seq_summary <- utils::read.delim(
    filePath,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_cols <- c(sample_col, length_col)
  missing_cols <- setdiff(required_cols, names(seq_summary))
  if (length(missing_cols) > 0L) {
    stop(
      "`filePath` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  seq_summary[[sample_col]] <- as.character(seq_summary[[sample_col]])
  seq_summary[[length_col]] <- suppressWarnings(as.numeric(seq_summary[[length_col]]))
  seq_summary <- seq_summary[!is.na(seq_summary[[sample_col]]), , drop = FALSE]
  if (!isTRUE(plot_unclassified)) {
    seq_summary <- seq_summary[seq_summary[[sample_col]] != unclassified, , drop = FALSE]
  }
  n_reads_raw <- nrow(seq_summary)

  keep_reads <- !is.na(seq_summary[[length_col]])
  if (!is.null(min_read_length)) {
    keep_reads <- keep_reads & seq_summary[[length_col]] >= min_read_length
  }
  if (!is.null(max_read_length)) {
    keep_reads <- keep_reads & seq_summary[[length_col]] <= max_read_length
  }
  seq_summary <- seq_summary[keep_reads, , drop = FALSE]
  n_reads_filtered <- nrow(seq_summary)

  sample_levels <- seqqc_sample_levels(
    labels = seq_summary[[sample_col]],
    unclassified = unclassified,
    plot_unclassified = plot_unclassified,
    barcodes = barcodes,
    barcode_digits = barcode_digits
  )
  seq_summary$sample <- factor(seq_summary[[sample_col]], levels = sample_levels)
  n_reads_unmatched_sample <- sum(is.na(seq_summary$sample))
  seq_summary <- seq_summary[!is.na(seq_summary$sample), , drop = FALSE]

  read_counts <- data.frame(
    sample = names(table(seq_summary$sample)),
    Freq = as.integer(table(seq_summary$sample)),
    stringsAsFactors = FALSE
  )

  num_reads <- seqqc_read_count_data(read_counts, unclassified = unclassified)

  g1 <- ggplot2::ggplot(num_reads, ggplot2::aes(.data$sample, .data$Freq)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$type), scales = "free_y", ncol = 1) +
    seqqc_theme(x_angle = 50, x_hjust = 1, x_vjust = 1) +
    ggplot2::ggtitle("Number of reads") +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL)

  seq_len_df <- data.frame(
    seq_len = seq_summary[[length_col]],
    sample = seq_summary$sample,
    stringsAsFactors = FALSE
  )
  seq_len_df <- seq_len_df[is.finite(seq_len_df$seq_len) & seq_len_df$seq_len > 0, ]
  len_read_size <- seqqc_len_read_plot_size(
    n_facets = length(sample_levels),
    ncol = len_read_ncol,
    width = len_read_width,
    height = len_read_height,
    unit_width = len_read_unit_width,
    unit_height = len_read_unit_height
  )

  g2 <- ggplot2::ggplot(seq_len_df, ggplot2::aes(.data$seq_len)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = .data$sample),
      show.legend = FALSE,
      bins = 30
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$sample),
      scales = "free",
      ncol = len_read_ncol,
      drop = FALSE
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::annotation_logticks(sides = "b") +
    seqqc_theme() +
    ggplot2::xlab("Sequence length") +
    ggplot2::ylab("Reads")

  files <- character()
  sample_len_files <- character()
  if (!is.null(device)) {
    plot_dir <- file.path(out_dir, device)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    files <- c(
      numRead = file.path(
        plot_dir,
        paste0(
          "Distribution_num_of_reads__",
          runName,
          filename_suffix,
          ".",
          device
        )
      ),
      lenRead = file.path(
        plot_dir,
        paste0(
          "Distribution_seqLength__",
          runName,
          filename_suffix,
          ".",
          device
        )
      )
    )

    ggplot2::ggsave(files[["numRead"]], plot = g1, width = 20, height = 6)
    ggplot2::ggsave(
      files[["lenRead"]],
      plot = g2,
      width = len_read_size$width,
      height = len_read_size$height, limitsize = FALSE
    )

    if (isTRUE(save_sample_len_plots)) {
      sample_len_files <- seqqc_save_sample_len_plots(
        seq_len_df = seq_len_df,
        sample_levels = sample_levels,
        runName = runName,
        filename_suffix = filename_suffix,
        device = device,
        plot_dir = plot_dir,
        sample_len_subdir = sample_len_subdir,
        sample_len_dir_mode = sample_len_dir_mode,
        fastq_pass_trim_dir = fastq_pass_trim_dir,
        width = sample_len_width,
        height = sample_len_height
      )
    }
  }

  invisible(list(
    numRead = g1,
    lenRead = g2,
    recovery = if (isTRUE(plot_unclassified)) {
      recoveryRate(read_counts, unclassified = unclassified)
    } else {
      NA_real_
    },
    read_counts = read_counts,
    n_reads_raw = n_reads_raw,
    n_reads_filtered = n_reads_filtered,
    n_reads_unmatched_sample = n_reads_unmatched_sample,
    filters = list(
      min_read_length = min_read_length,
      max_read_length = max_read_length
    ),
    filename_suffix = filename_suffix,
    len_read_size = len_read_size,
    sample_len_dir_mode = sample_len_dir_mode,
    files = files,
    sample_len_files = sample_len_files
  ))
}

seqqc_validate_length_filter <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop("`", name, "` must be a single non-negative numeric value or `NULL`.",
         call. = FALSE)
  }

  x
}

seqqc_validate_optional_positive_number <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }

  seqqc_validate_positive_number(x, name)
}

seqqc_validate_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a single positive numeric value.",
         call. = FALSE)
  }

  x
}

seqqc_validate_positive_integer <- function(x, name) {
  seqqc_validate_positive_number(x, name)
  if (x != as.integer(x)) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }

  as.integer(x)
}

seqqc_len_read_plot_size <- function(n_facets,
                                     ncol,
                                     width,
                                     height,
                                     unit_width,
                                     unit_height) {
  n_facets <- max(1L, n_facets)
  ncol_used <- min(ncol, n_facets)
  nrow_used <- ceiling(n_facets / ncol)

  list(
    width = if (is.null(width)) ncol_used * unit_width else width,
    height = if (is.null(height)) nrow_used * unit_height else height,
    ncol = ncol_used,
    nrow = nrow_used,
    unit_width = unit_width,
    unit_height = unit_height
  )
}

seqqc_save_sample_len_plots <- function(seq_len_df,
                                        sample_levels,
                                        runName,
                                        filename_suffix,
                                        device,
                                        plot_dir,
                                        sample_len_subdir,
                                        sample_len_dir_mode,
                                        fastq_pass_trim_dir,
                                        width,
                                        height) {
  sample_len_dirs <- seqqc_sample_len_dirs(
    sample_levels = sample_levels,
    plot_dir = plot_dir,
    sample_len_subdir = sample_len_subdir,
    sample_len_dir_mode = sample_len_dir_mode,
    fastq_pass_trim_dir = fastq_pass_trim_dir
  )
  sample_len_files <- stats::setNames(
    file.path(
      sample_len_dirs,
      paste0(
        "Distribution_seqLength__",
        runName,
        "__",
        seqqc_safe_filename(sample_levels),
        filename_suffix,
        ".",
        device
      )
    ),
    sample_levels
  )

  for (sample in sample_levels) {
    dir.create(dirname(sample_len_files[[sample]]), recursive = TRUE, showWarnings = FALSE)
    sample_df <- seq_len_df[seq_len_df$sample == sample, , drop = FALSE]
    sample_plot <- ggplot2::ggplot(sample_df, ggplot2::aes(.data$seq_len)) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::scale_x_log10() +
      ggplot2::annotation_logticks(sides = "b") +
      seqqc_theme() +
      ggplot2::ggtitle(sample) +
      ggplot2::xlab("Sequence length") +
      ggplot2::ylab("Reads")

    ggplot2::ggsave(
      sample_len_files[[sample]],
      plot = sample_plot,
      width = width,
      height = height, limitsize = FALSE
    )
  }

  sample_len_files
}

seqqc_sample_len_dirs <- function(sample_levels,
                                  plot_dir,
                                  sample_len_subdir,
                                  sample_len_dir_mode,
                                  fastq_pass_trim_dir) {
  if (identical(sample_len_dir_mode, "fastq_pass_trim")) {
    return(file.path(fastq_pass_trim_dir, sample_levels))
  }

  rep(file.path(plot_dir, sample_len_subdir), length(sample_levels))
}

seqqc_safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

seqqc_filter_filename_suffix <- function(min_read_length, max_read_length) {
  if (is.null(min_read_length) && is.null(max_read_length)) {
    return("")
  }

  min_label <- if (is.null(min_read_length)) {
    "min"
  } else {
    format(min_read_length, scientific = FALSE, trim = TRUE)
  }
  max_label <- if (is.null(max_read_length)) {
    "max"
  } else {
    format(max_read_length, scientific = FALSE, trim = TRUE)
  }

  paste0("__len", min_label, "-", max_label)
}

seqqc_sample_levels <- function(labels,
                                unclassified,
                                plot_unclassified,
                                barcodes,
                                barcode_digits) {
  if (is.null(barcodes)) {
    return(unique(labels))
  }

  if (!is.numeric(barcodes) || anyNA(barcodes) || any(barcodes < 0)) {
    stop("`barcodes` must be a non-negative numeric vector or `NULL`.",
         call. = FALSE)
  }

  width <- max(barcode_digits, nchar(as.character(max(barcodes))))
  levels <- paste0("barcode", stringr::str_pad(barcodes, width = width, pad = "0"))
  if (isTRUE(plot_unclassified)) {
    levels <- c(levels, unclassified)
  }

  levels
}

seqqc_read_count_data <- function(read_counts, unclassified) {
  num_reads <- read_counts
  num_reads$type <- "Freq"

  log_reads <- read_counts
  log_reads$Freq <- ifelse(log_reads$Freq > 0, log10(log_reads$Freq), 0)
  log_reads$type <- "Log10Freq"

  recovered_reads <- read_counts
  recovered_reads$Freq[recovered_reads$sample == unclassified] <- 0L
  recovered_reads$type <- "Recovered"

  rbind(num_reads, log_reads, recovered_reads)
}

seqqc_theme <- function(x_angle = 0, x_hjust = 0.5, x_vjust = 0.5) {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(
        angle = x_angle,
        hjust = x_hjust,
        vjust = x_vjust
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70")
    )
}

# Calculate the recovery rate of sample demultiplexing.
recoveryRate <- function(read_counts, unclassified = "unclassified") {
  total_reads <- sum(read_counts$Freq)
  if (total_reads == 0L) {
    return(NA_real_)
  }

  unclassified_reads <- sum(read_counts$Freq[read_counts$sample == unclassified])
  (1 - unclassified_reads / total_reads) * 100
}
