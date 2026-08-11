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
#' @param barcodes Integer barcode numbers to include and order. Use `NULL` to
#'   use sample labels as they appear in the sequencing summary file.
#' @param barcode_digits Minimum number of digits used for barcode labels.
#' @param out_dir Parent output directory for plot files.
#' @param sample_col Column containing barcode or sample labels. Defaults to
#'   `"alias"`, matching the existing ONTools workflow.
#' @param length_col Column containing read lengths. Defaults to
#'   `"sequence_length_template"`.
#'
#' @return Invisibly returns a list with the read-count plot (`numRead`), the
#'   read-length plot (`lenRead`), demultiplexing recovery percentage
#'   (`recovery`), per-sample read counts (`read_counts`), and paths written
#'   (`files`).
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
#' @importFrom rlang .data
#' @export
plot_seqQC <- function(filePath,
                       runName = NULL,
                       device = "pdf",
                       unclassified = "unclassified",
                       barcodes = 1:96,
                       barcode_digits = 2,
                       out_dir = "figs",
                       sample_col = "alias",
                       length_col = "sequence_length_template") {
  if (!is.character(filePath) || length(filePath) != 1L || is.na(filePath)) {
    stop("`filePath` must be a single file path.", call. = FALSE)
  }

  if (!file.exists(filePath)) {
    stop("`filePath` does not exist: ", filePath, call. = FALSE)
  }

  if (!is.null(device)) {
    device <- match.arg(device, c("pdf", "png", "svg", "jpg", "jpeg"))
  }

  if (is.null(runName)) {
    runName <- tools::file_path_sans_ext(basename(filePath))
  }

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

  sample_levels <- seqqc_sample_levels(
    labels = seq_summary[[sample_col]],
    unclassified = unclassified,
    barcodes = barcodes,
    barcode_digits = barcode_digits
  )
  seq_summary$sample <- factor(seq_summary[[sample_col]], levels = sample_levels)

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

  g2 <- ggplot2::ggplot(seq_len_df, ggplot2::aes(.data$seq_len)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = .data$sample),
      show.legend = FALSE,
      bins = 30
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$sample), scales = "free", ncol = 12, drop = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::annotation_logticks(sides = "b") +
    seqqc_theme() +
    ggplot2::xlab("Sequence length") +
    ggplot2::ylab("Reads")

  files <- character()
  if (!is.null(device)) {
    plot_dir <- file.path(out_dir, device)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    files <- c(
      numRead = file.path(
        plot_dir,
        paste0("Distribution_num_of_reads__", runName, ".", device)
      ),
      lenRead = file.path(
        plot_dir,
        paste0("Distribution_seqLength__", runName, ".", device)
      )
    )

    ggplot2::ggsave(files[["numRead"]], plot = g1, width = 10, height = 6)
    ggplot2::ggsave(
      files[["lenRead"]],
      plot = g2,
      width = 24,
      height = max(2, 2 * ceiling(length(sample_levels) / 12))
    )
  }

  invisible(list(
    numRead = g1,
    lenRead = g2,
    recovery = recoveryRate(read_counts, unclassified = unclassified),
    read_counts = read_counts,
    files = files
  ))
}

seqqc_sample_levels <- function(labels,
                                unclassified,
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
  c(
    paste0("barcode", stringr::str_pad(barcodes, width = width, pad = "0")),
    unclassified
  )
}

seqqc_read_count_data <- function(read_counts, unclassified) {
  num_reads <- read_counts
  num_reads$type <- "Freq"

  log_reads <- read_counts
  log_reads$Freq <- ifelse(log_reads$Freq > 0, log10(log_reads$Freq), NA_real_)
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
