#' Plot pairwise differences among FASTA sequences
#'
#' `plot_fasta_dist()` calculates the number of aligned differences for every
#' pair of DNA sequences and displays the resulting symmetric matrix as a
#' heatmap. Pairwise differences are calculated with [seq_diffs()], so the
#' forward and reverse-complement orientations are compared automatically.
#'
#' @param fasta A character vector or [Biostrings::DNAStringSet] containing at
#'   least two DNA sequences. Sequence names are used as matrix labels; missing
#'   or empty names are replaced by `seq1`, `seq2`, and so on.
#' @param type Alignment type passed to [seq_diffs()] and ultimately to
#'   [pwalign::pairwiseAlignment()]. Defaults to `"global"`.
#' @param cores Number of CPU cores used for pairwise comparisons. When `NULL`,
#'   the function uses at most one fewer than the number of detected cores and
#'   never more than the number of matrix rows being calculated. Parallel
#'   execution is unavailable on Windows, where the function runs sequentially.
#' @param cluster_rows,cluster_cols Logical values passed to
#'   [pheatmap::pheatmap()] indicating whether rows or columns should be
#'   clustered. Both default to `FALSE` so input order is retained.
#' @param display_numbers Logical; display the difference count in each heatmap
#'   cell. Defaults to `TRUE`.
#' @param fontsize_number Font size used for displayed difference counts.
#' @param ... Additional arguments passed to [pheatmap::pheatmap()], such as
#'   `color`, `main`, `border_color`, or `silent`.
#'
#' @return Invisibly returns a list with two components:
#'
#' * `D`: a symmetric integer matrix of pairwise difference counts;
#' * `plot`: the object returned by [pheatmap::pheatmap()].
#'
#' @details
#' Only the upper triangle is calculated. It is then copied to the lower
#' triangle, and the diagonal is set to zero. Consequently, each unordered
#' sequence pair is aligned exactly once.
#'
#' @examples
#' fasta <- c(
#'   reference = "AAAA",
#'   sample_1 = "AAAT",
#'   sample_2 = "AATT"
#' )
#'
#' result <- plot_fasta_dist(fasta, cores = 1, silent = TRUE)
#' result$D
#' result$plot
#'
#' @seealso [seq_diffs()], [pheatmap::pheatmap()]
#' @export
plot_fasta_dist <- function(fasta,
                            type = "global",
                            cores = NULL,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            display_numbers = TRUE,
                            fontsize_number = 13,
                            ...) {
  if (!is.character(fasta) && !inherits(fasta, "DNAStringSet")) {
    stop(
      "`fasta` must be a character vector or a DNAStringSet.",
      call. = FALSE
    )
  }

  if (length(fasta) < 2L) {
    stop("`fasta` must contain at least two sequences.", call. = FALSE)
  }

  sequences <- as.character(fasta)
  if (anyNA(sequences) || any(!nzchar(sequences))) {
    stop("`fasta` must not contain missing or empty sequences.", call. = FALSE)
  }

  sequence_names <- names(fasta)
  if (is.null(sequence_names)) {
    sequence_names <- rep("", length(sequences))
  }
  missing_names <- is.na(sequence_names) | !nzchar(sequence_names)
  sequence_names[missing_names] <- paste0("seq", which(missing_names))
  names(sequences) <- sequence_names

  n_sequences <- length(sequences)
  if (is.null(cores)) {
    detected_cores <- suppressWarnings(parallel::detectCores())
    if (is.na(detected_cores)) {
      detected_cores <- 1L
    }
    cores <- max(1L, min(detected_cores - 1L, n_sequences - 1L))
  }

  if (
    length(cores) != 1L || !is.numeric(cores) || !is.finite(cores) ||
      cores < 1 || cores != as.integer(cores)
  ) {
    stop("`cores` must be a single positive integer or NULL.", call. = FALSE)
  }
  cores <- as.integer(cores)

  calculate_row <- function(i) {
    row <- integer(n_sequences)
    for (j in seq.int(i + 1L, n_sequences)) {
      row[[j]] <- seq_diffs(sequences[c(i, j)], type = type)$n_diff
    }
    row
  }

  row_indices <- seq_len(n_sequences - 1L)
  if (cores > 1L && .Platform$OS.type != "windows") {
    rows <- parallel::mclapply(
      row_indices,
      calculate_row,
      mc.cores = min(cores, length(row_indices))
    )
  } else {
    rows <- lapply(row_indices, calculate_row)
  }
  rows[[n_sequences]] <- integer(n_sequences)

  distance_matrix <- do.call(rbind, rows)
  distance_matrix[lower.tri(distance_matrix)] <-
    t(distance_matrix)[lower.tri(distance_matrix)]
  dimnames(distance_matrix) <- list(sequence_names, sequence_names)

  heatmap <- pheatmap::pheatmap(
    distance_matrix,
    scale = "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    display_numbers = display_numbers,
    fontsize_number = fontsize_number,
    ...
  )

  invisible(list(D = distance_matrix, plot = heatmap))
}
