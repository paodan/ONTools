#' Compare two DNA sequences and report differences
#'
#' `seq_diffs()` aligns two DNA sequences in both forward and reverse-complement
#' orientations, chooses the orientation with fewer aligned differences, and
#' reports mismatch, insertion, and deletion positions.
#'
#' @param seqs A character vector containing exactly two DNA sequences. Names are
#'   optional but will be used in the output when present.
#' @param type Alignment type passed to [pwalign::pairwiseAlignment()]. Defaults
#'   to `"global"`.
#' @param ... Additional arguments passed to [pwalign::pairwiseAlignment()].
#'
#' @return A list with:
#' \describe{
#'   \item{orientation}{Best orientation for the second sequence: `"forward"` or
#'     `"reverse_complement"`.}
#'   \item{n_diff}{Number of aligned positions where the two sequences differ.}
#'   \item{differences}{A data frame containing difference positions, bases,
#'     sequence names, and orientation metadata.}
#'   \item{alignment}{The selected `PairwiseAlignmentsSingleSubject` object.}
#'   \item{aligned_seq1}{The aligned first sequence, including gap characters.}
#'   \item{aligned_seq2}{The aligned second sequence in the selected orientation,
#'     including gap characters.}
#' }
#'
#' @examples
#' seqs <- c(
#'   A = "ATGCCGTAAA",
#'   B = "TTACGGCATAA"
#' )
#'
#' res <- seq_diffs(seqs)
#' res$orientation
#' res$differences
#'
#' @export
seq_diffs <- function(seqs, type = "global", ...) {
  if (!is.character(seqs) || length(seqs) != 2) {
    stop("`seqs` must be a character vector containing exactly two sequences.",
         call. = FALSE)
  }

  if (anyNA(seqs) || any(!nzchar(seqs))) {
    stop("`seqs` must not contain missing or empty sequences.", call. = FALSE)
  }

  seq_names <- names(seqs)
  if (is.null(seq_names) || any(!nzchar(seq_names))) {
    seq_names <- c("seq1", "seq2")
  }

  seq1 <- Biostrings::DNAString(seqs[[1]])
  seq2 <- Biostrings::DNAString(seqs[[2]])
  seq2_rc <- Biostrings::reverseComplement(seq2)

  aln_fwd <- pwalign::pairwiseAlignment(seq1, seq2, type = type, ...)
  aln_rc <- pwalign::pairwiseAlignment(seq1, seq2_rc, type = type, ...)

  res_fwd <- make_seq_diff_result(
    aln = aln_fwd,
    orientation = "forward",
    seq_names = seq_names,
    seq2_length = length(seq2)
  )

  res_rc <- make_seq_diff_result(
    aln = aln_rc,
    orientation = "reverse_complement",
    seq_names = seq_names,
    seq2_length = length(seq2)
  )

  if (res_rc$n_diff < res_fwd$n_diff) {
    res_rc
  } else {
    res_fwd
  }
}

make_seq_diff_result <- function(aln, orientation, seq_names, seq2_length) {
  p <- strsplit(as.character(pwalign::alignedPattern(aln)), "")[[1]]
  s <- strsplit(as.character(pwalign::alignedSubject(aln)), "")[[1]]

  seq1_coord <- cumsum(p != "-")
  seq2_coord <- cumsum(s != "-")

  diff_idx <- which(p != s)

  diffs <- data.frame(
    aln_pos = diff_idx,
    seq1_name = rep(seq_names[[1]], length(diff_idx)),
    seq1_pos = ifelse(p[diff_idx] == "-", NA_integer_, seq1_coord[diff_idx]),
    seq1_base = p[diff_idx],
    seq2_name = rep(seq_names[[2]], length(diff_idx)),
    seq2_pos = ifelse(s[diff_idx] == "-", NA_integer_, seq2_coord[diff_idx]),
    seq2_base = s[diff_idx],
    orientation = rep(orientation, length(diff_idx)),
    seq2_aligned_as = rep(seq2_orientation_label(orientation), length(diff_idx)),
    stringsAsFactors = FALSE
  )

  if (orientation == "reverse_complement" && nrow(diffs) > 0) {
    diffs$seq2_pos <- ifelse(
      is.na(diffs$seq2_pos),
      NA_integer_,
      seq2_length - diffs$seq2_pos + 1L
    )
  }

  list(
    orientation = orientation,
    n_diff = length(diff_idx),
    differences = diffs,
    alignment = aln,
    aligned_seq1 = as.character(pwalign::alignedPattern(aln)),
    aligned_seq2 = as.character(pwalign::alignedSubject(aln))
  )
}

seq2_orientation_label <- function(orientation) {
  if (identical(orientation, "reverse_complement")) {
    "seq2_reverse_complement"
  } else {
    "seq2_original"
  }
}
