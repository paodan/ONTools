#' Align Sanger reads to a full-length reference
#'
#' `align_sanger_with_full_reference()` orients Sanger reads against a single
#' full-length reference sequence, then returns a multiple sequence alignment
#' containing the reference followed by the oriented Sanger reads.
#'
#' Each Sanger read is first aligned to the reference in both its original and
#' reverse-complement orientations using local pairwise alignment. The
#' orientation with the higher alignment score is retained. The reference and
#' oriented Sanger reads are then aligned together with
#' [DECIPHER::AlignSeqs()].
#'
#' @param reference Full-length reference sequence. Must be a FASTA file path or
#'   a [Biostrings::DNAStringSet()] containing exactly one sequence.
#' @param sanger Sanger read sequences. Must be a FASTA file path or a
#'   [Biostrings::DNAStringSet()] containing one or more sequences.
#' @param match_score,mismatch_score Numeric match and mismatch scores used to
#'   build the nucleotide substitution matrix passed to
#'   [pwalign::pairwiseAlignment()].
#' @param gap_opening,gap_extension Numeric gap opening and extension penalties
#'   passed to [pwalign::pairwiseAlignment()].
#' @param processors Positive integer passed to [DECIPHER::AlignSeqs()].
#'
#' @return A list with class `"SangerMultipleAlignment"` containing:
#' \describe{
#'   \item{alignment}{Multiple sequence alignment returned by
#'     [DECIPHER::AlignSeqs()]. The reference is the first sequence.}
#'   \item{orientation}{A data frame with one row per Sanger read and columns
#'     `sample`, `direction`, `forward_score`, `reverse_score`,
#'     `selected_score`, and `score_delta`.}
#'   \item{oriented_sanger}{A [Biostrings::DNAStringSet()] containing the Sanger
#'     reads in the orientation selected for alignment.}
#' }
#'
#' @examples
#' reference <- Biostrings::DNAStringSet(c(
#'   reference = "AAACCCGGGTTTAAACCCGGGTTT"
#' ))
#' sanger <- Biostrings::DNAStringSet(c(
#'   read_forward = "CCCGGGTTTAAA",
#'   read_reverse = "TTTAAACCCGGG"
#' ))
#'
#' res <- align_sanger_with_full_reference(reference, sanger)
#' res$orientation
#' names(res$alignment)
#'
#' @export
align_sanger_with_full_reference <- function(reference,
                                             sanger,
                                             match_score = 2,
                                             mismatch_score = -3,
                                             gap_opening = 5,
                                             gap_extension = 2,
                                             processors = 1) {
  reference <- read_dna_string_set_input(reference, "reference")
  sanger <- read_dna_string_set_input(sanger, "sanger")

  if (length(reference) != 1L) {
    stop("`reference` must contain exactly one sequence.", call. = FALSE)
  }

  if (length(sanger) < 1L) {
    stop("`sanger` must contain at least one sequence.", call. = FALSE)
  }

  match_score <- validate_scalar_number(match_score, "match_score")
  mismatch_score <- validate_scalar_number(mismatch_score, "mismatch_score")
  gap_opening <- validate_nonnegative_number(gap_opening, "gap_opening")
  gap_extension <- validate_nonnegative_number(gap_extension, "gap_extension")
  processors <- validate_positive_integer(processors, "processors")

  reference_names <- names(reference)
  reference_name <- if (is.null(reference_names)) "" else reference_names[[1]]
  if (is.null(reference_name) || is.na(reference_name) || !nzchar(reference_name)) {
    reference_name <- "reference"
  }
  names(reference) <- reference_name

  sanger <- ensure_dna_string_set_names(sanger, prefix = "sanger")
  reference <- clean_dna_string_set(reference)
  sanger <- clean_dna_string_set(sanger)

  substitution_matrix <- pwalign::nucleotideSubstitutionMatrix(
    match = match_score,
    mismatch = mismatch_score,
    baseOnly = FALSE
  )

  orientation_results <- lapply(
    seq_along(sanger),
    function(i) {
      forward_sequence <- sanger[[i]]
      reverse_sequence <- Biostrings::reverseComplement(forward_sequence)

      forward_alignment <- pwalign::pairwiseAlignment(
        pattern = forward_sequence,
        subject = reference[[1]],
        type = "local",
        substitutionMatrix = substitution_matrix,
        gapOpening = gap_opening,
        gapExtension = gap_extension
      )

      reverse_alignment <- pwalign::pairwiseAlignment(
        pattern = reverse_sequence,
        subject = reference[[1]],
        type = "local",
        substitutionMatrix = substitution_matrix,
        gapOpening = gap_opening,
        gapExtension = gap_extension
      )

      forward_score <- as.numeric(pwalign::score(forward_alignment))
      reverse_score <- as.numeric(pwalign::score(reverse_alignment))

      if (reverse_score > forward_score) {
        selected_sequence <- reverse_sequence
        direction <- "reverse_complement"
        selected_score <- reverse_score
      } else {
        selected_sequence <- forward_sequence
        direction <- "forward"
        selected_score <- forward_score
      }

      list(
        sequence = selected_sequence,
        direction = direction,
        forward_score = forward_score,
        reverse_score = reverse_score,
        selected_score = selected_score
      )
    }
  )

  oriented_sanger <- Biostrings::DNAStringSet(vapply(
    orientation_results,
    function(x) as.character(x$sequence),
    character(1)
  ))
  names(oriented_sanger) <- names(sanger)

  orientation_table <- data.frame(
    sample = names(sanger),
    direction = vapply(orientation_results, `[[`, character(1), "direction"),
    forward_score = vapply(orientation_results, `[[`, numeric(1), "forward_score"),
    reverse_score = vapply(orientation_results, `[[`, numeric(1), "reverse_score"),
    selected_score = vapply(orientation_results, `[[`, numeric(1), "selected_score"),
    stringsAsFactors = FALSE
  )
  orientation_table$score_delta <- abs(
    orientation_table$forward_score - orientation_table$reverse_score
  )

  all_sequences <- c(reference, oriented_sanger)
  alignment <- DECIPHER::AlignSeqs(
    all_sequences,
    processors = processors,
    verbose = FALSE
  )
  names(alignment) <- c(reference_name, names(oriented_sanger))

  structure(
    list(
      alignment = alignment,
      orientation = orientation_table,
      oriented_sanger = oriented_sanger
    ),
    class = c("SangerMultipleAlignment", "list")
  )
}

read_dna_string_set_input <- function(x, name) {
  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop("`", name, "` must be a single FASTA file path or a DNAStringSet.",
           call. = FALSE)
    }

    if (!file.exists(x)) {
      stop("`", name, "` does not exist: ", x, call. = FALSE)
    }

    x <- Biostrings::readDNAStringSet(x)
  }

  if (!inherits(x, "DNAStringSet")) {
    stop("`", name, "` must be a FASTA file path or a DNAStringSet.",
         call. = FALSE)
  }

  if (length(x) == 0L) {
    stop("`", name, "` must contain at least one sequence.", call. = FALSE)
  }

  x
}

ensure_dna_string_set_names <- function(x, prefix) {
  sequence_names <- names(x)
  if (is.null(sequence_names)) {
    sequence_names <- rep("", length(x))
  }

  missing_names <- is.na(sequence_names) | !nzchar(sequence_names)
  sequence_names[missing_names] <- paste0(prefix, "_", which(missing_names))
  names(x) <- sequence_names

  x
}

clean_dna_string_set <- function(x) {
  sequence_names <- names(x)
  x <- Biostrings::DNAStringSet(toupper(gsub("\\s+", "", as.character(x))))
  names(x) <- sequence_names
  x
}

validate_scalar_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be a single finite numeric value.", call. = FALSE)
  }

  x
}
