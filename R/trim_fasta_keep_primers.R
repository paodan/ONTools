#' Trim FASTA sequences while keeping primer sequences
#'
#' `trim_fasta_keep_primers()` trims extra bases outside a primer-bounded
#' amplicon region. The forward-orientation structure is expected to be
#' `extra + forward primer + insert + reverse-complement reverse primer + extra`.
#' The reverse-complement structure is also checked, and trimmed sequences can
#' optionally be returned in forward orientation.
#'
#' @param input_fasta Path to an input FASTA file.
#' @param output_fasta Path to the output FASTA file.
#' @param f_primer Forward primer sequence. IUPAC ambiguity codes are supported.
#' @param r_primer Reverse primer sequence. IUPAC ambiguity codes are supported.
#' @param max_mismatch Maximum mismatches allowed when matching each primer.
#' @param orient_to_forward Logical. If `TRUE`, sequences matched in
#'   reverse-complement orientation are reverse-complemented before writing.
#' @param width Number of bases per FASTA line in the output file.
#'
#' @return Invisibly returns a data frame with one row per input sequence:
#'   sequence name, whether trimming succeeded, detected orientation, trim start
#'   and end coordinates in the input sequence, input length, output length, and
#'   output name written to the FASTA file.
#'
#' @examples
#' f_primer <- "AGAGTTTGATCMTGGCTCAG"
#' r_primer <- "TACGGYTACCTTGTTACGACTT"
#' input <- tempfile(fileext = ".fasta")
#' output <- tempfile(fileext = ".fasta")
#' seq <- paste0(
#'   "AA",
#'   "AGAGTTTGATCATGGCTCAG",
#'   "ACGTACGT",
#'   "AAGTCGTAACAAGGTACCGTA",
#'   "TT"
#' )
#' Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(sample1 = seq)), input)
#'
#' trim_fasta_keep_primers(
#'   input_fasta = input,
#'   output_fasta = output,
#'   f_primer = f_primer,
#'   r_primer = r_primer
#' )
#'
#' @export
trim_fasta_keep_primers <- function(input_fasta,
                                    output_fasta,
                                    f_primer,
                                    r_primer,
                                    max_mismatch = 2,
                                    orient_to_forward = TRUE,
                                    width = 80) {
  check_file_arg(input_fasta, "input_fasta")
  check_scalar_character(output_fasta, "output_fasta")
  check_scalar_character(f_primer, "f_primer")
  check_scalar_character(r_primer, "r_primer")
  max_mismatch <- validate_nonnegative_integer(max_mismatch, "max_mismatch")
  width <- validate_positive_integer(width, "width")

  if (!is.logical(orient_to_forward) || length(orient_to_forward) != 1L ||
      is.na(orient_to_forward)) {
    stop("`orient_to_forward` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  seqs <- Biostrings::readDNAStringSet(input_fasta)
  if (length(seqs) == 0L) {
    stop("`input_fasta` must contain at least one sequence.", call. = FALSE)
  }

  if (is.null(names(seqs)) || any(!nzchar(names(seqs)))) {
    names(seqs) <- paste0("seq", seq_along(seqs))
  }

  trimmed <- lapply(
    seqs,
    trim_consensus_keep_primers,
    f_primer = f_primer,
    r_primer = r_primer,
    max_mismatch = max_mismatch,
    orient_to_forward = orient_to_forward
  )

  trimmed_sequences <- vapply(trimmed, `[[`, character(1), "sequence")
  out_seqs <- Biostrings::DNAStringSet(trimmed_sequences)

  summary <- data.frame(
    seq_name = names(seqs),
    trimmed = vapply(trimmed, `[[`, logical(1), "trimmed"),
    orientation = vapply(trimmed, `[[`, character(1), "orientation"),
    start = vapply(trimmed, `[[`, integer(1), "start"),
    end = vapply(trimmed, `[[`, integer(1), "end"),
    input_length = Biostrings::width(seqs),
    output_length = Biostrings::width(out_seqs),
    stringsAsFactors = FALSE
  )
  summary$output_name <- paste(
    summary$seq_name,
    paste0(
      "trimmed=", summary$trimmed,
      ";orientation=", summary$orientation,
      ";start=", summary$start,
      ";end=", summary$end
    )
  )

  names(out_seqs) <- summary$output_name
  Biostrings::writeXStringSet(out_seqs, output_fasta, width = width)

  invisible(summary)
}

trim_consensus_keep_primers <- function(seq,
                                        f_primer,
                                        r_primer,
                                        max_mismatch = 2,
                                        orient_to_forward = TRUE) {
  seq <- toupper(as.character(seq))
  f_primer <- toupper(f_primer)
  r_primer <- toupper(r_primer)

  rc_f_primer <- reverse_complement_string(f_primer)
  rc_r_primer <- reverse_complement_string(r_primer)

  candidates <- list()

  f_hit <- find_primer_hit(f_primer, seq, max_mismatch = max_mismatch)
  rc_r_hit <- find_primer_hit(rc_r_primer, seq, max_mismatch = max_mismatch)

  if (!is.null(f_hit) && !is.null(rc_r_hit) && f_hit$start < rc_r_hit$start) {
    candidates[["forward"]] <- list(
      orientation = "forward",
      start = f_hit$start,
      end = rc_r_hit$end
    )
  }

  r_hit <- find_primer_hit(r_primer, seq, max_mismatch = max_mismatch)
  rc_f_hit <- find_primer_hit(rc_f_primer, seq, max_mismatch = max_mismatch)

  if (!is.null(r_hit) && !is.null(rc_f_hit) && r_hit$start < rc_f_hit$start) {
    candidates[["reverse_complement"]] <- list(
      orientation = "reverse_complement",
      start = r_hit$start,
      end = rc_f_hit$end
    )
  }

  if (length(candidates) == 0L) {
    return(list(
      sequence = seq,
      trimmed = FALSE,
      orientation = NA_character_,
      start = NA_integer_,
      end = NA_integer_
    ))
  }

  candidate_widths <- vapply(
    candidates,
    function(x) x$end - x$start + 1L,
    numeric(1)
  )
  best <- candidates[[which.max(candidate_widths)]]
  trimmed_seq <- substr(seq, best$start, best$end)

  if (isTRUE(orient_to_forward) &&
      identical(best$orientation, "reverse_complement")) {
    trimmed_seq <- reverse_complement_string(trimmed_seq)
  }

  list(
    sequence = trimmed_seq,
    trimmed = TRUE,
    orientation = best$orientation,
    start = best$start,
    end = best$end
  )
}

find_primer_hit <- function(pattern, subject, max_mismatch) {
  hits <- Biostrings::matchPattern(
    pattern = Biostrings::DNAString(pattern),
    subject = Biostrings::DNAString(subject),
    max.mismatch = max_mismatch,
    fixed = FALSE
  )

  if (length(hits) == 0L) {
    return(NULL)
  }

  data.frame(
    start = BiocGenerics::start(hits),
    end = BiocGenerics::end(hits),
    width = Biostrings::width(hits)
  )[1, ]
}

reverse_complement_string <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}

validate_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x != as.integer(x)) {
    stop("`", name, "` must be a single non-negative integer.", call. = FALSE)
  }

  as.integer(x)
}

validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x != as.integer(x)) {
    stop("`", name, "` must be a single positive integer.", call. = FALSE)
  }

  as.integer(x)
}
