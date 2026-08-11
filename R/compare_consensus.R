#' Compare consensus sequences between two ONT runs
#'
#' `compare_consensus()` compares matched consensus sequences from two FASTA
#' files or two [Biostrings::DNAStringSet()] objects. Each matched pair is
#' aligned with [seq_diffs()], allowing the second sequence to be used in either
#' its original or reverse-complement orientation.
#'
#' @param consensus1,consensus2 Consensus sequences. Each input can be a path to
#'   a FASTA file or a `DNAStringSet` object.
#' @param sampleMap Optional data frame defining which sequences to compare. The
#'   first two columns are used as sequence names in `consensus1` and
#'   `consensus2`, respectively. Columns named `seq1_name` and `seq2_name` are
#'   preferred when present. When `NULL`, sequences with names present in both
#'   inputs are compared.
#' @param type Alignment type passed to [seq_diffs()]. Defaults to `"global"`.
#' @param ... Additional arguments passed to [seq_diffs()].
#'
#' @return A data frame with one row per matched consensus pair and columns:
#' \describe{
#'   \item{seq1_name, seq2_name}{Sequence names compared from `consensus1` and
#'     `consensus2`.}
#'   \item{distance}{Number of aligned differences returned by [seq_diffs()].}
#'   \item{seq1_len, seq2_len}{Ungapped sequence lengths.}
#'   \item{seq2_orientation}{Selected orientation for the second sequence.}
#'   \item{seq1, seq2}{Original sequence strings.}
#'   \item{aligned_seq1, aligned_seq2}{Aligned sequence strings including gaps.}
#' }
#'
#' @examples
#' consensus_a <- Biostrings::DNAStringSet(c(
#'   barcode01 = "ATGCCGTAAA",
#'   barcode02 = "AACCGGTT"
#' ))
#' consensus_b <- Biostrings::DNAStringSet(c(
#'   barcode01 = "ATGCCGGTAAAC",
#'   barcode02 = "AACCGGTT"
#' ))
#'
#' compare_consensus(consensus_a, consensus_b)
#'
#' sample_map <- data.frame(
#'   seq1_name = c("barcode01", "barcode02"),
#'   seq2_name = c("barcode02", "barcode01")
#' )
#' compare_consensus(consensus_a, consensus_b, sampleMap = sample_map)
#'
#' @export
compare_consensus <- function(consensus1,
                              consensus2,
                              sampleMap = NULL,
                              type = "global",
                              ...) {
  consensus1 <- read_consensus_input(consensus1, "consensus1")
  consensus2 <- read_consensus_input(consensus2, "consensus2")

  sample_map <- make_consensus_sample_map(
    consensus1 = consensus1,
    consensus2 = consensus2,
    sampleMap = sampleMap
  )

  if (nrow(sample_map) == 0L) {
    return(empty_compare_consensus_result())
  }

  results <- vector("list", nrow(sample_map))
  for (i in seq_len(nrow(sample_map))) {
    seq1_name <- sample_map$seq1_name[[i]]
    seq2_name <- sample_map$seq2_name[[i]]
    seq1 <- as.character(consensus1[[seq1_name]])
    seq2 <- as.character(consensus2[[seq2_name]])

    diff_res <- seq_diffs(
      c(seq1 = seq1, seq2 = seq2),
      type = type,
      ...
    )

    results[[i]] <- data.frame(
      seq1_name = seq1_name,
      seq2_name = seq2_name,
      distance = diff_res$n_diff,
      seq1_len = nchar(seq1),
      seq2_len = nchar(seq2),
      seq2_orientation = diff_res$orientation,
      seq1 = seq1,
      seq2 = seq2,
      aligned_seq1 = diff_res$aligned_seq1,
      aligned_seq2 = diff_res$aligned_seq2,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

read_consensus_input <- function(x, name) {
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

  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("`", name, "` sequences must all have non-empty names.", call. = FALSE)
  }

  x
}

make_consensus_sample_map <- function(consensus1, consensus2, sampleMap) {
  if (is.null(sampleMap)) {
    common_names <- intersect(names(consensus1), names(consensus2))
    return(data.frame(
      seq1_name = common_names,
      seq2_name = common_names,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(sampleMap) || ncol(sampleMap) < 2L) {
    stop("`sampleMap` must be a data frame with at least two columns.",
         call. = FALSE)
  }

  if (all(c("seq1_name", "seq2_name") %in% names(sampleMap))) {
    sample_map <- sampleMap[c("seq1_name", "seq2_name")]
  } else {
    sample_map <- sampleMap[seq_len(2)]
    names(sample_map) <- c("seq1_name", "seq2_name")
  }

  sample_map$seq1_name <- as.character(sample_map$seq1_name)
  sample_map$seq2_name <- as.character(sample_map$seq2_name)
  sample_map <- sample_map[
    !is.na(sample_map$seq1_name) &
      !is.na(sample_map$seq2_name) &
      nzchar(sample_map$seq1_name) &
      nzchar(sample_map$seq2_name),
    ,
    drop = FALSE
  ]

  sample_map <- sample_map[
    sample_map$seq1_name %in% names(consensus1) &
      sample_map$seq2_name %in% names(consensus2),
    ,
    drop = FALSE
  ]

  row.names(sample_map) <- NULL
  sample_map
}

empty_compare_consensus_result <- function() {
  data.frame(
    seq1_name = character(),
    seq2_name = character(),
    distance = integer(),
    seq1_len = integer(),
    seq2_len = integer(),
    seq2_orientation = character(),
    seq1 = character(),
    seq2 = character(),
    aligned_seq1 = character(),
    aligned_seq2 = character(),
    stringsAsFactors = FALSE
  )
}
