#' Calculate per-read quality metrics from a FASTQ file
#'
#' `calculate_fastq_read_q()` reads a FASTQ or FASTQ.GZ file with
#' `ShortRead::readFastq()` and returns one row per read. The calculations keep
#' the same quality-score convention used by standard Sanger/Illumina FASTQ
#' files: ASCII-encoded Phred scores with offset 33.
#'
#' @param fastq_file Path to an input FASTQ or FASTQ.GZ file.
#'
#' @return A data frame with one row per read. Row names are the read names
#'   parsed from the FASTQ identifiers. If duplicated read names are present,
#'   row names are made unique with [make.unique()] while the original read
#'   names are preserved in `read_id`.
#'
#' The returned columns are:
#' \describe{
#'   \item{read_id}{Read name parsed from the FASTQ identifier. Any text after
#'     the first space is removed, matching common FASTQ header conventions.}
#'   \item{length}{Read length in bases.}
#'   \item{arithmetic_mean_q}{Arithmetic mean of per-base Phred quality scores
#'     for the read. This is calculated as the sum of per-base Q scores divided
#'     by read length.}
#'   \item{mean_q}{Probability-derived mean Phred quality. Per-base Q scores
#'     are first converted to error probabilities with `10^(-Q / 10)`, the mean
#'     error probability is calculated, and that mean error probability is then
#'     converted back to a Phred score with `-10 * log10(mean_error)`.}
#'   \item{error_rate}{Mean per-base error probability implied by `mean_q`,
#'     calculated as `10^(-mean_q / 10)`.}
#'   \item{accuracy}{Mean per-base accuracy implied by `mean_q`, calculated as
#'     `1 - error_rate`.}
#' }
#'
#' @examples
#' \dontrun{
#' read_q <- calculate_fastq_read_q("reads.fastq.gz")
#' head(read_q)
#' }
#'
#' @export
calculate_fastq_read_q <- function(fastq_file) {
  check_file_arg(fastq_file, "fastq_file")
  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop(
      "Package `ShortRead` is required by `calculate_fastq_read_q()`. ",
      "Install it with BiocManager::install('ShortRead').",
      call. = FALSE
    )
  }

  fastq_file <- normalizePath(fastq_file, mustWork = TRUE)
  fq <- ShortRead::readFastq(fastq_file)
  n_reads <- length(fq)

  if (n_reads == 0L) {
    out <- data.frame(
      read_id = character(),
      length = integer(),
      arithmetic_mean_q = numeric(),
      mean_q = numeric(),
      error_rate = numeric(),
      accuracy = numeric(),
      stringsAsFactors = FALSE
    )
    rownames(out) <- character()
    return(out)
  }

  fq_quality <- ShortRead::quality(fq)
  qstrings <- rep("", length(fq_quality))
  for (mi in seq_along(qstrings)) {
    qstrings[mi] <- as.character(fq_quality[[mi]])
  }

  probability_mean_q <- vapply(
    qstrings,
    function(x) {
      base_q <- utf8ToInt(x) - 33
      mean_error <- mean(10^(-base_q / 10))
      -10 * log10(mean_error)
    },
    numeric(1)
  )

  read_id <- sub(" .*", "", as.character(ShortRead::id(fq)))
  read_length <- BiocGenerics::width(fq)
  out <- data.frame(
    read_id = read_id,
    length = read_length,
    arithmetic_mean_q = ShortRead::alphabetScore(fq_quality) / read_length,
    mean_q = probability_mean_q,
    error_rate = 10^(-probability_mean_q / 10),
    accuracy = 1 - 10^(-probability_mean_q / 10),
    stringsAsFactors = FALSE
  )
  rownames(out) <- make.unique(read_id)
  out
}
