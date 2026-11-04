#' Summarize genomic loci of BAM alignments
#'
#' `summarize_bam_alignment_loci()` identifies the reference sequences and
#' genomic intervals to which reads align. It is intended for diagnosing
#' off-target or non-specific amplicons in long-read data. Reference alignment
#' ends are calculated from reference-consuming CIGAR operations rather than
#' from query width, so soft clipping, insertions, and deletions are handled
#' correctly.
#'
#' @param bam Path to a BAM file. An index is not required because the BAM is
#'   scanned sequentially.
#' @param min_mapq Non-negative integer. Keep alignments with mapping quality
#'   greater than or equal to this value. Use `21` to reproduce a strict
#'   `MAPQ > 20` filter.
#' @param reference_name `NULL` or a character vector of reference sequence
#'   names to retain. `NULL` retains every reference represented after the
#'   alignment filters.
#' @param primary_only Logical. If `TRUE`, exclude secondary (`0x100`) and
#'   supplementary (`0x800`) alignments.
#' @param exclude_duplicates Logical. If `TRUE`, exclude records marked as
#'   duplicates (`0x400`).
#' @param exclude_qcfail Logical. If `TRUE`, exclude records marked as failing
#'   vendor quality control (`0x200`).
#' @param min_reference_reads Positive integer. A reference is included in the
#'   histogram when it contains at least this many retained alignment records.
#' @param bin_width Positive integer giving the histogram and genomic-bin width
#'   in reference bases.
#'
#' @return A list with:
#'   \itemize{
#'     \item `alignments`: one row per retained BAM alignment, including query
#'       width, CIGAR-derived reference width, and reference end.
#'     \item `mapq_distribution`: MAPQ counts after flag and optional reference
#'       filtering but before the `min_mapq` filter.
#'     \item `reference_summary`: alignment-record and unique-read counts by
#'       reference sequence.
#'     \item `peak_summary`: counts in fixed-width genomic bins.
#'     \item `plot`: a [ggplot2::ggplot()] histogram of alignment starts for
#'       references meeting `min_reference_reads`.
#'     \item `filter_summary`: counts at each filtering stage and the parameters
#'       used.
#'   }
#'
#' @details
#' Mapped records are required in all cases. When `primary_only = TRUE`, only
#' primary alignments are retained. The reference span is the sum of CIGAR
#' operations consuming reference bases (`M`, `D`, `N`, `=`, and `X`). Query
#' width is reported separately and is never used to calculate the reference
#' end.
#'
#' Counts in `alignments`, `reference_summary$alignment_records`, and
#' `peak_summary$alignment_records` refer to BAM alignment records. Unique read
#' counts are based on distinct query names and may be smaller when a query name
#' occurs more than once.
#'
#' @examples
#' \dontrun{
#' result <- summarize_bam_alignment_loci(
#'   "barcode150.wheat.bam",
#'   min_mapq = 21,
#'   min_reference_reads = 50,
#'   bin_width = 1000
#' )
#'
#' result$reference_summary
#' head(result$peak_summary)
#' result$plot
#' }
#'
#' @export
summarize_bam_alignment_loci <- function(
    bam,
    min_mapq = 20L,
    reference_name = NULL,
    primary_only = TRUE,
    exclude_duplicates = FALSE,
    exclude_qcfail = FALSE,
    min_reference_reads = 50L,
    bin_width = 1000L) {
  check_file_arg(bam, "bam")
  min_mapq <- validate_nonnegative_integer(min_mapq, "min_mapq")
  if (min_mapq > 255L) {
    stop("`min_mapq` cannot exceed 255.", call. = FALSE)
  }
  min_reference_reads <- validate_positive_integer(
    min_reference_reads,
    "min_reference_reads"
  )
  bin_width <- validate_positive_integer(bin_width, "bin_width")
  check_logical_scalar(primary_only, "primary_only")
  check_logical_scalar(exclude_duplicates, "exclude_duplicates")
  check_logical_scalar(exclude_qcfail, "exclude_qcfail")
  if (!is.null(reference_name)) {
    if (!is.character(reference_name) || length(reference_name) == 0L ||
        anyNA(reference_name) || any(!nzchar(trimws(reference_name)))) {
      stop(
        "`reference_name` must be `NULL` or contain non-empty names.",
        call. = FALSE
      )
    }
    reference_name <- unique(trimws(reference_name))
  }

  bam <- normalizePath(bam, mustWork = TRUE)
  total_records <- as.integer(Rsamtools::countBam(bam)$records[[1L]])
  flag_arguments <- list(isUnmappedQuery = FALSE)
  if (primary_only) {
    flag_arguments$isSecondaryAlignment <- FALSE
    flag_arguments$isSupplementaryAlignment <- FALSE
  }
  if (exclude_duplicates) flag_arguments$isDuplicate <- FALSE
  if (exclude_qcfail) flag_arguments$isNotPassingQualityControls <- FALSE
  bam_flag <- do.call(Rsamtools::scanBamFlag, flag_arguments)
  bam_param <- Rsamtools::ScanBamParam(
    flag = bam_flag,
    what = c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar")
  )
  records <- Rsamtools::scanBam(bam, param = bam_param)[[1L]]
  flag_pass_records <- length(records$qname)

  record_reference <- as.character(records$rname)
  reference_keep <- !is.na(record_reference)
  if (!is.null(reference_name)) {
    reference_keep <- reference_keep & record_reference %in% reference_name
  }
  mapq_input <- suppressWarnings(as.integer(records$mapq[reference_keep]))
  mapq_distribution <- bam_locus_mapq_distribution(mapq_input)
  mapq_keep <- !is.na(records$mapq) & records$mapq >= min_mapq
  cigar <- as.character(records$cigar)
  valid_cigar <- !is.na(cigar) & nzchar(cigar) & cigar != "*"
  keep <- reference_keep & mapq_keep & valid_cigar
  selected <- which(keep)

  alignments <- empty_bam_locus_alignments()
  if (length(selected) > 0L) {
    reference_width <- cigar_reference_lengths(cigar[selected])
    alignment_start <- as.integer(records$pos[selected])
    alignments <- data.frame(
      read_id = as.character(records$qname[selected]),
      ref_name = record_reference[selected],
      start = alignment_start,
      end = as.integer(alignment_start + reference_width - 1L),
      query_width = as.integer(records$qwidth[selected]),
      reference_width = as.integer(reference_width),
      mapq = as.integer(records$mapq[selected]),
      cigar = cigar[selected],
      flag = as.integer(records$flag[selected]),
      stringsAsFactors = FALSE
    )
  }

  reference_summary <- summarize_bam_locus_references(alignments)
  peak_summary <- summarize_bam_locus_bins(alignments, bin_width)
  plot_references <- reference_summary$ref_name[
    reference_summary$alignment_records >= min_reference_reads
  ]
  plot_data <- alignments[
    alignments$ref_name %in% plot_references,
    ,
    drop = FALSE
  ]
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$start)) +
    ggplot2::geom_histogram(
      binwidth = bin_width,
      boundary = 1,
      color = "white",
      fill = "#4E79A7"
    ) +
    ggplot2::labs(
      title = "BAM alignment start positions",
      subtitle = paste0("MAPQ >= ", min_mapq),
      x = "Reference start position",
      y = "Alignment records"
    ) +
    ggplot2::theme_bw()
  if (nrow(plot_data) > 0L) {
    plot <- plot + ggplot2::facet_wrap(~ref_name, scales = "free_x")
  }

  filter_summary <- list(
    bam = bam,
    total_records = total_records,
    flag_pass_records = as.integer(flag_pass_records),
    reference_pass_records = as.integer(sum(reference_keep)),
    mapq_pass_records = as.integer(sum(reference_keep & mapq_keep)),
    retained_records = as.integer(nrow(alignments)),
    retained_unique_reads = as.integer(length(unique(alignments$read_id))),
    min_mapq = min_mapq,
    reference_name = reference_name,
    primary_only = primary_only,
    exclude_duplicates = exclude_duplicates,
    exclude_qcfail = exclude_qcfail,
    min_reference_reads = min_reference_reads,
    bin_width = bin_width
  )

  list(
    alignments = alignments,
    mapq_distribution = mapq_distribution,
    reference_summary = reference_summary,
    peak_summary = peak_summary,
    plot = plot,
    filter_summary = filter_summary
  )
}

empty_bam_locus_alignments <- function() {
  data.frame(
    read_id = character(),
    ref_name = character(),
    start = integer(),
    end = integer(),
    query_width = integer(),
    reference_width = integer(),
    mapq = integer(),
    cigar = character(),
    flag = integer(),
    stringsAsFactors = FALSE
  )
}

bam_locus_mapq_distribution <- function(mapq) {
  mapq <- mapq[!is.na(mapq)]
  if (length(mapq) == 0L) {
    return(data.frame(mapq = integer(), alignment_records = integer()))
  }
  count <- table(mapq)
  data.frame(
    mapq = as.integer(names(count)),
    alignment_records = as.integer(count),
    row.names = NULL
  )
}

summarize_bam_locus_references <- function(alignments) {
  if (nrow(alignments) == 0L) {
    return(data.frame(
      ref_name = character(),
      alignment_records = integer(),
      unique_reads = integer(),
      median_mapq = numeric(),
      min_start = integer(),
      max_end = integer(),
      stringsAsFactors = FALSE
    ))
  }
  groups <- split(seq_len(nrow(alignments)), alignments$ref_name)
  rows <- lapply(names(groups), function(reference) {
    index <- groups[[reference]]
    data.frame(
      ref_name = reference,
      alignment_records = length(index),
      unique_reads = length(unique(alignments$read_id[index])),
      median_mapq = stats::median(alignments$mapq[index]),
      min_start = min(alignments$start[index]),
      max_end = max(alignments$end[index]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$unique_reads, -out$alignment_records, out$ref_name), ]
  rownames(out) <- NULL
  out
}

summarize_bam_locus_bins <- function(alignments, bin_width) {
  if (nrow(alignments) == 0L) {
    return(data.frame(
      ref_name = character(),
      bin_start = integer(),
      bin_end = integer(),
      alignment_records = integer(),
      unique_reads = integer(),
      median_mapq = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  bin_start <- ((alignments$start - 1L) %/% bin_width) * bin_width + 1L
  key <- paste(alignments$ref_name, bin_start, sep = "\034")
  groups <- split(seq_len(nrow(alignments)), key)
  rows <- lapply(groups, function(index) {
    current_start <- bin_start[index[[1L]]]
    data.frame(
      ref_name = alignments$ref_name[index[[1L]]],
      bin_start = as.integer(current_start),
      bin_end = as.integer(current_start + bin_width - 1L),
      alignment_records = length(index),
      unique_reads = length(unique(alignments$read_id[index])),
      median_mapq = stats::median(alignments$mapq[index]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$unique_reads, -out$alignment_records,
                   out$ref_name, out$bin_start), ]
  rownames(out) <- NULL
  out
}
