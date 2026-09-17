#' Summarize support for a deletion event from BAM CIGAR strings
#'
#' `summarize_deletion_event()` counts read support for a user-specified
#' deletion interval directly from BAM alignments. It is intended for gene
#' editing results where a deletion event is visually clear in IGV/BAM but may
#' be awkwardly represented in a VCF.
#'
#' @param bam Coordinate-sorted BAM file. A BAM index is recommended.
#' @param reference_name Reference/contig name containing the expected deletion.
#' @param deletion_start,deletion_end One-based inclusive reference coordinates
#'   of the expected deleted interval.
#' @param breakpoint_tolerance Maximum allowed distance in bp between observed
#'   and expected deletion boundaries for a read to count as `DEL_SUPPORT`.
#' @param min_mapq Minimum read mapping quality.
#' @param min_flank_coverage Number of aligned reference bases required on both
#'   sides of the deletion interval for a read to be informative. Near reference
#'   boundaries, the available flank length is used.
#' @param max_wt_deletion_bases Maximum number of deleted/skipped bases allowed
#'   inside the expected deletion interval for a read to count as `WT_SUPPORT`.
#' @param include_read_table Logical. If `TRUE`, include per-read classifications
#'   and observed deletion coordinates in the return value.
#'
#' @return A list with `summary`, `read_table`, and `parameters`. `summary`
#'   contains read counts and percentages for `DEL_SUPPORT`, `WT_SUPPORT`, and
#'   `AMBIGUOUS`.
#'
#' @examples
#' # A real indexed BAM is needed for execution.
#' # summarize_deletion_event(
#' #   bam = "sample.sort.bam",
#' #   reference_name = "PLA3_B_",
#' #   deletion_start = 100,
#' #   deletion_end = 160
#' # )
#'
#' @export
summarize_deletion_event <- function(bam,
                                     reference_name,
                                     deletion_start,
                                     deletion_end,
                                     breakpoint_tolerance = 5,
                                     min_mapq = 20,
                                     min_flank_coverage = 20,
                                     max_wt_deletion_bases = 0,
                                     include_read_table = FALSE) {
  check_file_arg(bam, "bam")
  check_scalar_character(reference_name, "reference_name")
  deletion_start <- validate_positive_integer(deletion_start, "deletion_start")
  deletion_end <- validate_positive_integer(deletion_end, "deletion_end")
  breakpoint_tolerance <- validate_nonnegative_integer(
    breakpoint_tolerance,
    "breakpoint_tolerance"
  )
  min_mapq <- validate_nonnegative_integer(min_mapq, "min_mapq")
  min_flank_coverage <- validate_nonnegative_integer(
    min_flank_coverage,
    "min_flank_coverage"
  )
  max_wt_deletion_bases <- validate_nonnegative_integer(
    max_wt_deletion_bases,
    "max_wt_deletion_bases"
  )
  check_logical_scalar(include_read_table, "include_read_table")
  if (deletion_start > deletion_end) {
    stop("`deletion_start` must be <= `deletion_end`.", call. = FALSE)
  }

  bam <- normalizePath(bam, mustWork = TRUE)
  if (!has_standard_bam_index(bam)) {
    warning(
      "No standard BAM index found next to `bam`; scanBam may be slow or fail ",
      "for some inputs. Create one with samtools index if needed.",
      call. = FALSE
    )
  }

  records <- read_deletion_event_bam_records(bam, min_mapq = min_mapq)
  if (length(records$qname) == 0L) {
    read_table <- empty_deletion_event_read_table()
  } else {
    rows <- lapply(seq_along(records$qname), function(i) {
      read <- deletion_event_bam_record(records, i)
      if (!identical(read$rname, reference_name)) return(NULL)
      classify_deletion_event_read(
        read = read,
        deletion_start = deletion_start,
        deletion_end = deletion_end,
        breakpoint_tolerance = breakpoint_tolerance,
        min_flank_coverage = min_flank_coverage,
        max_wt_deletion_bases = max_wt_deletion_bases
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    read_table <- if (length(rows) == 0L) {
      empty_deletion_event_read_table()
    } else {
      do.call(rbind, rows)
    }
  }

  summary <- summarize_deletion_event_table(read_table)
  parameters <- data.frame(
    bam = bam,
    reference_name = reference_name,
    deletion_start = deletion_start,
    deletion_end = deletion_end,
    breakpoint_tolerance = breakpoint_tolerance,
    min_mapq = min_mapq,
    min_flank_coverage = min_flank_coverage,
    max_wt_deletion_bases = max_wt_deletion_bases,
    total_primary_records_after_mapq = length(records$qname),
    reference_records = nrow(read_table),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    read_table = if (isTRUE(include_read_table)) read_table else NULL,
    parameters = parameters
  )
}

read_deletion_event_bam_records <- function(bam, min_mapq) {
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    ),
    what = c("qname", "rname", "pos", "cigar", "mapq", "strand")
  )
  x <- Rsamtools::scanBam(bam, param = param)[[1L]]
  keep <- !is.na(x$mapq) & x$mapq >= min_mapq
  lapply(x, function(col) col[keep])
}

deletion_event_bam_record <- function(x, i) {
  list(
    qname = x$qname[[i]],
    rname = as.character(x$rname[[i]]),
    pos = as.integer(x$pos[[i]]),
    cigar = x$cigar[[i]],
    mapq = x$mapq[[i]],
    strand = as.character(x$strand[[i]])
  )
}

classify_deletion_event_read <- function(read,
                                         deletion_start,
                                         deletion_end,
                                         breakpoint_tolerance,
                                         min_flank_coverage,
                                         max_wt_deletion_bases) {
  cigar_info <- parse_deletion_event_cigar(read$pos, read$cigar)
  left_interval <- c(
    max(1L, deletion_start - min_flank_coverage),
    deletion_start - 1L
  )
  right_interval <- c(
    deletion_end + 1L,
    deletion_end + min_flank_coverage
  )
  target_interval <- c(deletion_start, deletion_end)

  left_required <- deletion_start - left_interval[[1L]]
  right_required <- right_interval[[2L]] - deletion_end
  left_covered <- interval_covered_bases(cigar_info$aligned, left_interval)
  right_covered <- interval_covered_bases(cigar_info$aligned, right_interval)
  target_covered <- interval_covered_bases(cigar_info$aligned, target_interval)
  target_width <- deletion_end - deletion_start + 1L
  target_deleted <- interval_covered_bases(cigar_info$deletions, target_interval)
  best_deletion <- best_matching_deletion(
    cigar_info$deletions,
    deletion_start,
    deletion_end
  )

  left_ok <- left_covered >= left_required
  right_ok <- right_covered >= right_required
  flank_ok <- left_ok && right_ok
  matched <- !is.null(best_deletion) &&
    abs(best_deletion$start - deletion_start) <= breakpoint_tolerance &&
    abs(best_deletion$end - deletion_end) <= breakpoint_tolerance

  if (!flank_ok) {
    class <- "AMBIGUOUS"
    reason <- paste(
      c(if (!left_ok) "missing_left_flank", if (!right_ok) "missing_right_flank"),
      collapse = ";"
    )
  } else if (matched) {
    class <- "DEL_SUPPORT"
    reason <- "matched_deletion"
  } else if (!is.null(best_deletion) && best_deletion$overlap > 0L) {
    class <- "AMBIGUOUS"
    reason <- "deletion_boundary_mismatch"
  } else if (target_covered >= target_width - max_wt_deletion_bases &&
             target_deleted <= max_wt_deletion_bases) {
    class <- "WT_SUPPORT"
    reason <- "covers_region_no_target_deletion"
  } else {
    class <- "AMBIGUOUS"
    reason <- "target_region_not_fully_resolved"
  }

  data.frame(
    read_id = read$qname,
    reference = read$rname,
    read_start = read$pos,
    read_end = cigar_info$read_end,
    mapq = read$mapq,
    event_class = class,
    reason = reason,
    observed_deletion_start = if (is.null(best_deletion)) NA_integer_ else best_deletion$start,
    observed_deletion_end = if (is.null(best_deletion)) NA_integer_ else best_deletion$end,
    observed_deletion_length = if (is.null(best_deletion)) NA_integer_ else best_deletion$width,
    deletion_overlap_bases = if (is.null(best_deletion)) 0L else best_deletion$overlap,
    left_flank_covered = left_covered,
    right_flank_covered = right_covered,
    target_covered = target_covered,
    target_deleted = target_deleted,
    stringsAsFactors = FALSE
  )
}

parse_deletion_event_cigar <- function(pos, cigar) {
  lens <- as.integer(unlist(regmatches(cigar, gregexpr("[0-9]+", cigar))))
  ops <- unlist(regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar)))
  ref_pos <- pos
  aligned <- empty_interval_table()
  deletions <- empty_interval_table()

  for (i in seq_along(ops)) {
    op <- ops[[i]]
    n <- lens[[i]]
    if (op %in% c("M", "=", "X")) {
      aligned <- rbind(aligned, interval_row(ref_pos, ref_pos + n - 1L))
      ref_pos <- ref_pos + n
    } else if (op %in% c("D", "N")) {
      deletions <- rbind(deletions, interval_row(ref_pos, ref_pos + n - 1L))
      ref_pos <- ref_pos + n
    } else {
      next
    }
  }

  list(
    aligned = aligned,
    deletions = deletions,
    read_end = ref_pos - 1L
  )
}

empty_interval_table <- function() {
  data.frame(start = integer(), end = integer(), width = integer())
}

interval_row <- function(start, end) {
  data.frame(start = start, end = end, width = end - start + 1L)
}

interval_covered_bases <- function(intervals, query_interval) {
  if (nrow(intervals) == 0L) return(0L)
  start <- query_interval[[1L]]
  end <- query_interval[[2L]]
  if (is.na(start) || is.na(end) || start > end) return(0L)
  overlaps <- pmax(0L, pmin(intervals$end, end) - pmax(intervals$start, start) + 1L)
  sum(overlaps)
}

best_matching_deletion <- function(deletions, deletion_start, deletion_end) {
  if (nrow(deletions) == 0L) return(NULL)
  overlaps <- pmax(
    0L,
    pmin(deletions$end, deletion_end) - pmax(deletions$start, deletion_start) + 1L
  )
  boundary_distance <- abs(deletions$start - deletion_start) +
    abs(deletions$end - deletion_end)
  idx <- order(-overlaps, boundary_distance, -deletions$width)[[1L]]
  list(
    start = deletions$start[[idx]],
    end = deletions$end[[idx]],
    width = deletions$width[[idx]],
    overlap = overlaps[[idx]]
  )
}

summarize_deletion_event_table <- function(read_table) {
  levels <- c("DEL_SUPPORT", "WT_SUPPORT", "AMBIGUOUS")
  if (nrow(read_table) == 0L) {
    return(data.frame(
      event_class = levels,
      read_count = 0L,
      frequency = 0,
      frequency_percent = 0,
      stringsAsFactors = FALSE
    ))
  }
  counts <- table(factor(read_table$event_class, levels = levels))
  out <- data.frame(
    event_class = names(counts),
    read_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
  total <- sum(out$read_count)
  out$frequency <- if (total > 0L) out$read_count / total else 0
  out$frequency_percent <- out$frequency * 100
  out
}

empty_deletion_event_read_table <- function() {
  data.frame(
    read_id = character(),
    reference = character(),
    read_start = integer(),
    read_end = integer(),
    mapq = integer(),
    event_class = character(),
    reason = character(),
    observed_deletion_start = integer(),
    observed_deletion_end = integer(),
    observed_deletion_length = integer(),
    deletion_overlap_bases = integer(),
    left_flank_covered = integer(),
    right_flank_covered = integer(),
    target_covered = integer(),
    target_deleted = integer(),
    stringsAsFactors = FALSE
  )
}
