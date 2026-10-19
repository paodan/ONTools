parse_cigar_operations <- function(cigar) {
  if (!is.character(cigar) || length(cigar) != 1L || is.na(cigar) ||
      !nzchar(cigar) || identical(cigar, "*")) {
    stop("`cigar` must be one non-empty aligned CIGAR string.", call. = FALSE)
  }

  operations <- regmatches(
    cigar,
    gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl = TRUE)
  )[[1L]]
  if (length(operations) == 0L ||
      !identical(paste0(operations, collapse = ""), cigar)) {
    stop("Could not parse the complete CIGAR string: ", cigar,
         call. = FALSE)
  }

  operation_length <- suppressWarnings(
    as.integer(sub("[MIDNSHP=X]$", "", operations))
  )
  operation_code <- sub("^[0-9]+", "", operations)
  if (anyNA(operation_length) || any(operation_length <= 0L)) {
    stop("CIGAR operation lengths must be positive integers: ", cigar,
         call. = FALSE)
  }

  list(
    table = data.frame(
      length = operation_length,
      operation = operation_code,
      stringsAsFactors = FALSE
    ),
    reference_length = sum(
      operation_length[operation_code %in% c("M", "D", "N", "=", "X")]
    ),
    query_length = sum(
      operation_length[operation_code %in% c("M", "I", "S", "=", "X")]
    )
  )
}

validate_editing_region <- function(reference, region_start, region_end,
                                    min_anchor = 0L) {
  reference <- toupper(as.character(reference))
  if (length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
    stop("`reference` must contain exactly one non-empty sequence.",
         call. = FALSE)
  }
  region_start <- validate_positive_integer(region_start, "region_start")
  region_end <- validate_positive_integer(region_end, "region_end")
  min_anchor <- validate_nonnegative_integer(min_anchor, "min_anchor")
  if (region_end < region_start) {
    stop("`region_end` must be greater than or equal to `region_start`.",
         call. = FALSE)
  }
  if (region_end > nchar(reference, type = "bytes")) {
    stop("`region_end` exceeds the reference sequence length.",
         call. = FALSE)
  }
  list(
    reference = reference,
    region_start = region_start,
    region_end = region_end,
    min_anchor = min_anchor
  )
}

#' Extract an aligned read and reference region
#'
#' `extract_aligned_region()` applies a read's CIGAR string and returns aligned
#' reference and read strings for a closed reference interval. Insertions are
#' represented by gaps in the reference, while deletions and skipped reference
#' positions are represented by gaps in the read.
#'
#' @param reference A single reference sequence, supplied as a character string
#'   or an object coercible to one.
#' @param read_seq The read sequence in BAM alignment orientation.
#' @param cigar A single SAM/BAM CIGAR string.
#' @param ref_start Positive 1-based reference coordinate at which the alignment
#'   starts.
#' @param x,y Positive 1-based coordinates defining the closed reference
#'   interval `[x, y]`.
#' @param include_insertions Logical. Include insertions occurring immediately
#'   before the first interval base or between interval bases. Insertions after
#'   `y` are excluded.
#'
#' @return A list with `reference_aligned` and `read_aligned`. Both are `NA` if
#'   the alignment does not span the complete requested interval.
#'
#' @examples
#' extract_aligned_region(
#'   reference = "ACGTACGT",
#'   read_seq = "ACGTTACGT",
#'   cigar = "4M1I4M",
#'   ref_start = 1,
#'   x = 2,
#'   y = 7
#' )
#'
#' @export
extract_aligned_region <- function(reference,
                                   read_seq,
                                   cigar,
                                   ref_start,
                                   x,
                                   y,
                                   include_insertions = TRUE) {
  check_logical_scalar(include_insertions, "include_insertions")
  validated <- validate_editing_region(reference, x, y, min_anchor = 0L)
  parsed <- parse_cigar_operations(cigar)
  extract_aligned_region_parsed(
    reference = validated$reference,
    read_seq = read_seq,
    parsed_cigar = parsed,
    ref_start = ref_start,
    x = validated$region_start,
    y = validated$region_end,
    include_insertions = include_insertions
  )
}

extract_aligned_region_parsed <- function(reference,
                                          read_seq,
                                          parsed_cigar,
                                          ref_start,
                                          x,
                                          y,
                                          include_insertions) {
  read_seq <- toupper(as.character(read_seq))
  if (length(read_seq) != 1L || is.na(read_seq)) {
    stop("`read_seq` must contain exactly one sequence.", call. = FALSE)
  }
  ref_start <- validate_positive_integer(ref_start, "ref_start")
  if (nchar(read_seq, type = "bytes") != parsed_cigar$query_length) {
    stop(
      "Read length does not equal the query-consuming CIGAR length.",
      call. = FALSE
    )
  }

  ref_end <- ref_start + parsed_cigar$reference_length - 1L
  if (ref_end > nchar(reference, type = "bytes")) {
    stop("The CIGAR alignment extends beyond the reference sequence.",
         call. = FALSE)
  }
  if (x < ref_start || y > ref_end) {
    return(list(
      reference_aligned = NA_character_,
      read_aligned = NA_character_
    ))
  }

  cigar_table <- parsed_cigar$table
  reference_chunks <- character()
  read_chunks <- character()
  reference_position <- ref_start
  query_position <- 1L

  for (i in seq_len(nrow(cigar_table))) {
    n <- cigar_table$length[[i]]
    operation <- cigar_table$operation[[i]]

    if (operation %in% c("M", "=", "X")) {
      overlap_start <- max(reference_position, x)
      overlap_end <- min(reference_position + n - 1L, y)
      if (overlap_start <= overlap_end) {
        query_start <- query_position + overlap_start - reference_position
        query_end <- query_start + overlap_end - overlap_start
        reference_chunks <- c(
          reference_chunks,
          substr(reference, overlap_start, overlap_end)
        )
        read_chunks <- c(
          read_chunks,
          substr(read_seq, query_start, query_end)
        )
      }
      reference_position <- reference_position + n
      query_position <- query_position + n
    } else if (operation == "I") {
      position_after <- reference_position - 1L
      if (include_insertions && position_after >= x - 1L &&
          position_after < y) {
        reference_chunks <- c(reference_chunks, strrep("-", n))
        read_chunks <- c(
          read_chunks,
          substr(read_seq, query_position, query_position + n - 1L)
        )
      }
      query_position <- query_position + n
    } else if (operation %in% c("D", "N")) {
      overlap_start <- max(reference_position, x)
      overlap_end <- min(reference_position + n - 1L, y)
      if (overlap_start <= overlap_end) {
        reference_chunks <- c(
          reference_chunks,
          substr(reference, overlap_start, overlap_end)
        )
        read_chunks <- c(
          read_chunks,
          strrep("-", overlap_end - overlap_start + 1L)
        )
      }
      reference_position <- reference_position + n
    } else if (operation == "S") {
      query_position <- query_position + n
    }
  }

  list(
    reference_aligned = paste0(reference_chunks, collapse = ""),
    read_aligned = paste0(read_chunks, collapse = "")
  )
}

new_edit_event <- function(read_id,
                           type,
                           label,
                           region_start,
                           region_end,
                           reference_aligned = NA_character_,
                           read_aligned = NA_character_,
                           start = NA_integer_,
                           end = NA_integer_,
                           position_after = NA_integer_,
                           ref = NA_character_,
                           alt = NA_character_) {
  data.frame(
    read_id = read_id,
    type = type,
    start = as.integer(start),
    end = as.integer(end),
    position_after = as.integer(position_after),
    ref = ref,
    alt = alt,
    label = label,
    trim_region_start = region_start,
    trim_region_end = region_end,
    reference_aligned_trim = reference_aligned,
    read_aligned_trim = read_aligned,
    stringsAsFactors = FALSE
  )
}

#' Call editing events for one aligned read
#'
#' `call_read_edit_events()` interprets a read sequence and CIGAR against one
#' reference region. It reports wild type, SNP, adjacent multi-nucleotide
#' substitutions, insertions, deletions, complex indels, or non-callable reads.
#'
#' @param read_id Character scalar identifying the read.
#' @inheritParams extract_aligned_region
#' @param region_start,region_end Closed 1-based reference interval in which
#'   events are called.
#' @param min_anchor Non-negative number of aligned reference bases required on
#'   both sides of the calling interval.
#'
#' @return A list containing `status`, `read_id`, `combination`, and `events`.
#' `events` has one row per event. A read with no event receives one `WT` row.
#'
#' @examples
#' call <- call_read_edit_events(
#'   read_id = "read1",
#'   read_seq = "ACATACGT",
#'   cigar = "8M",
#'   ref_start = 1,
#'   reference = "ACGTACGT",
#'   region_start = 1,
#'   region_end = 8
#' )
#' call$events
#'
#' @export
call_read_edit_events <- function(read_id,
                                  read_seq,
                                  cigar,
                                  ref_start,
                                  reference,
                                  region_start = 1L,
                                  region_end = nchar(as.character(reference)),
                                  min_anchor = 0L) {
  check_scalar_character(read_id, "read_id")
  validated <- validate_editing_region(
    reference,
    region_start,
    region_end,
    min_anchor
  )
  reference <- validated$reference
  region_start <- validated$region_start
  region_end <- validated$region_end
  min_anchor <- validated$min_anchor
  read_seq <- toupper(as.character(read_seq))
  parsed <- parse_cigar_operations(cigar)
  ref_start <- validate_positive_integer(ref_start, "ref_start")

  if (length(read_seq) != 1L || is.na(read_seq) ||
      nchar(read_seq, type = "bytes") != parsed$query_length) {
    stop("`read_seq` length must equal the query-consuming CIGAR length.",
         call. = FALSE)
  }
  aligned_end <- ref_start + parsed$reference_length - 1L
  if (aligned_end > nchar(reference, type = "bytes")) {
    stop("The CIGAR alignment extends beyond the reference sequence.",
         call. = FALSE)
  }

  required_start <- max(1L, region_start - min_anchor)
  required_end <- region_end + min_anchor
  if (ref_start > required_start || aligned_end < required_end) {
    event <- new_edit_event(
      read_id = read_id,
      type = "NO_CALL",
      label = "NO_CALL",
      region_start = region_start,
      region_end = region_end
    )
    return(list(
      status = "NO_CALL",
      read_id = read_id,
      combination = "NO_CALL",
      events = event
    ))
  }

  aligned <- extract_aligned_region_parsed(
    reference = reference,
    read_seq = read_seq,
    parsed_cigar = parsed,
    ref_start = ref_start,
    x = region_start,
    y = region_end,
    include_insertions = TRUE
  )
  mismatch_position <- integer()
  mismatch_ref <- character()
  mismatch_alt <- character()
  indel_events <- list()
  reference_position <- ref_start
  query_position <- 1L

  for (i in seq_len(nrow(parsed$table))) {
    n <- parsed$table$length[[i]]
    operation <- parsed$table$operation[[i]]

    if (operation %in% c("M", "=", "X")) {
      overlap_start <- max(reference_position, region_start)
      overlap_end <- min(reference_position + n - 1L, region_end)
      if (overlap_start <= overlap_end) {
        query_start <- query_position + overlap_start - reference_position
        query_end <- query_start + overlap_end - overlap_start
        reference_base <- strsplit(
          substr(reference, overlap_start, overlap_end),
          "",
          fixed = TRUE
        )[[1L]]
        query_base <- strsplit(
          substr(read_seq, query_start, query_end),
          "",
          fixed = TRUE
        )[[1L]]
        mismatch <- which(reference_base != query_base)
        if (length(mismatch) > 0L) {
          mismatch_position <- c(
            mismatch_position,
            overlap_start + mismatch - 1L
          )
          mismatch_ref <- c(mismatch_ref, reference_base[mismatch])
          mismatch_alt <- c(mismatch_alt, query_base[mismatch])
        }
      }
      reference_position <- reference_position + n
      query_position <- query_position + n
    } else if (operation == "I") {
      position_after <- reference_position - 1L
      if (position_after >= region_start - 1L && position_after < region_end) {
        inserted_sequence <- substr(
          read_seq,
          query_position,
          query_position + n - 1L
        )
        indel_events[[length(indel_events) + 1L]] <- new_edit_event(
          read_id = read_id,
          type = "IN",
          position_after = position_after,
          ref = "-",
          alt = inserted_sequence,
          label = paste0("IN_after_", position_after, ":", inserted_sequence),
          region_start = region_start,
          region_end = region_end,
          reference_aligned = aligned$reference_aligned,
          read_aligned = aligned$read_aligned
        )
      }
      query_position <- query_position + n
    } else if (operation == "D") {
      deletion_start <- reference_position
      deletion_end <- reference_position + n - 1L
      if (deletion_end >= region_start && deletion_start <= region_end) {
        deleted_sequence <- substr(reference, deletion_start, deletion_end)
        indel_events[[length(indel_events) + 1L]] <- new_edit_event(
          read_id = read_id,
          type = "DEL",
          start = deletion_start,
          end = deletion_end,
          ref = deleted_sequence,
          alt = "-",
          label = paste0(
            "DEL_", deletion_start, "_", deletion_end, ":", deleted_sequence
          ),
          region_start = region_start,
          region_end = region_end,
          reference_aligned = aligned$reference_aligned,
          read_aligned = aligned$read_aligned
        )
      }
      reference_position <- reference_position + n
    } else if (operation == "N") {
      skipped_start <- reference_position
      skipped_end <- reference_position + n - 1L
      if (skipped_end >= region_start && skipped_start <= region_end) {
        event <- new_edit_event(
          read_id = read_id,
          type = "NO_CALL",
          label = "NO_CALL_N_skip",
          region_start = region_start,
          region_end = region_end,
          reference_aligned = aligned$reference_aligned,
          read_aligned = aligned$read_aligned
        )
        return(list(
          status = "NO_CALL",
          read_id = read_id,
          combination = "NO_CALL",
          events = event
        ))
      }
      reference_position <- reference_position + n
    } else if (operation == "S") {
      query_position <- query_position + n
    }
  }

  substitution_events <- list()
  if (length(mismatch_position) > 0L) {
    groups <- cumsum(c(TRUE, diff(mismatch_position) > 1L))
    for (group in unique(groups)) {
      index <- which(groups == group)
      event_start <- mismatch_position[index[[1L]]]
      event_end <- mismatch_position[index[[length(index)]]]
      ref_allele <- paste0(mismatch_ref[index], collapse = "")
      alt_allele <- paste0(mismatch_alt[index], collapse = "")
      event_type <- if (length(index) == 1L) "SNP" else "MNP"
      substitution_events[[length(substitution_events) + 1L]] <-
        new_edit_event(
          read_id = read_id,
          type = event_type,
          start = event_start,
          end = event_end,
          ref = ref_allele,
          alt = alt_allele,
          label = paste0(
            event_type, "_", event_start, "_", event_end,
            ":", ref_allele, ">", alt_allele
          ),
          region_start = region_start,
          region_end = region_end,
          reference_aligned = aligned$reference_aligned,
          read_aligned = aligned$read_aligned
        )
    }
  }

  all_events <- c(substitution_events, indel_events)
  if (length(all_events) == 0L) {
    event <- new_edit_event(
      read_id = read_id,
      type = "WT",
      label = "WT",
      region_start = region_start,
      region_end = region_end,
      reference_aligned = aligned$reference_aligned,
      read_aligned = aligned$read_aligned
    )
    return(list(
      status = "PASS",
      read_id = read_id,
      combination = "WT",
      events = event
    ))
  }

  event_table <- do.call(rbind, all_events)
  insertion_index <- which(event_table$type == "IN")
  deletion_index <- which(event_table$type == "DEL")
  complex_index <- integer()
  if (length(insertion_index) > 0L && length(deletion_index) > 0L) {
    for (insertion in insertion_index) {
      position <- event_table$position_after[[insertion]]
      overlaps <- deletion_index[
        position >= event_table$start[deletion_index] - 1L &
          position <= event_table$end[deletion_index]
      ]
      if (length(overlaps) > 0L) {
        complex_index <- c(complex_index, insertion, overlaps)
      }
    }
  }
  complex_index <- sort(unique(complex_index))
  if (length(complex_index) > 0L) {
    event_table$type[complex_index] <- "OTHER"
    event_table$label[complex_index] <- paste0(
      "OTHER{", event_table$label[complex_index], "}"
    )
  }

  sort_position <- ifelse(
    is.na(event_table$position_after),
    event_table$start,
    event_table$position_after + 0.5
  )
  event_table <- event_table[order(sort_position), , drop = FALSE]
  rownames(event_table) <- NULL
  list(
    status = "PASS",
    read_id = read_id,
    combination = paste(event_table$label, collapse = ";"),
    events = event_table
  )
}

#' Keep alignments spanning a reference interval
#'
#' @param aln Data frame containing at least `pos` and `cigar` columns.
#' @param region_start,region_end Closed positive 1-based reference interval.
#' @param min_anchor Non-negative number of reference bases required on both
#'   sides of the interval.
#'
#' @return `aln` restricted to alignments spanning the requested interval and
#' anchors.
#'
#' @examples
#' alignments <- data.frame(
#'   read = c("complete", "partial"),
#'   pos = c(1L, 4L),
#'   cigar = c("12M", "8M")
#' )
#' filter_reads_by_region(
#'   alignments,
#'   region_start = 4,
#'   region_end = 9,
#'   min_anchor = 2
#' )
#'
#' @export
filter_reads_by_region <- function(aln,
                                   region_start,
                                   region_end,
                                   min_anchor = 2L) {
  if (!is.data.frame(aln) || !all(c("pos", "cigar") %in% names(aln))) {
    stop("`aln` must be a data frame containing `pos` and `cigar`.",
         call. = FALSE)
  }
  region_start <- validate_positive_integer(region_start, "region_start")
  region_end <- validate_positive_integer(region_end, "region_end")
  min_anchor <- validate_nonnegative_integer(min_anchor, "min_anchor")
  if (region_end < region_start) {
    stop("`region_end` must be greater than or equal to `region_start`.",
         call. = FALSE)
  }
  if (nrow(aln) == 0L) return(aln)

  reference_length <- vapply(
    as.character(aln$cigar),
    function(cigar) parse_cigar_operations(cigar)$reference_length,
    numeric(1)
  )
  read_start <- suppressWarnings(as.integer(aln$pos))
  if (anyNA(read_start) || any(read_start < 1L)) {
    stop("`aln$pos` must contain positive integer coordinates.",
         call. = FALSE)
  }
  read_end <- read_start + reference_length - 1L
  required_start <- max(1L, region_start - min_anchor)
  required_end <- region_end + min_anchor
  keep <- read_start <= required_start & read_end >= required_end
  aln[keep, , drop = FALSE]
}

empty_variant_events <- function() {
  data.frame(
    read_id = character(),
    type = character(),
    start = integer(),
    end = integer(),
    position_after = integer(),
    ref = character(),
    alt = character(),
    label = character(),
    trim_region_start = integer(),
    trim_region_end = integer(),
    reference_aligned_trim = character(),
    read_aligned_trim = character(),
    sample_id = character(),
    Ref_name = character(),
    region_start = integer(),
    region_end = integer(),
    stringsAsFactors = FALSE
  )
}

normalize_editing_cores <- function(mc.cores) {
  if (is.null(mc.cores)) {
    detected <- parallel::detectCores(logical = TRUE)
    if (is.na(detected) || detected < 2L) return(1L)
    return(as.integer(detected - 1L))
  }
  cores <- validate_positive_integer(mc.cores, "mc.cores")
  if (.Platform$OS.type == "windows" && cores > 1L) {
    warning("`mclapply()` uses one core on Windows; setting `mc.cores = 1`.",
            call. = FALSE)
    cores <- 1L
  }
  cores
}

#' Generate per-read gene-editing events
#'
#' `make_variant_events()` filters alignments for complete interval coverage and
#' applies [call_read_edit_events()] to every retained read.
#'
#' @param aln_valid Data frame containing `qname`, `seq`, `cigar`, and `pos`, as
#'   produced from `Rsamtools::scanBam()`.
#' @param reference One named reference sequence.
#' @param region_start,region_end Closed positive 1-based calling interval.
#' @param min_anchor Non-negative number of aligned reference bases required on
#'   both sides of the interval.
#' @param sample_id Character scalar identifying the sample.
#' @param mc.cores Optional positive integer worker count. Defaults to one fewer
#'   than the detected cores. Windows automatically uses one core.
#'
#' @return A data frame with one row per called event, including sample,
#' reference, interval, and aligned-sequence columns. A typed zero-row data
#' frame is returned when no alignment spans the interval.
#'
#' @examples
#' aln <- data.frame(
#'   qname = c("wt", "snp"),
#'   seq = c("ACGTACGT", "ACATACGT"),
#'   cigar = c("8M", "8M"),
#'   pos = c(1L, 1L),
#'   stringsAsFactors = FALSE
#' )
#' events <- make_variant_events(
#'   aln,
#'   reference = c(amplicon = "ACGTACGT"),
#'   region_start = 1,
#'   region_end = 8,
#'   min_anchor = 0,
#'   mc.cores = 1
#' )
#' events[, c("read_id", "type", "label")]
#'
#' @export
make_variant_events <- function(aln_valid,
                                reference,
                                region_start,
                                region_end,
                                min_anchor = 2L,
                                sample_id = "sample",
                                mc.cores = NULL) {
  required_columns <- c("qname", "seq", "cigar", "pos")
  if (!is.data.frame(aln_valid) ||
      !all(required_columns %in% names(aln_valid))) {
    stop(
      "`aln_valid` must contain columns: ",
      paste(required_columns, collapse = ", "),
      call. = FALSE
    )
  }
  check_scalar_character(sample_id, "sample_id")
  reference_name <- names(reference)
  if (length(reference) != 1L || length(reference_name) != 1L ||
      is.na(reference_name) || !nzchar(reference_name)) {
    stop("`reference` must be exactly one named sequence.", call. = FALSE)
  }
  validated <- validate_editing_region(
    reference,
    region_start,
    region_end,
    min_anchor
  )
  cores <- normalize_editing_cores(mc.cores)
  aln_valid <- filter_reads_by_region(
    aln = aln_valid,
    region_start = validated$region_start,
    region_end = validated$region_end,
    min_anchor = validated$min_anchor
  )
  if (nrow(aln_valid) == 0L) return(empty_variant_events())

  worker <- function(i) {
    call_read_edit_events(
      read_id = as.character(aln_valid$qname[[i]]),
      read_seq = as.character(aln_valid$seq[[i]]),
      cigar = as.character(aln_valid$cigar[[i]]),
      ref_start = as.integer(aln_valid$pos[[i]]),
      reference = validated$reference,
      region_start = validated$region_start,
      region_end = validated$region_end,
      min_anchor = validated$min_anchor
    )$events
  }
  events <- if (.Platform$OS.type == "windows" || cores == 1L) {
    lapply(seq_len(nrow(aln_valid)), worker)
  } else {
    parallel::mclapply(seq_len(nrow(aln_valid)), worker, mc.cores = cores)
  }
  events_all <- do.call(rbind, events)
  events_all$sample_id <- sample_id
  events_all$Ref_name <- reference_name
  events_all$region_start <- validated$region_start
  events_all$region_end <- validated$region_end
  rownames(events_all) <- NULL
  events_all
}

empty_variant_summary <- function() {
  data.frame(
    sample_id = character(),
    Ref_name = character(),
    Index = integer(),
    Read_count = integer(),
    Ratio = numeric(),
    Variation_count = integer(),
    Variation_combination = character(),
    region_start = integer(),
    region_end = integer(),
    reference_aligned_trim = character(),
    read_aligned_trim = character(),
    stringsAsFactors = FALSE
  )
}

empty_event_frequency <- function() {
  data.frame(
    sample_id = character(),
    Ref_name = character(),
    label = character(),
    Supporting_reads = integer(),
    Callable_reads = integer(),
    Ratio = numeric(),
    Pass_filter = logical(),
    stringsAsFactors = FALSE
  )
}

validate_percentage_threshold <- function(x, name, warn_high = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x > 100) {
    stop("`", name, "` must be one finite number between 0 and 100.",
         call. = FALSE)
  }
  if (warn_high && x > 20) {
    warning("`", name, "` is ", x,
            "%; confirm that this high threshold is intended.",
            call. = FALSE)
  }
  as.numeric(x)
}

#' Summarize per-read gene-editing events
#'
#' `summarise_variant_events()` first calculates read support for each event
#' label. A non-WT label is retained only when it reaches both the absolute
#' support threshold `min_event_reads` and the relative threshold
#' `min_event_ratio`. Final variation combinations are filtered independently
#' with `min_comb_reads` and `min_comb_ratio`.
#' Reads with no retained non-WT label are classified as WT. Sample/reference
#' groups with fewer than `min_callable_reads` callable reads are retained in
#' `event_frequency` for inspection but omitted from the combination tables.
#'
#' @param events_all Event data frame returned by [make_variant_events()].
#' @param min_event_reads Positive integer. Minimum number of distinct callable
#'   reads that must support a non-WT event label.
#' @param min_event_ratio Number between 0 and 100. Minimum percentage of
#'   callable reads that must support a non-WT event label.
#' @param min_comb_reads Positive integer. Minimum number of reads that must
#'   support a final `Variation_combination`, including WT.
#' @param min_comb_ratio Number between 0 and 100. Minimum percentage of
#'   callable reads that must support a final `Variation_combination`.
#' @param min_callable_reads Positive integer. Minimum number of callable reads
#'   required for a sample/reference group to produce `var_info`. Groups below
#'   this depth remain visible in `event_frequency` with `Pass_filter = FALSE`.
#'
#' @return A list containing the original `events_all`, per-label
#' `event_frequency`, the complete combination table `var_info`, and its
#' thresholded form `var_info_filtered`. `event_frequency$Pass_filter` indicates
#' whether each label passed the callable-depth, support-count, and percentage
#' requirements. Reads with `type = "NO_CALL"` are excluded from all
#' denominators.
#'
#' @examples
#' aln <- data.frame(
#'   qname = c("wt", "snp"),
#'   seq = c("ACGTACGT", "ACATACGT"),
#'   cigar = c("8M", "8M"),
#'   pos = c(1L, 1L),
#'   stringsAsFactors = FALSE
#' )
#' events <- make_variant_events(
#'   aln,
#'   reference = c(amplicon = "ACGTACGT"),
#'   region_start = 1,
#'   region_end = 8,
#'   min_anchor = 0,
#'   mc.cores = 1
#' )
#' summarise_variant_events(
#'   events,
#'   min_event_reads = 1,
#'   min_event_ratio = 1,
#'   min_comb_reads = 1,
#'   min_comb_ratio = 1,
#'   min_callable_reads = 1
#' )$var_info
#'
#' @export
summarise_variant_events <- function(events_all,
                                     min_event_reads = 3L,
                                     min_event_ratio = 1,
                                     min_comb_reads = 3L,
                                     min_comb_ratio = 1,
                                     min_callable_reads = 10L) {
  min_event_reads <- validate_positive_integer(
    min_event_reads,
    "min_event_reads"
  )
  min_event_ratio <- validate_percentage_threshold(
    min_event_ratio,
    "min_event_ratio",
    warn_high = TRUE
  )
  min_comb_reads <- validate_positive_integer(
    min_comb_reads,
    "min_comb_reads"
  )
  min_comb_ratio <- validate_percentage_threshold(
    min_comb_ratio,
    "min_comb_ratio"
  )
  min_callable_reads <- validate_positive_integer(
    min_callable_reads,
    "min_callable_reads"
  )
  if (!is.data.frame(events_all)) {
    stop("`events_all` must be a data frame.", call. = FALSE)
  }
  empty_summary <- empty_variant_summary()
  empty_frequency <- empty_event_frequency()
  if (nrow(events_all) == 0L) {
    return(list(
      events_all = events_all,
      event_frequency = empty_frequency,
      var_info = empty_summary,
      var_info_filtered = empty_summary
    ))
  }

  required <- c("read_id", "type", "label", "Ref_name")
  missing_columns <- setdiff(required, names(events_all))
  if (length(missing_columns) > 0L) {
    stop(
      "`events_all` is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"sample_id" %in% names(events_all)) {
    events_all$sample_id <- "sample"
  }
  if (!"region_start" %in% names(events_all) &&
      "trim_region_start" %in% names(events_all)) {
    events_all$region_start <- events_all$trim_region_start
  }
  if (!"region_end" %in% names(events_all) &&
      "trim_region_end" %in% names(events_all)) {
    events_all$region_end <- events_all$trim_region_end
  }
  if (!all(c("region_start", "region_end") %in% names(events_all))) {
    stop("`events_all` must contain calling-region columns.", call. = FALSE)
  }

  callable_events <- events_all[
    !is.na(events_all$type) & events_all$type != "NO_CALL",
    ,
    drop = FALSE
  ]
  if (nrow(callable_events) == 0L) {
    return(list(
      events_all = events_all,
      event_frequency = empty_frequency,
      var_info = empty_summary,
      var_info_filtered = empty_summary
    ))
  }

  support_data <- unique(callable_events[
    , c("sample_id", "Ref_name", "read_id", "label"),
    drop = FALSE
  ])
  event_frequency <- stats::aggregate(
    read_id ~ sample_id + Ref_name + label,
    data = support_data,
    FUN = length
  )
  names(event_frequency)[names(event_frequency) == "read_id"] <-
    "Supporting_reads"
  callable_reads <- unique(callable_events[
    , c("sample_id", "Ref_name", "read_id"),
    drop = FALSE
  ])
  denominators <- stats::aggregate(
    read_id ~ sample_id + Ref_name,
    data = callable_reads,
    FUN = length
  )
  names(denominators)[names(denominators) == "read_id"] <- "Callable_reads"
  event_frequency <- merge(
    event_frequency,
    denominators,
    by = c("sample_id", "Ref_name"),
    all.x = TRUE,
    sort = FALSE
  )
  event_frequency$Ratio <- 100 * event_frequency$Supporting_reads /
    event_frequency$Callable_reads

  group_key <- function(sample_id, reference_name) {
    paste(sample_id, reference_name, sep = "\034")
  }
  eligible_groups <- denominators[
    denominators$Callable_reads >= min_callable_reads,
    ,
    drop = FALSE
  ]
  eligible_group_key <- group_key(
    eligible_groups$sample_id,
    eligible_groups$Ref_name
  )
  frequency_group_key <- group_key(
    event_frequency$sample_id,
    event_frequency$Ref_name
  )
  event_frequency$Pass_filter <-
    frequency_group_key %in% eligible_group_key &
    (
      event_frequency$label == "WT" |
        (
          event_frequency$Supporting_reads >= min_event_reads &
            event_frequency$Ratio >= min_event_ratio
        )
    )

  if (length(eligible_group_key) == 0L) {
    return(list(
      events_all = events_all,
      event_frequency = event_frequency,
      var_info = empty_summary,
      var_info_filtered = empty_summary
    ))
  }

  frequent_labels <- event_frequency[
    event_frequency$Pass_filter,
    c("sample_id", "Ref_name", "label"),
    drop = FALSE
  ]
  frequent_key <- paste(
    frequent_labels$sample_id,
    frequent_labels$Ref_name,
    frequent_labels$label,
    sep = "\034"
  )
  callable_group_key <- group_key(
    callable_events$sample_id,
    callable_events$Ref_name
  )
  analysis_events <- callable_events[
    callable_group_key %in% eligible_group_key,
    ,
    drop = FALSE
  ]
  event_key <- paste(
    analysis_events$sample_id,
    analysis_events$Ref_name,
    analysis_events$label,
    sep = "\034"
  )
  keep_event <- event_key %in% frequent_key
  read_key <- paste(
    analysis_events$sample_id,
    analysis_events$Ref_name,
    analysis_events$read_id,
    sep = "\034"
  )
  read_groups <- split(seq_len(nrow(analysis_events)), read_key)

  first_column <- function(data, column) {
    if (column %in% names(data)) as.character(data[[column]][[1L]]) else
      NA_character_
  }
  combinations <- lapply(read_groups, function(index) {
    events <- analysis_events[index, , drop = FALSE]
    retained <- unique(events$label[keep_event[index]])
    retained <- retained[retained != "WT"]
    if (length(retained) == 0L) retained <- "WT"
    data.frame(
      sample_id = as.character(events$sample_id[[1L]]),
      Ref_name = as.character(events$Ref_name[[1L]]),
      read_id = as.character(events$read_id[[1L]]),
      Variation_combination = paste(retained, collapse = ";"),
      Variation_count = sum(retained != "WT"),
      region_start = as.integer(events$region_start[[1L]]),
      region_end = as.integer(events$region_end[[1L]]),
      reference_aligned_trim = first_column(
        events,
        "reference_aligned_trim"
      ),
      read_aligned_trim = first_column(events, "read_aligned_trim"),
      stringsAsFactors = FALSE
    )
  })
  combinations <- do.call(rbind, combinations)
  rownames(combinations) <- NULL

  combination_key <- paste(
    combinations$sample_id,
    combinations$Ref_name,
    combinations$Variation_combination,
    sep = "\034"
  )
  combination_groups <- split(seq_len(nrow(combinations)), combination_key)
  mode_or_na <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(NA_character_)
    names(which.max(table(x)))
  }
  var_info <- lapply(combination_groups, function(index) {
    data <- combinations[index, , drop = FALSE]
    denominator <- denominators$Callable_reads[
      denominators$sample_id == data$sample_id[[1L]] &
        denominators$Ref_name == data$Ref_name[[1L]]
    ][[1L]]
    data.frame(
      sample_id = data$sample_id[[1L]],
      Ref_name = data$Ref_name[[1L]],
      Read_count = nrow(data),
      Ratio = round(100 * nrow(data) / denominator, 2),
      Variation_count = data$Variation_count[[1L]],
      Variation_combination = data$Variation_combination[[1L]],
      region_start = data$region_start[[1L]],
      region_end = data$region_end[[1L]],
      reference_aligned_trim = mode_or_na(data$reference_aligned_trim),
      read_aligned_trim = mode_or_na(data$read_aligned_trim),
      stringsAsFactors = FALSE
    )
  })
  var_info <- do.call(rbind, var_info)
  rownames(var_info) <- NULL
  var_info <- var_info[order(
    var_info$sample_id,
    var_info$Ref_name,
    -var_info$Read_count,
    var_info$Variation_combination
  ), , drop = FALSE]
  group_key <- paste(var_info$sample_id, var_info$Ref_name, sep = "\034")
  var_info$Index <- stats::ave(
    seq_len(nrow(var_info)),
    group_key,
    FUN = seq_along
  )
  var_info <- var_info[
    , c(
      "sample_id", "Ref_name", "Index", "Read_count", "Ratio",
      "Variation_count", "Variation_combination", "region_start",
      "region_end", "reference_aligned_trim", "read_aligned_trim"
    ),
    drop = FALSE
  ]
  var_info_filtered <- var_info[
    var_info$Read_count >= min_comb_reads &
      var_info$Ratio >= min_comb_ratio,
    ,
    drop = FALSE
  ]
  list(
    events_all = events_all,
    event_frequency = event_frequency,
    var_info = var_info,
    var_info_filtered = var_info_filtered
  )
}

#' Write a color-grouped gene-editing Excel workbook
#'
#' Rows are colored first by barcode and then by lightness for each reference
#' within a barcode. The data are otherwise written unchanged.
#'
#' @param data A data frame containing barcode and reference columns.
#' @param file Output `.xlsx` path. Parent directories are created.
#' @param barcode_col,ref_col Column names identifying barcode and reference.
#' @param sheet Worksheet name.
#' @param main_colors Non-empty character vector of valid R colors, recycled
#'   across barcodes.
#'
#' @return Invisibly returns the normalized output path.
#'
#' @examples
#' if (requireNamespace("openxlsx", quietly = TRUE)) {
#'   output <- tempfile(fileext = ".xlsx")
#'   color_excel_by_barcode_ref(
#'     data.frame(
#'       barcode = c("barcode01", "barcode01", "barcode02"),
#'       Ref_name = c("amp1", "amp2", "amp1"),
#'       Ratio = c(70, 20, 90)
#'     ),
#'     output
#'   )
#' }
#'
#' @export
color_excel_by_barcode_ref <- function(
    data,
    file,
    barcode_col = "barcode",
    ref_col = "Ref_name",
    sheet = "Results",
    main_colors = c("#4E79A7", "#59A14F", "#E15759")) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package `openxlsx` is required to write the colored workbook.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  check_scalar_character(file, "file")
  check_scalar_character(barcode_col, "barcode_col")
  check_scalar_character(ref_col, "ref_col")
  check_scalar_character(sheet, "sheet")
  missing_columns <- setdiff(c(barcode_col, ref_col), names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing coloring columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  if (!is.character(main_colors) || length(main_colors) == 0L ||
      anyNA(main_colors) || any(!nzchar(main_colors))) {
    stop("`main_colors` must be a non-empty character vector of colors.",
         call. = FALSE)
  }
  tryCatch(
    grDevices::col2rgb(main_colors),
    error = function(e) stop("`main_colors` contains an invalid color.",
                             call. = FALSE)
  )

  barcode <- as.character(data[[barcode_col]])
  reference_name <- as.character(data[[ref_col]])
  barcode[is.na(barcode) | !nzchar(barcode)] <- "unclassified"
  reference_name[is.na(reference_name) | !nzchar(reference_name)] <-
    "unclassified"
  barcodes <- unique(barcode)
  barcode_colors <- stats::setNames(
    rep(main_colors, length.out = length(barcodes)),
    barcodes
  )
  workbook <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeData(workbook, sheet, data)

  for (current_barcode in barcodes) {
    barcode_rows <- which(barcode == current_barcode)
    references <- unique(reference_name[barcode_rows])
    white_mix <- if (length(references) == 1L) 0.45 else
      seq(0.78, 0.25, length.out = length(references))
    base_rgb <- grDevices::col2rgb(
      barcode_colors[[current_barcode]]
    )[, 1L] / 255
    for (i in seq_along(references)) {
      rows <- which(
        barcode == current_barcode & reference_name == references[[i]]
      )
      fill_rgb <- base_rgb * (1 - white_mix[[i]]) + white_mix[[i]]
      fill_color <- grDevices::rgb(
        fill_rgb[[1L]], fill_rgb[[2L]], fill_rgb[[3L]],
        maxColorValue = 1
      )
      openxlsx::addStyle(
        workbook,
        sheet = sheet,
        style = openxlsx::createStyle(fgFill = fill_color),
        rows = rows + 1L,
        cols = seq_len(ncol(data)),
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }
  openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  file <- normalizePath(file, mustWork = FALSE)
  openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
  invisible(file)
}

editing_cache_signature <- function(bam_file,
                                    reference,
                                    reference_name,
                                    region_start,
                                    region_end,
                                    min_anchor,
                                    min_mapq,
                                    reference_name_map) {
  bam_info <- file.info(bam_file)
  list(
    schema_version = 2L,
    bam_path = normalizePath(bam_file, mustWork = TRUE),
    bam_size = as.numeric(bam_info$size[[1L]]),
    bam_mtime = as.numeric(bam_info$mtime[[1L]]),
    reference_name = reference_name,
    reference_sequence = toupper(as.character(reference)),
    region_start = as.integer(region_start),
    region_end = as.integer(region_end),
    min_anchor = as.integer(min_anchor),
    min_mapq = as.integer(min_mapq),
    reference_name_map = reference_name_map
  )
}

read_editing_cache <- function(path, signature) {
  if (!file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(cached) || !identical(cached$schema_version, 2L) ||
      !identical(cached$signature, signature) ||
      !is.data.frame(cached$events)) {
    return(NULL)
  }
  cached$events
}

write_editing_cache <- function(events, path, signature) {
  saveRDS(
    list(schema_version = 2L, signature = signature, events = events),
    path
  )
  invisible(path)
}

normalize_editing_sample_info <- function(sample_info) {
  if (is.character(sample_info) && length(sample_info) == 1L &&
      !is.na(sample_info)) {
    check_file_arg(sample_info, "sample_info")
    sample_info <- as.data.frame(readxl::read_excel(sample_info))
  } else if (is.data.frame(sample_info)) {
    sample_info <- as.data.frame(sample_info)
  } else {
    stop("`sample_info` must be a data frame or one Excel file path.",
         call. = FALSE)
  }
  aliases <- list(
    ref_name = "Reference_name",
    region_start = "Reference_Region_Start",
    region_end = "Reference_Region_End"
  )
  for (standard in names(aliases)) {
    if (!standard %in% names(sample_info)) {
      legacy <- aliases[[standard]]
      if (!legacy %in% names(sample_info)) {
        stop("`sample_info` is missing column `", standard,
             "` (or `", legacy, "`).", call. = FALSE)
      }
      sample_info[[standard]] <- sample_info[[legacy]]
    }
  }
  if (!"barcode" %in% names(sample_info)) {
    stop("`sample_info` is missing column `barcode`.", call. = FALSE)
  }
  if (nrow(sample_info) == 0L) {
    stop("`sample_info` must contain at least one row.", call. = FALSE)
  }
  barcode <- as.character(sample_info$barcode)
  ref_name <- as.character(sample_info$ref_name)
  if (anyNA(barcode) || any(!nzchar(trimws(barcode)))) {
    stop("`sample_info$barcode` contains missing or empty values.",
         call. = FALSE)
  }
  if (anyNA(ref_name) || any(!nzchar(trimws(ref_name)))) {
    stop("`sample_info$ref_name` contains missing or empty values.",
         call. = FALSE)
  }
  region_start <- suppressWarnings(as.integer(sample_info$region_start))
  region_end <- suppressWarnings(as.integer(sample_info$region_end))
  if (anyNA(region_start) || anyNA(region_end) || any(region_start < 1L) ||
      any(region_end < region_start)) {
    stop("Sample editing regions must contain positive integers with start <= end.",
         call. = FALSE)
  }
  sample_info$barcode <- barcode
  sample_info$ref_name <- ref_name
  sample_info$region_start <- region_start
  sample_info$region_end <- region_end
  sample_info
}

#' Run batch gene-editing event analysis from aligned BAM files
#'
#' `run_editing_batch()` reads primary alignments for each barcode, calls editing
#' events inside one or more configured reference intervals, caches per-read
#' events, summarizes edit combinations, and optionally writes a colored Excel
#' workbook.
#'
#' @param sample_info A data frame or Excel path. Required columns are
#'   `barcode`, `ref_name`, `region_start`, and `region_end`. Legacy names
#'   `Reference_name`, `Reference_Region_Start`, and `Reference_Region_End` are
#'   accepted. Multiple references in `ref_name` are separated with `|`.
#' @param ref_all Named reference sequences. Every configured reference must be
#'   present.
#' @param re_names Optional named character vector mapping BAM reference names
#'   to names in `ref_all`. Unmapped names are retained unchanged.
#' @param bam_root Directory containing
#'   `<barcode>/alignments/<barcode>.aligned.sorted.bam`.
#' @param cache_dir Directory for versioned event caches. Cache validity includes
#'   BAM path, size and modification time, complete reference sequence, region,
#'   anchor, MAPQ threshold, and reference-name mapping.
#' @param min_anchor Non-negative reference bases required on both sides of the
#'   calling interval.
#' @param min_event_reads Minimum number of distinct reads supporting a non-WT
#'   event label.
#' @param min_event_ratio Minimum percentage of callable reads supporting a
#'   non-WT event label.
#' @param min_comb_reads Minimum number of reads supporting a final
#'   `Variation_combination`, including WT.
#' @param min_comb_ratio Minimum percentage of callable reads supporting a final
#'   `Variation_combination`.
#' @param min_callable_reads Minimum callable read depth required to summarize a
#'   sample/reference group. Passed to [summarise_variant_events()].
#' @param min_mapq Integer between 0 and 255. Alignments below this mapping
#'   quality are excluded.
#' @param mc.cores Optional worker count passed to [make_variant_events()].
#' @param force_rebuild Logical. Ignore valid caches and rebuild all events.
#' @param colored_xlsx Optional `.xlsx` output path passed to
#'   [color_excel_by_barcode_ref()].
#'
#' @return A list with `table`, the combined filtered summary, and `per_sample`,
#' a nested list of summaries by input row and reference.
#'
#' @details
#' Only mapped primary, non-supplementary, non-duplicate, QC-passing alignments
#' are loaded. Reads must span the complete calling interval plus anchors. Cache
#' files from older implementations are automatically considered stale.
#'
#' @examples
#' \dontrun{
#' sample_sheet <- data.frame(
#'   barcode = c("barcode01", "barcode02"),
#'   ref_name = c("amplicon_A", "amplicon_B"),
#'   region_start = c(120L, 85L),
#'   region_end = c(150L, 115L)
#' )
#' references <- c(
#'   amplicon_A = "ACGT...",
#'   amplicon_B = "TGCA..."
#' )
#'
#' result <- run_editing_batch(
#'   sample_info = sample_sheet,
#'   ref_all = references,
#'   bam_root = "/path/to/barcode-directories",
#'   min_anchor = 10,
#'   min_event_reads = 3,
#'   min_event_ratio = 1,
#'   min_comb_reads = 3,
#'   min_comb_ratio = 1,
#'   min_callable_reads = 10,
#'   min_mapq = 20,
#'   colored_xlsx = "editing-summary.xlsx"
#' )
#' result$table
#' }
#'
#' @seealso [make_variant_events()], [summarise_variant_events()]
#' @export
run_editing_batch <- function(sample_info,
                              ref_all,
                              re_names = NULL,
                              bam_root,
                              cache_dir = file.path(bam_root, "event_cache"),
                              min_anchor = 0L,
                              min_event_reads = 3L,
                              min_event_ratio = 1,
                              min_comb_reads = 3L,
                              min_comb_ratio = 1,
                              min_callable_reads = 10L,
                              min_mapq = 20L,
                              mc.cores = NULL,
                              force_rebuild = FALSE,
                              colored_xlsx = NULL) {
  sample_info <- normalize_editing_sample_info(sample_info)
  if (length(ref_all) == 0L || is.null(names(ref_all)) ||
      anyNA(names(ref_all)) || any(!nzchar(names(ref_all))) ||
      anyDuplicated(names(ref_all))) {
    stop("`ref_all` must contain uniquely named reference sequences.",
         call. = FALSE)
  }
  if (!is.null(re_names)) {
    if (!is.character(re_names) || is.null(names(re_names)) ||
        anyNA(re_names) || anyNA(names(re_names)) ||
        any(!nzchar(names(re_names)))) {
      stop("`re_names` must be `NULL` or a named character vector.",
           call. = FALSE)
    }
  }
  check_dir_arg(bam_root, "bam_root")
  check_scalar_character(cache_dir, "cache_dir")
  min_anchor <- validate_nonnegative_integer(min_anchor, "min_anchor")
  min_event_reads <- validate_positive_integer(
    min_event_reads,
    "min_event_reads"
  )
  min_event_ratio <- validate_percentage_threshold(
    min_event_ratio,
    "min_event_ratio",
    warn_high = TRUE
  )
  min_comb_reads <- validate_positive_integer(
    min_comb_reads,
    "min_comb_reads"
  )
  min_comb_ratio <- validate_percentage_threshold(
    min_comb_ratio,
    "min_comb_ratio"
  )
  min_callable_reads <- validate_positive_integer(
    min_callable_reads,
    "min_callable_reads"
  )
  min_mapq <- validate_nonnegative_integer(min_mapq, "min_mapq")
  if (min_mapq > 255L) stop("`min_mapq` cannot exceed 255.", call. = FALSE)
  cores <- normalize_editing_cores(mc.cores)
  check_logical_scalar(force_rebuild, "force_rebuild")
  if (!is.null(colored_xlsx)) check_scalar_character(colored_xlsx, "colored_xlsx")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- normalizePath(cache_dir, mustWork = TRUE)
  bam_root <- normalizePath(bam_root, mustWork = TRUE)

  bam_param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE,
      isNotPassingQualityControls = FALSE
    ),
    what = c("qname", "rname", "pos", "cigar", "mapq", "seq", "flag")
  )
  output_rows <- list()
  per_sample_results <- list()
  output_index <- 1L
  safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

  for (i in seq_len(nrow(sample_info))) {
    barcode <- sample_info$barcode[[i]]
    reference_names <- trimws(strsplit(
      sample_info$ref_name[[i]],
      split = "\\|"
    )[[1L]])
    reference_names <- unique(reference_names[nzchar(reference_names)])
    region_start <- sample_info$region_start[[i]]
    region_end <- sample_info$region_end[[i]]
    message("Processing ", barcode, ": ", paste(reference_names, collapse = ", "))

    missing_references <- setdiff(reference_names, names(ref_all))
    if (length(missing_references) > 0L) {
      warning(
        "References not found in `ref_all` for ", barcode, ": ",
        paste(missing_references, collapse = ", "),
        call. = FALSE
      )
    }
    reference_names <- intersect(reference_names, names(ref_all))
    if (length(reference_names) == 0L) {
      per_sample_results[[paste0(barcode, "_", i)]] <- list()
      next
    }

    bam_file <- file.path(
      bam_root,
      barcode,
      "alignments",
      paste0(barcode, ".aligned.sorted.bam")
    )
    if (!file.exists(bam_file)) {
      warning("BAM file not found; skipping: ", bam_file, call. = FALSE)
      per_sample_results[[paste0(barcode, "_", i)]] <- list()
      next
    }
    cache_files <- stats::setNames(
      file.path(
        cache_dir,
        paste0(
          safe_name(barcode), "__", safe_name(reference_names), "__",
          region_start, "_", region_end, "__anchor", min_anchor,
          "__mapq", min_mapq, ".events.rds"
        )
      ),
      reference_names
    )
    signatures <- lapply(reference_names, function(reference_name) {
      editing_cache_signature(
        bam_file = bam_file,
        reference = ref_all[[reference_name]],
        reference_name = reference_name,
        region_start = region_start,
        region_end = region_end,
        min_anchor = min_anchor,
        min_mapq = min_mapq,
        reference_name_map = re_names
      )
    })
    names(signatures) <- reference_names
    cached_events <- lapply(reference_names, function(reference_name) {
      if (force_rebuild) return(NULL)
      read_editing_cache(
        cache_files[[reference_name]],
        signatures[[reference_name]]
      )
    })
    names(cached_events) <- reference_names
    needs_rebuild <- vapply(cached_events, is.null, logical(1))

    if (any(needs_rebuild)) {
      alignment <- as.data.frame(
        Rsamtools::scanBam(bam_file, param = bam_param)[[1L]]
      )
      if (nrow(alignment) > 0L) {
        alignment <- alignment[
          !is.na(alignment$mapq) & alignment$mapq >= min_mapq,
          ,
          drop = FALSE
        ]
      }
      original_name <- as.character(alignment$rname)
      mapped_name <- original_name
      if (!is.null(re_names) && length(original_name) > 0L) {
        candidate <- unname(re_names[original_name])
        mapped_name[!is.na(candidate)] <- candidate[!is.na(candidate)]
      }
      alignment$rname <- mapped_name

      for (reference_name in reference_names[needs_rebuild]) {
        alignment_reference <- alignment[
          !is.na(alignment$rname) & alignment$rname == reference_name,
          ,
          drop = FALSE
        ]
        events <- make_variant_events(
          aln_valid = alignment_reference,
          reference = stats::setNames(
            as.character(ref_all[[reference_name]]),
            reference_name
          ),
          region_start = region_start,
          region_end = region_end,
          min_anchor = min_anchor,
          sample_id = barcode,
          mc.cores = cores
        )
        write_editing_cache(
          events,
          cache_files[[reference_name]],
          signatures[[reference_name]]
        )
        cached_events[[reference_name]] <- events
      }
    }

    barcode_result <- list()
    for (reference_name in reference_names) {
      events <- cached_events[[reference_name]]
      if (is.null(events) || nrow(events) == 0L) next
      summary_result <- summarise_variant_events(
        events,
        min_event_reads = min_event_reads,
        min_event_ratio = min_event_ratio,
        min_comb_reads = min_comb_reads,
        min_comb_ratio = min_comb_ratio,
        min_callable_reads = min_callable_reads
      )
      table <- summary_result$var_info_filtered
      if (nrow(table) == 0L) next
      table$barcode <- barcode
      if ("Sample_HiTOM" %in% names(sample_info)) {
        table$Sample_HiTOM <- as.character(sample_info$Sample_HiTOM[[i]])
      }
      barcode_result[[reference_name]] <- table
      output_rows[[output_index]] <- table
      output_index <- output_index + 1L
    }
    per_sample_results[[paste0(barcode, "_", i)]] <- barcode_result
  }

  result_table <- if (length(output_rows) == 0L) {
    data.frame()
  } else {
    output <- do.call(rbind, output_rows)
    rownames(output) <- NULL
    front_columns <- intersect(c("Sample_HiTOM", "barcode"), names(output))
    output[, c(front_columns, setdiff(names(output), front_columns)), drop = FALSE]
  }
  if (!is.null(colored_xlsx)) {
    if (nrow(result_table) == 0L) {
      warning("No summarized events were available; Excel output was skipped.",
              call. = FALSE)
    } else {
      color_excel_by_barcode_ref(
        data = result_table,
        file = colored_xlsx,
        barcode_col = "barcode",
        ref_col = "Ref_name",
        sheet = "Results"
      )
    }
  }
  list(table = result_table, per_sample = per_sample_results)
}
