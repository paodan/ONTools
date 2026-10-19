#' Subset and sample reads from a BAM file
#'
#' `subset_bam_reads()` selects alignments by reference sequence, genomic
#' interval, and mapping quality, optionally removes secondary and
#' supplementary alignments, and writes at most `max_reads` records to a new
#' BAM file. The input BAM is scanned sequentially, so an index is not required.
#'
#' @param bam Input BAM file.
#' @param output_bam Output BAM file. It must differ from `bam`.
#' @param reference_name Non-empty character vector naming one or more BAM
#'   reference sequences. Duplicate names are removed while preserving order.
#' @param region_start,region_end Positive 1-based coordinates defining a
#'   closed reference interval.
#' @param min_mapq Integer from 0 to 255. Only alignments with mapping quality
#'   greater than or equal to this value are eligible.
#' @param max_reads Positive integer giving the maximum number of alignment
#'   records to write, either overall or per reference according to
#'   `max_reads_per_reference`.
#' @param max_reads_per_reference Logical. If `FALSE` (the default), all
#'   eligible records across `reference_name` are pooled and at most
#'   `max_reads` are sampled. If `TRUE`, each reference is sampled separately
#'   and contributes at most `max_reads` records.
#' @param coverage Character scalar. `"complete"` requires an alignment to span
#'   the complete interval. `"overlap"` keeps any alignment whose reference
#'   span overlaps the interval.
#' @param primary_only Logical. If `TRUE`, secondary and supplementary
#'   alignments are excluded.
#' @param seed `NULL` or one non-negative integer used for reproducible random
#'   sampling when more than `max_reads` records are eligible. The caller's
#'   random-number state is restored before the function returns.
#' @param index_output Logical. Create a BAI index for the output BAM. This
#'   requires the retained records to remain coordinate sorted.
#' @param overwrite Logical. Overwrite an existing output BAM and its standard
#'   `.bai` index.
#'
#' @details
#' Reference spans are calculated from reference-consuming CIGAR operations
#' (`M`, `D`, `N`, `=`, and `X`). Consequently, a deletion or skipped region
#' contributes to whether an alignment spans the requested interval. The BAM
#' header and complete selected records, including sequences, qualities, and
#' optional tags, are retained in the output.
#'
#' The function first identifies all eligible records and then makes a random
#' sample only when their number exceeds the applicable limit. With
#' `max_reads_per_reference = FALSE`, the output contains at most `max_reads`
#' records in total. With `max_reads_per_reference = TRUE`, each requested
#' reference contributes at most `max_reads` records, so the total can exceed
#' `max_reads`. Sampling is performed on BAM records, not merely on query
#' names, so repeated query names do not increase the output beyond the
#' applicable maximum.
#'
#' @return Invisibly returns a list with `output_bam`, `eligible_reads`,
#'   `selected_reads`, `sampled`, `per_reference`, and the filtering
#'   parameters.
#'
#' @examples
#' dontrun{
#' result <- subset_bam_reads(
#'   bam = "barcode01.aligned.sorted.bam",
#'   output_bam = "barcode01.region_1000_1500.bam",
#'   reference_name = "amplicon_A",
#'   region_start = 1000,
#'   region_end = 1500,
#'   min_mapq = 20,
#'   max_reads = 2000,
#'   max_reads_per_reference = FALSE,
#'   seed = 2026
#' )
#'
#' result$eligible_reads
#' result$selected_reads
#'
#' # Sample at most 2,000 records from each requested reference.
#' per_reference <- subset_bam_reads(
#'   bam = "barcode01.aligned.sorted.bam",
#'   output_bam = "barcode01.references.bam",
#'   reference_name = c("amplicon_A", "amplicon_B"),
#'   region_start = 1000,
#'   region_end = 1500,
#'   max_reads = 2000,
#'   max_reads_per_reference = TRUE,
#'   seed = 2026
#' )
#' per_reference$per_reference
#' }
#'
#' @export
subset_bam_reads <- function(bam,
                             output_bam,
                             reference_name,
                             region_start,
                             region_end,
                             min_mapq = 20L,
                             max_reads = 2000L,
                             max_reads_per_reference = FALSE,
                             coverage = c("complete", "overlap"),
                             primary_only = TRUE,
                             seed = 1L,
                             index_output = FALSE,
                             overwrite = FALSE) {
  check_file_arg(bam, "bam")
  check_scalar_character(output_bam, "output_bam")
  if (!is.character(reference_name) || length(reference_name) == 0L ||
      anyNA(reference_name) || any(!nzchar(trimws(reference_name)))) {
    stop("`reference_name` must contain one or more non-empty names.",
         call. = FALSE)
  }
  reference_name <- unique(trimws(reference_name))
  check_logical_scalar(max_reads_per_reference, "max_reads_per_reference")
  check_logical_scalar(primary_only, "primary_only")
  check_logical_scalar(index_output, "index_output")
  check_logical_scalar(overwrite, "overwrite")
  coverage <- match.arg(coverage)
  region_start <- validate_positive_integer(region_start, "region_start")
  region_end <- validate_positive_integer(region_end, "region_end")
  if (region_end < region_start) {
    stop("`region_end` must be greater than or equal to `region_start`.",
         call. = FALSE)
  }
  min_mapq <- validate_nonnegative_integer(min_mapq, "min_mapq")
  if (min_mapq > 255L) {
    stop("`min_mapq` cannot exceed 255.", call. = FALSE)
  }
  max_reads <- validate_positive_integer(max_reads, "max_reads")
  if (!is.null(seed)) {
    integer_seed <- suppressWarnings(as.integer(seed))
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed < 0 || is.na(integer_seed) ||
        seed != integer_seed) {
      stop("`seed` must be `NULL` or one non-negative integer.",
           call. = FALSE)
    }
    seed <- integer_seed
  }
  if (!grepl("\\.bam$", output_bam, ignore.case = TRUE)) {
    stop("`output_bam` must end in `.bam`.", call. = FALSE)
  }

  bam <- normalizePath(bam, mustWork = TRUE)
  output_bam <- normalizePath(output_bam, mustWork = FALSE)
  if (identical(bam, output_bam)) {
    stop("`output_bam` must differ from the input `bam`.", call. = FALSE)
  }
  output_indexes <- unique(c(
    paste0(output_bam, ".bai"),
    sub("\\.bam$", ".bai", output_bam, ignore.case = TRUE)
  ))
  existing <- c(output_bam, output_indexes)
  existing <- existing[file.exists(existing)]
  if (length(existing) > 0L && !overwrite) {
    stop("Output already exists; use `overwrite = TRUE` to replace it: ",
         existing[[1L]], call. = FALSE)
  }
  dir.create(dirname(output_bam), recursive = TRUE, showWarnings = FALSE)
  if (overwrite && length(existing) > 0L) unlink(existing)

  bam_flag <- if (primary_only) {
    Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    )
  } else {
    Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
  }
  record_fields <- c("qname", "rname", "pos", "cigar", "mapq", "flag")
  bam_param <- Rsamtools::ScanBamParam(
    flag = bam_flag,
    mapqFilter = min_mapq,
    what = record_fields
  )
  records <- Rsamtools::scanBam(bam, param = bam_param)[[1L]]
  record_count <- length(records$qname)

  if (record_count > 0L) {
    cigar <- as.character(records$cigar)
    valid_cigar <- !is.na(cigar) & nzchar(cigar) & cigar != "*"
    reference_length <- rep(NA_real_, record_count)
    reference_length[valid_cigar] <- cigar_reference_lengths(
      cigar[valid_cigar]
    )
    alignment_start <- as.integer(records$pos)
    alignment_end <- alignment_start + reference_length - 1
    reference_match <- !is.na(records$rname) &
      as.character(records$rname) %in% reference_name
    region_match <- if (coverage == "complete") {
      alignment_start <= region_start & alignment_end >= region_end
    } else {
      alignment_start <= region_end & alignment_end >= region_start
    }
    eligible <- which(
      valid_cigar & reference_match & !is.na(region_match) & region_match
    )
  } else {
    eligible <- integer()
  }

  eligible_reference <- as.character(records$rname[eligible])
  selected <- sample_bam_record_indices(
    eligible = eligible,
    eligible_reference = eligible_reference,
    max_reads = max_reads,
    max_reads_per_reference = max_reads_per_reference,
    seed = seed
  )
  sampled <- length(selected) < length(eligible)
  eligible_count <- tabulate(
    match(eligible_reference, reference_name),
    nbins = length(reference_name)
  )
  selected_count <- tabulate(
    match(as.character(records$rname[selected]), reference_name),
    nbins = length(reference_name)
  )
  per_reference <- data.frame(
    reference_name = reference_name,
    eligible_reads = as.integer(eligible_count),
    selected_reads = as.integer(selected_count),
    sampled = selected_count < eligible_count,
    stringsAsFactors = FALSE
  )

  record_key <- function(x) {
    paste(
      as.character(x$qname),
      as.character(x$rname),
      as.integer(x$pos),
      as.character(x$cigar),
      as.integer(x$mapq),
      as.integer(x$flag),
      sep = "\034"
    )
  }
  all_keys <- record_key(records)
  occurrence <- integer(length(all_keys))
  if (length(all_keys) > 0L) {
    occurrence <- stats::ave(seq_along(all_keys), all_keys, FUN = seq_along)
  }
  selected_ids <- paste(all_keys[selected], occurrence[selected], sep = "\035")
  seen <- new.env(hash = TRUE, parent = emptyenv())
  keep_selected <- function(x) {
    key <- record_key(x)
    current_occurrence <- integer(length(key))
    key_groups <- split(seq_along(key), key)
    for (current_key in names(key_groups)) {
      index <- key_groups[[current_key]]
      previous <- if (exists(current_key, envir = seen, inherits = FALSE)) {
        get(current_key, envir = seen, inherits = FALSE)
      } else {
        0L
      }
      current_occurrence[index] <- previous + seq_along(index)
      assign(
        current_key,
        previous + length(index),
        envir = seen
      )
    }
    paste(key, current_occurrence, sep = "\035") %in% selected_ids
  }
  filter_rules <- S4Vectors::FilterRules(list(selected = keep_selected))
  Rsamtools::filterBam(
    file = bam,
    destination = output_bam,
    index = character(),
    filter = filter_rules,
    indexDestination = index_output,
    param = bam_param
  )

  invisible(list(
    output_bam = output_bam,
    eligible_reads = length(eligible),
    selected_reads = length(selected),
    sampled = sampled,
    per_reference = per_reference,
    reference_name = reference_name,
    region_start = region_start,
    region_end = region_end,
    coverage = coverage,
    min_mapq = min_mapq,
    max_reads = max_reads,
    max_reads_per_reference = max_reads_per_reference
  ))
}

sample_bam_record_indices <- function(eligible,
                                      eligible_reference,
                                      max_reads,
                                      max_reads_per_reference,
                                      seed) {
  if (length(eligible) == 0L) return(integer())
  sampler <- function() {
    if (!max_reads_per_reference) {
      if (length(eligible) <= max_reads) return(eligible)
      return(eligible[sample.int(length(eligible), max_reads)])
    }
    groups <- split(eligible, eligible_reference)
    unname(unlist(lapply(groups, function(index) {
      if (length(index) <= max_reads) return(index)
      index[sample.int(length(index), max_reads)]
    }), use.names = FALSE))
  }
  if (is.null(seed)) return(sampler())
  random_seed_exists <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (random_seed_exists) {
    old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (random_seed_exists) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  sampler()
}
