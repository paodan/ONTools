#' Summarize linked mutation frequencies from a BAM file
#'
#' `summarize_linked_mutations()` counts read-level mutation patterns across a
#' set of candidate variants. Candidate variants define the sites to inspect;
#' the BAM determines which variants co-occur on the same read.
#'
#' @param bam Coordinate-sorted BAM file. A BAM index is recommended.
#' @param reference_fasta Optional reference FASTA. It is currently used only
#'   for input provenance in the return value; candidate `REF` alleles are taken
#'   from `variants`.
#' @param region Optional region string such as `"chr1"` or `"chr1:100-250"`.
#'   Reads outside this region are ignored. Variants are also restricted to this
#'   region when supplied.
#' @param variants Candidate variants as a data frame or VCF path. Data frames
#'   must contain `CHROM`, `POS`, `REF`, and `ALT` columns. VCF paths are read
#'   with [read_vcf()]. Multi-allelic `ALT` values are split into separate
#'   candidate variants.
#' @param min_mapq Minimum read mapping quality.
#' @param min_baseq Minimum base quality required for bases used to classify
#'   SNP/MNP and anchor bases for INDELs.
#' @param require_complete Logical. If `TRUE`, haplotype frequencies use only
#'   reads that can be classified at every candidate variant. If `FALSE`, reads
#'   with `NO_CALL` at one or more sites are retained in the denominator.
#' @param include_read_table Logical. If `TRUE`, include one row per read in the
#'   return value.
#' @param no_call_label Label for variants not covered or not classifiable in a
#'   read.
#'
#' @return A list with:
#'   \itemize{
#'     \item `haplotypes`: haplotype counts and frequencies.
#'     \item `read_table`: optional read-level classifications.
#'     \item `variants`: normalized candidate variants.
#'     \item `summary`: input and read-count summary.
#'   }
#'
#' @examples
#' variants <- data.frame(
#'   CHROM = "chr1",
#'   POS = c(10, 20),
#'   REF = c("A", "C"),
#'   ALT = c("G", "T")
#' )
#' # A real indexed BAM is needed for execution.
#' # summarize_linked_mutations("sample.bam", variants = variants)
#'
#' @export
summarize_linked_mutations <- function(bam,
                                       reference_fasta = NULL,
                                       region = NULL,
                                       variants,
                                       min_mapq = 20,
                                       min_baseq = 10,
                                       require_complete = TRUE,
                                       include_read_table = FALSE,
                                       no_call_label = "NO_CALL") {
  check_file_arg(bam, "bam")
  if (!is.null(reference_fasta)) check_file_arg(reference_fasta, "reference_fasta")
  if (!is.null(region)) check_scalar_character(region, "region")
  check_logical_scalar(require_complete, "require_complete")
  check_logical_scalar(include_read_table, "include_read_table")
  check_scalar_character(no_call_label, "no_call_label")
  min_mapq <- validate_nonnegative_integer(min_mapq, "min_mapq")
  min_baseq <- validate_nonnegative_integer(min_baseq, "min_baseq")

  bam <- normalizePath(bam, mustWork = TRUE)
  if (!is.null(reference_fasta)) {
    reference_fasta <- normalizePath(reference_fasta, mustWork = TRUE)
  }
  if (!has_standard_bam_index(bam)) {
    warning(
      "No standard BAM index found next to `bam`; scanBam may be slow or fail ",
      "for some inputs. Create one with samtools index if needed.",
      call. = FALSE
    )
  }

  region_info <- parse_linked_region(region)
  variants <- normalize_linked_variants(variants, region_info)
  if (nrow(variants) == 0L) {
    stop("No candidate variants remain after normalization/filtering.",
         call. = FALSE)
  }

  bam_data <- read_linked_bam_records(bam, min_mapq = min_mapq)
  total_records <- length(bam_data$qname)
  if (total_records == 0L) {
    haplotypes <- empty_linked_haplotype_table()
    return(list(
      haplotypes = haplotypes,
      read_table = if (isTRUE(include_read_table)) empty_linked_read_table() else NULL,
      variants = variants,
      summary = linked_summary(
        bam, reference_fasta, region, min_mapq, min_baseq,
        total_records = 0L, region_records = 0L, classified_reads = 0L,
        informative_reads = 0L
      )
    ))
  }

  read_rows <- vector("list", total_records)
  keep_read <- logical(total_records)
  for (i in seq_len(total_records)) {
    read <- linked_bam_record(bam_data, i)
    if (!linked_read_in_region(read, region_info)) {
      next
    }
    calls <- classify_read_variants(read, variants, min_baseq, no_call_label)
    has_no_call <- any(calls$call == no_call_label)
    if (isTRUE(require_complete) && has_no_call) {
      read_rows[[i]] <- linked_read_row(read, calls, no_call_label, informative = FALSE)
      keep_read[[i]] <- isTRUE(include_read_table)
    } else {
      read_rows[[i]] <- linked_read_row(read, calls, no_call_label, informative = !has_no_call)
      keep_read[[i]] <- TRUE
    }
  }

  read_table <- do.call(
    rbind,
    read_rows[keep_read & !vapply(read_rows, is.null, logical(1))]
  )
  if (is.null(read_table)) read_table <- empty_linked_read_table()

  haplo_input <- if (isTRUE(require_complete)) {
    read_table[!grepl(no_call_label, read_table$haplotype, fixed = TRUE), , drop = FALSE]
  } else {
    read_table
  }
  haplotypes <- summarize_linked_haplotype_table(haplo_input)

  list(
    haplotypes = haplotypes,
    read_table = if (isTRUE(include_read_table)) read_table else NULL,
    variants = variants,
    summary = linked_summary(
      bam, reference_fasta, region, min_mapq, min_baseq,
      total_records = total_records,
      region_records = sum(vapply(read_rows, Negate(is.null), logical(1))),
      classified_reads = nrow(read_table),
      informative_reads = nrow(haplo_input)
    )
  )
}

normalize_linked_variants <- function(variants, region_info = NULL) {
  if (is.character(variants) && length(variants) == 1L) {
    check_file_arg(variants, "variants")
    variants <- read_vcf(variants, parse_genotypes = FALSE)
  }
  if (is.null(variants)) return(empty_linked_variants())
  if (!is.data.frame(variants)) {
    stop("`variants` must be a data frame or a VCF file path.", call. = FALSE)
  }
  required <- c("CHROM", "POS", "REF", "ALT")
  missing_cols <- setdiff(required, names(variants))
  if (length(missing_cols) > 0L) {
    stop("`variants` is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  variants <- variants[, required, drop = FALSE]
  variants$CHROM <- as.character(variants$CHROM)
  variants$POS <- as.integer(variants$POS)
  variants$REF <- toupper(as.character(variants$REF))
  variants$ALT <- toupper(as.character(variants$ALT))
  variants <- variants[!is.na(variants$POS) & nzchar(variants$ALT) & variants$ALT != ".", , drop = FALSE]

  rows <- lapply(seq_len(nrow(variants)), function(i) {
    alts <- strsplit(variants$ALT[[i]], ",", fixed = TRUE)[[1L]]
    data.frame(
      CHROM = variants$CHROM[[i]],
      POS = variants$POS[[i]],
      REF = variants$REF[[i]],
      ALT = alts,
      stringsAsFactors = FALSE
    )
  })
  out <- if (length(rows) == 0L) empty_linked_variants() else do.call(rbind, rows)
  out$variant_id <- paste0(out$CHROM, ":", out$POS, ":", out$REF, ">", out$ALT)
  out$type <- vcf_variant_type(out$REF, out$ALT)
  out$ref_end <- out$POS + nchar(out$REF) - 1L

  if (!is.null(region_info)) {
    out <- out[out$CHROM == region_info$reference, , drop = FALSE]
    if (!is.na(region_info$start)) {
      out <- out[out$ref_end >= region_info$start, , drop = FALSE]
    }
    if (!is.na(region_info$end)) {
      out <- out[out$POS <= region_info$end, , drop = FALSE]
    }
  }
  rownames(out) <- NULL
  out
}

empty_linked_variants <- function() {
  data.frame(
    CHROM = character(), POS = integer(), REF = character(), ALT = character(),
    variant_id = character(), type = character(), ref_end = integer(),
    stringsAsFactors = FALSE
  )
}

parse_linked_region <- function(region) {
  if (is.null(region)) return(NULL)
  m <- regexec("^([^:]+)(?::([0-9,]+)-([0-9,]+))?$", region)
  parts <- regmatches(region, m)[[1L]]
  if (length(parts) == 0L) {
    stop("`region` must look like 'chr1' or 'chr1:100-250'.", call. = FALSE)
  }
  start <- if (length(parts) >= 3L && nzchar(parts[[3L]])) as.integer(gsub(",", "", parts[[3L]])) else NA_integer_
  end <- if (length(parts) >= 4L && nzchar(parts[[4L]])) as.integer(gsub(",", "", parts[[4L]])) else NA_integer_
  if (!is.na(start) && !is.na(end) && start > end) {
    stop("`region` start must be <= end.", call. = FALSE)
  }
  list(reference = parts[[2L]], start = start, end = end)
}

read_linked_bam_records <- function(bam, min_mapq) {
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    ),
    what = c("qname", "rname", "pos", "cigar", "seq", "qual", "mapq", "strand")
  )
  x <- Rsamtools::scanBam(bam, param = param)[[1L]]
  keep <- !is.na(x$mapq) & x$mapq >= min_mapq
  lapply(x, function(col) col[keep])
}

linked_bam_record <- function(x, i) {
  seq <- as.character(x$seq[[i]])
  qual <- as.integer(x$qual[[i]])
  list(
    qname = x$qname[[i]],
    rname = as.character(x$rname[[i]]),
    pos = as.integer(x$pos[[i]]),
    cigar = x$cigar[[i]],
    seq = toupper(seq),
    qual = qual,
    mapq = x$mapq[[i]],
    strand = as.character(x$strand[[i]])
  )
}

linked_read_in_region <- function(read, region_info) {
  if (is.null(region_info)) return(TRUE)
  if (!identical(read$rname, region_info$reference)) return(FALSE)
  end <- linked_cigar_ref_end(read$pos, read$cigar)
  if (!is.na(region_info$start) && end < region_info$start) return(FALSE)
  if (!is.na(region_info$end) && read$pos > region_info$end) return(FALSE)
  TRUE
}

linked_cigar_ref_end <- function(pos, cigar) {
  pos + sum_cigar_ops(cigar, "[MDN=X]") - 1L
}

sum_cigar_ops <- function(cigar, ops_pattern) {
  parts <- regmatches(cigar, gregexpr(paste0("[0-9]+", ops_pattern), cigar))[[1L]]
  if (length(parts) == 0L || identical(parts, character(0))) return(0L)
  sum(as.integer(sub("[A-Z=]$", "", parts)))
}

classify_read_variants <- function(read, variants, min_baseq, no_call_label) {
  map <- build_read_reference_map(read)
  calls <- vapply(seq_len(nrow(variants)), function(i) {
    classify_read_variant(read, map, variants[i, , drop = FALSE], min_baseq, no_call_label)
  }, character(1))
  data.frame(
    variant_id = variants$variant_id,
    call = calls,
    stringsAsFactors = FALSE
  )
}

build_read_reference_map <- function(read) {
  cigar <- read$cigar
  lengths <- as.integer(unlist(regmatches(cigar, gregexpr("[0-9]+", cigar))))
  ops <- unlist(regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar)))
  ref_pos <- read$pos
  query_pos <- 1L
  ref_to_query <- integer()
  ref_to_query_qual <- integer()
  deletions <- integer()
  insertions <- list()
  last_ref <- NA_integer_

  for (i in seq_along(ops)) {
    op <- ops[[i]]
    n <- lengths[[i]]
    if (op %in% c("M", "=", "X")) {
      refs <- seq.int(ref_pos, length.out = n)
      queries <- seq.int(query_pos, length.out = n)
      ref_to_query[as.character(refs)] <- queries
      ref_to_query_qual[as.character(refs)] <- read$qual[queries]
      last_ref <- refs[[length(refs)]]
      ref_pos <- ref_pos + n
      query_pos <- query_pos + n
    } else if (op == "I") {
      ins_seq <- substr(read$seq, query_pos, query_pos + n - 1L)
      ins_qual <- read$qual[seq.int(query_pos, length.out = n)]
      anchor <- if (is.na(last_ref)) ref_pos - 1L else last_ref
      insertions[[as.character(anchor)]] <- list(seq = ins_seq, qual = ins_qual)
      query_pos <- query_pos + n
    } else if (op %in% c("D", "N")) {
      refs <- seq.int(ref_pos, length.out = n)
      deletions <- c(deletions, refs)
      ref_pos <- ref_pos + n
      last_ref <- refs[[length(refs)]]
    } else if (op %in% c("S")) {
      query_pos <- query_pos + n
    } else if (op %in% c("H", "P")) {
      next
    }
  }

  list(
    ref_to_query = ref_to_query,
    ref_to_query_qual = ref_to_query_qual,
    deletions = deletions,
    insertions = insertions
  )
}

classify_read_variant <- function(read, map, variant, min_baseq, no_call_label) {
  ref <- variant$REF[[1L]]
  alt <- variant$ALT[[1L]]
  pos <- variant$POS[[1L]]
  if (nchar(ref) == nchar(alt)) {
    return(classify_read_substitution(read, map, pos, ref, alt, min_baseq, no_call_label))
  }
  if (startsWith(ref, substr(alt, 1L, 1L)) && nchar(ref) > nchar(alt)) {
    return(classify_read_deletion(read, map, pos, ref, alt, min_baseq, no_call_label))
  }
  if (startsWith(alt, ref) && nchar(alt) > nchar(ref)) {
    return(classify_read_insertion(read, map, pos, ref, alt, min_baseq, no_call_label))
  }
  no_call_label
}

classify_read_substitution <- function(read, map, pos, ref, alt, min_baseq, no_call_label) {
  refs <- seq.int(pos, length.out = nchar(ref))
  qpos <- map$ref_to_query[as.character(refs)]
  qqual <- map$ref_to_query_qual[as.character(refs)]
  if (any(is.na(qpos)) || any(is.na(qqual)) || any(qqual < min_baseq)) return(no_call_label)
  allele <- paste0(substr(read$seq, qpos, qpos), collapse = "")
  if (identical(allele, ref)) return("REF")
  if (identical(allele, alt)) return("ALT")
  paste0("OTHER:", allele)
}

classify_read_insertion <- function(read, map, pos, ref, alt, min_baseq, no_call_label) {
  anchor_qpos <- unname(map$ref_to_query[as.character(pos)])
  anchor_qual <- unname(map$ref_to_query_qual[as.character(pos)])
  if (length(anchor_qpos) == 0L || is.na(anchor_qpos) ||
      length(anchor_qual) == 0L || is.na(anchor_qual) || anchor_qual < min_baseq) {
    return(no_call_label)
  }
  anchor_base <- substr(read$seq, anchor_qpos, anchor_qpos)
  if (!identical(anchor_base, ref)) return(paste0("OTHER:", anchor_base))

  ins <- map$insertions[[as.character(pos)]]
  expected_insert <- substr(alt, nchar(ref) + 1L, nchar(alt))
  if (is.null(ins)) return("REF")
  if (any(ins$qual < min_baseq)) return(no_call_label)
  if (identical(ins$seq, expected_insert)) return("ALT")
  paste0("OTHER:", paste0(ref, ins$seq))
}

classify_read_deletion <- function(read, map, pos, ref, alt, min_baseq, no_call_label) {
  anchor_qpos <- unname(map$ref_to_query[as.character(pos)])
  anchor_qual <- unname(map$ref_to_query_qual[as.character(pos)])
  if (length(anchor_qpos) == 0L || is.na(anchor_qpos) ||
      length(anchor_qual) == 0L || is.na(anchor_qual) || anchor_qual < min_baseq) {
    return(no_call_label)
  }
  anchor_base <- substr(read$seq, anchor_qpos, anchor_qpos)
  if (!identical(anchor_base, substr(ref, 1L, 1L))) return(paste0("OTHER:", anchor_base))

  deleted_refs <- seq.int(pos + nchar(alt), pos + nchar(ref) - 1L)
  if (all(deleted_refs %in% map$deletions)) return("ALT")

  refs <- seq.int(pos, length.out = nchar(ref))
  qpos <- map$ref_to_query[as.character(refs)]
  qqual <- map$ref_to_query_qual[as.character(refs)]
  if (!any(is.na(qpos)) && !any(is.na(qqual)) && all(qqual >= min_baseq)) {
    allele <- paste0(substr(read$seq, qpos, qpos), collapse = "")
    if (identical(allele, ref)) return("REF")
    return(paste0("OTHER:", allele))
  }
  no_call_label
}

linked_read_row <- function(read, calls, no_call_label, informative) {
  alt_ids <- calls$variant_id[calls$call == "ALT"]
  other_ids <- paste0(calls$variant_id[calls$call != "REF" & calls$call != "ALT"], "=", calls$call[calls$call != "REF" & calls$call != "ALT"])
  mutation_pattern <- c(alt_ids, other_ids)
  if (length(mutation_pattern) == 0L) mutation_pattern <- "WT"
  haplotype <- paste(mutation_pattern, collapse = ";")
  data.frame(
    read_id = read$qname,
    reference = read$rname,
    read_start = read$pos,
    read_end = linked_cigar_ref_end(read$pos, read$cigar),
    mapq = read$mapq,
    haplotype = haplotype,
    calls = paste(paste0(calls$variant_id, "=", calls$call), collapse = ";"),
    informative = informative,
    stringsAsFactors = FALSE
  )
}

summarize_linked_haplotype_table <- function(read_table) {
  if (nrow(read_table) == 0L) return(empty_linked_haplotype_table())
  tab <- as.data.frame(table(read_table$haplotype), stringsAsFactors = FALSE)
  names(tab) <- c("haplotype", "read_count")
  tab <- tab[order(-tab$read_count, tab$haplotype), , drop = FALSE]
  total <- sum(tab$read_count)
  tab$frequency <- tab$read_count / total
  tab$frequency_percent <- tab$frequency * 100
  rownames(tab) <- NULL
  tab
}

empty_linked_haplotype_table <- function() {
  data.frame(
    haplotype = character(), read_count = integer(),
    frequency = numeric(), frequency_percent = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_linked_read_table <- function() {
  data.frame(
    read_id = character(), reference = character(), read_start = integer(),
    read_end = integer(), mapq = integer(), haplotype = character(),
    calls = character(), informative = logical(),
    stringsAsFactors = FALSE
  )
}

linked_summary <- function(bam, reference_fasta, region, min_mapq, min_baseq,
                           total_records, region_records, classified_reads,
                           informative_reads) {
  data.frame(
    bam = bam,
    reference_fasta = if (is.null(reference_fasta)) NA_character_ else reference_fasta,
    region = if (is.null(region)) NA_character_ else region,
    min_mapq = min_mapq,
    min_baseq = min_baseq,
    total_primary_records_after_mapq = total_records,
    region_records = region_records,
    classified_reads = classified_reads,
    informative_reads = informative_reads,
    stringsAsFactors = FALSE
  )
}
