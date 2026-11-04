#' Summarize genomic loci of BAM alignments
#'
#' `summarize_bam_alignment_loci()` identifies the reference sequences and
#' genomic intervals to which reads align. It is intended for diagnosing
#' off-target or non-specific amplicons in long-read data. Reference alignment
#' ends are calculated from reference-consuming CIGAR operations rather than
#' from query width, so soft clipping, insertions, and deletions are handled
#' correctly.
#'
#' @param bam Path to a BAM file containing reads aligned to a genomic reference.
#'   The file must contain a valid SAM header and aligned records with query
#'   name, reference name, leftmost reference position, CIGAR, flag, and MAPQ.
#'   Coordinate sorting and a BAM index are recommended for interoperability
#'   with `samtools` and IGV, but an index is not required by this function
#'   because the BAM is scanned sequentially. For finding non-specific wheat
#'   amplicons, align the reads to the complete wheat genome rather than only to
#'   the intended amplicon sequence.
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
#' The input BAM may contain unmapped, secondary, supplementary, duplicate, and
#' QC-fail records. Keeping an unfiltered source BAM is useful because filtering
#' choices can then be changed without repeating alignment. This function always
#' removes unmapped records; the remaining flag filters are controlled by
#' `primary_only`, `exclude_duplicates`, and `exclude_qcfail`. The MAPQ filter is
#' inclusive: `min_mapq = 20` means `MAPQ >= 20`, whereas the original R
#' expression `mapq > 20` is reproduced with `min_mapq = 21`.
#'
#' @section Preparing a wheat-genome BAM from ONT reads:
#'
#' To determine where non-specific PCR products originate, map the demultiplexed
#' ONT reads against a complete wheat reference genome. For example, the Chinese
#' Spring bread-wheat assembly IWGSC CS RefSeq v2.1 is available from NCBI as
#' assembly `GCF_018294505.1`. Assume the downloaded and decompressed FASTA is:
#'
#' ```text
#' GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
#' ```
#'
#' and the sample reads are:
#'
#' ```text
#' barcode150.fastq.gz
#' ```
#'
#' First build a reusable minimap2 index. This is done once for the wheat
#' reference and can require substantial memory and disk space because the
#' assembly is approximately 14 Gb:
#'
#' ```bash
#' minimap2 -d wheat_CS_RefSeq_v2.1.mmi \
#'   GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
#' ```
#'
#' Minimap2 normally indexes a limited number of target bases in one batch
#' (`-I`, commonly 8G by default). The wheat assembly is larger than one default
#' batch, so its `.mmi` can be a multi-part index. When mapping against such an
#' index, provide `--split-prefix` so minimap2 has a writable prefix for its
#' temporary intermediate files and can combine the index parts correctly.
#'
#' Align the ONT reads with the `map-ont` preset, convert the SAM stream to a
#' coordinate-sorted BAM, and create its index:
#'
#' ```bash
#' minimap2 -ax map-ont -t 32 \
#'   --split-prefix=barcode150.wheat.minimap2_split \
#'   -R '@RG\tID:barcode150\tSM:barcode150\tPL:ONT' \
#'   wheat_CS_RefSeq_v2.1.mmi barcode150.fastq.gz \
#' | samtools sort -@ 8 -m 2G -o barcode150.wheat.bam -
#'
#' samtools index -@ 8 barcode150.wheat.bam
#' ```
#'
#' `-t 32` controls minimap2 threads, while `samtools sort -@ 8 -m 2G` uses up
#' to 2 GB per sorting thread. Adjust both settings to the available computer.
#' `--split-prefix` is a filename prefix, not an output BAM name. Put it on a
#' fast local disk with sufficient free space when possible, for example
#' `/scratch/project/barcode150.minimap2_split`. It must be writable and must be
#' unique for every simultaneously running sample; reusing one prefix in
#' parallel jobs can corrupt or overwrite temporary files. Interrupted jobs may
#' leave files with this prefix, which can be removed after confirming that no
#' minimap2 process is still using them.
#'
#' The `-R` read-group field is optional for this function but records the sample
#' identity in the BAM. `samtools sort` is essential before indexing; the final
#' `-` tells it to read SAM/BAM data from the pipe.
#'
#' On a machine with enough RAM, an alternative is to rebuild a single-part
#' wheat index by setting `-I` above the reference size, for example:
#'
#' ```bash
#' minimap2 -I 20G -d wheat_CS_RefSeq_v2.1.single.mmi \
#'   GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
#' ```
#'
#' This can require substantially more memory. For typical servers, retaining
#' the multi-part index and using a unique `--split-prefix` is safer.
#'
#' Do not apply a high MAPQ filter while creating the archival BAM unless disk
#' space requires it. Bread wheat is hexaploid and contains many highly similar
#' A-, B-, and D-subgenome loci. Genuine reads from homologous or repetitive
#' regions can therefore receive MAPQ 0 or have secondary alignments. Retaining
#' these records allows comparison of `min_mapq = 0`, `20`, or `21` later. With
#' `primary_only = TRUE`, the function reports only the primary placement of
#' each alignment record; use `primary_only = FALSE` when investigating all
#' alternative placements.
#'
#' Validate and summarize the BAM before analysis:
#'
#' ```bash
#' samtools quickcheck -v barcode150.wheat.bam
#' samtools flagstat barcode150.wheat.bam
#' samtools idxstats barcode150.wheat.bam \
#'   | sort -k3,3nr \
#'   | head -30
#' ```
#'
#' `quickcheck` prints nothing when the BAM passes its basic checks. `flagstat`
#' reports mapped, secondary, supplementary, duplicate, and other flag counts.
#' In `idxstats`, the third column is the number of mapped alignment records for
#' each reference sequence. NCBI FASTA files commonly use accessions such as
#' `NC_...` as BAM reference names. Use those exact names in `reference_name`;
#' they can be inspected with `samtools idxstats` or in
#' `result$reference_summary`.
#'
#' If only a smaller diagnostic BAM is desired, mapped primary records with
#' strict `MAPQ > 20` can be written as follows:
#'
#' ```bash
#' samtools view -b -F 2308 -q 21 barcode150.wheat.bam \
#' | samtools sort -@ 8 -o barcode150.wheat.primary.mapq21.bam -
#' samtools index barcode150.wheat.primary.mapq21.bam
#' ```
#'
#' Here `-F 2308` excludes unmapped (`0x4`), secondary (`0x100`), and
#' supplementary (`0x800`) records, and `-q 21` means MAPQ at least 21. This
#' smaller file cannot later be used to inspect the excluded placements, so the
#' unfiltered `barcode150.wheat.bam` should normally be retained.
#'
#' @section Interpreting the returned wheat results:
#'
#' `reference_summary` ranks wheat chromosomes or unplaced scaffolds by retained
#' unique reads and alignment records. `peak_summary` then divides each reference
#' into `bin_width`-base intervals and ranks occupied intervals. A strong peak at
#' the expected target locus supports the intended amplicon; strong peaks on
#' other chromosomes or subgenomes suggest non-specific amplification or
#' homeologous targets. MAPQ should be interpreted together with read count:
#' low-MAPQ peaks may represent ambiguous homeologous placements rather than a
#' unique off-target locus.
#'
#' For a wheat-scale genome, avoid using very small histogram bins across whole
#' chromosomes. A 1-kb bin across hundreds of megabases can create hundreds of
#' thousands of empty plotting bins. `peak_summary` contains only occupied bins
#' and is therefore preferable for ranking loci; use a larger `bin_width` such
#' as `100000` or `1000000` for chromosome-wide plots, then rerun on selected
#' references or regions at higher resolution.
#'
#' @examples
#' \dontrun{
#' result <- summarize_bam_alignment_loci(
#'   bam = "barcode150.wheat.bam",
#'   min_mapq = 21,
#'   primary_only = TRUE,
#'   min_reference_reads = 51,
#'   bin_width = 100000
#' )
#'
#' result$filter_summary
#' head(result$mapq_distribution)
#' result$reference_summary
#' head(result$peak_summary, 30)
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
