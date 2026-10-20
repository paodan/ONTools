#' Create a synthetic AB1 trace from a consensus FASTA and BAM
#'
#' `synthetic_ab1_from_bam()` writes a minimal ABIF/AB1 file with simulated
#' Sanger-style four-channel peak data. It uses `samtools mpileup` to estimate
#' A/C/G/T support at each consensus position, then turns those counts into
#' Gaussian peaks using a selectable peak-value mode.
#'
#' This function does not reconstruct real capillary-sequencing fluorescence
#' signal from the BAM. Use the result for visualization or compatibility with
#' tools that require an AB1-like chromatogram, not as raw Sanger evidence.
#'
#' @param consensus Consensus FASTA used as the BAM reference. The sequence
#'   names must match the reference names in `bam`. If the FASTA has no `.fai`
#'   index, the function copies it to a temporary directory and runs
#'   `samtools faidx` there. A single generated AB1 file can contain one
#'   consensus sequence; the current writer supports sequences up to 16,000 bp.
#' @param bam Coordinate-sorted BAM aligned to `consensus`. A standard BAM
#'   index (`.bai` or `.csi`) should exist next to the BAM; the function warns
#'   if no index is found, because `samtools mpileup` may fail without it.
#' @param output_ab1 Output `.ab1` file. Parent directories are created
#'   automatically. Existing files are protected unless `overwrite = TRUE`.
#' @param reference_name FASTA/BAM reference name. Required when `consensus`
#'   contains more than one sequence. Default `NULL` means the only FASTA record
#'   is used automatically.
#' @param sample Sample name stored in the AB1 `SMPL1` tag. Defaults to the
#'   output filename without extension.
#' @param samtools Command name or full path used to launch `samtools`. Default
#'   is `"samtools"`.
#' @param conda_env Optional conda environment name. If supplied, `samtools` is
#'   run with `conda run -n <conda_env>`. Default `NULL` calls `samtools` from
#'   the current environment.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#'   Default is `"conda"`.
#' @param min_base_quality Minimum base quality passed to `samtools mpileup -Q`.
#'   Default `0` keeps all bases reported by mpileup.
#' @param spacing Number of trace points between adjacent called bases. Default
#'   `12` gives separated Sanger-like peaks.
#' @param sigma Gaussian peak width in trace points. Default is `2`; larger
#'   values make peaks wider and more overlapped.
#' @param peak_height Peak-height scale used by `peak_value_mode = "relative"`
#'   and `"fraction"`. It is ignored by `"count"`, where peak amplitude above
#'   `baseline` equals the A/C/G/T supporting-read count. Default is `900`.
#' @param baseline Baseline signal added to every A/C/G/T channel. Default is
#'   `15`.
#' @param noise_fraction Minimum off-channel peak fraction used only by
#'   `peak_value_mode = "relative"`. Default `0.03` keeps non-major channels
#'   visible at low height; must be between 0 and 1. Count and fraction modes
#'   do not add artificial off-channel peaks.
#' @param peak_value_mode Meaning of peak amplitude. `"relative"` (default)
#'   divides each A/C/G/T count by the largest count at that position, so the
#'   major peak has amplitude `peak_height`. `"count"` uses the raw supporting
#'   read count as amplitude. `"fraction"` divides each count by the sum of
#'   A/C/G/T counts and multiplies by `peak_height`. The baseline is added after
#'   calculating these amplitudes.
#' @param overwrite Logical. If `FALSE` (default), stop when `output_ab1`
#'   already exists. Set `TRUE` to replace it.
#' @param echo Logical. If `TRUE` (default), print the `samtools mpileup`
#'   command.
#' @param stderr Passed to [system2()]. Default `""` streams errors to the R
#'   console.
#'
#' @return Invisibly returns a list with `status`, `paths`, `reference_name`,
#'   `called_sequence`, `sequence_length`, `trace_length`, `peak_value_mode`,
#'   and `command`.
#'
#' @details
#' The function runs `samtools mpileup` with `-aa`, so positions with zero read
#' coverage are still represented. For each consensus position, pileup symbols
#' are reduced to A/C/G/T counts. At covered positions, the most-supported
#' pileup base determines the `PBAS` call. Peak amplitudes are calculated using
#' `peak_value_mode`. At tied positions, the FASTA base is preferred when it is
#' among the tied bases. The FASTA base is used as a `PBAS` fallback only when
#' no A/C/G/T pileup observation is available. In relative mode, that fallback
#' also receives a full-height peak. In count and fraction modes, an uncovered
#' position contains baseline only, because there are no supporting reads.
#'
#' Count and fraction modes currently use the sum of observed A/C/G/T bases.
#' Insertions and deletions are not represented by separate AB1 channels and
#' are therefore not included in those four peak amplitudes. In count mode,
#' values above the signed 16-bit ABIF trace limit are clipped to `32767` with
#' a warning. Use `baseline = 0` when the stored trace value itself should be
#' as close as possible to the supporting-read count, or `baseline = 1` to
#' retain a visible minimum line in regions without A/C/G/T observations.
#'
#' The AB1 file contains common ABIF tags used by many viewers, including
#' sequence calls (`PBAS`), peak locations (`PLOC`), base qualities (`PCON`),
#' sample name (`SMPL`), and four trace channels (`DATA9`-`DATA12`). The quality
#' values are synthetic and are estimated from the pileup support for the
#' final `PBAS` call; they are not original sequencer quality values.
#'
#' Typical use in ONTools delivery workflows is to create a Sanger-style visual
#' companion for `barcode*/alignments/*.bam` plus the matching consensus or
#' reference FASTA. The output is useful for review and customer-facing trace
#' visualization, but it should be labeled as synthetic when delivered.
#'
#' @examples
#' consensus <- tempfile(fileext = ".fasta")
#' bam <- tempfile(fileext = ".bam")
#' ab1 <- tempfile(fileext = ".ab1")
#' writeLines(c(">contig1", "ACGTACGT"), consensus)
#' # A real coordinate-sorted BAM and its index are needed for execution.
#' # synthetic_ab1_from_bam(consensus, bam, ab1)
#' # synthetic_ab1_from_bam(
#' #   consensus, bam, ab1,
#' #   peak_value_mode = "count", baseline = 1
#' # )
#'
#' @export
synthetic_ab1_from_bam <- function(consensus,
                                   bam,
                                   output_ab1,
                                   reference_name = NULL,
                                   sample = NULL,
                                   samtools = "samtools",
                                   conda_env = NULL,
                                   conda = "conda",
                                   min_base_quality = 0,
                                   spacing = 12,
                                   sigma = 2,
                                   peak_height = 900,
                                   baseline = 15,
                                   noise_fraction = 0.03,
                                   peak_value_mode = c(
                                     "relative", "count", "fraction"
                                   ),
                                   overwrite = FALSE,
                                   echo = TRUE,
                                   stderr = "") {
  check_file_arg(consensus, "consensus")
  check_file_arg(bam, "bam")
  check_scalar_character(output_ab1, "output_ab1")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  if (!is.null(reference_name)) check_scalar_character(reference_name, "reference_name")
  if (!is.null(sample)) check_scalar_character(sample, "sample")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(echo, "echo")
  peak_value_mode <- match.arg(peak_value_mode)

  min_base_quality <- validate_nonnegative_integer(min_base_quality, "min_base_quality")
  spacing <- validate_positive_integer(spacing, "spacing")
  sigma <- validate_positive_number(sigma, "sigma")
  peak_height <- validate_positive_integer(peak_height, "peak_height")
  baseline <- validate_nonnegative_integer(baseline, "baseline")
  noise_fraction <- validate_nonnegative_number(noise_fraction, "noise_fraction")
  if (noise_fraction > 1) {
    stop("`noise_fraction` must be between 0 and 1.", call. = FALSE)
  }

  consensus <- normalizePath(consensus, mustWork = TRUE)
  bam <- normalizePath(bam, mustWork = TRUE)
  if (!has_standard_bam_index(bam)) {
    warning(
      "No standard BAM index found next to `bam`; samtools mpileup may fail. ",
      "Create one with samtools index if needed.",
      call. = FALSE
    )
  }
  dir.create(dirname(output_ab1), recursive = TRUE, showWarnings = FALSE)
  output_ab1 <- normalizePath(output_ab1, mustWork = FALSE)
  if (file.exists(output_ab1) && !isTRUE(overwrite)) {
    stop("`output_ab1` already exists and `overwrite` is FALSE: ", output_ab1,
         call. = FALSE)
  }

  fasta <- read_single_consensus_sequence(consensus, reference_name)
  reference_name <- fasta$name
  sequence <- fasta$sequence
  if (!nzchar(sequence)) {
    stop("Consensus sequence is empty.", call. = FALSE)
  }
  if (nchar(sequence) > 16000L) {
    stop("This minimal AB1 writer supports consensus sequences up to 16,000 bases.",
         call. = FALSE)
  }

  indexed <- prepare_indexed_consensus_fasta(consensus, samtools, conda_env, conda, stderr)
  on.exit(if (!is.null(indexed$tmpdir)) unlink(indexed$tmpdir, recursive = TRUE), add = TRUE)

  pileup_call <- dehost_fastq_external_call(
    samtools,
    c(
      "mpileup",
      "-aa",
      "-Q", as.character(min_base_quality),
      "-f", indexed$fasta,
      "-r", reference_name,
      bam
    ),
    conda_env,
    conda
  )
  command <- paste(c(shQuote(pileup_call$command), shQuote(pileup_call$args)), collapse = " ")
  if (isTRUE(echo)) message(command)

  if (is.null(conda_env)) {
    require_external_command(samtools)
  } else {
    require_external_command(conda)
  }

  pileup_text <- system2(
    pileup_call$command,
    args = pileup_call$args,
    stdout = TRUE,
    stderr = stderr
  )
  status <- attr(pileup_text, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("samtools mpileup failed with exit status: ", status, call. = FALSE)
  }

  counts <- parse_mpileup_counts(pileup_text)
  trace <- synthetic_ab1_trace(
    sequence = sequence,
    pileup_counts = counts,
    spacing = spacing,
    sigma = sigma,
    peak_height = peak_height,
    baseline = baseline,
    noise_fraction = noise_fraction,
    peak_value_mode = peak_value_mode
  )
  if (is.null(sample)) {
    sample <- sub("[.][^.]+$", "", basename(output_ab1))
  }
  write_abif_file(
    output_ab1,
    sequence = trace$sequence,
    channels = trace$channels,
    peak_locations = trace$peak_locations,
    qualities = trace$qualities,
    sample = sample
  )

  invisible(list(
    status = 0L,
    paths = list(consensus = consensus, bam = bam, ab1 = output_ab1),
    reference_name = reference_name,
    called_sequence = trace$sequence,
    sequence_length = nchar(sequence),
    trace_length = length(trace$channels$A),
    peak_value_mode = peak_value_mode,
    command = command
  ))
}

read_single_consensus_sequence <- function(consensus, reference_name) {
  records <- Biostrings::readDNAStringSet(consensus)
  if (length(records) < 1L) {
    stop("No FASTA records found in `consensus`.", call. = FALSE)
  }
  record_names <- names(records)
  if (is.null(record_names)) {
    record_names <- rep("", length(records))
  }

  if (!is.null(reference_name)) {
    idx <- match(reference_name, record_names)
    if (is.na(idx)) {
      stop("`reference_name` was not found in `consensus`: ", reference_name,
           call. = FALSE)
    }
  } else {
    if (length(records) > 1L) {
      stop("`consensus` contains multiple records; pass `reference_name`.",
           call. = FALSE)
    }
    idx <- 1L
    reference_name <- record_names[[idx]]
    if (!nzchar(reference_name)) reference_name <- "reference"
  }

  list(name = reference_name, sequence = toupper(as.character(records[[idx]])))
}

prepare_indexed_consensus_fasta <- function(consensus, samtools, conda_env, conda, stderr) {
  if (file.exists(paste0(consensus, ".fai"))) {
    return(list(fasta = consensus, tmpdir = NULL))
  }

  tmpdir <- tempfile("synthetic-ab1-faidx-")
  dir.create(tmpdir)
  tmp_fasta <- file.path(tmpdir, basename(consensus))
  file.copy(consensus, tmp_fasta, overwrite = TRUE)

  faidx_call <- dehost_fastq_external_call(
    samtools,
    c("faidx", tmp_fasta),
    conda_env,
    conda
  )
  faidx_status <- system2(
    faidx_call$command,
    args = faidx_call$args,
    stdout = FALSE,
    stderr = stderr
  )
  if (!identical(faidx_status, 0L)) {
    unlink(tmpdir, recursive = TRUE)
    stop("samtools faidx failed with exit status: ", faidx_status, call. = FALSE)
  }

  list(fasta = tmp_fasta, tmpdir = tmpdir)
}

parse_mpileup_counts <- function(lines) {
  counts <- list()
  if (length(lines) == 0L) return(counts)

  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 5L) next
    pos <- as.integer(fields[[2]])
    ref <- toupper(fields[[3]])
    counts[[as.character(pos)]] <- parse_pileup_bases(fields[[5]], ref)
  }
  counts
}

parse_pileup_bases <- function(base_string, ref_base) {
  counts <- stats::setNames(integer(4L), c("A", "C", "G", "T"))
  chars <- strsplit(base_string, "", fixed = TRUE)[[1]]
  i <- 1L
  while (i <= length(chars)) {
    char <- chars[[i]]
    if (identical(char, "^")) {
      i <- i + 2L
      next
    }
    if (identical(char, "$")) {
      i <- i + 1L
      next
    }
    if (char %in% c("+", "-")) {
      i <- i + 1L
      digits <- character()
      while (i <= length(chars) && grepl("^[0-9]$", chars[[i]])) {
        digits <- c(digits, chars[[i]])
        i <- i + 1L
      }
      indel_len <- if (length(digits)) as.integer(paste0(digits, collapse = "")) else 0L
      i <- i + indel_len
      next
    }
    if (char %in% c(".", ",")) {
      if (ref_base %in% names(counts)) counts[[ref_base]] <- counts[[ref_base]] + 1L
      i <- i + 1L
      next
    }
    upper <- toupper(char)
    if (upper %in% names(counts)) counts[[upper]] <- counts[[upper]] + 1L
    i <- i + 1L
  }
  counts
}

synthetic_ab1_trace <- function(sequence,
                                pileup_counts,
                                spacing,
                                sigma,
                                peak_height,
                                baseline,
                                noise_fraction,
                                peak_value_mode = c(
                                  "relative", "count", "fraction"
                                )) {
  peak_value_mode <- match.arg(peak_value_mode)
  bases <- c("A", "C", "G", "T")
  n_bases <- nchar(sequence)
  seq_chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  called_chars <- pileup_called_bases(seq_chars, pileup_counts)
  left_pad <- spacing * 2L
  right_pad <- spacing * 2L
  peak_locations <- left_pad + ((seq_len(n_bases) - 1L) * spacing)
  trace_length <- peak_locations[[n_bases]] + right_pad + 1L
  channels <- stats::setNames(
    replicate(4L, rep.int(baseline, trace_length), simplify = FALSE),
    bases
  )

  radius <- max(2L, ceiling(sigma * 4))
  for (idx in seq_len(n_bases)) {
    values <- ab1_base_values(
      pileup_counts[[as.character(idx)]],
      seq_chars[[idx]],
      noise_fraction,
      peak_value_mode
    )
    center <- peak_locations[[idx]]
    for (base in bases) {
      height <- if (peak_value_mode == "count") {
        values[[base]]
      } else {
        peak_height * values[[base]]
      }
      positions <- seq.int(max(1L, center - radius + 1L), min(trace_length, center + radius + 1L))
      distances <- (positions - 1L) - center
      channels[[base]][positions] <- channels[[base]][positions] +
        round(height * exp(-(distances * distances) / (2 * sigma * sigma)))
    }
  }

  clipped <- any(vapply(channels, max, numeric(1), na.rm = TRUE) > 32767)
  if (clipped) {
    warning(
      "One or more synthetic AB1 trace values exceeded 32767 and were clipped.",
      call. = FALSE
    )
  }
  channels <- lapply(
    channels,
    function(x) pmin(32767L, pmax(0L, as.integer(x)))
  )
  list(
    sequence = paste0(called_chars, collapse = ""),
    channels = channels,
    peak_locations = as.integer(peak_locations),
    qualities = pileup_ab1_qualities(called_chars, pileup_counts)
  )
}

ab1_base_values <- function(counts,
                            consensus_base,
                            noise_fraction,
                            peak_value_mode = c(
                              "relative", "count", "fraction"
                            )) {
  peak_value_mode <- match.arg(peak_value_mode)
  bases <- c("A", "C", "G", "T")
  values <- stats::setNames(rep(0, 4L), bases)
  if (!is.null(counts)) {
    values[names(counts)] <- as.numeric(counts)
  }
  total <- sum(values)
  if (total > 0) {
    if (peak_value_mode == "count") return(values)
    if (peak_value_mode == "fraction") return(values / total)
    relative <- values / max(values)
    relative <- pmax(relative, noise_fraction)
    return(relative / max(relative))
  }
  if (peak_value_mode != "relative") return(values)
  consensus_base <- toupper(consensus_base)
  if (consensus_base %in% bases) {
    values[[consensus_base]] <- 1
  } else {
    values[] <- noise_fraction
  }
  pmax(values, noise_fraction)
}

pileup_called_bases <- function(fallback_bases, pileup_counts) {
  bases <- c("A", "C", "G", "T")
  vapply(seq_along(fallback_bases), function(idx) {
    fallback <- toupper(fallback_bases[[idx]])
    counts <- pileup_counts[[as.character(idx)]]
    if (is.null(counts)) return(fallback)
    values <- stats::setNames(rep(0, length(bases)), bases)
    values[names(counts)] <- as.numeric(counts)
    if (sum(values) <= 0) return(fallback)
    candidates <- names(values)[values == max(values)]
    if (fallback %in% candidates) return(fallback)
    candidates[[1L]]
  }, character(1))
}

pileup_ab1_qualities <- function(seq_chars, pileup_counts) {
  qualities <- integer(length(seq_chars))
  bases <- c("A", "C", "G", "T")
  for (idx in seq_along(seq_chars)) {
    base <- toupper(seq_chars[[idx]])
    counts <- pileup_counts[[as.character(idx)]]
    if (is.null(counts) || !(base %in% bases) || sum(counts) == 0) {
      qualities[[idx]] <- 20L
      next
    }
    support <- counts[[base]] / sum(counts)
    error_probability <- max(0.001, 1 - support)
    qualities[[idx]] <- as.integer(min(60, max(2, round(-10 * log10(error_probability)))))
  }
  qualities
}

write_abif_file <- function(path, sequence, channels, peak_locations, qualities, sample) {
  entries <- list(
    abif_entry("FWO_", 1L, "char", "GATC"),
    abif_entry("PBAS", 1L, "char", sequence),
    abif_entry("PBAS", 2L, "char", sequence),
    abif_entry("PLOC", 1L, "short", peak_locations),
    abif_entry("PLOC", 2L, "short", peak_locations),
    abif_entry("PCON", 1L, "byte", qualities),
    abif_entry("PCON", 2L, "byte", qualities),
    abif_entry("SMPL", 1L, "char", substr(sample, 1L, 255L)),
    abif_entry("DATA", 9L, "short", channels$G),
    abif_entry("DATA", 10L, "short", channels$A),
    abif_entry("DATA", 11L, "short", channels$T),
    abif_entry("DATA", 12L, "short", channels$C)
  )
  order_idx <- order(vapply(entries, `[[`, character(1), "tag"),
                     vapply(entries, `[[`, integer(1), "number"))
  entries <- entries[order_idx]

  header_size <- 34L
  directory_offset <- 128L
  directory_size <- length(entries) * 28L
  data_offset <- directory_offset + directory_size
  payload <- raw()

  for (i in seq_along(entries)) {
    data <- entries[[i]]$data
    if (length(data) <= 4L) {
      entries[[i]]$offset_or_data <- c(data, raw(4L - length(data)))
    } else {
      if ((data_offset + length(payload)) %% 2L == 1L) payload <- c(payload, as.raw(0))
      entries[[i]]$offset_or_data <- write_uint32_be(data_offset + length(payload))
      payload <- c(payload, data)
    }
  }

  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("ABIF"), con)
  writeBin(as.integer(101L), con, size = 2L, endian = "big")
  writeBin(pack_abif_directory_entry("tdir", 1L, 1023L, 28L, length(entries),
                                     directory_size, directory_offset, 0L), con)
  writeBin(raw(directory_offset - header_size), con)
  for (entry in entries) {
    writeBin(pack_abif_directory_entry(
      entry$tag,
      entry$number,
      entry$type,
      entry$element_size,
      entry$count,
      length(entry$data),
      entry$offset_or_data,
      0L
    ), con)
  }
  writeBin(payload, con)
}

abif_entry <- function(tag, number, kind, values) {
  packed <- pack_abif_value(kind, values)
  list(
    tag = tag,
    number = as.integer(number),
    type = packed$type,
    element_size = packed$element_size,
    count = packed$count,
    data = packed$data
  )
}

pack_abif_value <- function(kind, values) {
  if (identical(kind, "char")) {
    data <- charToRaw(paste0(values, collapse = ""))
    return(list(type = 2L, element_size = 1L, count = length(data), data = data))
  }
  if (identical(kind, "byte")) {
    values <- as.integer(values)
    return(list(type = 1L, element_size = 1L, count = length(values),
                data = as.raw(pmin(255L, pmax(0L, values)))))
  }
  if (identical(kind, "short")) {
    values <- as.integer(values)
    return(list(type = 4L, element_size = 2L, count = length(values),
                data = write_int_be(values, size = 2L)))
  }
  stop("Unsupported ABIF value kind: ", kind, call. = FALSE)
}

pack_abif_directory_entry <- function(tag, number, type, element_size, count,
                                      data_size, offset_or_data, handle) {
  c(
    charToRaw(tag),
    write_uint32_be(number),
    write_uint16_be(type),
    write_uint16_be(element_size),
    write_uint32_be(count),
    write_uint32_be(data_size),
    if (is.raw(offset_or_data)) offset_or_data else write_uint32_be(offset_or_data),
    write_uint32_be(handle)
  )
}

write_int_be <- function(values, size) {
  con <- rawConnection(raw(), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.integer(values), con, size = size, endian = "big")
  rawConnectionValue(con)
}

write_uint16_be <- function(value) {
  value <- as.integer(value)
  as.raw(c(bitwAnd(bitwShiftR(value, 8L), 255L), bitwAnd(value, 255L)))
}

write_uint32_be <- function(value) {
  value <- as.numeric(value)
  as.raw(c(
    floor(value / 16777216) %% 256,
    floor(value / 65536) %% 256,
    floor(value / 256) %% 256,
    value %% 256
  ))
}
