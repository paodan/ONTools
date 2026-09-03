#' Create a synthetic AB1 trace from a consensus FASTA and BAM
#'
#' `synthetic_ab1_from_bam()` writes a minimal ABIF/AB1 file with simulated
#' Sanger-style four-channel peak data. It uses `samtools mpileup` to estimate
#' A/C/G/T support at each consensus position, then turns those proportions into
#' Gaussian peaks.
#'
#' This function does not reconstruct real capillary-sequencing fluorescence
#' signal from the BAM. Use the result for visualization or compatibility with
#' tools that require an AB1-like chromatogram, not as raw Sanger evidence.
#'
#' @param consensus Consensus FASTA used as the BAM reference.
#' @param bam Coordinate-sorted BAM aligned to `consensus`.
#' @param output_ab1 Output `.ab1` file.
#' @param reference_name FASTA/BAM reference name. Required when `consensus`
#'   contains more than one sequence.
#' @param sample Sample name stored in the AB1 `SMPL1` tag. Defaults to the
#'   output filename without extension.
#' @param samtools Command name or full path used to launch `samtools`.
#' @param conda_env Optional conda environment name. If supplied, `samtools` is
#'   run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param min_base_quality Minimum base quality passed to `samtools mpileup -Q`.
#' @param spacing Number of trace points between adjacent called bases.
#' @param sigma Gaussian peak width in trace points.
#' @param peak_height Main peak height.
#' @param baseline Baseline signal added to every channel.
#' @param noise_fraction Minimum off-channel peak fraction.
#' @param overwrite Logical. If `FALSE`, stop when `output_ab1` already exists.
#' @param echo Logical. If `TRUE`, print the `samtools mpileup` command.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `paths`, `reference_name`,
#'   `sequence_length`, `trace_length`, and `command`.
#'
#' @examples
#' consensus <- tempfile(fileext = ".fasta")
#' bam <- tempfile(fileext = ".bam")
#' ab1 <- tempfile(fileext = ".ab1")
#' writeLines(c(">contig1", "ACGTACGT"), consensus)
#' # A real coordinate-sorted BAM and its index are needed for execution.
#' # synthetic_ab1_from_bam(consensus, bam, ab1)
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
    noise_fraction = noise_fraction
  )
  if (is.null(sample)) {
    sample <- sub("[.][^.]+$", "", basename(output_ab1))
  }
  write_abif_file(
    output_ab1,
    sequence = sequence,
    channels = trace$channels,
    peak_locations = trace$peak_locations,
    qualities = trace$qualities,
    sample = sample
  )

  invisible(list(
    status = 0L,
    paths = list(consensus = consensus, bam = bam, ab1 = output_ab1),
    reference_name = reference_name,
    sequence_length = nchar(sequence),
    trace_length = length(trace$channels$A),
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
                                noise_fraction) {
  bases <- c("A", "C", "G", "T")
  n_bases <- nchar(sequence)
  seq_chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
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
    fractions <- ab1_base_fractions(
      pileup_counts[[as.character(idx)]],
      seq_chars[[idx]],
      noise_fraction
    )
    center <- peak_locations[[idx]]
    for (base in bases) {
      height <- peak_height * fractions[[base]]
      positions <- seq.int(max(1L, center - radius + 1L), min(trace_length, center + radius + 1L))
      distances <- (positions - 1L) - center
      channels[[base]][positions] <- channels[[base]][positions] +
        round(height * exp(-(distances * distances) / (2 * sigma * sigma)))
    }
  }

  channels <- lapply(channels, function(x) pmin(32767L, pmax(0L, as.integer(x))))
  list(
    channels = channels,
    peak_locations = as.integer(peak_locations),
    qualities = consensus_ab1_qualities(seq_chars, pileup_counts)
  )
}

ab1_base_fractions <- function(counts, consensus_base, noise_fraction) {
  bases <- c("A", "C", "G", "T")
  values <- stats::setNames(rep(0, 4L), bases)
  if (!is.null(counts)) {
    values[names(counts)] <- as.numeric(counts)
  }
  total <- sum(values)
  if (total > 0) values <- values / total
  consensus_base <- toupper(consensus_base)
  if (consensus_base %in% bases) {
    values[[consensus_base]] <- max(values[[consensus_base]], 1)
  } else {
    values[] <- 0.25
  }
  values <- pmax(values, noise_fraction)
  values / max(values)
}

consensus_ab1_qualities <- function(seq_chars, pileup_counts) {
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
