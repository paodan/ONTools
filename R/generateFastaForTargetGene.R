#' Generate target-gene and replacement FASTA fragments
#'
#' `generateFastaForTargetGene()` maps one gene sequence to a genome and writes
#' FASTA files for the gene, homology arms, flanking regions, and replacement
#' constructs.
#'
#' @param output_path Output directory. One subdirectory is created for each
#'   accepted mapping hit.
#' @param genome_fasta Genome/reference FASTA file.
#' @param gene_seq_fasta FASTA file containing exactly one target gene sequence.
#' @param replace_seq_fastq FASTA file containing exactly one replacement
#'   sequence. The argument name is kept for backward compatibility.
#' @param h1_len,h2_len Homology-arm lengths upstream and downstream of the
#'   target gene.
#' @param h1up_len,h2down_len Additional flank lengths outside H1 and H2.
#' @param preset Minimap2 preset used when aligning the gene to the genome.
#' @param min_mapq Optional minimum mapping quality for accepted PAF hits. If
#'   `NULL`, mapping quality is not used for filtering.
#' @param max_unaligned Optional maximum allowed number of unaligned query bases
#'   across both ends. If `NULL`, query end coverage is not used for filtering.
#' @param paf_file Optional existing PAF file. If supplied, minimap2 is not run.
#' @param minimap2 Command name or executable path.
#' @param conda_env Optional conda environment name. If supplied, minimap2 is run
#'   with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param output_ext Extension appended to generated FASTA file names.
#' @param combined_fasta_name File name for the combined FASTA written in each
#'   hit directory. Set to `NULL` to skip writing the combined FASTA.
#' @param trim_to_bounds Logical. If `FALSE`, stop when requested flank regions
#'   extend outside the contig. If `TRUE`, trim regions to contig bounds.
#' @param dry_run Logical. If `TRUE`, return planned commands without running
#'   minimap2 or writing FASTA files.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stderr Passed to [system2()].
#'
#' @return Invisibly returns a list with `paf`, `seq`, `output_path`,
#'   `written_files`, `commands`, and `paths`.
#'
#' @examples
#' \dontrun{
#' res <- generateFastaForTargetGene(
#'   output_path = "./fasta/v2",
#'   genome_fasta = "genome/GCA_019096115.1_ASM1909611v1_genomic.fna",
#'   gene_seq_fasta = "fasta/v1/3_GME292_g.fasta",
#'   replace_seq_fastq = "fasta/v1/3_NTC.fa",
#'   h1_len = 1000,
#'   h2_len = 1000,
#'   h1up_len = 1000,
#'   h2down_len = 1000
#' )
#' }
#'
#' @export
generateFastaForTargetGene <- function(output_path,
                                       genome_fasta = "genome/GCA_019096115.1_ASM1909611v1_genomic.fna",
                                       gene_seq_fasta = "fasta/v1/3_GME292_g.fasta",
                                       replace_seq_fastq = "fasta/v1/3_NTC.fa",
                                       h1_len = 1000,
                                       h2_len = 1000,
                                       h1up_len = 1000,
                                       h2down_len = 1000,
                                       preset = "asm5",
                                       min_mapq = NULL,
                                       max_unaligned = NULL,
                                       paf_file = NULL,
                                       minimap2 = "minimap2",
                                       conda_env = NULL,
                                       conda = "conda",
                                       output_ext = ".fasta",
                                       combined_fasta_name = "all_sequences.fasta",
                                       trim_to_bounds = FALSE,
                                       dry_run = FALSE,
                                       echo = TRUE,
                                       stderr = "") {
  if (missing(output_path)) stop("`output_path` is required.", call. = FALSE)
  check_scalar_character(output_path, "output_path")
  check_file_arg(genome_fasta, "genome_fasta")
  check_file_arg(gene_seq_fasta, "gene_seq_fasta")
  check_file_arg(replace_seq_fastq, "replace_seq_fastq")
  check_scalar_character(preset, "preset")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(conda, "conda")
  check_scalar_character(output_ext, "output_ext")
  check_logical_scalar(trim_to_bounds, "trim_to_bounds")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  if (!is.null(paf_file)) check_file_arg(paf_file, "paf_file")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  if (!is.null(combined_fasta_name)) check_scalar_character(combined_fasta_name, "combined_fasta_name")

  h1_len <- validate_nonnegative_integer(h1_len, "h1_len")
  h2_len <- validate_nonnegative_integer(h2_len, "h2_len")
  h1up_len <- validate_nonnegative_integer(h1up_len, "h1up_len")
  h2down_len <- validate_nonnegative_integer(h2down_len, "h2down_len")
  if (!is.null(min_mapq)) min_mapq <- validate_nonnegative_number(min_mapq, "min_mapq")
  if (!is.null(max_unaligned)) max_unaligned <- validate_nonnegative_integer(max_unaligned, "max_unaligned")

  genome_fasta <- normalizePath(genome_fasta, mustWork = TRUE)
  gene_seq_fasta <- normalizePath(gene_seq_fasta, mustWork = TRUE)
  replace_seq_fastq <- normalizePath(replace_seq_fastq, mustWork = TRUE)

  genome <- Biostrings::readDNAStringSet(genome_fasta)
  names(genome) <- sub(" .*", "", names(genome))

  gene_seq <- Biostrings::readDNAStringSet(gene_seq_fasta)
  gene_name <- sub(" .*", "", names(gene_seq))
  if (length(gene_seq) != 1L) {
    stop("`gene_seq_fasta` must contain exactly one sequence.", call. = FALSE)
  }

  replace_seq <- Biostrings::readDNAStringSet(replace_seq_fastq)
  names(replace_seq) <- sub(" .*", "", names(replace_seq))
  if (length(replace_seq) != 1L) {
    stop("`replace_seq_fastq` must contain exactly one sequence.", call. = FALSE)
  }

  temp_paf <- is.null(paf_file)
  if (temp_paf) {
    paf_file <- tempfile(fileext = ".paf")
  } else {
    paf_file <- normalizePath(paf_file, mustWork = TRUE)
  }

  minimap2_args <- c("-x", preset, genome_fasta, gene_seq_fasta)
  minimap2_call <- dehost_fastq_external_call(minimap2, minimap2_args, conda_env, conda)
  commands <- list(
    minimap2 = paste(c(shQuote(minimap2_call$command), shQuote(minimap2_call$args), ">", shQuote(paf_file)), collapse = " ")
  )

  if (temp_paf && isTRUE(echo)) {
    message(commands$minimap2)
  }

  paths <- list(
    output_path = output_path,
    genome_fasta = genome_fasta,
    gene_seq_fasta = gene_seq_fasta,
    replace_seq_fasta = replace_seq_fastq,
    paf = paf_file
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      paf = NULL,
      seq = list(),
      output_path = output_path,
      written_files = character(),
      commands = commands,
      paths = paths,
      conda_env = conda_env
    )))
  }

  if (temp_paf) {
    if (is.null(conda_env)) {
      require_external_command(minimap2)
    } else {
      require_external_command(conda)
    }
    on.exit(unlink(paf_file), add = TRUE)

    paf_status <- system2(
      minimap2_call$command,
      args = minimap2_call$args,
      stdout = paf_file,
      stderr = stderr
    )
    if (!identical(paf_status, 0L)) {
      stop("minimap2 failed with exit status: ", paf_status, call. = FALSE)
    }
  }

  paf <- read_paf(paf_file)
  if (is.null(paf) || nrow(paf) == 0L) {
    stop("No target-gene alignment was found in the genome.", call. = FALSE)
  }

  paf <- generate_target_gene_filter_paf(
    paf = paf,
    min_mapq = min_mapq,
    max_unaligned = max_unaligned
  )
  if (nrow(paf) == 0L) {
    stop("No target-gene alignment passed the filtering criteria.", call. = FALSE)
  }

  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  output_path <- normalizePath(output_path, mustWork = TRUE)

  seq <- list()
  written_files <- character()
  contig_lengths <- stats::setNames(BiocGenerics::width(genome), names(genome))

  for (i in seq_len(nrow(paf))) {
    hit <- generate_target_gene_hit_interval(paf[i, , drop = FALSE])
    contig <- hit$target_name
    if (!contig %in% names(genome)) {
      stop("Aligned contig not found in `genome_fasta`: ", contig, call. = FALSE)
    }

    regions <- generate_target_gene_regions(
      contig = contig,
      gene_name = gene_name,
      gene_start = hit$genome_start,
      gene_end = hit$genome_end,
      h1_len = h1_len,
      h2_len = h2_len,
      h1up_len = h1up_len,
      h2down_len = h2down_len
    )
    regions <- generate_target_gene_check_regions(
      regions = regions,
      contig_length = contig_lengths[[contig]],
      trim_to_bounds = trim_to_bounds
    )

    seq_i <- Biostrings::subseq(
      genome[regions$target_name],
      start = regions$genome_start,
      end = regions$genome_end
    )
    names(seq_i) <- regions$region

    h1_replace_h2 <- Biostrings::DNAStringSet(c(seq_i[["h1"]], replace_seq[[1]], seq_i[["h2"]]))
    names(h1_replace_h2) <- paste(
      names(replace_seq),
      "h1_replace_h2",
      contig,
      hit$genome_start,
      hit$genome_end,
      sep = "_"
    )

    h1up_h1_replace_h2_h2down <- Biostrings::DNAStringSet(c(
      seq_i[["h1up"]],
      seq_i[["h1"]],
      replace_seq[[1]],
      seq_i[["h2"]],
      seq_i[["h2down"]]
    ))
    names(h1up_h1_replace_h2_h2down) <- paste(
      names(replace_seq),
      "h1up_h1_replace_h2_h2down",
      contig,
      hit$genome_start,
      hit$genome_end,
      sep = "_"
    )

    names(seq_i) <- regions$full_name
    seq[[i]] <- c(seq_i, replace_seq, h1_replace_h2, h1up_h1_replace_h2_h2down)

    hit_dir <- file.path(output_path, as.character(i))
    dir.create(hit_dir, showWarnings = FALSE, recursive = TRUE)
    if (!is.null(combined_fasta_name)) {
      combined_file <- file.path(hit_dir, combined_fasta_name)
      Biostrings::writeXStringSet(seq[[i]], filepath = combined_file, format = "fasta")
      written_files <- c(written_files, combined_file)
    }
    for (sequence_name in names(seq[[i]])) {
      file_name <- paste0(generate_target_gene_safe_name(sequence_name), output_ext)
      out_file <- file.path(hit_dir, file_name)
      Biostrings::writeXStringSet(seq[[i]][sequence_name], filepath = out_file, format = "fasta")
      written_files <- c(written_files, out_file)
    }
  }

  invisible(list(
    paf = paf,
    seq = seq,
    output_path = output_path,
    written_files = written_files,
    commands = commands,
    paths = paths,
    conda_env = conda_env
  ))
}

generate_target_gene_filter_paf <- function(paf, min_mapq, max_unaligned) {
  if (is.null(min_mapq) && is.null(max_unaligned)) {
    return(paf)
  }

  keep <- rep(TRUE, nrow(paf))
  if (!is.null(min_mapq)) {
    keep <- keep &
      !is.na(paf$mapping_quality) &
      paf$mapping_quality >= min_mapq
  }

  unaligned <- paf$query_start + (paf$query_length - paf$query_end)
  if (!is.null(max_unaligned)) {
    keep <- keep &
      !is.na(unaligned) &
      unaligned <= max_unaligned
  }

  paf[keep, , drop = FALSE]
}

generate_target_gene_hit_interval <- function(paf_row) {
  if (identical(paf_row$strand, "-")) {
    genome_start <- paf_row$target_start - (paf_row$query_length - paf_row$query_end) + 1
    genome_end <- paf_row$target_end + paf_row$query_start
  } else {
    genome_start <- paf_row$target_start - paf_row$query_start + 1
    genome_end <- paf_row$target_end + (paf_row$query_length - paf_row$query_end)
  }

  list(
    target_name = paf_row$target_name,
    genome_start = as.integer(genome_start),
    genome_end = as.integer(genome_end)
  )
}

generate_target_gene_regions <- function(contig,
                                         gene_name,
                                         gene_start,
                                         gene_end,
                                         h1_len,
                                         h2_len,
                                         h1up_len,
                                         h2down_len) {
  regions <- data.frame(
    region = c(
      "h1up",
      "h1",
      "genes",
      "h2",
      "h2down",
      "h1_gene_h2",
      "h1up_h1_gene_h2_h2down"
    ),
    target_name = contig,
    genome_start = c(
      gene_start - h1_len - h1up_len,
      gene_start - h1_len,
      gene_start,
      gene_end + 1L,
      gene_end + h2_len + 1L,
      gene_start - h1_len,
      gene_start - h1_len - h1up_len
    ),
    genome_end = c(
      gene_start - h1_len - 1L,
      gene_start - 1L,
      gene_end,
      gene_end + h2_len,
      gene_end + h2_len + h2down_len,
      gene_end + h2_len,
      gene_end + h2_len + h2down_len
    ),
    name = c(
      paste0(gene_name, "_h1up"),
      paste0(gene_name, "_h1"),
      gene_name,
      paste0(gene_name, "_h2"),
      paste0(gene_name, "_h2down"),
      paste0(gene_name, "_h1_gene_h2"),
      paste0(gene_name, "_h1up_h1_gene_h2_h2down")
    ),
    stringsAsFactors = FALSE
  )
  regions$full_name <- paste0(
    regions$name,
    "_",
    regions$target_name,
    "_",
    regions$genome_start,
    "_",
    regions$genome_end
  )
  regions
}

generate_target_gene_check_regions <- function(regions, contig_length, trim_to_bounds) {
  out_of_bounds <- regions$genome_start < 1L | regions$genome_end > contig_length
  if (any(out_of_bounds) && !isTRUE(trim_to_bounds)) {
    bad <- regions[out_of_bounds, , drop = FALSE]
    stop(
      "Requested region extends outside contig bounds: ",
      paste(bad$full_name, collapse = ", "),
      call. = FALSE
    )
  }

  regions$genome_start <- pmax(1L, regions$genome_start)
  regions$genome_end <- pmin(contig_length, regions$genome_end)
  invalid <- regions$genome_start > regions$genome_end
  if (any(invalid)) {
    stop(
      "Requested region has non-positive length after boundary handling: ",
      paste(regions$full_name[invalid], collapse = ", "),
      call. = FALSE
    )
  }

  regions$full_name <- paste0(
    regions$name,
    "_",
    regions$target_name,
    "_",
    regions$genome_start,
    "_",
    regions$genome_end
  )
  row.names(regions) <- regions$region
  regions
}

generate_target_gene_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("^_+", "", sub("_+$", "", x))
}

validate_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x != as.integer(x)) {
    stop("`", name, "` must be a single non-negative integer.", call. = FALSE)
  }

  as.integer(x)
}
