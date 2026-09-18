#' Annotate consensus sequences with BLASTN
#'
#' `annotate_consensus_blast()` runs `blastn` against a local nucleotide BLAST
#' database, parses the tabular result, calculates query/reference coverage,
#' extracts common UNITE-style taxonomy fields when present in subject titles,
#' and assigns a conservative annotation level for each hit.
#'
#' This helper is intended for independent review of fungal/eukaryotic ITS
#' consensus sequences, especially by searching consensus or representative
#' sequences against a local UNITE Species Hypothesis (SH) FASTA database.
#'
#' @param consensus_fasta Query consensus/representative sequence FASTA file.
#'   This may contain one or multiple consensus sequences. Each FASTA record is
#'   reported separately in the BLAST table through `qseqid`.
#' @param db BLAST database prefix passed to `blastn -db`. If `NULL`,
#'   `db_fasta` must be supplied and a database is built under `out_dir`.
#'   Use this when the UNITE database has already been prepared with
#'   `makeblastdb`.
#' @param db_fasta Optional FASTA file used to build a BLAST nucleotide database
#'   with `makeblastdb`. This is useful when you have downloaded a UNITE FASTA
#'   file but have not yet built local BLAST index files.
#' @param out_dir Output directory for BLAST results and optional database
#'   files. Default is `"work/consensus_annotation"`.
#' @param output_tsv Output BLAST TSV path. If `NULL`, writes
#'   `<out_dir>/<prefix>.blast.tsv`.
#' @param prefix Prefix used when constructing default output paths. If `NULL`,
#'   it is inferred from `consensus_fasta`.
#' @param threads Positive integer thread count passed to `blastn`. Default is
#'   `10`.
#' @param max_target_seqs Maximum number of target sequences reported per query.
#'   Default is `20`.
#' @param evalue E-value cutoff passed to `blastn`. Default is `"1e-20"`.
#' @param task Optional BLASTN task, for example `"blastn"`, `"megablast"`, or
#'   `"dc-megablast"`. Default `NULL` lets BLAST use its own default.
#' @param word_size Optional word size passed to `blastn -word_size`. Default is
#'   `NULL`.
#' @param strand Optional strand passed to `blastn -strand`, for example
#'   `"both"`, `"plus"`, or `"minus"`. Default is `NULL`.
#' @param dust Optional dust setting passed to `blastn -dust`, for example
#'   `"yes"` or `"no"`. Default is `NULL`.
#' @param perc_identity Optional minimum percent identity passed to
#'   `blastn -perc_identity`. Default is `NULL`.
#' @param extra_args Optional additional `blastn` arguments as a character
#'   vector, for example `c("-ungapped")`. Default is `NULL`.
#' @param blastn,makeblastdb Command names or executable paths. Defaults are
#'   `"blastn"` and `"makeblastdb"`.
#' @param conda_env Optional conda environment name. If supplied, external
#'   commands are run with `conda run -n <conda_env>`. Default `NULL` calls
#'   BLAST tools from the current environment.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#'   Default is `"conda"`.
#' @param taxdump_dir Optional NCBI taxonomy dump directory containing
#'   `nodes.dmp` and `names.dmp`. When supplied, missing taxonomy fields are
#'   filled from BLAST `staxids`, which is useful for NCBI 16S ribosomal RNA
#'   BLAST databases whose FASTA titles do not contain full lineage strings.
#' @param dry_run Logical. If `TRUE`, return planned commands without running.
#'   Default is `FALSE`.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#'   Default is `TRUE`.
#' @param stdout,stderr Passed to [system2()]. Defaults stream command output
#'   and errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `commands`, `paths`, `blast`,
#'   `top_hits`, and `conda_env`.
#'
#' @details
#' `annotate_consensus_blast()` runs BLAST in tabular format, then rewrites the
#' output TSV with column headers. The BLAST-derived fields are: `qseqid`,
#' `qlen`, `qstart`, `qend`, `sseqid`, `slen`, `sstart`, `send`, `length`,
#' `pident`, `qcovs`, `mismatch`, `gapopen`, `evalue`, `bitscore`,
#' `salltitles`, and `staxids`.
#'
#' After reading the BLAST output, the function adds:
#'
#' * `query_coverage`: percent of the query sequence spanned by the alignment,
#'   calculated from `qstart`, `qend`, and `qlen`.
#' * `reference_coverage`: percent of the subject/reference sequence spanned by
#'   the alignment, calculated from `sstart`, `send`, and `slen`.
#' * `taxonomy_path`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`,
#'   and `species`: parsed from UNITE/SILVA-style `k__...;p__...;...;s__...`
#'   taxonomy strings in `salltitles` when present, and optionally filled from
#'   NCBI taxdump lineage using `staxids` when `taxdump_dir` is supplied.
#' * `rank`: hit rank within each query, ordered by higher `bitscore`, lower
#'   `evalue`, higher `pident`, then higher `query_coverage`.
#'
#' The function intentionally reports the BLAST metrics and parsed taxonomy
#' fields without assigning species/genus/family calls or novel-candidate flags.
#' Users should interpret the result from identity, query/reference coverage,
#' hit separation, and database metadata.
#'
#' For UNITE databases, ONTools expects FASTA headers that contain an SH
#' identifier and taxonomy string. A typical header looks like:
#'
#' `>Claroideoglomus_sp|AM076567|SH1229972.10FU|reps|k__Fungi;p__Glomeromycota;c__Glomeromycetes;o__Entrophosporales;f__Entrophosporaceae;g__Claroideoglomus;s__Claroideoglomus_sp`
#'
#' `reps` indicates representative sequences chosen automatically, whereas
#' `refs` indicates reference sequences chosen, overridden, or confirmed by
#' users with taxonomic expertise. UNITE SH FASTA files are appropriate for
#' local BLAST review of ITS consensus sequences.
#'
#' Example local UNITE database record used by the authors: downloaded from
#' <https://unite.ut.ee/repository.php> on 2026-09-12. The all-eukaryotes file
#' was `UNITE_eukaryotes_all.fasta`, corresponding to
#' `sh_general_release_s_all_19.02.2025`, version 10.0, release date
#' 2025-02-19, taxon group all eukaryotes, 20,802 RefS, 245,787 RepS, DOI
#' <https://doi.org/10.15156/BIO/3301232>. The fungi-only file was
#' `UNITE_fungi.fasta`, corresponding to `sh_general_release_s_19.02.2025`,
#' version 10.0, release date 2025-02-19, taxon group Fungi, 20,295 RefS,
#' 147,735 RepS, DOI <https://doi.org/10.15156/BIO/3301230>.
#'
#' To build a BLAST database manually, run for example:
#'
#' `makeblastdb -in UNITE_eukaryotes_all.fasta -dbtype nucl -out unite_eukaryotes`
#'
#' Then call this function with `db = "unite_eukaryotes"`. Alternatively, pass
#' the FASTA directly with `db_fasta = "UNITE_eukaryotes_all.fasta"` and the
#' function will run `makeblastdb` under `out_dir`.
#'
#' @examples
#' query <- tempfile(fileext = ".fasta")
#' writeLines(c(">consensus1", "ACGTACGT"), query)
#' res <- annotate_consensus_blast(
#'   consensus_fasta = query,
#'   db = "unite_eukaryotes",
#'   dry_run = TRUE
#' )
#' res$commands$blastn
#'
#' # Build a BLAST database from a downloaded UNITE FASTA, then annotate.
#' # res <- annotate_consensus_blast(
#' #   consensus_fasta = "all-consensus-seqs_trimmed.fasta",
#' #   db_fasta = "UNITE_eukaryotes_all.fasta",
#' #   out_dir = "work/unite_consensus_annotation"
#' # )
#'
#' # Use an already-built local UNITE database.
#' # res <- annotate_consensus_blast(
#' #   consensus_fasta = "all-consensus-seqs_trimmed.fasta",
#' #   db = "/data/reference/UNITE/unite_eukaryotes",
#' #   output_tsv = "unite_consensus_annotation/consensus.blast.tsv",
#' #   threads = 20
#' # )
#'
#' @export
annotate_consensus_blast <- function(consensus_fasta,
                                     db = NULL,
                                     db_fasta = NULL,
                                     out_dir = "work/consensus_annotation",
                                     output_tsv = NULL,
                                     prefix = NULL,
                                     threads = 10,
                                     max_target_seqs = 20,
                                     evalue = "1e-20",
                                     task = NULL,
                                     word_size = NULL,
                                     strand = NULL,
                                     dust = NULL,
                                     perc_identity = NULL,
                                     extra_args = NULL,
                                     blastn = "blastn",
                                     makeblastdb = "makeblastdb",
                                     conda_env = NULL,
                                     conda = "conda",
                                     taxdump_dir = NULL,
                                     dry_run = FALSE,
                                     echo = TRUE,
                                     stdout = "",
                                     stderr = "") {
  check_scalar_character(consensus_fasta, "consensus_fasta")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(evalue, "evalue")
  check_scalar_character(blastn, "blastn")
  check_scalar_character(makeblastdb, "makeblastdb")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  if (!is.null(db)) check_scalar_character(db, "db")
  if (!is.null(db_fasta)) check_file_arg(db_fasta, "db_fasta")
  if (is.null(db) && is.null(db_fasta)) {
    stop("Supply either `db` or `db_fasta`.", call. = FALSE)
  }
  if (!is.null(output_tsv)) check_scalar_character(output_tsv, "output_tsv")
  if (!is.null(prefix)) check_scalar_character(prefix, "prefix")
  if (!is.null(task)) check_scalar_character(task, "task")
  if (!is.null(strand)) check_scalar_character(strand, "strand")
  if (!is.null(dust)) check_scalar_character(dust, "dust")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  if (!is.null(taxdump_dir)) check_dir_arg(taxdump_dir, "taxdump_dir")
  if (!is.null(extra_args) &&
      (!is.character(extra_args) || anyNA(extra_args))) {
    stop("`extra_args` must be a character vector without missing values.",
         call. = FALSE)
  }

  threads <- validate_positive_integer(threads, "threads")
  max_target_seqs <- validate_positive_integer(max_target_seqs, "max_target_seqs")
  if (!is.null(word_size)) word_size <- validate_positive_integer(word_size, "word_size")
  if (!is.null(perc_identity)) {
    perc_identity <- validate_nonnegative_number(perc_identity, "perc_identity")
  }
  if (file.exists(consensus_fasta)) {
    consensus_fasta <- normalizePath(consensus_fasta, mustWork = TRUE)
  } else if (!isTRUE(dry_run)) {
    stop("`consensus_fasta` does not exist: ", consensus_fasta, call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  if (!is.null(db_fasta)) db_fasta <- normalizePath(db_fasta, mustWork = TRUE)
  if (is.null(prefix)) prefix <- consensus_blast_prefix(consensus_fasta)
  if (is.null(db)) {
    db <- file.path(out_dir, consensus_blast_prefix(db_fasta))
  }
  if (is.null(output_tsv)) {
    output_tsv <- file.path(out_dir, paste0(prefix, ".blast.tsv"))
  } else {
    dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
  }
  output_tsv <- normalizePath(output_tsv, mustWork = FALSE)

  outfmt_fields <- c(
    "qseqid", "qlen", "qstart", "qend", "sseqid", "slen", "sstart",
    "send", "length", "pident", "qcovs", "mismatch", "gapopen",
    "evalue", "bitscore", "salltitles", "staxids"
  )
  blast_args <- c(
    "-query", consensus_fasta,
    "-db", db,
    "-num_threads", as.character(threads),
    "-max_target_seqs", as.character(max_target_seqs),
    "-evalue", evalue,
    "-outfmt", paste("6", paste(outfmt_fields, collapse = " ")),
    "-out", output_tsv
  )
  if (!is.null(task)) blast_args <- c(blast_args, "-task", task)
  if (!is.null(word_size)) blast_args <- c(blast_args, "-word_size", as.character(word_size))
  if (!is.null(strand)) blast_args <- c(blast_args, "-strand", strand)
  if (!is.null(dust)) blast_args <- c(blast_args, "-dust", dust)
  if (!is.null(perc_identity)) {
    blast_args <- c(blast_args, "-perc_identity", as.character(perc_identity))
  }
  if (!is.null(extra_args)) blast_args <- c(blast_args, extra_args)

  makeblastdb_args <- NULL
  if (!is.null(db_fasta)) {
    makeblastdb_args <- c("-in", db_fasta, "-dbtype", "nucl", "-out", db)
  }

  blast_call <- dehost_fastq_external_call(blastn, blast_args, conda_env, conda)
  makeblastdb_call <- NULL
  if (!is.null(makeblastdb_args)) {
    makeblastdb_call <- dehost_fastq_external_call(
      makeblastdb,
      makeblastdb_args,
      conda_env,
      conda
    )
  }

  commands <- list()
  if (!is.null(makeblastdb_call)) {
    commands$makeblastdb <- paste(
      c(shQuote(makeblastdb_call$command), shQuote(makeblastdb_call$args)),
      collapse = " "
    )
  }
  commands$blastn <- paste(
    c(shQuote(blast_call$command), shQuote(blast_call$args)),
    collapse = " "
  )

  if (isTRUE(echo)) {
    if (!is.null(commands$makeblastdb)) message(commands$makeblastdb)
    message(commands$blastn)
  }

  paths <- list(
    consensus_fasta = consensus_fasta,
    db = db,
    db_fasta = db_fasta,
    output_tsv = output_tsv,
    out_dir = out_dir
  )
  if (!is.null(taxdump_dir)) {
    paths$taxdump_dir <- normalizePath(taxdump_dir, mustWork = TRUE)
  }
  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      commands = commands,
      paths = paths,
      blast = NULL,
      top_hits = NULL,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(blastn)
    if (!is.null(makeblastdb_call)) require_external_command(makeblastdb)
  } else {
    require_external_command(conda)
  }

  if (!is.null(makeblastdb_call)) {
    make_status <- system2(
      makeblastdb_call$command,
      args = quote_system2_args(makeblastdb_call$args),
      stdout = stdout,
      stderr = stderr
    )
    if (!identical(make_status, 0L)) {
      stop("makeblastdb failed with exit status: ", make_status, call. = FALSE)
    }
  }

  blast_status <- system2(
    blast_call$command,
    args = quote_system2_args(blast_call$args),
    stdout = stdout,
    stderr = stderr
  )
  if (!identical(blast_status, 0L)) {
    stop("blastn failed with exit status: ", blast_status, call. = FALSE)
  }

  blast <- read_consensus_blast_table(
    output_tsv,
    outfmt_fields,
    taxdump_dir = taxdump_dir
  )
  top_hits <- blast[blast$rank == 1L, , drop = FALSE]
  utils::write.table(
    blast,
    output_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(list(
    status = 0L,
    commands = commands,
    paths = paths,
    blast = blast,
    top_hits = top_hits,
    conda_env = conda_env
  ))
}

quote_system2_args <- function(args) {
  if (length(args) == 0L) return(args)
  needs_quote <- grepl("[[:space:]'\"\\\\$`!#&;<>*?(){}\\[\\]|]", args)
  args[needs_quote] <- shQuote(args[needs_quote])
  args
}

read_consensus_blast_table <- function(path, fields, taxdump_dir = NULL) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    empty <- data.frame(matrix(ncol = length(fields), nrow = 0L))
    names(empty) <- fields
    return(empty)
  }

  blast <- utils::read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  names(blast) <- fields[seq_len(ncol(blast))]

  numeric_cols <- intersect(
    c(
      "qlen", "qstart", "qend", "slen", "sstart", "send", "length",
      "pident", "qcovs", "mismatch", "gapopen", "evalue", "bitscore"
    ),
    names(blast)
  )
  for (col in numeric_cols) {
    blast[[col]] <- suppressWarnings(as.numeric(blast[[col]]))
  }

  blast$query_coverage <- ifelse(
    !is.na(blast$qlen) & blast$qlen > 0,
    abs(blast$qend - blast$qstart) + 1L,
    NA_real_
  )
  blast$query_coverage <- blast$query_coverage / blast$qlen * 100
  blast$reference_coverage <- ifelse(
    !is.na(blast$slen) & blast$slen > 0,
    (abs(blast$send - blast$sstart) + 1L) / blast$slen * 100,
    NA_real_
  )

  taxonomy <- parse_consensus_blast_taxonomy(blast$salltitles)
  if (!is.null(taxdump_dir)) {
    taxonomy <- fill_consensus_taxonomy_from_taxdump(
      taxonomy = taxonomy,
      staxids = blast$staxids,
      taxdump_dir = taxdump_dir
    )
  }
  blast <- cbind(blast, taxonomy)
  blast <- rank_consensus_blast_hits(blast)
  blast
}

rank_consensus_blast_hits <- function(blast) {
  if (nrow(blast) == 0L) {
    blast$rank <- integer()
    return(blast)
  }

  order_index <- order(
    blast$qseqid,
    -blast$bitscore,
    blast$evalue,
    -blast$pident,
    -blast$query_coverage,
    na.last = TRUE
  )
  blast <- blast[order_index, , drop = FALSE]
  blast$rank <- ave(
    seq_len(nrow(blast)),
    blast$qseqid,
    FUN = seq_along
  )
  rownames(blast) <- NULL
  blast
}

parse_consensus_blast_taxonomy <- function(titles) {
  rows <- lapply(titles, parse_one_consensus_blast_title)
  do.call(rbind, rows)
}

parse_one_consensus_blast_title <- function(title) {
  out <- data.frame(
    taxonomy_path = NA_character_,
    kingdom = NA_character_,
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    species = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.na(title) || !nzchar(title)) return(out)

  match <- regexpr(
    "([dkpcofgs]__[^;|[:space:]]+[;|]?){2,}",
    title,
    perl = TRUE
  )
  if (match[1L] < 0L) return(out)

  tax <- regmatches(title, match)
  tax <- sub("[|]$", "", tax)
  parts <- strsplit(tax, ";", fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  values <- stats::setNames(rep(NA_character_, 8L), c(
    "d", "k", "p", "c", "o", "f", "g", "s"
  ))
  for (part in parts) {
    key <- sub("^([dkpcofgs])__.*$", "\\1", part)
    value <- sub("^[dkpcofgs]__", "", part)
    value <- gsub("_", " ", value, fixed = TRUE)
    if (key %in% names(values) && nzchar(value)) values[[key]] <- value
  }

  out$taxonomy_path <- paste(stats::na.omit(unname(values[c("k", "p", "c", "o", "f", "g", "s")])), collapse = ";")
  if (!nzchar(out$taxonomy_path)) out$taxonomy_path <- NA_character_
  out$kingdom <- values[["k"]]
  out$phylum <- values[["p"]]
  out$class <- values[["c"]]
  out$order <- values[["o"]]
  out$family <- values[["f"]]
  out$genus <- values[["g"]]
  out$species <- values[["s"]]
  out
}

fill_consensus_taxonomy_from_taxdump <- function(taxonomy, staxids, taxdump_dir) {
  missing <- is.na(taxonomy$taxonomy_path) | !nzchar(taxonomy$taxonomy_path)
  if (!any(missing)) return(taxonomy)

  taxdump <- read_ncbi_taxdump(taxdump_dir)
  for (i in which(missing)) {
    taxid <- first_blast_taxid(staxids[[i]])
    if (is.na(taxid)) next
    lineage <- ncbi_taxid_lineage(taxid, taxdump)
    if (is.null(lineage)) next
    for (col in intersect(names(lineage), names(taxonomy))) {
      taxonomy[[col]][[i]] <- lineage[[col]]
    }
  }
  taxonomy
}

first_blast_taxid <- function(staxid) {
  if (is.na(staxid) || !nzchar(as.character(staxid))) return(NA_character_)
  ids <- unlist(strsplit(as.character(staxid), "[;,[:space:]]+", perl = TRUE), use.names = FALSE)
  ids <- ids[nzchar(ids) & ids != "0"]
  if (length(ids) == 0L) return(NA_character_)
  ids[[1L]]
}

read_ncbi_taxdump <- function(taxdump_dir) {
  nodes_file <- file.path(taxdump_dir, "nodes.dmp")
  names_file <- file.path(taxdump_dir, "names.dmp")
  check_file_arg(nodes_file, "nodes.dmp")
  check_file_arg(names_file, "names.dmp")

  nodes <- utils::read.delim(
    nodes_file,
    sep = "|",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  names <- utils::read.delim(
    names_file,
    sep = "|",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  nodes <- data.frame(
    taxid = trimws(nodes[[1L]]),
    parent = trimws(nodes[[2L]]),
    rank = trimws(nodes[[3L]]),
    stringsAsFactors = FALSE
  )
  names <- data.frame(
    taxid = trimws(names[[1L]]),
    name = trimws(names[[2L]]),
    class = trimws(names[[4L]]),
    stringsAsFactors = FALSE
  )
  sci <- names[names$class == "scientific name", c("taxid", "name"), drop = FALSE]
  node_index <- stats::setNames(seq_len(nrow(nodes)), nodes$taxid)
  name_index <- stats::setNames(sci$name, sci$taxid)
  list(nodes = nodes, node_index = node_index, names = name_index)
}

ncbi_taxid_lineage <- function(taxid, taxdump) {
  rank_map <- c(
    superkingdom = "kingdom",
    kingdom = "kingdom",
    phylum = "phylum",
    class = "class",
    order = "order",
    family = "family",
    genus = "genus",
    species = "species"
  )
  values <- stats::setNames(rep(NA_character_, 7L), unique(unname(rank_map)))
  seen <- character()
  current <- as.character(taxid)
  while (!is.na(current) && nzchar(current) && !current %in% seen) {
    seen <- c(seen, current)
    row_id <- taxdump$node_index[[current]]
    if (is.null(row_id) || is.na(row_id)) break
    node <- taxdump$nodes[row_id, , drop = FALSE]
    rank <- node$rank[[1L]]
    out_col <- if (rank %in% names(rank_map)) rank_map[[rank]] else NULL
    if (!is.null(out_col) && is.na(values[[out_col]])) {
      values[[out_col]] <- taxdump$names[[current]]
    }
    parent <- node$parent[[1L]]
    if (identical(parent, current)) break
    current <- parent
  }

  if (all(is.na(values))) return(NULL)
  taxonomy_path <- paste(stats::na.omit(unname(values)), collapse = ";")
  if (!nzchar(taxonomy_path)) taxonomy_path <- NA_character_
  as.list(c(taxonomy_path = taxonomy_path, values))
}

consensus_blast_prefix <- function(path) {
  prefix <- basename(path)
  prefix <- sub("[.]gz$", "", prefix)
  prefix <- sub("[.]fasta$", "", prefix, ignore.case = TRUE)
  prefix <- sub("[.]fa$", "", prefix, ignore.case = TRUE)
  prefix <- sub("[.]fna$", "", prefix, ignore.case = TRUE)
  prefix
}
