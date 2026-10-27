#' Build a maximum-likelihood phylogeny from multiple loci
#'
#' `build_multilocus_phylogeny()` reads one FASTA file per sample and locus,
#' optionally trims primer-bounded regions, aligns each locus with MAFFT,
#' concatenates the aligned loci by sample, and builds a maximum-likelihood
#' phylogeny with optional bootstrap support.
#'
#' File-name parsing is deliberately kept outside this function. Supply a
#' `file_info` data frame containing the resolved file path, locus name, and
#' sample identifier. Every selected locus must contain exactly one sequence
#' for every selected sample.
#'
#' @param file_info Data frame with columns `file`, `gene`, and `sample`.
#'   Each row represents one FASTA file, and every FASTA file must contain
#'   exactly one DNA sequence.
#' @param primers Optional data frame with columns `gene`, `f`, and `r`, giving
#'   the forward and reverse primers for each locus. When `NULL`, primer
#'   trimming is skipped.
#' @param out_dir Output directory. It is created recursively when needed.
#' @param gene_order Optional character vector defining which loci to use and
#'   their concatenation order. When `NULL`, loci are used in their first
#'   appearance order in `file_info`.
#' @param sample_order Optional character vector defining sample order. When
#'   `NULL`, samples are used in their first appearance order in `file_info`.
#' @param tip_prefix Character prefix added to sample identifiers in the
#'   concatenated alignment and tree tip labels.
#' @param output_prefix Non-empty filename prefix for concatenated alignment,
#'   tree, PDF, and difference-matrix outputs.
#' @param mafft Path or command name for the MAFFT executable.
#' @param mafft_args Character vector of arguments passed to MAFFT before the
#'   input FASTA path.
#' @param primer_max_mismatch Non-negative integer passed to
#'   [trim_fasta_keep_primers()].
#' @param orient_to_forward Logical. If `TRUE`, primer-trimmed sequences are
#'   oriented from the forward primer to the reverse primer.
#' @param trim_failure How to handle sequences for which primers cannot be
#'   identified: issue a warning, stop with an error, or continue silently.
#' @param distance_model Substitution model passed to
#'   [phangorn::dist.ml()] when constructing the initial NJ tree.
#' @param ml_model Maximum-likelihood substitution model passed to
#'   [phangorn::optim.pml()].
#' @param opt_inv Logical. Optimize the proportion of invariant sites.
#' @param opt_gamma Logical. Optimize gamma-distributed rate variation.
#' @param rearrangement Tree rearrangement method passed to
#'   [phangorn::optim.pml()].
#' @param bootstrap_replicates Non-negative integer number of bootstrap
#'   replicates. Use `0` to skip bootstrap analysis.
#' @param bootstrap_opt_nni Logical. Optimize NNI rearrangements for each
#'   bootstrap replicate. This is slower but performs a more thorough search.
#' @param bootstrap_cutoff Numeric support threshold displayed on the tree.
#' @param seed Optional integer random seed for bootstrap sampling. The caller's
#'   random-number state is restored before the function returns.
#' @param write_alignment_pdfs Logical. Write colored alignment PDFs for each
#'   locus and for the concatenated alignment.
#' @param alignment_pdf_args Named list of additional arguments passed to
#'   [write_colored_alignment_pdf()]. The `alignment`, `pdf_file`, `overwrite`,
#'   and `quiet` arguments are managed internally and must not be supplied.
#' @param write_difference_matrix Logical. Write a symmetric matrix containing
#'   the number of unequal aligned columns between every pair of samples. Gaps
#'   count as differences.
#' @param tree_pdf_width,tree_pdf_height Positive numeric dimensions of the
#'   tree PDF in inches.
#' @param tree_title Character title displayed above the tree.
#' @param scale_bar_length Optional positive numeric length for the phylogram
#'   scale bar. When `NULL`, [ape::add.scale.bar()] chooses a suitable length.
#' @param overwrite Logical. If `FALSE`, stop before replacing any planned
#'   output file.
#' @param quiet Logical. If `FALSE`, report progress and output paths.
#'
#' @return Invisibly returns a list containing `aligned_by_gene`,
#'   `trim_summary_by_gene`, `concatenated_alignment`, `phy_data`, `nj_tree`,
#'   `fit`, `bootstrap_trees`, `tree`, `difference_matrix`, `sample_order`,
#'   `gene_order`, and `files`.
#'
#' @examples
#' \dontrun{
#' fasta_files <- list.files(
#'   "0_RawData/project/Sequence",
#'   pattern = "\\.fasta$",
#'   full.names = TRUE
#' )
#'
#' # Parse project-specific filenames outside the reusable function.
#' file_info <- data.frame(
#'   file = fasta_files,
#'   gene = sub("^.*-([^-]+)-[0-9]+\\.consensus\\.fasta$", "\\1",
#'              basename(fasta_files)),
#'   sample = sub("^.*-([0-9]+)\\.consensus\\.fasta$", "\\1",
#'                basename(fasta_files))
#' )
#'
#' primers <- data.frame(
#'   gene = c("ITS", "rbcL"),
#'   f = c("TCCGTAGGTGAACCTGCGG", "TGTCACCACAAACAGAAAC"),
#'   r = c("TCCTCCGCTTATTGATATG", "CGCATGATCCTGCAGTAG")
#' )
#'
#' result <- build_multilocus_phylogeny(
#'   file_info = file_info,
#'   primers = primers,
#'   out_dir = "phylogeny",
#'   gene_order = c("ITS", "rbcL"),
#'   tip_prefix = "Sample_",
#'   bootstrap_replicates = 100
#' )
#' result$tree
#' }
#'
#' @export
build_multilocus_phylogeny <- function(
    file_info,
    primers = NULL,
    out_dir = "phylogeny",
    gene_order = NULL,
    sample_order = NULL,
    tip_prefix = "",
    output_prefix = "multilocus",
    mafft = "mafft",
    mafft_args = "--auto",
    primer_max_mismatch = 2L,
    orient_to_forward = TRUE,
    trim_failure = c("warn", "error", "ignore"),
    distance_model = "JC69",
    ml_model = "GTR",
    opt_inv = TRUE,
    opt_gamma = TRUE,
    rearrangement = "NNI",
    bootstrap_replicates = 100L,
    bootstrap_opt_nni = FALSE,
    bootstrap_cutoff = 50,
    seed = 123L,
    write_alignment_pdfs = TRUE,
    alignment_pdf_args = list(
      bases_per_block = 100L,
      blocks_per_page = 3L,
      font_size = 8,
      character_spacing = 0.8,
      line_spacing = 0.9
    ),
    write_difference_matrix = TRUE,
    tree_pdf_width = 8,
    tree_pdf_height = 6,
    tree_title = "Maximum-likelihood phylogeny",
    scale_bar_length = NULL,
    overwrite = FALSE,
    quiet = FALSE) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop(
      "Package `ape` is required. Install it with `install.packages(\"ape\")`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    stop(
      "Package `phangorn` is required. Install it with ",
      "`install.packages(\"phangorn\")`.",
      call. = FALSE
    )
  }

  validated <- validate_multilocus_file_info(
    file_info = file_info,
    gene_order = gene_order,
    sample_order = sample_order
  )
  file_info <- validated$file_info
  gene_order <- validated$gene_order
  sample_order <- validated$sample_order

  check_scalar_character(out_dir, "out_dir")
  if (!is.character(tip_prefix) || length(tip_prefix) != 1L ||
      is.na(tip_prefix)) {
    stop("`tip_prefix` must be a single character string.", call. = FALSE)
  }
  check_scalar_character(output_prefix, "output_prefix")
  if (grepl("[/\\\\]", output_prefix)) {
    stop("`output_prefix` must not contain path separators.", call. = FALSE)
  }
  check_scalar_character(mafft, "mafft")
  check_scalar_character(distance_model, "distance_model")
  check_scalar_character(ml_model, "ml_model")
  check_scalar_character(rearrangement, "rearrangement")
  check_scalar_character(tree_title, "tree_title")
  if (!is.character(mafft_args) || anyNA(mafft_args)) {
    stop("`mafft_args` must be a character vector without missing values.",
         call. = FALSE)
  }

  primer_max_mismatch <- validate_nonnegative_integer(
    primer_max_mismatch,
    "primer_max_mismatch"
  )
  bootstrap_replicates <- validate_nonnegative_integer(
    bootstrap_replicates,
    "bootstrap_replicates"
  )
  if (bootstrap_replicates == 1L) {
    stop("`bootstrap_replicates` must be 0 or at least 2.", call. = FALSE)
  }
  bootstrap_cutoff <- validate_multilocus_number(
    bootstrap_cutoff,
    "bootstrap_cutoff",
    lower = 0,
    upper = 100
  )
  tree_pdf_width <- validate_multilocus_number(
    tree_pdf_width,
    "tree_pdf_width",
    lower = 0,
    lower_inclusive = FALSE
  )
  tree_pdf_height <- validate_multilocus_number(
    tree_pdf_height,
    "tree_pdf_height",
    lower = 0,
    lower_inclusive = FALSE
  )
  if (!is.null(scale_bar_length)) {
    scale_bar_length <- validate_multilocus_number(
      scale_bar_length,
      "scale_bar_length",
      lower = 0,
      lower_inclusive = FALSE
    )
  }
  if (!is.null(seed)) {
    seed <- validate_nonnegative_integer(seed, "seed")
  }

  check_logical_scalar(orient_to_forward, "orient_to_forward")
  check_logical_scalar(opt_inv, "opt_inv")
  check_logical_scalar(opt_gamma, "opt_gamma")
  check_logical_scalar(bootstrap_opt_nni, "bootstrap_opt_nni")
  check_logical_scalar(write_alignment_pdfs, "write_alignment_pdfs")
  check_logical_scalar(write_difference_matrix, "write_difference_matrix")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(quiet, "quiet")
  trim_failure <- match.arg(trim_failure)

  if (!is.list(alignment_pdf_args) ||
      (length(alignment_pdf_args) > 0L &&
       (is.null(names(alignment_pdf_args)) ||
        any(!nzchar(names(alignment_pdf_args)))))) {
    stop("`alignment_pdf_args` must be a fully named list.", call. = FALSE)
  }
  reserved_pdf_args <- c("alignment", "pdf_file", "overwrite", "quiet")
  duplicated_pdf_args <- intersect(names(alignment_pdf_args), reserved_pdf_args)
  if (length(duplicated_pdf_args) > 0L) {
    stop(
      "Do not supply internally managed arguments in `alignment_pdf_args`: ",
      paste(duplicated_pdf_args, collapse = ", "),
      call. = FALSE
    )
  }

  primers <- validate_multilocus_primers(primers, gene_order)
  mafft_path <- Sys.which(mafft)
  if (!nzchar(mafft_path)) {
    stop(
      "MAFFT executable was not found: ", mafft,
      ". Install MAFFT or provide its full path with `mafft`.",
      call. = FALSE
    )
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  safe_gene_names <- sanitize_multilocus_names(gene_order)
  if (anyDuplicated(safe_gene_names)) {
    stop(
      "Locus names are not unique after filename sanitization: ",
      paste(gene_order, collapse = ", "),
      call. = FALSE
    )
  }
  names(safe_gene_names) <- gene_order

  files <- make_multilocus_output_paths(
    out_dir = out_dir,
    output_prefix = output_prefix,
    gene_order = gene_order,
    safe_gene_names = safe_gene_names,
    trim_primers = !is.null(primers),
    write_alignment_pdfs = write_alignment_pdfs,
    write_difference_matrix = write_difference_matrix
  )
  planned_outputs <- unlist(files, use.names = FALSE)
  existing_outputs <- planned_outputs[file.exists(planned_outputs)]
  if (length(existing_outputs) > 0L && !isTRUE(overwrite)) {
    stop(
      "Output file already exists. Set `overwrite = TRUE` to replace: ",
      paste(existing_outputs, collapse = ", "),
      call. = FALSE
    )
  }

  internal_ids <- sprintf("OTSEQ%06d", seq_along(sample_order))
  names(internal_ids) <- sample_order
  aligned_by_gene <- structure(
    vector("list", length(gene_order)),
    names = gene_order
  )
  trim_summary_by_gene <- structure(
    vector("list", length(gene_order)),
    names = gene_order
  )

  for (gene_name in gene_order) {
    if (!isTRUE(quiet)) {
      message("Processing locus: ", gene_name)
    }
    rows <- file_info[file_info$gene == gene_name, , drop = FALSE]
    rows <- rows[match(sample_order, rows$sample), , drop = FALSE]
    sequences <- do.call(
      c,
      lapply(seq_len(nrow(rows)), function(i) {
        sequence <- Biostrings::readDNAStringSet(rows$file[[i]])
        if (length(sequence) != 1L) {
          stop(
            "Each FASTA file must contain exactly one sequence: ",
            rows$file[[i]],
            call. = FALSE
          )
        }
        names(sequence) <- internal_ids[[rows$sample[[i]]]]
        sequence
      })
    )

    gene_files <- files$by_gene[[gene_name]]
    Biostrings::writeXStringSet(sequences, gene_files$unaligned_fasta)
    mafft_input <- gene_files$unaligned_fasta

    if (!is.null(primers)) {
      primer_row <- primers[primers$gene == gene_name, , drop = FALSE]
      trim_summary <- trim_fasta_keep_primers(
        input_fasta = gene_files$unaligned_fasta,
        output_fasta = gene_files$trimmed_fasta,
        f_primer = primer_row$f[[1L]],
        r_primer = primer_row$r[[1L]],
        max_mismatch = primer_max_mismatch,
        orient_to_forward = orient_to_forward
      )
      trim_summary$sample <- sample_order[
        match(trim_summary$seq_name, internal_ids)
      ]
      trim_summary_by_gene[[gene_name]] <- trim_summary
      failed_samples <- trim_summary$sample[!trim_summary$trimmed]
      if (length(failed_samples) > 0L) {
        trim_message <- paste0(
          "Primer trimming failed for locus ", gene_name, ": ",
          paste(failed_samples, collapse = ", "),
          ". Untrimmed sequences were retained."
        )
        if (identical(trim_failure, "error")) {
          stop(trim_message, call. = FALSE)
        }
        if (identical(trim_failure, "warn")) {
          warning(trim_message, call. = FALSE)
        }
      }
      mafft_input <- gene_files$trimmed_fasta
    }

    mafft_log <- tempfile(
      pattern = paste0("mafft-", safe_gene_names[[gene_name]], "-"),
      tmpdir = out_dir,
      fileext = ".log"
    )
    status <- system2(
      mafft_path,
      args = c(mafft_args, shQuote(mafft_input)),
      stdout = gene_files$aligned_fasta,
      stderr = mafft_log
    )
    if (!identical(status, 0L)) {
      log_text <- if (file.exists(mafft_log)) {
        paste(utils::tail(readLines(mafft_log, warn = FALSE), 10L), collapse = "\n")
      } else {
        "No MAFFT diagnostic log was produced."
      }
      stop(
        "MAFFT alignment failed for locus `", gene_name,
        "` with exit status ", status, ".\n", log_text,
        call. = FALSE
      )
    }
    unlink(mafft_log)

    alignment <- Biostrings::readDNAStringSet(gene_files$aligned_fasta)
    alignment_ids <- sub("[[:space:]].*$", "", names(alignment))
    if (!setequal(alignment_ids, internal_ids)) {
      stop(
        "MAFFT output identifiers do not match the expected samples for locus `",
        gene_name, "`.",
        call. = FALSE
      )
    }
    names(alignment) <- sample_order[match(alignment_ids, internal_ids)]
    alignment <- alignment[match(sample_order, names(alignment))]
    if (length(unique(Biostrings::width(alignment))) != 1L) {
      stop("Aligned sequences have unequal widths for locus `", gene_name,
           "`.", call. = FALSE)
    }
    aligned_by_gene[[gene_name]] <- alignment

    if (isTRUE(write_alignment_pdfs)) {
      do.call(
        write_colored_alignment_pdf,
        c(
          list(
            alignment = alignment,
            pdf_file = gene_files$alignment_pdf,
            overwrite = overwrite,
            quiet = quiet
          ),
          alignment_pdf_args
        )
      )
    }
  }

  concatenated_sequences <- vapply(sample_order, function(sample_id) {
    paste0(
      vapply(aligned_by_gene, function(alignment) {
        as.character(alignment[match(sample_id, names(alignment))][[1L]])
      }, character(1)),
      collapse = ""
    )
  }, character(1))
  concatenated_alignment <- Biostrings::DNAStringSet(concatenated_sequences)
  names(concatenated_alignment) <- paste0(tip_prefix, sample_order)
  Biostrings::writeXStringSet(
    concatenated_alignment,
    files$concatenated_fasta
  )

  if (isTRUE(write_alignment_pdfs)) {
    do.call(
      write_colored_alignment_pdf,
      c(
        list(
          alignment = concatenated_alignment,
          pdf_file = files$concatenated_alignment_pdf,
          overwrite = overwrite,
          quiet = quiet
        ),
        alignment_pdf_args
      )
    )
  }

  difference_matrix <- multilocus_difference_matrix(concatenated_alignment)
  if (isTRUE(write_difference_matrix)) {
    utils::write.table(
      difference_matrix,
      file = files$difference_matrix,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }

  phy_data <- phangorn::read.phyDat(
    file = files$concatenated_fasta,
    format = "fasta",
    type = "DNA"
  )
  distance_matrix <- phangorn::dist.ml(
    phy_data,
    model = distance_model
  )
  nj_tree <- phangorn::NJ(distance_matrix)
  fit <- phangorn::pml(nj_tree, data = phy_data)
  fit <- phangorn::optim.pml(
    fit,
    model = ml_model,
    optInv = opt_inv,
    optGamma = opt_gamma,
    rearrangement = rearrangement,
    control = phangorn::pml.control(trace = 0)
  )

  bootstrap_trees <- NULL
  tree <- fit$tree
  if (bootstrap_replicates > 0L) {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                              inherits = FALSE)
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    if (!isTRUE(quiet)) {
      message("Running ", bootstrap_replicates, " bootstrap replicates")
    }
    bootstrap_trees <- tryCatch(
      phangorn::bootstrap.pml(
        fit,
        bs = bootstrap_replicates,
        optNni = bootstrap_opt_nni,
        control = phangorn::pml.control(trace = 0)
      ),
      error = function(error) {
        stop(
          "Bootstrap analysis failed. The alignment may contain too little ",
          "phylogenetic information; try a simpler model or set ",
          "`bootstrap_replicates = 0` to inspect the ML tree first. ",
          "Original error: ", conditionMessage(error),
          call. = FALSE
        )
      }
    )
    tree <- phangorn::plotBS(
      fit$tree,
      bootstrap_trees,
      p = bootstrap_cutoff,
      plot = FALSE
    )
  }

  ape::write.tree(tree, file = files$tree)
  grDevices::pdf(
    files$tree_pdf,
    width = tree_pdf_width,
    height = tree_pdf_height,
    useDingbats = FALSE
  )
  tree_device_open <- TRUE
  on.exit({
    if (isTRUE(tree_device_open)) {
      grDevices::dev.off()
    }
  }, add = TRUE)
  if (is.null(bootstrap_trees)) {
    ape::plot.phylo(fit$tree, type = "phylogram", main = tree_title)
  } else {
    phangorn::plotBS(
      fit$tree,
      bootstrap_trees,
      p = bootstrap_cutoff,
      type = "phylogram",
      main = tree_title
    )
  }
  scale_args <- list(x = "bottomleft", lwd = 2, cex = 0.8)
  if (!is.null(scale_bar_length)) {
    scale_args$length <- scale_bar_length
  }
  do.call(ape::add.scale.bar, scale_args)
  grDevices::dev.off()
  tree_device_open <- FALSE

  if (!isTRUE(quiet)) {
    message("Concatenated alignment: ", files$concatenated_fasta)
    message("Tree file: ", files$tree)
    message("Tree PDF: ", files$tree_pdf)
  }

  invisible(list(
    aligned_by_gene = aligned_by_gene,
    trim_summary_by_gene = trim_summary_by_gene,
    concatenated_alignment = concatenated_alignment,
    phy_data = phy_data,
    nj_tree = nj_tree,
    fit = fit,
    bootstrap_trees = bootstrap_trees,
    tree = tree,
    difference_matrix = difference_matrix,
    sample_order = sample_order,
    gene_order = gene_order,
    files = files
  ))
}

validate_multilocus_file_info <- function(file_info,
                                          gene_order = NULL,
                                          sample_order = NULL) {
  if (!is.data.frame(file_info)) {
    stop("`file_info` must be a data frame.", call. = FALSE)
  }
  required_columns <- c("file", "gene", "sample")
  missing_columns <- setdiff(required_columns, names(file_info))
  if (length(missing_columns) > 0L) {
    stop(
      "`file_info` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(file_info) == 0L) {
    stop("`file_info` must contain at least one row.", call. = FALSE)
  }
  file_info <- file_info[, required_columns, drop = FALSE]
  file_info[] <- lapply(file_info, as.character)
  invalid <- vapply(file_info, function(column) {
    anyNA(column) || any(!nzchar(column))
  }, logical(1))
  if (any(invalid)) {
    stop(
      "`file_info` columns must not contain missing or empty values: ",
      paste(names(invalid)[invalid], collapse = ", "),
      call. = FALSE
    )
  }
  missing_files <- unique(file_info$file[!file.exists(file_info$file)])
  if (length(missing_files) > 0L) {
    stop(
      "FASTA files do not exist: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  file_info$file <- normalizePath(file_info$file, mustWork = TRUE)

  if (is.null(gene_order)) {
    gene_order <- unique(file_info$gene)
  } else {
    gene_order <- validate_multilocus_id_order(gene_order, "gene_order")
    unknown_genes <- setdiff(gene_order, file_info$gene)
    if (length(unknown_genes) > 0L) {
      stop("Unknown loci in `gene_order`: ",
           paste(unknown_genes, collapse = ", "), call. = FALSE)
    }
    file_info <- file_info[file_info$gene %in% gene_order, , drop = FALSE]
  }
  if (length(gene_order) < 1L) {
    stop("At least one locus is required.", call. = FALSE)
  }

  available_samples <- unique(file_info$sample)
  if (is.null(sample_order)) {
    sample_order <- available_samples
  } else {
    sample_order <- validate_multilocus_id_order(sample_order, "sample_order")
    if (!setequal(sample_order, available_samples)) {
      stop(
        "`sample_order` must contain exactly the samples in the selected loci.",
        call. = FALSE
      )
    }
  }
  if (length(sample_order) < 3L) {
    stop("At least three samples are required to build a tree.", call. = FALSE)
  }

  pair_keys <- paste(file_info$gene, file_info$sample, sep = "\r")
  if (anyDuplicated(pair_keys)) {
    duplicated_rows <- unique(pair_keys[duplicated(pair_keys)])
    stop(
      "Each locus/sample combination must occur exactly once. Duplicates: ",
      paste(gsub("\r", "/", duplicated_rows, fixed = TRUE), collapse = ", "),
      call. = FALSE
    )
  }
  expected_keys <- as.vector(outer(gene_order, sample_order, paste, sep = "\r"))
  missing_keys <- setdiff(expected_keys, pair_keys)
  if (length(missing_keys) > 0L) {
    stop(
      "Every locus must contain every sample. Missing combinations: ",
      paste(gsub("\r", "/", missing_keys, fixed = TRUE), collapse = ", "),
      call. = FALSE
    )
  }

  list(
    file_info = file_info,
    gene_order = gene_order,
    sample_order = sample_order
  )
}

validate_multilocus_primers <- function(primers, gene_order) {
  if (is.null(primers)) {
    return(NULL)
  }
  if (!is.data.frame(primers)) {
    stop("`primers` must be `NULL` or a data frame.", call. = FALSE)
  }
  required_columns <- c("gene", "f", "r")
  missing_columns <- setdiff(required_columns, names(primers))
  if (length(missing_columns) > 0L) {
    stop("`primers` is missing required columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  primers <- primers[, required_columns, drop = FALSE]
  primers[] <- lapply(primers, as.character)
  if (anyNA(primers) || any(!nzchar(as.matrix(primers)))) {
    stop("`primers` must not contain missing or empty values.", call. = FALSE)
  }
  primer_counts <- table(factor(primers$gene, levels = gene_order))
  if (any(primer_counts != 1L)) {
    stop(
      "`primers` must contain exactly one row for every selected locus.",
      call. = FALSE
    )
  }
  primers
}

validate_multilocus_id_order <- function(x, name) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x)) ||
      anyDuplicated(x)) {
    stop("`", name,
         "` must be a non-empty character vector of unique identifiers.",
         call. = FALSE)
  }
  x
}

validate_multilocus_number <- function(x,
                                       name,
                                       lower = -Inf,
                                       upper = Inf,
                                       lower_inclusive = TRUE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x <= upper
  if (isTRUE(lower_inclusive)) {
    valid <- valid && x >= lower
  } else {
    valid <- valid && x > lower
  }
  if (!valid) {
    stop("`", name, "` is outside its allowed numeric range.", call. = FALSE)
  }
  as.numeric(x)
}

sanitize_multilocus_names <- function(x) {
  sanitized <- gsub("[^A-Za-z0-9._-]+", "_", x)
  sanitized[!nzchar(sanitized)] <- "locus"
  sanitized
}

make_multilocus_output_paths <- function(out_dir,
                                         output_prefix,
                                         gene_order,
                                         safe_gene_names,
                                         trim_primers,
                                         write_alignment_pdfs,
                                         write_difference_matrix) {
  by_gene <- structure(lapply(gene_order, function(gene_name) {
    safe_gene <- safe_gene_names[[gene_name]]
    paths <- list(
      unaligned_fasta = file.path(
        out_dir,
        paste0(safe_gene, ".unaligned.fasta")
      ),
      aligned_fasta = file.path(
        out_dir,
        paste0(safe_gene, ".aligned.fasta")
      )
    )
    if (isTRUE(trim_primers)) {
      paths$trimmed_fasta <- file.path(
        out_dir,
        paste0(safe_gene, ".unaligned_trim.fasta")
      )
    }
    if (isTRUE(write_alignment_pdfs)) {
      paths$alignment_pdf <- file.path(
        out_dir,
        paste0("DNA_alignment_", safe_gene, ".pdf")
      )
    }
    paths
  }), names = gene_order)

  paths <- list(
    by_gene = by_gene,
    concatenated_fasta = file.path(
      out_dir,
      paste0(output_prefix, "_concatenated.aligned.fasta")
    ),
    tree = file.path(
      out_dir,
      paste0(output_prefix, "_ML_bootstrap.tree")
    ),
    tree_pdf = file.path(
      out_dir,
      paste0(output_prefix, "_ML_bootstrap.pdf")
    )
  )
  if (isTRUE(write_alignment_pdfs)) {
    paths$concatenated_alignment_pdf <- file.path(
      out_dir,
      paste0("DNA_alignment_", output_prefix, "_concatenated.pdf")
    )
  }
  if (isTRUE(write_difference_matrix)) {
    paths$difference_matrix <- file.path(
      out_dir,
      paste0(output_prefix, "_pairwise_difference_counts.tsv")
    )
  }
  paths
}

multilocus_difference_matrix <- function(alignment) {
  sequence_matrix <- do.call(
    rbind,
    strsplit(as.character(alignment), split = "")
  )
  number_of_sequences <- nrow(sequence_matrix)
  result <- matrix(
    0L,
    nrow = number_of_sequences,
    ncol = number_of_sequences,
    dimnames = list(names(alignment), names(alignment))
  )
  if (number_of_sequences > 1L) {
    for (i in seq_len(number_of_sequences - 1L)) {
      for (j in seq.int(i + 1L, number_of_sequences)) {
        difference_count <- sum(sequence_matrix[i, ] != sequence_matrix[j, ])
        result[i, j] <- difference_count
        result[j, i] <- difference_count
      }
    }
  }
  result
}
