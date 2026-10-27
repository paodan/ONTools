test_that("build_multilocus_phylogeny validates the experimental design", {
  fasta <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet("ACGT"), fasta)

  expect_error(
    build_multilocus_phylogeny(
      data.frame(file = fasta, gene = "gene1"),
      write_alignment_pdfs = FALSE,
      bootstrap_replicates = 0
    ),
    "missing required columns"
  )

  incomplete <- data.frame(
    file = rep(fasta, 5),
    gene = c("gene1", "gene1", "gene1", "gene2", "gene2"),
    sample = c("s1", "s2", "s3", "s1", "s2")
  )
  expect_error(
    build_multilocus_phylogeny(
      incomplete,
      write_alignment_pdfs = FALSE,
      bootstrap_replicates = 0
    ),
    "Missing combinations"
  )
})

test_that("build_multilocus_phylogeny runs a small MAFFT analysis", {
  skip_if(Sys.which("mafft") == "", "MAFFT is not installed")
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")

  out_dir <- tempfile("multilocus-output-")
  input_dir <- tempfile("multilocus-input-")
  dir.create(input_dir)

  sample_ids <- paste0("s", 1:4)
  sequences <- list(
    gene1 = c(
      s1 = "ACGTACGTACGTACGTACGTACGTACGT",
      s2 = "ACGTACGTACGTACGTACGTACGTACGA",
      s3 = "ACGTACGTACGTTCGTACGTACGTACGA",
      s4 = "ACGTACGTTCGTTCGTACGTACGTACGA"
    ),
    gene2 = c(
      s1 = "TTGCAATTGCAATTGCAATTGCAATTGC",
      s2 = "TTGCAATTGCAATTGCAATTGCAATTGA",
      s3 = "TTGCAATTGCAACTGCAATTGCAATTGA",
      s4 = "TTGCAACTGCAACTGCAATTGCAATTGA"
    )
  )

  rows <- list()
  row_index <- 1L
  for (gene_name in names(sequences)) {
    for (sample_id in sample_ids) {
      fasta <- file.path(input_dir, paste0(gene_name, "-", sample_id, ".fasta"))
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(sequences[[gene_name]][[sample_id]]),
        fasta
      )
      rows[[row_index]] <- data.frame(
        file = fasta,
        gene = gene_name,
        sample = sample_id
      )
      row_index <- row_index + 1L
    }
  }
  file_info <- do.call(rbind, rows)

  result <- suppressWarnings(build_multilocus_phylogeny(
    file_info = file_info,
    primers = NULL,
    out_dir = out_dir,
    gene_order = c("gene1", "gene2"),
    tip_prefix = "Sample_",
    output_prefix = "test",
    ml_model = "JC",
    opt_inv = FALSE,
    opt_gamma = FALSE,
    bootstrap_replicates = 0,
    write_alignment_pdfs = FALSE,
    quiet = TRUE
  ))

  expect_s3_class(result$tree, "phylo")
  expect_equal(result$gene_order, c("gene1", "gene2"))
  expect_equal(result$sample_order, sample_ids)
  expect_equal(length(result$concatenated_alignment), 4)
  expect_true(all(file.exists(c(
    result$files$concatenated_fasta,
    result$files$tree,
    result$files$tree_pdf,
    result$files$difference_matrix
  ))))
  expect_equal(result$difference_matrix, t(result$difference_matrix))
  expect_equal(unname(diag(result$difference_matrix)), rep(0L, 4))
})
