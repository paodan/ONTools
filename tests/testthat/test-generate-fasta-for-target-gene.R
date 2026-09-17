test_that("generateFastaForTargetGene writes target and replacement FASTA files from PAF", {
  genome <- tempfile(fileext = ".fasta")
  gene <- tempfile(fileext = ".fasta")
  replacement <- tempfile(fileext = ".fasta")
  paf <- tempfile(fileext = ".paf")
  out_dir <- tempfile()

  writeLines(c(">chr1 some description", "AAAACCCCGGGGTTTTAAAA"), genome)
  writeLines(c(">gene1 target", "CCCG"), gene)
  writeLines(c(">replacement1", "TT"), replacement)
  writeLines(
    paste(
      "gene1", 4, 0, 4, "+",
      "chr1", 20, 5, 9,
      4, 4, 60,
      sep = "\t"
    ),
    paf
  )

  res <- generateFastaForTargetGene(
    output_path = out_dir,
    genome_fasta = genome,
    gene_seq_fasta = gene,
    replace_seq_fastq = replacement,
    h1_len = 2,
    h2_len = 2,
    h1up_len = 1,
    h2down_len = 1,
    paf_file = paf,
    echo = FALSE
  )

  expect_s3_class(res$paf, "data.frame")
  expect_length(res$seq, 1)
  expect_length(res$seq[[1]], 10)
  expect_length(res$written_files, 10)
  expect_true(all(file.exists(res$written_files)))
  expect_true(all(grepl("[.]fasta$", basename(res$written_files))))

  seq_widths <- stats::setNames(as.integer(BiocGenerics::width(res$seq[[1]])), names(res$seq[[1]]))
  expect_equal(seq_widths["gene1_h1_chr1_4_5"], c(gene1_h1_chr1_4_5 = 2L))
  expect_equal(seq_widths["gene1_h2_chr1_10_11"], c(gene1_h2_chr1_10_11 = 2L))
  expect_equal(seq_widths["replacement1_h1_replace_h2_chr1_6_9"], c(replacement1_h1_replace_h2_chr1_6_9 = 6L))
  expect_equal(seq_widths["replacement1_h1up_h1_replace_h2_h2down_chr1_6_9"], c(replacement1_h1up_h1_replace_h2_h2down_chr1_6_9 = 8L))
})

test_that("generateFastaForTargetGene can return a dry-run plan", {
  genome <- tempfile(fileext = ".fasta")
  gene <- tempfile(fileext = ".fasta")
  replacement <- tempfile(fileext = ".fasta")
  out_dir <- tempfile()

  writeLines(c(">chr1", "AAAACCCCGGGG"), genome)
  writeLines(c(">gene1", "CCCC"), gene)
  writeLines(c(">replacement1", "TT"), replacement)

  res <- generateFastaForTargetGene(
    output_path = out_dir,
    genome_fasta = genome,
    gene_seq_fasta = gene,
    replace_seq_fastq = replacement,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_null(res$paf)
  expect_length(res$written_files, 0)
  expect_match(res$commands$minimap2, "minimap2")
})
