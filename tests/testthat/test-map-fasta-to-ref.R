test_that("map_fasta_to_ref plans minimap2 and samtools commands", {
  query <- tempfile(fileext = ".fasta")
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">query1", "ACGTACGT"), query)
  writeLines(c(">chr1", "ACGTACGTACGT"), ref)

  res <- map_fasta_to_ref(
    fastaFile = query,
    ref = ref,
    output_file = bam,
    threads = 2,
    secondary = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_identical(res$status, NA_integer_)
  expect_equal(res$output_file, normalizePath(bam, mustWork = FALSE))
  expect_match(res$commands$minimap2_paf, "-x")
  expect_match(res$commands$minimap2_paf, "asm5")
  expect_match(res$commands$minimap2_paf, "--secondary=no", fixed = TRUE)
  expect_match(res$commands$samtools_sort, "sort")
  expect_match(res$commands$samtools_index, "index")
})
