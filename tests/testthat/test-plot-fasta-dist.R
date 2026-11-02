test_that("plot_fasta_dist returns a symmetric difference matrix", {
  fasta <- c(reference = "AAAA", sample_1 = "AAAT", sample_2 = "AATT")

  result <- plot_fasta_dist(fasta, cores = 1, silent = TRUE)

  expect_type(result$D, "integer")
  expect_equal(result$D, t(result$D))
  expect_equal(unname(diag(result$D)), rep(0L, 3L))
  expect_equal(
    result$D,
    matrix(
      c(0L, 1L, 2L, 1L, 0L, 1L, 2L, 1L, 0L),
      nrow = 3L,
      dimnames = list(names(fasta), names(fasta))
    )
  )
  expect_s3_class(result$plot, "pheatmap")
})

test_that("plot_fasta_dist accepts DNAStringSet and fills missing names", {
  fasta <- Biostrings::DNAStringSet(c("AAAA", "AAAT"))
  names(fasta) <- c("reference", "")

  result <- plot_fasta_dist(fasta, cores = 1, silent = TRUE)

  expect_equal(rownames(result$D), c("reference", "seq2"))
  expect_equal(colnames(result$D), c("reference", "seq2"))
  expect_equal(result$D[1, 2], 1L)
})

test_that("plot_fasta_dist validates its inputs", {
  expect_error(plot_fasta_dist("AAAA"), "at least two")
  expect_error(plot_fasta_dist(c("AAAA", NA_character_)), "missing or empty")
  expect_error(plot_fasta_dist(c("AAAA", "AAAT"), cores = 0), "positive integer")
  expect_error(plot_fasta_dist(list("AAAA", "AAAT")), "character vector")
})
