test_that("compare_consensus compares consensus sequences with shared names", {
  consensus1 <- Biostrings::DNAStringSet(c(
    barcode01 = "ATGCCGTAAA",
    barcode02 = "AACCGGTT"
  ))
  consensus2 <- Biostrings::DNAStringSet(c(
    barcode01 = "ATGCCGGTAAAC",
    barcode02 = "AACCGGTT"
  ))

  res <- compare_consensus(consensus1, consensus2)

  expect_s3_class(res, "data.frame")
  expect_equal(res$seq1_name, c("barcode01", "barcode02"))
  expect_equal(res$seq2_name, c("barcode01", "barcode02"))
  expect_true(res$distance[[1]] > 0)
  expect_equal(res$distance[[2]], 0L)
  expect_equal(res$seq1_len, c(10L, 8L))
  expect_equal(res$seq2_len, c(12L, 8L))
  expect_true(all(nzchar(res$aligned_seq1)))
  expect_true(all(nzchar(res$aligned_seq2)))
})

test_that("compare_consensus uses an explicit sample map", {
  consensus1 <- Biostrings::DNAStringSet(c(A = "ACGT", B = "AAAA"))
  consensus2 <- Biostrings::DNAStringSet(c(X = "ACGT", Y = "TTTT"))
  sample_map <- data.frame(
    seq1_name = c("A", "B", "missing"),
    seq2_name = c("X", "Y", "X")
  )

  res <- compare_consensus(consensus1, consensus2, sampleMap = sample_map)

  expect_equal(nrow(res), 2L)
  expect_equal(res$seq1_name, c("A", "B"))
  expect_equal(res$seq2_name, c("X", "Y"))
  expect_equal(res$distance[[1]], 0L)
})

test_that("compare_consensus reads FASTA file inputs", {
  fasta1 <- tempfile(fileext = ".fasta")
  fasta2 <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(barcode01 = "ACGT", barcode02 = "AAAA")),
    fasta1
  )
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(barcode01 = "ACGT", barcode02 = "AAAT")),
    fasta2
  )

  res <- compare_consensus(fasta1, fasta2)

  expect_equal(res$seq1_name, c("barcode01", "barcode02"))
  expect_equal(res$distance, c(0L, 1L))
})

test_that("compare_consensus returns an empty result when no pairs match", {
  consensus1 <- Biostrings::DNAStringSet(c(A = "ACGT"))
  consensus2 <- Biostrings::DNAStringSet(c(B = "ACGT"))

  res <- compare_consensus(consensus1, consensus2)

  expect_equal(nrow(res), 0L)
  expect_named(
    res,
    c(
      "seq1_name", "seq2_name", "distance", "seq1_len", "seq2_len",
      "seq2_orientation", "seq1", "seq2", "aligned_seq1", "aligned_seq2"
    )
  )
})

test_that("compare_consensus validates inputs", {
  consensus1 <- Biostrings::DNAStringSet(c(A = "ACGT"))
  consensus2 <- Biostrings::DNAStringSet(c(A = "ACGT"))

  expect_error(
    compare_consensus("missing.fasta", consensus2),
    "consensus1"
  )
  expect_error(
    compare_consensus(consensus1, consensus2, sampleMap = data.frame(A = "A")),
    "sampleMap"
  )
})
