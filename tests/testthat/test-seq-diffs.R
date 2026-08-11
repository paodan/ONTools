test_that("seq_diffs chooses reverse complement and reports differences", {
  seqs <- c(
    A = "ATGCCGTAAA",
    B = "TTACGGCATAA"
  )

  res <- seq_diffs(seqs)

  expect_equal(res$orientation, "reverse_complement")
  expect_equal(res$n_diff, 3L)
  expect_equal(nrow(res$differences), 3L)
  expect_true(all(res$differences$seq1_name == "A"))
  expect_true(all(res$differences$seq2_name == "B"))
  expect_equal(res$differences$seq1_pos, c(NA_integer_, NA_integer_, 10L))
  expect_equal(res$differences$seq2_pos, c(11L, 10L, NA_integer_))
  expect_equal(res$differences$seq1_base, c("-", "-", "A"))
  expect_equal(res$differences$seq2_base, c("T", "T", "-"))
  expect_true(all(res$differences$seq2_aligned_as == "seq2_reverse_complement"))
})

test_that("seq_diffs reports forward orientation differences", {
  seqs <- c(
    A = "ATGCCGTAAA",
    B = "ATGCCGGTAAAC"
  )

  res <- seq_diffs(seqs)

  expect_equal(res$orientation, "forward")
  expect_true(res$n_diff >= 1L)
  expect_true(all(res$differences$orientation == "forward"))
  expect_true(all(res$differences$seq2_aligned_as == "seq2_original"))
})

test_that("seq_diffs handles identical sequences with an empty differences table", {
  seqs <- c(
    A = "ATGCCGTAAA",
    B = "ATGCCGTAAA"
  )

  res <- seq_diffs(seqs)

  expect_equal(res$orientation, "forward")
  expect_equal(res$n_diff, 0L)
  expect_equal(nrow(res$differences), 0L)
  expect_named(
    res$differences,
    c(
      "aln_pos", "seq1_name", "seq1_pos", "seq1_base", "seq2_name",
      "seq2_pos", "seq2_base", "orientation", "seq2_aligned_as"
    )
  )
})
