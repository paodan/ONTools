test_that("align_sanger_with_full_reference orients Sanger reads", {
  reference <- Biostrings::DNAStringSet(c(
    reference = "ATGCGTACCAAAGGGTTTCCCA"
  ))
  forward_read <- "GTACCAAAGG"
  reverse_read <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAString("AAAGGGTTTC")
  ))
  sanger <- Biostrings::DNAStringSet(c(
    read_forward = forward_read,
    read_reverse = reverse_read
  ))

  res <- align_sanger_with_full_reference(reference, sanger)

  expect_s3_class(res, "SangerMultipleAlignment")
  expect_equal(names(res$alignment), c("reference", "read_forward", "read_reverse"))
  expect_equal(res$orientation$sample, c("read_forward", "read_reverse"))
  expect_equal(res$orientation$direction[[1]], "forward")
  expect_equal(res$orientation$direction[[2]], "reverse_complement")
  expect_equal(
    as.character(res$oriented_sanger[["read_reverse"]]),
    "AAAGGGTTTC"
  )
  expect_true(all(res$orientation$selected_score >= 0))
})

test_that("align_sanger_with_full_reference reads FASTA paths and fills names", {
  reference_fasta <- tempfile(fileext = ".fasta")
  sanger_fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">refA", "ATGCGTACCAAAGGGTTTCCCA"), reference_fasta)
  writeLines(c(">readA", "GTACCAAAGG"), sanger_fasta)

  res <- align_sanger_with_full_reference(reference_fasta, sanger_fasta)

  expect_equal(names(res$alignment), c("refA", "readA"))
  expect_equal(res$orientation$sample, "readA")

  unnamed_reference <- Biostrings::DNAStringSet("ATGCGTACCAAAGGGTTTCCCA")
  unnamed_sanger <- Biostrings::DNAStringSet("GTACCAAAGG")
  res2 <- align_sanger_with_full_reference(unnamed_reference, unnamed_sanger)

  expect_equal(names(res2$alignment), c("reference", "sanger_1"))
})

test_that("align_sanger_with_full_reference validates inputs", {
  reference <- Biostrings::DNAStringSet(c(ref = "ATGCGTACCAAAGGGTTTCCCA"))
  sanger <- Biostrings::DNAStringSet(c(read = "GTACCAAAGG"))

  expect_error(
    align_sanger_with_full_reference(Biostrings::DNAStringSet(c("AAA", "CCC")), sanger),
    "exactly one"
  )
  expect_error(
    align_sanger_with_full_reference(reference, Biostrings::DNAStringSet()),
    "at least one"
  )
  expect_error(
    align_sanger_with_full_reference(reference, sanger, processors = 0),
    "processors"
  )
  expect_error(
    align_sanger_with_full_reference(reference, sanger, gap_opening = -1),
    "gap_opening"
  )
  expect_error(
    align_sanger_with_full_reference("missing.fasta", sanger),
    "does not exist"
  )
})
