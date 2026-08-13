test_that("trim_fasta_keep_primers trims extra bases while keeping primers", {
  f_primer <- "AGAGTTTGATCMTGGCTCAG"
  r_primer <- "TACGGYTACCTTGTTACGACTT"
  insert <- "ACGTACGT"
  expected <- paste0(
    "AGAGTTTGATCATGGCTCAG",
    insert,
    reverse_complement_string("TACGGTTACCTTGTTACGACTT")
  )
  input_seq <- paste0("AA", expected, "TT")
  input <- tempfile(fileext = ".fasta")
  output <- tempfile(fileext = ".fasta")

  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(sample1 = input_seq)),
    input
  )

  summary <- trim_fasta_keep_primers(
    input_fasta = input,
    output_fasta = output,
    f_primer = f_primer,
    r_primer = r_primer
  )
  out <- Biostrings::readDNAStringSet(output)

  expect_true(file.exists(output))
  expect_equal(as.character(out[[1]]), expected)
  expect_true(summary$trimmed[[1]])
  expect_equal(summary$orientation[[1]], "forward")
  expect_equal(summary$start[[1]], 3L)
  expect_equal(summary$end[[1]], nchar(input_seq) - 2L)
  expect_equal(summary$output_length[[1]], nchar(expected))
  expect_match(names(out)[[1]], "trimmed=TRUE")
})

test_that("trim_fasta_keep_primers orients reverse-complement matches forward", {
  f_primer <- "AGAGTTTGATCATGGCTCAG"
  r_primer <- "TACGGTTACCTTGTTACGACTT"
  forward_amplicon <- paste0(
    f_primer,
    "ACGTACGT",
    reverse_complement_string(r_primer)
  )
  reverse_input <- paste0("GG", reverse_complement_string(forward_amplicon), "CC")
  input <- tempfile(fileext = ".fasta")
  output <- tempfile(fileext = ".fasta")

  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(sample1 = reverse_input)),
    input
  )

  summary <- trim_fasta_keep_primers(
    input_fasta = input,
    output_fasta = output,
    f_primer = f_primer,
    r_primer = r_primer
  )
  out <- Biostrings::readDNAStringSet(output)

  expect_equal(as.character(out[[1]]), forward_amplicon)
  expect_equal(summary$orientation[[1]], "reverse_complement")
})

test_that("trim_fasta_keep_primers preserves untrimmed sequences when primers are absent", {
  input <- tempfile(fileext = ".fasta")
  output <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(sample1 = "ACGTACGT")),
    input
  )

  summary <- trim_fasta_keep_primers(
    input_fasta = input,
    output_fasta = output,
    f_primer = "AAAA",
    r_primer = "TTTT",
    max_mismatch = 0
  )
  out <- Biostrings::readDNAStringSet(output)

  expect_false(summary$trimmed[[1]])
  expect_true(is.na(summary$orientation[[1]]))
  expect_equal(as.character(out[[1]]), "ACGTACGT")
})

test_that("trim_fasta_keep_primers validates inputs", {
  input <- tempfile(fileext = ".fasta")
  output <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(sample1 = "ACGT")),
    input
  )

  expect_error(
    trim_fasta_keep_primers(
      input,
      output,
      f_primer = "AAAA",
      r_primer = "TTTT",
      max_mismatch = -1
    ),
    "max_mismatch"
  )
  expect_error(
    trim_fasta_keep_primers(
      "missing.fasta",
      output,
      f_primer = "AAAA",
      r_primer = "TTTT"
    ),
    "input_fasta"
  )
})
