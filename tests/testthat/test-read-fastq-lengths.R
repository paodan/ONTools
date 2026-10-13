test_that("read_fastq_lengths reads plain FASTQ in chunks", {
  fastq <- tempfile(fileext = ".fastq")
  writeLines(
    c(
      "@read1", "ACGT", "+", "!!!!",
      "@read2", "ACGTAC", "+", "!!!!!!",
      "@read3", "AC", "+", "!!"
    ),
    fastq
  )

  expect_equal(
    read_fastq_lengths(fastq, chunk_records = 1L),
    c(4L, 6L, 2L)
  )
})

test_that("read_fastq_lengths reads gzip and literal quote qualities", {
  fastq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fastq, open = "wt")
  writeLines(
    c(
      "@read1", "ACGT", "+", "\"\"\"\"",
      "@read2", "ACGTA", "+", "!!!!!"
    ),
    con
  )
  close(con)

  expect_equal(read_fastq_lengths(fastq), c(4L, 5L))
})

test_that("read_fastq_lengths returns integer zero for an empty FASTQ", {
  fastq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fastq, open = "wt")
  close(con)

  lengths <- read_fastq_lengths(fastq)
  expect_type(lengths, "integer")
  expect_length(lengths, 0L)
})

test_that("read_fastq_lengths validates arguments and records", {
  wrong_suffix <- tempfile(fileext = ".fq")
  file.create(wrong_suffix)
  expect_error(read_fastq_lengths(wrong_suffix), "must end")
  expect_error(read_fastq_lengths("missing.fastq"), "does not exist")

  valid <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), valid)
  expect_error(read_fastq_lengths(valid, chunk_records = 0), "chunk_records")

  incomplete <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+"), incomplete)
  expect_error(read_fastq_lengths(incomplete), "Incomplete FASTQ")

  invalid_header <- tempfile(fileext = ".fastq")
  writeLines(c("read1", "ACGT", "+", "!!!!"), invalid_header)
  expect_error(read_fastq_lengths(invalid_header), "Invalid four-line")

  unequal <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!"), unequal)
  expect_error(read_fastq_lengths(unequal), "lengths differ")
})
