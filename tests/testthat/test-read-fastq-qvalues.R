test_that("read_fastq_qvalues calculates Phred+33 means in chunks", {
  fastq <- tempfile(fileext = ".fastq")
  writeLines(
    c(
      "@read1", "ACGT", "+", "IIII",
      "@read2", "ACGT", "+", "!+5?",
      "@read3", "AC", "+", "\"\""
    ),
    fastq
  )

  expect_equal(
    read_fastq_qvalues(fastq, chunk_records = 1L),
    c(40, -10 * log10(mean(10^(-c(0, 10, 20, 30) / 10))), 1)
  )
})

test_that("read_fastq_qvalues follows Dorado leading-base trimming", {
  fastq <- tempfile(fileext = ".fastq")
  writeLines(
    c(
      "@read1",
      paste0(rep("A", 61L), collapse = ""),
      "+",
      paste0(c(rep("!", 60L), "I"), collapse = "")
    ),
    fastq
  )

  expect_equal(read_fastq_qvalues(fastq), 40)
})

test_that("read_fastq_qvalues reads gzip and returns numeric zero", {
  fastq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fastq, open = "wt")
  writeLines(c("@read1", "ACGTA", "+", "55555"), con)
  close(con)
  expect_equal(read_fastq_qvalues(fastq), 20)

  empty <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(empty, open = "wt")
  close(con)
  qvalues <- read_fastq_qvalues(empty)
  expect_type(qvalues, "double")
  expect_length(qvalues, 0L)
})

test_that("read_fastq_qvalues validates arguments and records", {
  wrong_suffix <- tempfile(fileext = ".fq")
  file.create(wrong_suffix)
  expect_error(read_fastq_qvalues(wrong_suffix), "must end")
  expect_error(read_fastq_qvalues("missing.fastq"), "does not exist")

  valid <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "IIII"), valid)
  expect_error(read_fastq_qvalues(valid, chunk_records = 0), "chunk_records")

  incomplete <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+"), incomplete)
  expect_error(read_fastq_qvalues(incomplete), "Incomplete FASTQ")

  invalid_header <- tempfile(fileext = ".fastq")
  writeLines(c("read1", "ACGT", "+", "IIII"), invalid_header)
  expect_error(read_fastq_qvalues(invalid_header), "Invalid four-line")

  unequal <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "III"), unequal)
  expect_error(read_fastq_qvalues(unequal), "lengths differ")
})
