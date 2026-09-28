test_that("calculate_fastq_read_q validates input files", {
  expect_error(
    calculate_fastq_read_q(tempfile(fileext = ".fastq")),
    "does not exist"
  )
})

test_that("calculate_fastq_read_q returns per-read quality metrics", {
  skip_if_not_installed("ShortRead")

  fastq <- tempfile(fileext = ".fastq")
  writeLines(c(
    "@read1 comment",
    "ACGT",
    "+",
    "IIII",
    "@read2",
    "ACGT",
    "+",
    "!!!!"
  ), fastq)

  res <- calculate_fastq_read_q(fastq)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_equal(rownames(res), c("read1", "read2"))
  expect_equal(res$read_id, c("read1", "read2"))
  expect_equal(res$length, c(4L, 4L))
  expect_equal(res$arithmetic_mean_q, c(40, 0))
  expect_equal(res$mean_q, c(40, 0), tolerance = 1e-8)
  expect_equal(res$error_rate, c(1e-4, 1), tolerance = 1e-12)
  expect_equal(res$accuracy, c(0.9999, 0), tolerance = 1e-12)
})

test_that("calculate_fastq_read_q makes duplicate read-name row names unique", {
  skip_if_not_installed("ShortRead")

  fastq <- tempfile(fileext = ".fastq")
  writeLines(c(
    "@read1",
    "ACGT",
    "+",
    "IIII",
    "@read1",
    "ACGT",
    "+",
    "IIII"
  ), fastq)

  res <- calculate_fastq_read_q(fastq)

  expect_equal(res$read_id, c("read1", "read1"))
  expect_equal(rownames(res), c("read1", "read1.1"))
})
