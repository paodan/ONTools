test_that("fastq_read_info parses plain and compressed FASTQ headers", {
  fastq_dir <- tempfile("fastq-info-")
  dir.create(fastq_dir)
  plain <- file.path(fastq_dir, "a.fastq")
  compressed <- file.path(fastq_dir, "b.fastq.gz")

  writeLines(c(
    "@read1 runid=run-a start_time=2025-01-01T00:00:00Z barcode=barcode01",
    "ACGT", "+", "!!!!",
    "@read2 runid=run-a start_time=2025-01-01T00:00:01.5Z barcode=barcode01 extra=a=b",
    "TGCA", "+", "####"
  ), plain)
  con <- gzfile(compressed, open = "wt")
  writeLines(c(
    "@read3 runid=run-b start_time=2025-01-01T08:00:02+08:00 barcode=barcode02",
    "AAAA", "+", "$$$$"
  ), con)
  close(con)

  plain_info <- fastq_read_info(
    fastq_dir,
    file_format = ".fastq",
    plot_time = FALSE
  )
  gzip_info <- fastq_read_info(
    fastq_dir,
    file_format = ".fastq.gz",
    plot_time = FALSE
  )

  expect_equal(plain_info$read, c("read1", "read2"))
  expect_equal(plain_info$runid, c("run-a", "run-a"))
  expect_equal(plain_info$extra, c(NA_character_, "a=b"))
  expect_s3_class(plain_info$start_time, "POSIXct")
  expect_equal(as.numeric(diff(plain_info$start_time)), 1.5)
  expect_equal(gzip_info$read, "read3")
  expect_equal(
    gzip_info$start_time,
    as.POSIXct("2025-01-01 00:00:02", tz = "UTC")
  )
})

test_that("fastq_read_info supports legacy argument names", {
  fastq_dir <- tempfile("fastq-info-")
  dir.create(fastq_dir)
  writeLines(
    c("@read1 barcode=barcode01", "ACGT", "+", "!!!!"),
    file.path(fastq_dir, "reads.fq")
  )

  info <- fastq_read_info(
    fastq_dir,
    fileFormat = ".fq",
    plotTime = FALSE
  )
  expect_equal(info$read, "read1")
})

test_that("fastq_read_info handles no matches and validates FASTQ records", {
  fastq_dir <- tempfile("fastq-info-")
  dir.create(fastq_dir)

  empty <- fastq_read_info(fastq_dir, plot_time = FALSE)
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("read", "start_time"))

  bad_fastq <- file.path(fastq_dir, "bad.fastq")
  writeLines(c("@read1", "ACGT", "+"), bad_fastq)
  expect_error(
    fastq_read_info(fastq_dir, file_format = ".fastq", plot_time = FALSE),
    "Incomplete FASTQ record"
  )
  expect_error(
    fastq_read_info(tempfile(), plot_time = FALSE),
    "existing directory"
  )
  expect_error(
    fastq_read_info(fastq_dir, plot_time = NA),
    "plot_time"
  )
})
