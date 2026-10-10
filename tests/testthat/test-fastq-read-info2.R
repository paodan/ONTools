test_that("fastq_read_info2 parses and standardizes SAM-style headers", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  writeLines(c(
    paste0(
      "@new-read qs:f:15.715393 mx:i:4 ch:i:1803 rn:i:438 ",
      "st:Z:2026-09-20T14:53:20.912722+08:00 ns:i:3694 ts:i:10 ",
      "sm:f:93.99999 sd:f:24 sv:Z:pa du:f:0.7388 dx:i:0 ",
      "RG:Z:run_model_barcode86 DT:Z:2026-09-20T14:27:18+08:00 ",
      "PU:Z:PBO11046 LB:Z:20260920pcr SM:Z:barcode86 al:Z:sample86"
    ),
    "ACGT", "+", "!!!!"
  ), file.path(fastq_dir, "reads.fastq"))

  info <- fastq_read_info2(
    fastq_dir,
    file_format = ".fastq",
    plot_time = FALSE
  )

  expect_equal(info$read, "new-read")
  expect_equal(info$header_format, "sam_tag")
  expect_equal(info$mean_qscore, 15.715393)
  expect_equal(info$mux, 4L)
  expect_equal(info$ch, 1803L)
  expect_equal(info$read_number, 438L)
  expect_equal(info$signal_samples, 3684L)
  expect_equal(info$duration, 0.7388)
  expect_equal(info$barcode, "barcode86")
  expect_equal(info$barcode_alias, "sample86")
  expect_equal(info$flow_cell_id, "PBO11046")
  expect_s3_class(info$start_time, "POSIXct")
  expect_equal(
    info$start_time,
    as.POSIXct("2026-09-20 06:53:20.912722", tz = "UTC")
  )
})

test_that("fastq_read_info2 recognizes legacy and mixed header formats", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  writeLines(c(
    "@old-read runid=RUN1 read=12 ch=7 start_time=2026-01-01T00:00:00Z barcode=barcode01",
    "ACGT", "+", "!!!!",
    "@mixed-read ch=8 qs:f:12.5 st:Z:2026-01-01T00:00:01Z SM:Z:barcode02",
    "TGCA", "+", "####"
  ), file.path(fastq_dir, "reads.fastq"))

  info <- fastq_read_info2(
    fastq_dir,
    file_format = ".fastq",
    plot_time = FALSE
  )

  expect_equal(info$read, c("old-read", "mixed-read"))
  expect_equal(info$header_format, c("key_value", "mixed"))
  expect_equal(info$read_number, c(12L, NA_integer_))
  expect_equal(info$ch, c(7L, 8L))
  expect_equal(info$barcode, c("barcode01", "barcode02"))
  expect_equal(info$run_id, c("RUN1", NA_character_))
})

test_that("fastq_read_info2 handles aliases, unknown tags, and empty input", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  writeLines(
    c("@read1 al:Z:sample-A XY:Z:value", "ACGT", "+", "!!!!"),
    file.path(fastq_dir, "reads.fastq")
  )

  info <- fastq_read_info2(
    fastq_dir,
    file_format = ".fastq",
    plot_time = FALSE
  )
  expect_equal(info$barcode, "sample-A")
  expect_equal(info$tag_XY, "value")

  empty <- fastq_read_info2(fastq_dir, file_format = ".fq", plot_time = FALSE)
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("read", "source_file", "header_format", "start_time"))
})

test_that("fastq_read_info2 validates arguments and FASTQ records", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  expect_error(fastq_read_info2(tempfile(), plot_time = FALSE), "existing directory")
  expect_error(fastq_read_info2(fastq_dir, plot_time = NA), "plot_time")

  writeLines(
    c("@read1 ch:i:1", "ACGT", "+"),
    file.path(fastq_dir, "bad.fastq")
  )
  expect_error(
    fastq_read_info2(fastq_dir, file_format = ".fastq", plot_time = FALSE),
    "Incomplete FASTQ record"
  )
})
