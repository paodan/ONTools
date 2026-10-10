test_that("fastq_read_info2 vectorizes MinKNOW SAM-style headers", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  writeLines(c(
    paste0(
      "@read-1 qs:f:15.715393 mx:i:4 ch:i:1803 rn:i:438 ",
      "st:Z:2026-09-20T14:53:20.912722+08:00 ns:i:3694 ts:i:10 ",
      "sm:f:93.99999 sd:f:24 sv:Z:pa du:f:0.7388 dx:i:0 ",
      "RG:Z:run_model_barcode86 DT:Z:2026-09-20T14:27:18+08:00 ",
      "PU:Z:PBO11046 LB:Z:20260920pcr SM:Z:barcode86 al:Z:sample86"
    ),
    "ACGT", "+", "!!!!",
    paste0(
      "@read-2 ch:i:17 qs:f:12.5 rn:i:439 ",
      "st:Z:2026-09-20T14:53:21+08:00 al:Z:sample87"
    ),
    "TGCA", "+", "####"
  ), file.path(fastq_dir, "reads.fastq"))

  info <- fastq_read_info2(
    fastq_dir,
    file_format = ".fastq",
    plot_time = FALSE
  )

  expect_equal(info$read, c("read-1", "read-2"))
  expect_equal(info$header_format, rep("sam_tag", 2))
  expect_equal(info$mean_qscore, c(15.715393, 12.5))
  expect_equal(info$mux, c(4L, NA_integer_))
  expect_equal(info$ch, c(1803L, 17L))
  expect_equal(info$read_number, c(438L, 439L))
  expect_equal(info$signal_samples, c(3684L, NA_integer_))
  expect_equal(info$duration, c(0.7388, NA_real_))
  expect_equal(info$barcode, c("barcode86", "sample87"))
  expect_equal(info$barcode_alias, c("sample86", "sample87"))
  expect_equal(info$flow_cell_id, c("PBO11046", NA_character_))
  expect_s3_class(info$start_time, "POSIXct")
  expect_equal(
    info$start_time[[1]],
    as.POSIXct("2026-09-20 06:53:20.912722", tz = "UTC")
  )
})

test_that("fastq_read_info2 returns stable typed optional columns", {
  info <- parse_minknow_fastq_headers(
    c("read-1 ch:i:1", "read-2 ch:i:2"),
    "reads.fastq"
  )

  expect_type(info$ch, "integer")
  expect_type(info$mean_qscore, "double")
  expect_type(info$barcode, "character")
  expect_s3_class(info$start_time, "POSIXct")
  expect_true(all(is.na(info$mean_qscore)))
  expect_true(all(is.na(info$barcode)))
})

test_that("fastq_read_info2 rejects legacy headers", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  writeLines(
    c("@read1 ch=1 barcode=barcode01", "ACGT", "+", "!!!!"),
    file.path(fastq_dir, "reads.fastq")
  )
  expect_error(
    fastq_read_info2(fastq_dir, file_format = ".fastq", plot_time = FALSE),
    "SAM-style"
  )
})

test_that("fastq_read_info2 handles empty input and validates arguments", {
  fastq_dir <- tempfile("fastq-info2-")
  dir.create(fastq_dir)
  empty <- fastq_read_info2(fastq_dir, file_format = ".fq", plot_time = FALSE)
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("read", "source_file", "header_format", "start_time"))

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
