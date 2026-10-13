write_auto_fastq <- function(path, header, sequence = "ACGT", quality = "!!!!") {
  writeLines(c(header, sequence, "+", quality), path)
  invisible(path)
}

test_that("fastq_read_info_auto selects the legacy parser", {
  fastq <- tempfile(fileext = ".fastq")
  write_auto_fastq(
    fastq,
    "@legacy-read ch=42 start_time=2025-01-01T00:00:00Z barcode=barcode01"
  )

  info <- fastq_read_info_auto(fastqs = fastq, plot_time = FALSE)

  expect_equal(info$read, "legacy-read")
  expect_equal(info$ch, "42")
  expect_equal(info$barcode, "barcode01")
  expect_identical(attr(info, "header_format"), "legacy_key_value")
  expect_identical(attr(info, "parser"), "fastq_read_info")
})

test_that("fastq_read_info_auto selects the MinKNOW/Dorado parser", {
  fastq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fastq, open = "wt")
  writeLines(
    c(
      "@dorado-read qs:f:15.7 ch:i:42 SM:Z:barcode01",
      "ACGT", "+", "!!!!"
    ),
    con
  )
  close(con)

  info <- fastq_read_info_auto(fastqs = fastq, plot_time = FALSE)

  expect_equal(info$read, "dorado-read")
  expect_equal(info$mean_qscore, 15.7)
  expect_equal(info$ch, 42L)
  expect_equal(info$barcode, "barcode01")
  expect_identical(attr(info, "header_format"), "sam_tag")
  expect_identical(attr(info, "parser"), "fastq_read_info2")
})

test_that("fastq_read_info_auto detects from the first selected file", {
  fastq_dir <- tempfile("fastq-auto-")
  dir.create(fastq_dir)
  write_auto_fastq(
    file.path(fastq_dir, "a.fastq"),
    "@read-a ch:i:1 SM:Z:a"
  )
  write_auto_fastq(
    file.path(fastq_dir, "b.fastq"),
    "@read-b ch:i:2 SM:Z:b"
  )

  expect_message(
    info <- fastq_read_info_auto(
      path = fastq_dir,
      fastqs = c("b.fastq", "a.fastq"),
      plot_time = FALSE,
      verbose = TRUE
    ),
    "using `fastq_read_info2"
  )
  expect_equal(info$read, c("read-b", "read-a"))
})

test_that("fastq_read_info_auto rejects unclassifiable or ambiguous headers", {
  unknown <- tempfile(fileext = ".fastq")
  write_auto_fastq(unknown, "@read-without-metadata")
  expect_error(
    fastq_read_info_auto(fastqs = unknown, plot_time = FALSE),
    "neither recognizable"
  )

  ambiguous <- tempfile(fileext = ".fastq")
  write_auto_fastq(ambiguous, "@read ch=1 ch:i:1")
  expect_error(
    fastq_read_info_auto(fastqs = ambiguous, plot_time = FALSE),
    "ambiguous"
  )
})

test_that("fastq_read_info_auto validates the first record and arguments", {
  empty <- tempfile(fileext = ".fastq")
  file.create(empty)
  expect_error(
    fastq_read_info_auto(fastqs = empty, plot_time = FALSE),
    "is empty"
  )

  malformed <- tempfile(fileext = ".fastq")
  writeLines(c("@read ch=1", "ACGT", "+", "!!!"), malformed)
  expect_error(
    fastq_read_info_auto(fastqs = malformed, plot_time = FALSE),
    "lengths differ"
  )

  expect_error(fastq_read_info_auto(), "No FASTQ files|Supply `path`")
  expect_error(
    fastq_read_info_auto(fastqs = malformed, verbose = NA),
    "verbose"
  )
  expect_error(
    fastq_read_info_auto(fastqs = malformed, barcode_column = NA_character_),
    "barcode_column"
  )
})
