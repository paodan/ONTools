test_that("plot_seqQC returns plots and recovery without writing files", {
  summary_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(
      alias = c("barcode01", "barcode01", "barcode02", "unclassified"),
      sequence_length_template = c(1000, 1500, 800, 500)
    ),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  qc <- plot_seqQC(summary_file, runName = "run1", device = NULL, barcodes = 1:2)

  expect_s3_class(qc$numRead, "ggplot")
  expect_s3_class(qc$lenRead, "ggplot")
  expect_equal(qc$recovery, 75)
  expect_equal(qc$files, character())
  expect_equal(qc$n_reads_raw, 4L)
  expect_equal(qc$n_reads_filtered, 4L)
  expect_null(qc$filters$min_read_length)
  expect_null(qc$filters$max_read_length)
  expect_equal(qc$filename_suffix, "")
  expect_equal(qc$read_counts$sample, c("barcode01", "barcode02", "unclassified"))
  expect_equal(qc$read_counts$Freq, c(2L, 1L, 1L))
})

test_that("plot_seqQC filters reads by length before QC summaries", {
  summary_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(
      alias = c("barcode01", "barcode01", "barcode02", "unclassified"),
      sequence_length_template = c(1000, 1200, 1800, 2000)
    ),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  qc <- plot_seqQC(
    summary_file,
    runName = "run-filtered",
    device = NULL,
    barcodes = 1:2,
    min_read_length = 1200,
    max_read_length = 1800
  )

  expect_equal(qc$n_reads_raw, 4L)
  expect_equal(qc$n_reads_filtered, 2L)
  expect_equal(qc$filters$min_read_length, 1200)
  expect_equal(qc$filters$max_read_length, 1800)
  expect_equal(qc$filename_suffix, "__len1200-1800")
  expect_equal(qc$read_counts$sample, c("barcode01", "barcode02", "unclassified"))
  expect_equal(qc$read_counts$Freq, c(1L, 1L, 0L))
  expect_equal(qc$recovery, 100)
})

test_that("plot_seqQC adds filter information to output filenames", {
  summary_file <- tempfile(fileext = ".txt")
  out_dir <- tempfile("seqqc-figs-")
  write.table(
    data.frame(
      alias = c("barcode01", "barcode01"),
      sequence_length_template = c(1200, 1800)
    ),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  qc <- plot_seqQC(
    summary_file,
    runName = "run-filtered",
    device = "png",
    barcodes = 1,
    out_dir = out_dir,
    min_read_length = 1200,
    max_read_length = 1800
  )

  expect_true(all(grepl("__len1200-1800\\.png$", qc$files)))
  expect_true(all(file.exists(qc$files)))
})

test_that("plot_seqQC validates length filters", {
  summary_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(
      alias = "barcode01",
      sequence_length_template = 1000
    ),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    plot_seqQC(summary_file, device = NULL, min_read_length = -1),
    "min_read_length"
  )
  expect_error(
    plot_seqQC(
      summary_file,
      device = NULL,
      min_read_length = 1800,
      max_read_length = 1200
    ),
    "less than or equal"
  )
})

test_that("plot_seqQC validates required columns", {
  summary_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(alias = "barcode01"),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    plot_seqQC(summary_file, device = NULL),
    "missing required column"
  )
})

test_that("plot_seqQC can write plot files to a selected directory", {
  skip_if_not(
    capabilities("png"),
    "PNG graphics device is not available on this platform."
  )

  summary_file <- tempfile(fileext = ".txt")
  out_dir <- tempfile("seqqc-figs-")
  write.table(
    data.frame(
      alias = c("barcode01", "unclassified"),
      sequence_length_template = c(1000, 500)
    ),
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  qc <- plot_seqQC(
    summary_file,
    runName = "run2",
    device = "png",
    barcodes = 1,
    out_dir = out_dir
  )

  expect_true(all(file.exists(qc$files)))
})
