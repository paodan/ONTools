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
  expect_equal(qc$read_counts$sample, c("barcode01", "barcode02", "unclassified"))
  expect_equal(qc$read_counts$Freq, c(2L, 1L, 1L))
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
