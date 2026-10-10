test_that("plot_fastq_read_distribution returns a composed QC plot", {
  info <- data.frame(
    barcode = c("barcode01", "barcode01", "barcode02", NA),
    ch = c(1, 1, 2, 3),
    stringsAsFactors = FALSE
  )

  plot <- plot_fastq_read_distribution(info, title = "Run QC")

  expect_s3_class(plot, "patchwork")
  expect_equal(
    attr(plot, "sample_depth"),
    data.frame(
      sample = c("barcode01", "barcode02", "unclassified"),
      depth = c(2L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(nrow(attr(plot, "quantile_data")), 202L)
  expect_equal(unique(attr(plot, "quantile_data")$percentile), 0:100)
})

test_that("plot_fastq_read_distribution supports another barcode column", {
  info <- data.frame(sample_id = c("A", "A", "B"))

  plot <- plot_fastq_read_distribution(
    info,
    barcode_column = "sample_id"
  )

  expect_s3_class(plot, "patchwork")
  expect_equal(attr(plot, "sample_depth")$depth, c(2L, 1L))
})

test_that("plot_fastq_read_distribution validates its inputs", {
  expect_error(plot_fastq_read_distribution(1:3), "data frame")
  expect_error(
    plot_fastq_read_distribution(data.frame(barcode = character())),
    "at least one read"
  )
  expect_error(
    plot_fastq_read_distribution(data.frame(sample = "A")),
    "was not found"
  )
  expect_error(
    plot_fastq_read_distribution(data.frame(barcode = "A"), title = 1),
    "title"
  )
})
