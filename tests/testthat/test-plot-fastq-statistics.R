write_statistics_fastq <- function(path, lengths, with_barcode = TRUE) {
  records <- unlist(lapply(seq_along(lengths), function(i) {
    sequence <- paste(rep("A", lengths[[i]]), collapse = "")
    quality <- paste(rep("!", lengths[[i]]), collapse = "")
    barcode <- if (with_barcode) " SM:Z:barcode01" else ""
    c(
      paste0(
        "@read", i, " qs:f:15.5 ch:i:", (i %% 4L) + 1L,
        barcode
      ),
      sequence,
      "+",
      quality
    )
  }), use.names = FALSE)
  writeLines(records, path)
  invisible(path)
}

test_that("plot_fastq_statistics returns an integrated patchwork", {
  fastq_dir <- tempfile("fastq-statistics-")
  dir.create(fastq_dir)
  write_statistics_fastq(
    file.path(fastq_dir, "sample.fastq"),
    c(95:105, 295:305)
  )

  plot <- plot_fastq_statistics(
    fastqs = "sample.fastq",
    path = fastq_dir,
    title = "Example run",
    min_distance = 100,
    density_n = 1024,
    title_hjust = 0
  )

  expect_s3_class(plot, "patchwork")
  expect_s3_class(plot, "ggplot")
  expect_identical(attr(plot, "header_format"), "sam_tag")
  expect_identical(attr(plot, "parser"), "fastq_read_info2")
  expect_identical(attr(plot, "barcode_column"), "barcode")
  expect_s3_class(attr(plot, "read_distribution"), "patchwork")
  expect_s3_class(attr(plot, "peak_statistics"), "patchwork")
  expect_silent(ggplot2::ggplotGrob(plot))
})

test_that("plot_fastq_statistics falls back to source_file grouping", {
  fastq <- tempfile(fileext = ".fastq")
  write_statistics_fastq(fastq, c(95:105, 295:305), with_barcode = FALSE)

  plot <- plot_fastq_statistics(
    fastqs = fastq,
    min_distance = 100,
    density_n = 1024
  )

  expect_identical(attr(plot, "barcode_column"), "source_file")
})

test_that("plot_fastq_statistics validates layout and grouping arguments", {
  fastq <- tempfile(fileext = ".fastq")
  write_statistics_fastq(fastq, c(95:105, 295:305))

  expect_error(
    plot_fastq_statistics(fastq, heights = 1),
    "heights"
  )
  expect_error(
    plot_fastq_statistics(fastq, heights = c(1, 0)),
    "heights"
  )
  expect_error(
    plot_fastq_statistics(fastq, title = NA_character_),
    "title"
  )
  expect_error(
    plot_fastq_statistics(fastq, barcode_column = "missing"),
    "barcode_column"
  )
})

test_that("plot_peaks_statistics validates and supports uncombined output", {
  fastq <- tempfile(fileext = ".fastq")
  write_statistics_fastq(fastq, c(95:105, 295:305))

  plots <- plot_peaks_statistics(
    fastqs = fastq,
    min_distance = 100,
    density_n = 1024,
    combine_plot = FALSE
  )
  expect_named(
    plots,
    c("peak_num_hist", "peak_pos_hist", "peak_pos_hist_major")
  )
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_error(
    plot_peaks_statistics(fastqs = fastq, combine_plot = NA),
    "combine_plot"
  )
})
