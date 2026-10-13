write_peaks_info_fastq <- function(path, lengths) {
  records <- unlist(lapply(seq_along(lengths), function(i) {
    sequence <- paste(rep("A", lengths[[i]]), collapse = "")
    quality <- paste(rep("!", lengths[[i]]), collapse = "")
    c(paste0("@read", i), sequence, "+", quality)
  }), use.names = FALSE)
  writeLines(records, path)
  invisible(path)
}

test_that("get_peaks_info combines and annotates per-file peaks", {
  fastq_dir <- tempfile("peak-info-")
  dir.create(fastq_dir)
  write_peaks_info_fastq(
    file.path(fastq_dir, "sample1.fastq"),
    c(95:105, 295:305)
  )
  write_peaks_info_fastq(
    file.path(fastq_dir, "sample2.fastq"),
    c(145:155, 395:405)
  )

  peaks <- get_peaks_info(
    fastqs = c("sample1.fastq", "sample2.fastq"),
    path = fastq_dir,
    min_distance = 100,
    density_n = 1024
  )

  expect_named(
    peaks,
    c(
      "file", "peak", "major", "position", "read_length", "density",
      "window_start", "window_end", "read_count", "read_percent"
    )
  )
  expect_setequal(unique(peaks$file), c("sample1.fastq", "sample2.fastq"))
  expect_true(all(peaks$peak >= 2L))
  expect_equal(
    as.vector(tapply(peaks$major, peaks$file, sum)),
    c(1L, 1L)
  )
  expect_true(all(vapply(
    split(peaks$read_length, peaks$file),
    function(x) identical(x, sort(x, decreasing = TRUE)),
    logical(1)
  )))
})

test_that("get_peaks_info preserves detection-scale and base-pair positions", {
  fastq <- tempfile(fileext = ".fastq")
  write_peaks_info_fastq(fastq, c(95:105, 295:305))

  peaks <- get_peaks_info(
    fastq,
    log_scale = TRUE,
    min_distance = 100,
    density_n = 1024
  )

  expect_equal(10^peaks$position, peaks$read_length, tolerance = 1e-8)
})

test_that("get_peaks_info returns a typed empty data frame", {
  fastq <- tempfile(fileext = ".fastq")
  file.create(fastq)

  expect_warning(
    peaks <- get_peaks_info(fastq),
    "Skipping"
  )
  expect_s3_class(peaks, "data.frame")
  expect_equal(nrow(peaks), 0L)
  expect_type(peaks$peak, "integer")
  expect_type(peaks$major, "logical")
})

test_that("plot_peaks_statistics returns a three-panel patchwork", {
  fastq <- tempfile(fileext = ".fastq")
  write_peaks_info_fastq(fastq, c(95:105, 295:305))

  plots <- plot_peaks_statistics(
    fastq,
    min_distance = 100,
    density_n = 1024,
    title_hjust = 0
  )

  expect_s3_class(plots, "patchwork")
  expect_s3_class(plots, "ggplot")
  expect_length(plots$patches$plots, 2L)
})

test_that("plot_peaks_statistics errors clearly when no peaks are available", {
  fastq <- tempfile(fileext = ".fastq")
  file.create(fastq)

  expect_warning(
    expect_error(
      plot_peaks_statistics(fastq),
      "No read-length peaks"
    ),
    "Skipping"
  )
})
