write_peak_test_fastq <- function(path, lengths, quality_character = "!") {
  records <- unlist(lapply(seq_along(lengths), function(i) {
    sequence <- paste(rep("A", lengths[[i]]), collapse = "")
    quality <- paste(rep(quality_character, lengths[[i]]), collapse = "")
    c(paste0("@read", i), sequence, "+", quality)
  }), use.names = FALSE)
  writeLines(records, path)
  invisible(records)
}

test_that("plot_fastq_peaks plots multiple files using chunked FASTQ reading", {
  fastq_dir <- tempfile("fastq-peaks-")
  dir.create(fastq_dir)
  write_peak_test_fastq(
    file.path(fastq_dir, "sample1.fastq"),
    c(95:105, 295:305),
    quality_character = "\""
  )
  write_peak_test_fastq(
    file.path(fastq_dir, "sample2.fastq"),
    c(145:155, 395:405)
  )

  plots <- plot_fastq_peaks(
    fastqs = c("sample1.fastq", "sample2.fastq"),
    path = fastq_dir,
    min_distance = 100,
    density_n = 1024,
    x_breaks = NULL,
    print_plots = FALSE
  )

  expect_named(plots, c("sample1.fastq", "sample2.fastq"))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_equal(plots[[1L]]$labels$subtitle, "sample1.fastq")
  expect_equal(plots[[2L]]$labels$subtitle, "sample2.fastq")
  expect_true(is.data.frame(attr(plots[[1L]], "peak_counts")))
})

test_that("plot_fastq_peaks passes titles and per-file subtitles", {
  fastq_dir <- tempfile("fastq-peaks-")
  dir.create(fastq_dir)
  files <- file.path(fastq_dir, c("a.fastq", "b.fastq"))
  write_peak_test_fastq(files[[1L]], c(95:105, 295:305))
  write_peak_test_fastq(files[[2L]], c(145:155, 395:405))

  plots <- plot_fastq_peaks(
    fastqs = files,
    title = "Amplicon lengths",
    subtitle = c("Sample A", "Sample B"),
    min_distance = 100,
    x_breaks = NULL,
    print_plots = FALSE
  )

  expect_equal(plots[[1L]]$labels$title, "Amplicon lengths")
  expect_equal(plots[[2L]]$labels$title, "Amplicon lengths")
  expect_equal(plots[[1L]]$labels$subtitle, "Sample A")
  expect_equal(plots[[2L]]$labels$subtitle, "Sample B")
})

test_that("plot_fastq_peaks warns and skips inputs without enough lengths", {
  fastq_dir <- tempfile("fastq-peaks-")
  dir.create(fastq_dir)
  empty <- file.path(fastq_dir, "empty.fastq")
  short <- file.path(fastq_dir, "short.fastq")
  file.create(empty)
  write_peak_test_fastq(short, c(100L, 100L, 100L))

  expect_warning(
    empty_plot <- plot_fastq_peaks(empty, print_plots = FALSE),
    "Skipping `empty[.]fastq`"
  )
  expect_warning(
    short_plot <- plot_fastq_peaks(short, print_plots = FALSE),
    "Skipping `short[.]fastq`"
  )

  expect_length(empty_plot, 1L)
  expect_length(short_plot, 1L)
  expect_null(empty_plot[[1L]])
  expect_null(short_plot[[1L]])
})

test_that("plot_fastq_peaks validates paths, subtitles, and FASTQ records", {
  fastq <- tempfile(fileext = ".fastq")
  write_peak_test_fastq(fastq, c(10:12))

  expect_error(
    plot_fastq_peaks(fastq, subtitle = c("a", "b"), print_plots = FALSE),
    "subtitle"
  )
  expect_error(
    plot_fastq_peaks("missing.fastq", print_plots = FALSE),
    "does not exist"
  )

  malformed <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!"), malformed)
  expect_error(
    plot_fastq_peaks(malformed, print_plots = FALSE),
    "lengths differ"
  )
})
