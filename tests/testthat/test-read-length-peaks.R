test_that("detect_length_peaks finds separated simulated modes", {
  set.seed(101)
  lengths <- c(
    rnorm(500, 250, 15),
    rnorm(1200, 1500, 35),
    100000, 105000,
    NA, -1
  )

  result <- detect_length_peaks(
    lengths,
    log_scale = TRUE,
    adjust = 0.8,
    min_relative_height = 0.05,
    min_distance = 300,
    n = 2048,
    detection_range = c(100, 2500)
  )

  expect_equal(nrow(result$peaks), 2L)
  expect_equal(result$peaks$read_length, c(250, 1500), tolerance = 30)
  expect_equal(result$n_reads_total, 1702L)
  expect_equal(result$n_reads_detection, 1700L)
  expect_equal(result$detection_range, c(100, 2500))
})

test_that("detection_range excludes modes outside the search interval", {
  set.seed(102)
  lengths <- c(rnorm(500, 200, 12), rnorm(800, 1500, 30))

  result <- detect_length_peaks(
    lengths,
    adjust = 0.8,
    min_distance = 100,
    detection_range = c(1000, 2000)
  )

  expect_equal(nrow(result$peaks), 1L)
  expect_equal(result$peaks$read_length, 1500, tolerance = 30)
  expect_equal(result$n_reads_detection, 800L)
})

test_that("plot_read_length_peaks returns plot and attached data", {
  set.seed(103)
  lengths <- c(rnorm(400, 250, 15), rnorm(1000, 1500, 35), 100000)

  plot <- plot_read_length_peaks(
    lengths,
    log_scale = TRUE,
    adjust = 0.8,
    min_relative_height = 0.05,
    min_distance = 300,
    density_n = 2048,
    detection_range = c(100, 2500),
    left_offset = 50,
    right_offset = 100,
    binwidth = 25,
    x_limits = c(100, 2500)
  )

  expect_s3_class(plot, "ggplot")
  expect_true(is.list(attr(plot, "peak_result")))
  expect_named(
    attr(plot, "peak_counts"),
    c("read_length", "window_start", "window_end", "read_count", "read_percent")
  )
  expect_true(all(attr(plot, "peak_counts")$read_count > 0L))
  expect_equal(sum(attr(plot, "histogram_data")$count), 1401L)
  expect_equal(attr(plot, "parameters")$peak_window$left_offset, 50)
})

test_that("plot_read_length_peaks supports linear axes and custom breaks", {
  set.seed(104)
  lengths <- c(rnorm(300, 300, 10), rnorm(500, 1500, 25))
  plot <- plot_read_length_peaks(
    lengths,
    log_scale = FALSE,
    x_scale = "linear",
    x_breaks = seq(0, 2000, 500),
    x_limits = c(0, 2000),
    show_peak_labels = FALSE,
    show_peak_windows = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(attr(plot, "parameters")$display$x_scale, "linear")
  expect_equal(attr(plot, "parameters")$display$x_limits, c(0, 2000))
})

test_that("read-length peak functions validate inputs", {
  expect_error(detect_length_peaks("100"), "read_length")
  expect_error(detect_length_peaks(c(100, 100, 100)), "variation")
  expect_error(detect_length_peaks(1:10, adjust = 0), "adjust")
  expect_error(detect_length_peaks(1:10, min_relative_height = 2), "between 0 and 1")
  expect_error(detect_length_peaks(1:10, n = 2), "at least 3")
  expect_error(detect_length_peaks(1:10, detection_range = c(20, 30)),
               "Fewer than three")
  expect_error(plot_read_length_peaks(1:10, binwidth = 0), "binwidth")
  expect_error(plot_read_length_peaks(1:10, x_scale = "log10", x_breaks = 0:2),
               "positive")
  expect_error(plot_read_length_peaks(1:10, x_limits = c(10, 1)), "x_limits")
})
