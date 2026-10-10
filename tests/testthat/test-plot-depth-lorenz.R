test_that("gini_coefficient describes equality and inequality", {
  expect_equal(gini_coefficient(rep(10, 5)), 0)
  expect_equal(gini_coefficient(c(0, 0, 0, 10)), 1)
  expect_equal(gini_coefficient(1:5), gini_coefficient((1:5) * 10))
  expect_true(is.na(gini_coefficient(0)))
  expect_error(gini_coefficient(c(1, NA)), "NA")
  expect_error(gini_coefficient(c(1, -1)), "non-negative")
})

test_that("summarize_depth_uniformity returns documented metrics", {
  result <- summarize_depth_uniformity(c(0, 10, 10, 20, NA))

  expect_s3_class(result, "data.frame")
  expect_equal(result$samples, 4)
  expect_equal(result$total_depth, 40)
  expect_equal(result$median_depth, 10)
  expect_equal(result$zero_depth_percent, 25)
  expect_named(
    result,
    c(
      "samples", "total_depth", "mean_depth", "median_depth",
      "minimum_depth", "maximum_depth", "cv", "robust_cv",
      "iqr_over_median", "gini", "p95_p5", "zero_depth_percent",
      "within_half_to_two_fold_percent", "within_20_percent"
    )
  )
})

test_that("plot_depth_lorenz returns plot data and summary", {
  plot <- plot_depth_lorenz(c(0, 10, 20, 30), title = "Test")

  expect_s3_class(plot, "ggplot")
  expect_s3_class(attr(plot, "depth_summary"), "data.frame")
  expect_equal(
    attr(plot, "lorenz_data")$sample_fraction,
    seq(0, 1, by = 0.25)
  )
  expect_equal(tail(attr(plot, "lorenz_data")$depth_fraction, 1), 1)

  expect_error(plot_depth_lorenz(c(0, 0)), "greater than zero")
  expect_error(plot_depth_lorenz(c(1, -1)), "negative")
})
