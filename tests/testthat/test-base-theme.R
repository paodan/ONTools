test_that("base_theme returns the documented ggplot2 theme", {
  theme <- base_theme()

  expect_s3_class(theme, "theme")
  expect_equal(theme$panel.grid.major$colour, "grey85")
  expect_equal(theme$panel.grid.major$linewidth, 0.2)
  expect_s3_class(theme$panel.grid.minor, "element_blank")
  expect_equal(theme$plot.background$fill, "transparent")
  expect_true(is.na(theme$plot.background$colour))
})

test_that("base_theme forwards base font settings", {
  theme <- base_theme(base_size = 14, base_family = "sans")

  expect_equal(theme$text$size, 14)
  expect_equal(theme$text$family, "sans")
})

test_that("base_theme can be added to a ggplot", {
  plot <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 3:1),
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point() +
    base_theme()

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$theme$panel.grid.major$colour, "grey85")
})

test_that("base_theme validates font arguments", {
  expect_error(base_theme(base_size = 0), "base_size")
  expect_error(base_theme(base_size = c(10, 12)), "base_size")
  expect_error(base_theme(base_family = NA_character_), "base_family")
})
