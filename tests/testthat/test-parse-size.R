test_that("parse_size parses base-pair size strings", {
  expect_equal(parse_size("1000"), 1000)
  expect_equal(parse_size("1000 bp"), 1000)
  expect_equal(parse_size("180k"), 180000)
  expect_equal(parse_size("180kbp"), 180000)
  expect_equal(parse_size("2.5Mbp"), 2500000)
  expect_equal(parse_size("0.5Gb"), 5e8)
  expect_equal(parse_size(c("1bp", "2kb", "3Mb")), c(1, 2000, 3e6))
})

test_that("parseSize remains available as a compatibility alias", {
  expect_equal(parseSize(c("180kbp", "2Mbp")), c(180000, 2e6))
})

test_that("parse_size supports custom units", {
  expect_equal(parse_size("3foo", units = c(bp = 1, foo = 7)), 21)
})

test_that("parse_size validates inputs", {
  expect_error(parse_size(1000), "`x`")
  expect_error(parse_size("abc"), "Invalid size string")
  expect_error(parse_size("10tb"), "Unknown size unit")
  expect_error(parse_size("10kb", units = c(kb = NA_real_)), "`units`")
  expect_error(parse_size("10kb", units = c(1000)), "`units`")
})
