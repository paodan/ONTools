test_that("split_sample_info splits sample info and fills read length ranges", {
  sample_info <- data.frame(
    Sample_ID = c("S1", "S2", "S3"),
    Project_ID = c("P1", "P1", "P2"),
    Barcode_ID = c("NB01", "NB02", "NB03"),
    Expected_Size_bp = c("1000bp", "1kb", "bad"),
    stringsAsFactors = FALSE
  )

  res <- split_sample_info(sample_info)

  expect_named(res, c("data", "files"))
  expect_named(res$data, c("P1_1000bp", "P1_1kb", "P2_bad"))
  expect_equal(res$data[["P1_1000bp"]]$Min_Read_Length, 700)
  expect_equal(res$data[["P1_1000bp"]]$Max_Read_Length, 1200)
  expect_equal(res$data[["P1_1kb"]]$Min_Read_Length, 700)
  expect_equal(res$data[["P1_1kb"]]$Max_Read_Length, 1200)
  expect_true(is.na(res$data[["P2_bad"]]$Min_Read_Length))
  expect_true(is.na(res$data[["P2_bad"]]$Max_Read_Length))
  expect_equal(res$files, list())
})

test_that("split_sample_info writes one CSV per folder when output_dir is supplied", {
  sample_info <- data.frame(
    Sample_ID = c("S1", "S2"),
    Project_ID = c("P1", "P1"),
    Barcode_ID = c("NB01", "NB02"),
    Expected_Size_bp = c("1000bp", "1500bp"),
    stringsAsFactors = FALSE
  )
  output_dir <- tempfile("split-sample-info-")

  res <- split_sample_info(sample_info, output_dir = output_dir)

  expect_named(res$files, c("P1_1000bp", "P1_1500bp"))
  expect_true(file.exists(res$files[["P1_1000bp"]]))
  expect_true(file.exists(res$files[["P1_1500bp"]]))
  written <- utils::read.csv(res$files[["P1_1500bp"]])
  expect_equal(written$Min_Read_Length, 1200)
  expect_equal(written$Max_Read_Length, 1700)
})

test_that("split_sample_info supports custom absolute read length deltas", {
  sample_info <- data.frame(
    Sample_ID = "S1",
    Project_ID = "P1",
    Barcode_ID = "NB01",
    Expected_Size_bp = "1000bp",
    stringsAsFactors = FALSE
  )

  res <- split_sample_info(
    sample_info,
    min_read_delta = 100,
    max_read_delta = 50
  )

  expect_equal(res$data[["P1_1000bp"]]$Min_Read_Length, 900)
  expect_equal(res$data[["P1_1000bp"]]$Max_Read_Length, 1050)
})
