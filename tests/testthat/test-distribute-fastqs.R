make_distribute_fastq_files <- function(names) {
  root <- tempfile("distribute-fastqs-source-")
  dir.create(root)
  paths <- file.path(root, names)
  invisible(file.create(paths))
  paths
}

test_that("distribute_fastqs separates peak 1 and bins remaining files", {
  fastqs <- make_distribute_fastq_files(c(
    "sample01_peak_1_1600bp.fastq.gz",
    "sample01_peak_2_720bp.fastq.gz",
    "sample02_peak_3_1300bp.fastq.gz"
  ))
  output <- tempfile("distribute-fastqs-output-")

  result <- distribute_fastqs(
    fastqs,
    out_path = output,
    save_diff_size_folder = TRUE
  )

  expect_equal(
    result$group,
    c(
      "peak1_target_range",
      "read_len_range__500_1000",
      "read_len_range__1000_2000"
    )
  )
  expect_true(all(result$copied))
  expect_true(all(file.exists(result$destination)))
  expect_match(
    result$destination[[2L]],
    file.path(
      "read_len_range__500_1000",
      "sample01_peak_2_720bp",
      "sample01_peak_2_720bp.fastq.gz"
    ),
    fixed = TRUE
  )
})

test_that("distribute_fastqs can bin peak 1 by length", {
  fastq <- make_distribute_fastq_files(
    "sample01_peak_1_1600bp.fastq.gz"
  )
  result <- distribute_fastqs(
    fastq,
    out_path = tempfile("distribute-fastqs-output-"),
    separate_peak1 = FALSE
  )

  expect_equal(result$group, "read_len_range__1000_2000")
  expect_true(result$copied)
})

test_that("distribute_fastqs includes the lowest boundary", {
  fastqs <- make_distribute_fastq_files(c(
    "sample01_peak_2_500bp.fastq.gz",
    "sample02_peak_2_1000bp.fastq.gz",
    "sample03_peak_2_2000bp.fastq.gz"
  ))
  result <- distribute_fastqs(
    fastqs,
    out_path = tempfile("distribute-fastqs-output-")
  )

  expect_equal(
    result$group,
    c(
      "read_len_range__500_1000",
      "read_len_range__500_1000",
      "read_len_range__1000_2000"
    )
  )
})

test_that("distribute_fastqs reports out-of-range files", {
  fastq <- make_distribute_fastq_files(
    "sample01_peak_2_2500bp.fastq.gz"
  )

  expect_warning(
    result <- distribute_fastqs(
      fastq,
      out_path = tempfile("distribute-fastqs-output-")
    ),
    "outside"
  )
  expect_equal(result$status, "out_of_range")
  expect_false(result$copied)
  expect_true(is.na(result$destination))

  expect_error(
    distribute_fastqs(
      fastq,
      out_path = tempfile("distribute-fastqs-output-"),
      out_of_range = "error"
    ),
    "outside"
  )
})

test_that("distribute_fastqs validates filenames and arguments", {
  malformed <- make_distribute_fastq_files("sample01.fastq.gz")
  expect_error(
    distribute_fastqs(
      malformed,
      out_path = tempfile("distribute-fastqs-output-")
    ),
    "do not match"
  )
  expect_error(
    distribute_fastqs(
      malformed,
      out_path = tempfile("distribute-fastqs-output-"),
      size_breaks = c(1000, 500)
    ),
    "strictly increasing"
  )
  expect_error(
    distribute_fastqs(
      c(malformed, malformed),
      out_path = tempfile("distribute-fastqs-output-")
    ),
    "duplicated"
  )
})

test_that("distribute_fastqs honors overwrite and preserves timestamps", {
  fastq <- make_distribute_fastq_files(
    "sample01_peak_2_720bp.fastq.gz"
  )
  expected_time <- as.POSIXct("2020-01-02 03:04:05", tz = "UTC")
  Sys.setFileTime(fastq, expected_time)
  output <- tempfile("distribute-fastqs-output-")

  first <- distribute_fastqs(fastq, out_path = output)
  expect_equal(first$status, "copied")
  expect_equal(
    as.numeric(file.info(first$destination)$mtime),
    as.numeric(file.info(fastq)$mtime)
  )

  second <- distribute_fastqs(
    fastq,
    out_path = output,
    overwrite = FALSE
  )
  expect_equal(second$status, "destination_exists")
  expect_false(second$copied)
})
