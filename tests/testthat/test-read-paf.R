test_that("read_paf reads standard PAF columns", {
  paf <- tempfile(fileext = ".paf")
  writeLines(c(
    "read1\t1000\t0\t900\t+\tchr1\t2000\t100\t1000\t890\t900\t60",
    "read2\t500\t10\t450\t-\tchr2\t1500\t20\t460\t420\t440\t20"
  ), paf)

  res <- read_paf(paf, parse_tags = FALSE)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(
    names(res)[1:12],
    c(
      "query_name", "query_length", "query_start", "query_end", "strand",
      "target_name", "target_length", "target_start", "target_end",
      "residue_matches", "alignment_block_length", "mapping_quality"
    )
  )
  expect_equal(res$query_name, c("read1", "read2"))
  expect_equal(res$query_length, c(1000, 500))
  expect_equal(res$strand, c("+", "-"))
  expect_equal(res$mapping_quality, c(60, 20))
})

test_that("read_paf parses optional tags", {
  paf <- tempfile(fileext = ".paf")
  writeLines(c(
    "read1\t1000\t0\t900\t+\tchr1\t2000\t100\t1000\t890\t900\t60\tNM:i:10\tcg:Z:900M",
    "read2\t500\t10\t450\t-\tchr2\t1500\t20\t460\t420\t440\t20\tNM:i:3\tdv:f:0.01"
  ), paf)

  res <- read_paf(paf)

  expect_equal(res$NM, c(10, 3))
  expect_equal(res$cg, c("900M", NA))
  expect_equal(res$dv, c(NA, 0.01))
  expect_false(any(grepl("^tag_", names(res))))
})

test_that("read_paf can keep raw optional tag columns", {
  paf <- tempfile(fileext = ".paf")
  writeLines(
    "read1\t1000\t0\t900\t+\tchr1\t2000\t100\t1000\t890\t900\t60\tNM:i:10",
    paf
  )

  res <- read_paf(paf, parse_tags = FALSE)

  expect_equal(res$tag_1, "NM:i:10")
})

test_that("read_paf reads gzipped PAF files", {
  paf <- tempfile(fileext = ".paf.gz")
  con <- gzfile(paf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(
    "read1\t1000\t0\t900\t+\tchr1\t2000\t100\t1000\t890\t900\t60",
    con
  )
  close(con)
  on.exit(NULL)

  res <- read_paf(paf)

  expect_equal(nrow(res), 1)
  expect_equal(res$query_name, "read1")
  expect_equal(res$mapping_quality, 60)
})

test_that("read_paf handles empty files and validates inputs", {
  paf <- tempfile(fileext = ".paf")
  file.create(paf)

  expect_null(read_paf(paf))
  expect_error(read_paf("missing.paf"), "paf_file")

  bad <- tempfile(fileext = ".paf")
  writeLines("read1\t1000", bad)
  expect_error(read_paf(bad), "at least 12")
})
