make_subset_test_bam <- function(directory) {
  sam <- file.path(directory, "input.sam")
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      "@SQ\tSN:refA\tLN:100",
      "@SQ\tSN:refB\tLN:100",
      paste(
        "r1", 0, "refA", 1, 30, "20M", "*", 0, 0,
        "ACGTACGTACGTACGTACGT", "FFFFFFFFFFFFFFFFFFFF",
        sep = "\t"
      ),
      paste(
        "r2", 0, "refA", 5, 20, "20M", "*", 0, 0,
        "ACGTACGTACGTACGTACGT", "FFFFFFFFFFFFFFFFFFFF",
        sep = "\t"
      ),
      paste(
        "r3", 0, "refA", 8, 40, "20M", "*", 0, 0,
        "ACGTACGTACGTACGTACGT", "FFFFFFFFFFFFFFFFFFFF",
        sep = "\t"
      ),
      paste(
        "r4", 0, "refA", 1, 10, "20M", "*", 0, 0,
        "ACGTACGTACGTACGTACGT", "FFFFFFFFFFFFFFFFFFFF",
        sep = "\t"
      ),
      paste(
        "r5", 0, "refB", 1, 50, "20M", "*", 0, 0,
        "ACGTACGTACGTACGTACGT", "FFFFFFFFFFFFFFFFFFFF",
        sep = "\t"
      )
    ),
    sam
  )
  Rsamtools::asBam(
    sam,
    destination = file.path(directory, "input"),
    overwrite = TRUE,
    indexDestination = FALSE
  )
}

test_that("subset_bam_reads works without an input index and limits records", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "sampled.bam")

  set.seed(99)
  random_seed <- .Random.seed
  result <- subset_bam_reads(
    bam = bam,
    output_bam = output,
    reference_name = "refA",
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 2,
    seed = 2026
  )

  expect_false(file.exists(paste0(bam, ".bai")))
  expect_equal(result$eligible_reads, 3L)
  expect_equal(result$selected_reads, 2L)
  expect_true(result$sampled)
  expect_equal(.Random.seed, random_seed)
  expect_equal(Rsamtools::countBam(output)$records, 2L)
  records <- Rsamtools::scanBam(
    output,
    param = Rsamtools::ScanBamParam(
      what = c("rname", "pos", "cigar", "mapq")
    )
  )[[1L]]
  expect_true(all(as.character(records$rname) == "refA"))
  expect_true(all(records$mapq >= 20L))
})

test_that("subset_bam_reads retains all eligible records below the limit", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "all.bam")

  result <- subset_bam_reads(
    bam,
    output,
    reference_name = "refA",
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 2000
  )

  expect_equal(result$eligible_reads, 3L)
  expect_equal(result$selected_reads, 3L)
  expect_false(result$sampled)
  expect_equal(Rsamtools::countBam(output)$records, 3L)
})

test_that("subset_bam_reads writes an empty valid BAM when nothing qualifies", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "empty.bam")

  result <- subset_bam_reads(
    bam,
    output,
    reference_name = "missing-reference",
    region_start = 1,
    region_end = 10
  )

  expect_equal(result$eligible_reads, 0L)
  expect_equal(result$selected_reads, 0L)
  expect_equal(Rsamtools::countBam(output)$records, 0L)
  header <- Rsamtools::scanBamHeader(output)[[1L]]
  expect_true("refA" %in% names(header$targets))
})

test_that("multiple references can be sampled globally", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "multiple-global.bam")

  result <- subset_bam_reads(
    bam,
    output,
    reference_name = c("refA", "refB"),
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 2,
    seed = 2026
  )

  expect_equal(result$eligible_reads, 4L)
  expect_equal(result$selected_reads, 2L)
  expect_equal(sum(result$per_reference$selected_reads), 2L)
  expect_false(result$max_reads_per_reference)
  expect_equal(Rsamtools::countBam(output)$records, 2L)
})

test_that("multiple references can be sampled independently", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "multiple-per-reference.bam")

  result <- subset_bam_reads(
    bam,
    output,
    reference_name = c("refA", "refB"),
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 2,
    max_reads_per_reference = TRUE,
    seed = 2026
  )

  expect_equal(result$eligible_reads, 4L)
  expect_equal(result$selected_reads, 3L)
  expect_equal(result$per_reference$eligible_reads, c(3L, 1L))
  expect_equal(result$per_reference$selected_reads, c(2L, 1L))
  expect_true(result$max_reads_per_reference)
  expect_equal(Rsamtools::countBam(output)$records, 3L)
})

test_that("NULL reference_name samples across every reference", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "all-references.bam")

  result <- subset_bam_reads(
    bam = bam,
    output_bam = output,
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 3,
    seed = 2026
  )

  expect_null(result$reference_name)
  expect_equal(result$eligible_reads, 4L)
  expect_equal(result$selected_reads, 3L)
  expect_equal(result$per_reference$reference_name, c("refA", "refB"))
  expect_equal(result$per_reference$eligible_reads, c(3L, 1L))
  expect_equal(Rsamtools::countBam(output)$records, 3L)
})

test_that("NULL reference_name supports per-reference limits", {
  directory <- tempfile("subset-bam-")
  dir.create(directory)
  bam <- make_subset_test_bam(directory)
  output <- file.path(directory, "all-references-separated.bam")

  result <- subset_bam_reads(
    bam = bam,
    output_bam = output,
    region_start = 10,
    region_end = 15,
    min_mapq = 20,
    max_reads = 2,
    max_reads_per_reference = TRUE,
    seed = 2026
  )

  expect_equal(result$selected_reads, 3L)
  expect_equal(result$per_reference$selected_reads, c(2L, 1L))
  expect_equal(Rsamtools::countBam(output)$records, 3L)
})
