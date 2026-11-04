make_bam_locus_test_bam <- function(directory) {
  sam <- file.path(directory, "loci.sam")
  sequence_10 <- "ACGTACGTAC"
  sequence_15 <- "ACGTACGTACGTACG"
  quality_10 <- strrep("F", 10)
  quality_15 <- strrep("F", 15)
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      "@SQ\tSN:refA\tLN:200",
      "@SQ\tSN:refB\tLN:200",
      paste("softclip", 0, "refA", 10, 20, "5S10M", "*", 0, 0,
            sequence_15, quality_15, sep = "\t"),
      paste("deletion", 0, "refA", 20, 30, "5M3D5M", "*", 0, 0,
            sequence_10, quality_10, sep = "\t"),
      paste("duplicate", 1024, "refA", 40, 40, "10M", "*", 0, 0,
            sequence_10, quality_10, sep = "\t"),
      paste("qcfail", 512, "refB", 50, 50, "10M", "*", 0, 0,
            sequence_10, quality_10, sep = "\t"),
      paste("secondary", 256, "refB", 60, 60, "10M", "*", 0, 0,
            sequence_10, quality_10, sep = "\t"),
      paste("supplementary", 2048, "refB", 70, 60, "10M", "*", 0, 0,
            sequence_10, quality_10, sep = "\t"),
      paste("unmapped", 4, "*", 0, 0, "*", "*", 0, 0,
            sequence_10, quality_10, sep = "\t")
    ),
    sam
  )
  Rsamtools::asBam(
    sam,
    destination = file.path(directory, "loci"),
    overwrite = TRUE,
    indexDestination = FALSE
  )
}

test_that("summarize_bam_alignment_loci calculates reference spans from CIGAR", {
  directory <- tempfile("bam-loci-")
  dir.create(directory)
  bam <- make_bam_locus_test_bam(directory)

  result <- summarize_bam_alignment_loci(
    bam,
    min_mapq = 20,
    min_reference_reads = 1,
    bin_width = 10
  )

  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$filter_summary$total_records, 7L)
  expect_equal(result$filter_summary$flag_pass_records, 4L)
  expect_equal(result$filter_summary$retained_records, 4L)
  expect_equal(result$alignments$end[result$alignments$read_id == "softclip"], 19L)
  expect_equal(result$alignments$reference_width[
    result$alignments$read_id == "softclip"
  ], 10L)
  expect_equal(result$alignments$end[result$alignments$read_id == "deletion"], 32L)
  expect_equal(result$alignments$reference_width[
    result$alignments$read_id == "deletion"
  ], 13L)
  expect_equal(result$reference_summary$alignment_records, c(3L, 1L))
  expect_equal(result$peak_summary$bin_start[1], 1L)
})

test_that("summarize_bam_alignment_loci applies optional filters", {
  directory <- tempfile("bam-loci-")
  dir.create(directory)
  bam <- make_bam_locus_test_bam(directory)

  result <- summarize_bam_alignment_loci(
    bam,
    min_mapq = 21,
    reference_name = "refA",
    exclude_duplicates = TRUE,
    exclude_qcfail = TRUE,
    min_reference_reads = 1
  )

  expect_equal(result$alignments$read_id, "deletion")
  expect_equal(result$filter_summary$reference_pass_records, 2L)
  expect_equal(result$filter_summary$retained_unique_reads, 1L)
  expect_equal(result$mapq_distribution$mapq, c(20L, 30L))
})

test_that("summarize_bam_alignment_loci validates reference names", {
  directory <- tempfile("bam-loci-")
  dir.create(directory)
  bam <- make_bam_locus_test_bam(directory)

  expect_error(
    summarize_bam_alignment_loci(bam, reference_name = ""),
    "non-empty"
  )
})

test_that("summarize_bam_alignment_loci returns typed empty results", {
  directory <- tempfile("bam-loci-")
  dir.create(directory)
  bam <- make_bam_locus_test_bam(directory)

  result <- summarize_bam_alignment_loci(
    bam,
    reference_name = "missing-reference",
    min_reference_reads = 1
  )

  expect_equal(nrow(result$alignments), 0L)
  expect_equal(nrow(result$reference_summary), 0L)
  expect_equal(nrow(result$peak_summary), 0L)
  expect_equal(result$filter_summary$retained_records, 0L)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})
