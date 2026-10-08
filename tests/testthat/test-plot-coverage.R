make_fake_coverage_samtools <- function(lines, status = 0L) {
  fake_bin <- tempfile("coverage-samtools-")
  dir.create(fake_bin)
  script <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "[[ \"${1:-}\" == 'depth' ]] || exit 2",
    sprintf("exit_status=%d", status),
    "[[ \"$exit_status\" -eq 0 ]] || exit \"$exit_status\"",
    sprintf("printf '%%s\\n' %s", paste(shQuote(lines), collapse = " "))
  )
  path <- file.path(fake_bin, "samtools")
  writeLines(script, path)
  Sys.chmod(path, mode = "0755")
  path
}

test_that("plot_coverage returns coverage plot and threshold summary", {
  bam <- tempfile(fileext = ".bam")
  file.create(bam, paste0(bam, ".bai"))
  samtools <- make_fake_coverage_samtools(c(
    "#CHROM\tPOS\tsample.bam",
    "contig1\t1\t0",
    "contig1\t2\t10",
    "contig1\t3\t20",
    "contig2\t1\t50"
  ))

  plot <- plot_coverage(
    bam_file = bam,
    min_base_quality = 5,
    min_mapping_quality = 10,
    depth_thresholds = c(10, 20, 50),
    major_breaks = 100,
    minor_breaks = 10,
    samtools = samtools,
    echo = FALSE,
    stderr = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_named(plot$data, c("reference", "position", "depth"))
  expect_equal(plot$data$depth, c(0L, 10L, 20L, 50L))
  expect_equal(
    attr(plot, "depth_summary"),
    data.frame(
      threshold = c(10L, 20L, 50L),
      positions = c(3L, 2L, 1L),
      percent = c(75, 50, 25)
    )
  )
  expect_match(attr(plot, "command"), "'depth' '-aa' '-H'", fixed = TRUE)
  expect_match(attr(plot, "command"), "'-q' '5' '-Q' '10'", fixed = TRUE)
  expect_true(inherits(plot$facet, "FacetWrap"))
})

test_that("plot_coverage supports original argument names", {
  bam <- tempfile(fileext = ".bam")
  file.create(bam, paste0(bam, ".bai"))
  samtools <- make_fake_coverage_samtools(c(
    "#CHROM\tPOS\tsample.bam",
    "contig1\t1\t10"
  ))

  plot <- plot_coverage(
    bamFile = bam,
    minSeqQ = 1,
    minMapQ = 2,
    depth = c(10, 20),
    breaks = 200,
    minor_breaks = 20,
    samtools = samtools,
    echo = FALSE,
    stderr = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(attr(plot, "depth_summary")$positions, c(1L, 0L))
  expect_match(attr(plot, "command"), "'-q' '1' '-Q' '2'", fixed = TRUE)
})

test_that("plot_coverage handles empty output and validates arguments", {
  bam <- tempfile(fileext = ".bam")
  file.create(bam, paste0(bam, ".bai"))
  samtools <- make_fake_coverage_samtools("#CHROM\tPOS\tsample.bam")

  expect_warning(
    expect_null(plot_coverage(bam, samtools = samtools, echo = FALSE, stderr = FALSE)),
    "No coverage data"
  )
  expect_error(plot_coverage(), "bam_file")
  expect_error(
    plot_coverage(bam, depth_thresholds = c(10, -1), samtools = samtools,
                  echo = FALSE),
    "depth_thresholds"
  )
  expect_error(
    plot_coverage(bam, major_breaks = 0, samtools = samtools, echo = FALSE),
    "major_breaks"
  )
  expect_error(
    plot_coverage(bam_file = bam, bamFile = bam, samtools = samtools,
                  echo = FALSE),
    "only one"
  )
})

test_that("plot_coverage reports samtools failures", {
  bam <- tempfile(fileext = ".bam")
  file.create(bam, paste0(bam, ".bai"))
  samtools <- make_fake_coverage_samtools(character(), status = 7L)

  expect_error(
    plot_coverage(bam, samtools = samtools, echo = FALSE, stderr = FALSE),
    "exit status: 7"
  )
})
