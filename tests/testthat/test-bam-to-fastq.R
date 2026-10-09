make_fake_bam_to_fastq_samtools <- function(status = 0L, write_output = TRUE) {
  fake_bin <- tempfile("bam-to-fastq-bin-")
  dir.create(fake_bin)
  path <- file.path(fake_bin, "samtools")
  script <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "[[ \"${1:-}\" == 'fastq' ]] || exit 2",
    sprintf("status=%d", status),
    "[[ \"$status\" -eq 0 ]] || exit \"$status\"",
    "shift",
    "output=''",
    "while [[ $# -gt 0 ]]; do",
    "  case \"$1\" in",
    "    -0|-1|-2) [[ -n \"$output\" ]] || output=\"$2\"; shift 2 ;;",
    "    -@|-F) shift 2 ;;",
    "    -n|-t|-O) shift ;;",
    "    *) shift ;;",
    "  esac",
    "done"
  )
  if (isTRUE(write_output)) {
    script <- c(
      script,
      "printf '@read1\\nACGT\\n+\\n!!!!\\n' > \"$output\""
    )
  }
  writeLines(script, path)
  Sys.chmod(path, mode = "0755")
  path
}

test_that("bam_to_fastq builds a Nanopore single-end command", {
  bam <- tempfile(fileext = ".bam")
  output <- tempfile(fileext = ".fastq.gz")
  file.create(bam)
  result <- bam_to_fastq(
    bam, output, threads = 8, dry_run = TRUE, echo = FALSE
  )
  expect_equal(result$status, NA_integer_)
  expect_equal(result$read_category, "unpaired")
  expect_match(result$command_string, "'samtools' 'fastq' '-@' '8'", fixed = TRUE)
  expect_match(result$command_string, "'-n' '-F' '0x900' '-0'", fixed = TRUE)
  expect_false(file.exists(output))
})

test_that("bam_to_fastq converts to a non-empty output", {
  bam <- tempfile(fileext = ".bam")
  output <- tempfile(fileext = ".fastq")
  file.create(bam)
  samtools <- make_fake_bam_to_fastq_samtools()
  result <- bam_to_fastq(
    bam, output, samtools = samtools, echo = FALSE, stderr = FALSE
  )
  expect_equal(result$status, 0L)
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_equal(readLines(output), c("@read1", "ACGT", "+", "!!!!"))
})

test_that("bam_to_fastq supports all read categories and optional tags", {
  bam <- tempfile(fileext = ".bam")
  output <- tempfile(fileext = ".fastq")
  file.create(bam)
  result <- bam_to_fastq(
    bam,
    output,
    read_category = "all",
    copy_tags = TRUE,
    use_original_quality = TRUE,
    keep_read_names = FALSE,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )
  expect_match(
    result$command_string,
    "'conda' 'run' '-n' 'ont-tools' 'samtools'",
    fixed = TRUE
  )
  expect_match(result$command_string, "'-t' '-O' '-F' '0x900'", fixed = TRUE)
  expect_equal(sum(result$args == "-0"), 1L)
  expect_equal(sum(result$args == "-1"), 1L)
  expect_equal(sum(result$args == "-2"), 1L)
  expect_equal(sum(result$args == "-n"), 1L)
})

test_that("bam_to_fastq protects outputs and reports failures", {
  bam <- tempfile(fileext = ".bam")
  output <- tempfile(fileext = ".fastq")
  file.create(bam, output)
  expect_error(
    bam_to_fastq(bam, output, dry_run = TRUE, echo = FALSE),
    "already exists"
  )
  expect_error(
    bam_to_fastq(bam, tempfile(), threads = 0, dry_run = TRUE, echo = FALSE),
    "threads"
  )
  failed_output <- tempfile(fileext = ".fastq")
  failed_samtools <- make_fake_bam_to_fastq_samtools(status = 7L)
  expect_error(
    bam_to_fastq(
      bam,
      failed_output,
      samtools = failed_samtools,
      echo = FALSE,
      stderr = FALSE
    ),
    "exit status: 7"
  )
  expect_false(file.exists(failed_output))
})

test_that("bam_to_fastq rejects empty samtools output", {
  bam <- tempfile(fileext = ".bam")
  output <- tempfile(fileext = ".fastq")
  file.create(bam)
  samtools <- make_fake_bam_to_fastq_samtools(write_output = FALSE)
  expect_error(
    bam_to_fastq(
      bam, output, samtools = samtools, echo = FALSE, stderr = FALSE
    ),
    "non-empty FASTQ"
  )
  expect_false(file.exists(output))
})
