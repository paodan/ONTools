test_that("filter_fastq builds the expected dry-run command", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- filter_fastq(
    fastq = fastq,
    output_fastq = output_fastq,
    min_length = 5000,
    max_length = 20000,
    min_quality = 10,
    reverse = TRUE,
    upper_case = TRUE,
    threads = 8,
    extra_args = c("--seq-type", "dna"),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$command, "seqkit")
  expect_equal(res$args[1:5], c("seq", "-m", "5000", "-M", "20000"))
  expect_true("-Q" %in% res$args)
  expect_true("--reverse" %in% res$args)
  expect_true("--upper-case" %in% res$args)
  expect_true("--threads" %in% res$args)
  expect_true("--seq-type" %in% res$args)
  expect_equal(res$paths$fastq, normalizePath(fastq))
  expect_equal(res$paths$output_fastq, normalizePath(output_fastq, mustWork = FALSE))
  expect_match(res$command_string, "seqkit", fixed = TRUE)
})

test_that("filter_fastq dry-run supports conda_env", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- filter_fastq(
    fastq = fastq,
    output_fastq = output_fastq,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$command, "conda")
  expect_equal(res$args[1:5], c("run", "-n", "ont-tools", "seqkit", "seq"))
  expect_match(res$command_string, "'conda' 'run' '-n' 'ont-tools' 'seqkit' 'seq'", fixed = TRUE)
})

test_that("filter_fastq runs seqkit", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_fastq <- tempfile(fileext = ".fastq")
  fake_bin <- tempfile("seqkit-bin-")
  dir.create(fake_bin)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "[[ \"$1\" == 'seq' ]]",
      "[[ \"$2\" == '-m' ]]",
      "[[ \"$3\" == '5000' ]]",
      "cat \"$4\""
    ),
    file.path(fake_bin, "seqkit")
  )
  Sys.chmod(file.path(fake_bin, "seqkit"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- filter_fastq(
    fastq = fastq,
    output_fastq = output_fastq,
    min_length = 5000,
    echo = FALSE
  )

  expect_equal(res$status, 0L)
  expect_equal(readLines(output_fastq), readLines(fastq))
})

test_that("filter_fastq_by_length remains available as a compatibility wrapper", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- filter_fastq_by_length(
    fastq = fastq,
    output_fastq = output_fastq,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$args[1:3], c("seq", "-m", "5000"))
})

test_that("filter_fastq validates arguments", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  expect_error(
    filter_fastq(fastq, output_fastq, min_length = 0, dry_run = TRUE, echo = FALSE),
    "min_length"
  )
  expect_error(
    filter_fastq("missing.fastq.gz", output_fastq, dry_run = TRUE, echo = FALSE),
    "fastq"
  )
  expect_error(
    filter_fastq(fastq, output_fastq, conda_env = "", dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
  expect_error(
    filter_fastq(fastq, output_fastq, max_quality = -1, dry_run = TRUE, echo = FALSE),
    "max_quality"
  )
})
