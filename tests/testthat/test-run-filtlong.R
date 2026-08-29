test_that("run_filtlong builds expected dry-run command", {
  reads <- tempfile(fileext = ".fastq")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- run_filtlong(
    reads = reads,
    output_fastq = output_fastq,
    target_bases = 80000000,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$command, "filtlong")
  expect_equal(res$args, c("--target_bases", "80000000", normalizePath(reads)))
  expect_equal(res$paths$output_fastq, normalizePath(output_fastq, mustWork = FALSE))
  expect_match(res$command_string, "--target_bases", fixed = TRUE)
  expect_match(res$command_string, ">", fixed = TRUE)
})

test_that("run_filtlong supports common options and conda_env", {
  reads <- tempfile(fileext = ".fastq")
  illumina_1 <- tempfile(fileext = ".fastq")
  illumina_2 <- tempfile(fileext = ".fastq")
  assembly <- tempfile(fileext = ".fasta")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
  writeLines(c("@short1", "ACGT", "+", "!!!!"), illumina_1)
  writeLines(c("@short2", "ACGT", "+", "!!!!"), illumina_2)
  writeLines(c(">contig1", "ACGT"), assembly)

  res <- run_filtlong(
    reads = reads,
    output_fastq = output_fastq,
    target_bases = 80000000,
    keep_percent = 90,
    min_length = 5000,
    min_mean_q = 9,
    min_window_q = 7,
    length_weight = 1,
    mean_q_weight = 10,
    window_q_weight = 0,
    window_size = 250,
    trim = TRUE,
    split = 500,
    illumina_1 = illumina_1,
    illumina_2 = illumina_2,
    assembly = assembly,
    verbose = TRUE,
    extra_args = c("--some-future-option"),
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$command, "conda")
  expect_equal(res$args[1:4], c("run", "-n", "ont-tools", "filtlong"))
  expect_true(all(c(
    "--target_bases", "80000000",
    "--keep_percent", "90",
    "--min_length", "5000",
    "--min_mean_q", "9",
    "--min_window_q", "7",
    "--length_weight", "1",
    "--mean_q_weight", "10",
    "--window_q_weight", "0",
    "--window_size", "250",
    "--trim",
    "--split", "500",
    "--illumina_1", normalizePath(illumina_1),
    "--illumina_2", normalizePath(illumina_2),
    "--assembly", normalizePath(assembly),
    "--verbose",
    "--some-future-option",
    normalizePath(reads)
  ) %in% res$args))
})

test_that("run_filtlong runs filtlong command", {
  reads <- tempfile(fileext = ".fastq")
  output_fastq <- tempfile(fileext = ".fastq")
  fake_bin <- tempfile("filtlong-bin-")
  dir.create(fake_bin)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "input=\"${@: -1}\"",
      "cat \"$input\""
    ),
    file.path(fake_bin, "filtlong")
  )
  Sys.chmod(file.path(fake_bin, "filtlong"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- run_filtlong(
    reads = reads,
    output_fastq = output_fastq,
    target_bases = 4,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_equal(readLines(output_fastq), readLines(reads))
})

test_that("run_filtlong validates arguments", {
  reads <- tempfile(fileext = ".fastq")
  output_fastq <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  expect_error(
    run_filtlong(reads, output_fastq, target_bases = 0,
                 dry_run = TRUE, echo = FALSE),
    "target_bases"
  )
  expect_error(
    run_filtlong(reads, output_fastq, keep_percent = 101,
                 dry_run = TRUE, echo = FALSE),
    "keep_percent"
  )
  expect_error(
    run_filtlong(reads, output_fastq, min_length = 0,
                 dry_run = TRUE, echo = FALSE),
    "min_length"
  )
  expect_error(
    run_filtlong(reads, output_fastq, conda_env = "",
                 dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
  expect_error(
    run_filtlong("missing.fastq", output_fastq, dry_run = TRUE, echo = FALSE),
    "reads"
  )
})
