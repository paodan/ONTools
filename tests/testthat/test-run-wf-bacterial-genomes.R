test_that("run_wf_bacterial_genomes builds filter and Nextflow commands", {
  fastq <- tempfile("fastq-")
  dir.create(file.path(fastq, "barcode009"), recursive = TRUE)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq, "barcode009", "barcode009.fastq.gz"))
  sample_sheet <- tempfile(fileext = ".csv")
  writeLines("sample,barcode", sample_sheet)

  res <- run_wf_bacterial_genomes(
    fastq = fastq,
    out_dir = "results",
    work_dir = "work",
    filtered_fastq_dir = file.path(tempdir(), "filtered-fastq"),
    threads = 20,
    min_len = 3000,
    max_len = 12000,
    sample_sheet = sample_sheet,
    flye_genome_size = "5m",
    flye_asm_coverage = 300,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$command, "nextflow")
  expect_true(res$filter$enabled)
  expect_length(res$filter$command_strings, 1)
  expect_match(res$filter$command_strings, "seqkit")
  expect_match(res$filter$command_strings, "gzip")
  expect_match(res$filter$command_strings, "-m '3000' -M '12000'", fixed = TRUE)
  expect_true(any(res$args == "epi2me-labs/wf-bacterial-genomes"))
  expect_true(any(res$args == "--fastq"))
  expect_true(any(res$args == res$paths$workflow_fastq))
  expect_true(any(res$args == "--min_read_length"))
  expect_true(any(res$args == "3000"))
  expect_true(any(res$args == "--threads"))
  expect_true(any(res$args == "20"))
  expect_true(any(res$args == "-work-dir"))
  expect_true(any(res$args == "work"))
  expect_true(any(res$args == "--sample_sheet"))
  expect_true(any(res$args == normalizePath(sample_sheet)))
  expect_true(any(res$args == "--flye_genome_size"))
  expect_true(any(res$args == "5m"))
  expect_true(any(res$args == "--flye_asm_coverage"))
  expect_true(any(res$args == "300"))
  expect_true(any(res$args == "-resume"))
  expect_true(res$uses_shell)
  expect_equal(res$execution_command, "sh")
  expect_match(res$command_string, "-offline", fixed = TRUE)
})

test_that("run_wf_bacterial_genomes can skip pre-filtering", {
  res <- run_wf_bacterial_genomes(
    fastq = "fastq",
    out_dir = "results",
    work_dir = "work",
    max_len = NULL,
    resume = FALSE,
    extra_args = NULL,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_false(res$filter$enabled)
  expect_equal(res$filter$command_strings, character())
  expect_equal(res$paths$workflow_fastq, "fastq")
  expect_equal(res$execution_command, "nextflow")
  expect_false(res$uses_shell)
  expect_false(any(res$args == "-resume"))
})

test_that("run_wf_bacterial_genomes supports isolates and quiet mode", {
  res <- run_wf_bacterial_genomes(
    fastq = "fastq",
    max_len = NULL,
    isolates = TRUE,
    quiet = TRUE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$args[1:3], c("-q", "run", "epi2me-labs/wf-bacterial-genomes"))
  expect_true(any(res$args == "--isolates"))
})

test_that("run_wf_bacterial_genomes validates filter inputs", {
  expect_error(
    run_wf_bacterial_genomes(min_len = 12000, max_len = 3000,
                             dry_run = TRUE, echo = FALSE),
    "less than or equal"
  )
  expect_error(
    run_wf_bacterial_genomes(fastq = tempfile("missing-"),
                             dry_run = TRUE, echo = FALSE),
    "fastq"
  )

  fastq <- tempfile("fastq-")
  dir.create(file.path(fastq, "barcode001"), recursive = TRUE)
  expect_error(
    run_wf_bacterial_genomes(fastq = fastq, dry_run = TRUE, echo = FALSE),
    "not found"
  )
})
