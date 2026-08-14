test_that("make_ont_fastq_delivery builds the default command", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  sample_sheet <- tempfile(fileext = ".csv")
  dir.create(fastq_dir)
  dir.create(output_dir)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq_dir, "sample.fastq"))
  writeLines("barcode,sample\nbarcode01,sampleA", sample_sheet)

  res <- make_ont_fastq_delivery(
    input = fastq_dir,
    output = output_dir,
    project = "PROJECT001",
    sample_sheet = sample_sheet,
    threads = 4,
    run_nanoplot = FALSE,
    run_multiqc = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$execution_command, "bash")
  expect_equal(res$status, NA_integer_)
  expect_true(any(res$args == "--input"))
  expect_true(any(res$args == normalizePath(fastq_dir)))
  expect_true(any(res$args == "--output"))
  expect_true(any(res$args == normalizePath(output_dir)))
  expect_true(any(res$args == "--project"))
  expect_true(any(res$args == "PROJECT001"))
  expect_true(any(res$args == "--threads"))
  expect_true(any(res$args == "4"))
  expect_true(any(res$args == "--sample-sheet"))
  expect_true(any(res$args == normalizePath(sample_sheet)))
  expect_true(any(res$args == "--skip-nanoplot"))
  expect_true(any(res$args == "--skip-multiqc"))
  expect_match(res$command_string, "make_ont_fastq_delivery.sh", fixed = TRUE)
  expect_equal(res$paths$delivery_dir, file.path(normalizePath(output_dir), "PROJECT001_delivery"))
  expect_equal(res$paths$archive, file.path(normalizePath(output_dir), "PROJECT001_delivery.tar.gz"))
})

test_that("make_ont_fastq_delivery validates arguments", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  dir.create(fastq_dir)

  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      threads = 0,
      dry_run = TRUE,
      echo = FALSE
    ),
    "threads"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      run_nanoplot = "no",
      dry_run = TRUE,
      echo = FALSE
    ),
    "run_nanoplot"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      sample_sheet = file.path(output_dir, "missing.csv"),
      dry_run = TRUE,
      echo = FALSE
    ),
    "sample_sheet"
  )
})
