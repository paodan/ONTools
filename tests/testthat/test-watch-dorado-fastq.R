fake_realtime_dorado <- function() {
  path <- tempfile("fake-dorado-")
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "subcommand=\"${1:-}\"",
      "shift",
      "if [[ \"$subcommand\" == 'basecaller' ]]; then",
      "  printf 'fake-bam\\n'",
      "  exit 0",
      "fi",
      "if [[ \"$subcommand\" == 'demux' ]]; then",
      "  out=''",
      "  while [[ $# -gt 0 ]]; do",
      "    if [[ \"$1\" == '--output-dir' ]]; then out=\"$2\"; shift 2; else shift; fi",
      "  done",
      "  mkdir -p \"$out/run/fastq_pass/sampleA\" \"$out/run/fastq_pass/sampleB\"",
      "  printf '@a\\nAC\\n+\\nII\\n' > \"$out/run/fastq_pass/sampleA/chunk.fastq\"",
      "  printf '@b\\nGT\\n+\\nII\\n' > \"$out/run/fastq_pass/sampleB/chunk.fastq\"",
      "  exit 0",
      "fi",
      "exit 2"
    ),
    path
  )
  Sys.chmod(path, "0755")
  path
}

test_that("basecall_one adds barcode and sample sheet arguments", {
  pod5 <- tempfile(fileext = ".pod5")
  sheet <- tempfile(fileext = ".csv")
  writeLines("pod5", pod5)
  writeLines("experiment_id,flow_cell_id,kit", sheet)

  result <- basecall_one(
    pod5,
    tempfile(fileext = ".bam"),
    kit_name = "SQK-NBD114-96",
    sample_sheet = sheet,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true(all(c("--kit-name", "SQK-NBD114-96", "--sample-sheet") %in% result$args))
  expect_equal(result$args[[length(result$args)]], normalizePath(sheet))
})

test_that("process_one_pod5 publishes per-alias FASTQ files and a marker", {
  root <- tempfile("realtime-")
  pod5_dir <- file.path(root, "pod5")
  bam_dir <- file.path(root, "bam")
  fastq_dir <- file.path(root, "fastq")
  state_dir <- file.path(root, "state")
  dir.create(pod5_dir, recursive = TRUE)
  pod5 <- file.path(pod5_dir, "batch01.pod5")
  sheet <- file.path(root, "sample_sheet.csv")
  writeLines("pod5", pod5)
  writeLines("experiment_id,flow_cell_id,kit", sheet)

  result <- process_one_pod5(
    pod5,
    bam_dir,
    fastq_dir,
    state_dir,
    kit_name = "SQK-NBD114-96",
    sample_sheet = sheet,
    dorado = fake_realtime_dorado(),
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(result$status, "completed")
  expect_true(file.exists(file.path(bam_dir, "batch01.bam")))
  expect_true(file.exists(file.path(fastq_dir, "sampleA", "batch01__sampleA.fastq")))
  expect_true(file.exists(file.path(fastq_dir, "sampleB", "batch01__sampleB.fastq")))
  expect_true(file.exists(file.path(state_dir, "batch01.done")))

  repeated <- process_one_pod5(
    pod5,
    bam_dir,
    fastq_dir,
    state_dir,
    kit_name = "SQK-NBD114-96",
    sample_sheet = sheet,
    dorado = fake_realtime_dorado(),
    echo = FALSE,
    stderr = FALSE
  )
  expect_equal(repeated$status, "skipped")
})

test_that("finalize_run_fastqs concatenates batches for each alias", {
  root <- tempfile("finalize-")
  fastq_dir <- file.path(root, "fastq")
  final_dir <- file.path(root, "final")
  dir.create(file.path(fastq_dir, "sampleA"), recursive = TRUE)
  dir.create(file.path(fastq_dir, "sampleB"), recursive = TRUE)
  writeLines(c("@a", "AC", "+", "II"), file.path(fastq_dir, "sampleA", "batch01__sampleA.fastq"))
  writeLines(c("@b", "GT", "+", "II"), file.path(fastq_dir, "sampleA", "batch02__sampleA.fastq"))
  writeLines(c("@c", "AA", "+", "II"), file.path(fastq_dir, "sampleB", "batch01__sampleB.fastq"))

  result <- finalize_run_fastqs(fastq_dir, final_dir)

  expect_setequal(result$alias, c("sampleA", "sampleB"))
  expect_equal(readLines(file.path(final_dir, "sampleA.fastq")),
               c("@a", "AC", "+", "II", "@b", "GT", "+", "II"))
  expect_equal(result$n_fragments[result$alias == "sampleA"], 2L)
})

test_that("finalize_run_fastqs preserves alias subdirectories and writes MD5", {
  root <- tempfile("finalize-layout-")
  fastq_dir <- file.path(root, "fastq")
  final_dir <- file.path(root, "final")
  dir.create(file.path(fastq_dir, "sampleA"), recursive = TRUE)
  writeLines(c("@a", "AC", "+", "II"),
             file.path(fastq_dir, "sampleA", "batch01__sampleA.fastq"))
  writeLines(c("@b", "GT", "+", "II"),
             file.path(fastq_dir, "sampleA", "batch02__sampleA.fastq"))

  result <- finalize_run_fastqs(
    fastq_path = fastq_dir,
    output_path = final_dir,
    preserve_subdirs = TRUE,
    write_md5 = TRUE
  )

  expected_fastq <- file.path(final_dir, "sampleA", "sampleA.fastq")
  expected_md5 <- paste0(expected_fastq, ".md5")
  expect_equal(result$output_fastq, normalizePath(expected_fastq))
  expect_equal(result$md5_file, normalizePath(expected_md5))
  expect_true(file.exists(expected_md5))
  expect_equal(
    readLines(expected_md5),
    paste(unname(tools::md5sum(expected_fastq)), basename(expected_fastq))
  )
  expect_true(dir.exists(file.path(fastq_dir, "sampleA")))
})

test_that("watch_dorado_fastq processes stable POD5 files once", {
  root <- tempfile("watch-realtime-")
  pod5_dir <- file.path(root, "pod5")
  dir.create(pod5_dir, recursive = TRUE)
  writeLines("pod5", file.path(pod5_dir, "batch01.pod5"))
  sheet <- file.path(root, "sample_sheet.csv")
  writeLines("experiment_id,flow_cell_id,kit", sheet)

  result <- watch_dorado_fastq(
    pod5_path = pod5_dir,
    bam_path = file.path(root, "bam"),
    fastq_path = file.path(root, "fastq"),
    state_path = file.path(root, "state"),
    kit_name = "SQK-NBD114-96",
    sample_sheet = sheet,
    max_run_time = 0.5,
    check_interval = 0.01,
    idle_interval = 0.01,
    dorado = fake_realtime_dorado(),
    echo = FALSE,
    stderr = FALSE
  )

  expect_named(result$batches, "batch01")
  expect_equal(result$batches$batch01$status, "completed")
  expect_true(file.exists(file.path(root, "state", "batch01.done")))
})
