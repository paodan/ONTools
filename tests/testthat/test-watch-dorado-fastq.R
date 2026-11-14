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

test_that("create_dorado_run_structure creates the standardized layout", {
  root <- tempfile("standard-run-")
  layout <- create_dorado_run_structure(root)

  expect_s3_class(layout, "dorado_run_structure")
  expect_equal(layout$run_path, normalizePath(root))
  expect_true(all(vapply(
    layout[c(
      "bam_path", "fastq_path", "final_fastq_path",
      "qc_path", "qc_cache_path", "state_path", "log_path"
    )],
    dir.exists,
    logical(1)
  )))
  expect_false("pod5_path" %in% names(layout))
  expect_false("sample_sheet" %in% names(layout))
})

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
    update_qc = FALSE,
    dorado = fake_realtime_dorado(),
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(result$status, "completed")
  expect_true(file.exists(file.path(bam_dir, "batch01.bam")))
  expect_true(file.exists(file.path(fastq_dir, "sampleA", "batch01__sampleA.fastq.gz")))
  expect_true(file.exists(file.path(fastq_dir, "sampleB", "batch01__sampleB.fastq.gz")))
  expect_true(file.exists(file.path(state_dir, "batch01.done")))

  repeated <- process_one_pod5(
    pod5,
    bam_dir,
    fastq_dir,
    state_dir,
    kit_name = "SQK-NBD114-96",
    sample_sheet = sheet,
    update_qc = FALSE,
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

  result <- finalize_run_fastqs(fastq_dir, final_dir, compress = FALSE)

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

  expected_fastq <- file.path(final_dir, "sampleA", "sampleA.fastq.gz")
  expected_md5 <- paste0(expected_fastq, ".md5")
  expect_equal(result$output_fastq, normalizePath(expected_fastq))
  expect_equal(result$md5_file, normalizePath(expected_md5))
  expect_true(file.exists(expected_md5))
  con <- gzfile(expected_fastq, open = "rt")
  expect_equal(readLines(con), c("@a", "AC", "+", "II", "@b", "GT", "+", "II"))
  close(con)
  expect_equal(
    readLines(expected_md5),
    paste(unname(tools::md5sum(expected_fastq)), basename(expected_fastq))
  )
  expect_true(dir.exists(file.path(fastq_dir, "sampleA")))
})

test_that("finalize_run_fastqs recompresses FASTQ.GZ fragments", {
  root <- tempfile("finalize-gzip-")
  alias_dir <- file.path(root, "fastq", "sampleA")
  dir.create(alias_dir, recursive = TRUE)
  first <- gzfile(file.path(alias_dir, "batch01__sampleA.fastq.gz"), open = "wt")
  writeLines(c("@a", "AC", "+", "II"), first)
  close(first)
  second <- gzfile(file.path(alias_dir, "batch02__sampleA.fastq.gz"), open = "wt")
  writeLines(c("@b", "GT", "+", "II"), second)
  close(second)

  result <- finalize_run_fastqs(
    fastq_path = file.path(root, "fastq"),
    output_path = file.path(root, "final")
  )

  expect_match(result$output_fastq, "[.]fastq[.]gz$")
  merged <- gzfile(result$output_fastq, open = "rt")
  expect_equal(readLines(merged), c("@a", "AC", "+", "II", "@b", "GT", "+", "II"))
  close(merged)
})

test_that("finalize_run_fastqs resolves standardized run paths", {
  root <- tempfile("finalize-standard-run-")
  layout <- create_dorado_run_structure(root)
  alias_dir <- file.path(layout$fastq_path, "sampleA")
  dir.create(alias_dir)
  writeLines(
    c("@a", "AC", "+", "II"),
    file.path(alias_dir, "batch01__sampleA.fastq")
  )

  result <- finalize_run_fastqs(run_path = root, compress = FALSE)

  expect_equal(
    result$output_fastq,
    normalizePath(file.path(layout$final_fastq_path, "sampleA.fastq"))
  )
})

test_that("directory QC merges each alias and preserves the last good plot", {
  root <- tempfile("live-qc-")
  alias_a <- file.path(root, "fastq", "sampleA")
  alias_b <- file.path(root, "fastq", "sampleB")
  unclassified <- file.path(root, "fastq", "unclassified")
  qc_dir <- file.path(root, "qc")
  dir.create(alias_a, recursive = TRUE)
  dir.create(alias_b, recursive = TRUE)
  dir.create(unclassified, recursive = TRUE)

  write_qc_fastq <- function(path, start, lengths) {
    records <- unlist(lapply(seq_along(lengths), function(i) {
      sequence <- paste(rep("A", lengths[[i]]), collapse = "")
      quality <- paste(rep("I", lengths[[i]]), collapse = "")
      c(
        paste0("@read", start + i, " qs:f:15.5 ch:i:1 SM:Z:barcode01"),
        sequence,
        "+",
        quality
      )
    }), use.names = FALSE)
    connection <- gzfile(path, open = "wt")
    writeLines(records, connection)
    close(connection)
  }
  write_qc_fastq(
    file.path(alias_a, "batch01__sampleA.fastq.gz"),
    0L,
    95:105
  )
  write_qc_fastq(
    file.path(alias_a, "batch02__sampleA.fastq.gz"),
    100L,
    295:305
  )
  write_qc_fastq(
    file.path(alias_b, "batch01__sampleB.fastq.gz"),
    200L,
    c(145:155, 395:405)
  )
  write_qc_fastq(
    file.path(unclassified, "batch01__unclassified.fastq.gz"),
    300L,
    c(195:205, 495:505)
  )

  result <- update_fastq_directory_qc(
    fastq_path = file.path(root, "fastq"),
    qc_path = qc_dir,
    include_unclassified = FALSE,
    plot_args = list(min_distance = 100, density_n = 1024),
    error_action = "stop"
  )

  expect_equal(result$status, "updated")
  expect_equal(result$n_aliases, 2L)
  expect_equal(result$n_fragments, 3L)
  expect_true(file.exists(result$plot_file))
  first_md5 <- unname(tools::md5sum(result$plot_file))

  failed <- update_fastq_directory_qc(
    fastq_path = file.path(root, "fastq"),
    qc_path = qc_dir,
    include_unclassified = FALSE,
    plot_args = list(heights = 1),
    error_action = "silent"
  )
  expect_equal(failed$status, "error")
  expect_true(file.exists(failed$error_file))
  expect_equal(unname(tools::md5sum(failed$plot_file)), first_md5)
})

test_that("watch_dorado_fastq processes stable POD5 files once", {
  root <- tempfile("watch-realtime-")
  layout <- create_dorado_run_structure(root)
  input_dir <- tempfile("watch-input-")
  dir.create(input_dir)
  pod5_dir <- file.path(input_dir, "pod5")
  dir.create(pod5_dir)
  pod5_file <- file.path(pod5_dir, "batch01.pod5")
  sample_sheet <- file.path(input_dir, "sample_sheet.csv")
  writeLines("pod5", pod5_file)
  writeLines("experiment_id,flow_cell_id,kit", sample_sheet)

  result <- watch_dorado_fastq(
    pod5_path = pod5_dir,
    sample_sheet = sample_sheet,
    run_path = root,
    kit_name = "SQK-NBD114-96",
    max_run_time = 0.5,
    check_interval = 0.01,
    idle_interval = 0.01,
    update_qc = FALSE,
    dorado = fake_realtime_dorado(),
    echo = FALSE,
    stderr = FALSE
  )

  expect_named(result$batches, "batch01")
  expect_equal(result$batches$batch01$status, "completed")
  expect_true(file.exists(file.path(layout$state_path, "batch01.done")))
  expect_equal(result$paths$run_path, normalizePath(root))
  expect_true(file.exists(pod5_file))
  expect_true(file.exists(sample_sheet))
})
