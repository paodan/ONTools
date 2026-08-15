test_that("run_dorado_demux_to_fastq exposes barcode-both-ends control", {
  proj <- tempfile("dorado-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  default_cmd <- run_dorado_demux_to_fastq(
    proj = proj,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_true("--barcode-both-ends" %in% default_cmd$args)

  relaxed_cmd <- run_dorado_demux_to_fastq(
    proj = proj,
    barcode_both_ends = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_false("--barcode-both-ends" %in% relaxed_cmd$args)
})
