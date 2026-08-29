test_that("igv_snapshot writes an IGV batch file in dry-run mode", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  bai <- paste0(bam, ".bai")
  out_dir <- tempfile("igv-out-")
  batch_file <- tempfile(fileext = ".igv")

  writeLines(c(">barcode09", "ACGTACGTACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(bai)

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "barcode09",
    start = 1,
    end = 12,
    out_dir = out_dir,
    snapshot_name = "barcode09.png",
    batch_file = batch_file,
    dry_run = TRUE
  )

  expect_equal(res$command, "igv.sh")
  expect_equal(res$args, c("-b", batch_file))
  expect_equal(res$status, 0L)
  expect_equal(res$format, "png")
  expect_true(res$quiet)
  expect_false(res$stdout)
  expect_false(res$stderr)
  expect_true(file.exists(batch_file))
  expect_equal(readLines(batch_file), res$commands)
  expect_true(any(grepl("^genome ", res$commands)))
  expect_true(any(grepl("^load ", res$commands)))
  expect_equal(sum(res$commands == "goto barcode09:1-12"), 2L)
  expect_true(any(res$commands == "sort base"))
  expect_true(any(res$commands == "collapse"))
  expect_true(any(res$commands == "sleep 1"))
  expect_true(any(res$commands == "snapshot barcode09.png"))
  expect_equal(res$snapshot, file.path(normalizePath(out_dir), "barcode09.png"))
})

test_that("igv_snapshot can leave command output unsuppressed for debugging", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    quiet = FALSE,
    dry_run = TRUE
  )

  expect_false(res$quiet)
  expect_equal(res$stdout, "")
  expect_equal(res$stderr, "")
})

test_that("igv_snapshot can generate PDF snapshot names", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    format = "pdf",
    dry_run = TRUE
  )

  expect_equal(res$format, "pdf")
  expect_true(any(res$commands == "snapshot contig1_1-4.pdf"))
  expect_match(res$snapshot, "\\.pdf$")
})

test_that("igv_snapshot can show a whole contig when start and end are omitted", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    dry_run = TRUE
  )

  expect_true(any(res$commands == "goto contig1"))
  expect_true(any(res$commands == "snapshot contig1.png"))
})

test_that("igv_snapshot can skip snapshot delay", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    snapshot_delay = 0,
    dry_run = TRUE
  )

  expect_false(any(grepl("^sleep ", res$commands)))
  expect_equal(sum(res$commands == "goto contig1:1-4"), 2L)
})

test_that("igv_snapshot validates snapshot delay", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_error(
    igv_snapshot(ref, bam, chr = "contig1", start = 1, end = 4,
                 snapshot_delay = -1, dry_run = TRUE),
    "snapshot_delay"
  )
})

test_that("igv_snapshot requires start and end together", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_error(
    igv_snapshot(
      genome_fasta = ref,
      bam = bam,
      chr = "contig1",
      start = 1,
      dry_run = TRUE
    ),
    "both be supplied"
  )
})

test_that("igv_snapshot checks snapshot extension against format", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_error(
    igv_snapshot(
      genome_fasta = ref,
      bam = bam,
      chr = "contig1",
      start = 1,
      end = 4,
      format = "pdf",
      snapshot_name = "contig1.png",
      dry_run = TRUE
    ),
    "extension must match"
  )
})

test_that("igv_snapshot supports xvfb command wrapping", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    igv = "/opt/IGV/igv.sh",
    use_xvfb = TRUE,
    dry_run = TRUE
  )

  expect_equal(res$command, "xvfb-run")
  expect_equal(res$args[1:2], c("-a", "/opt/IGV/igv.sh"))
})

test_that("igv_snapshot creates a missing FASTA index with Rsamtools", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), ref)
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- igv_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 8,
    dry_run = TRUE
  )

  expect_true(file.exists(paste0(ref, ".fai")))
  expect_equal(res$fasta_index, normalizePath(paste0(ref, ".fai")))
  expect_true(res$fasta_index_created)
})

test_that("igv_snapshot can leave missing FASTA index as a warning", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), ref)
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_warning(
    res <- igv_snapshot(
      genome_fasta = ref,
      bam = bam,
      chr = "contig1",
      start = 1,
      end = 8,
      auto_index_fasta = FALSE,
      dry_run = TRUE
    ),
    "No FASTA index"
  )

  expect_false(file.exists(paste0(ref, ".fai")))
  expect_false(res$fasta_index_created)
})

test_that("igv_snapshot validates inputs before writing a batch file", {
  expect_error(
    igv_snapshot(
      "missing.fasta",
      "missing.bam",
      chr = "chr1",
      start = 1,
      end = 10,
      dry_run = TRUE
    ),
    "genome_fasta"
  )
})

test_that("igv_snapshot validates region coordinates", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">chr1", "ACGT"), ref)
  file.create(bam)

  expect_error(
    igv_snapshot(ref, bam, chr = "chr1", start = 10, end = 1, dry_run = TRUE),
    "start"
  )
})
