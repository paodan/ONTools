test_that("bamsnap_snapshot builds a BamSnap command in dry-run mode", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  out_dir <- tempfile("bamsnap-out-")

  writeLines(c(">barcode09", "ACGTACGTACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- bamsnap_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "barcode09",
    start = 1,
    end = 12,
    out_dir = out_dir,
    snapshot_name = "barcode09.png",
    dry_run = TRUE
  )

  expect_equal(res$command, "bamsnap")
  expect_equal(res$status, 0L)
  expect_equal(res$region, "barcode09:1-12")
  expect_equal(res$format, "png")
  expect_true(res$quiet)
  expect_false(res$stdout)
  expect_false(res$stderr)
  expect_equal(res$snapshot, file.path(normalizePath(out_dir), "barcode09.png"))
  expect_true(all(c("-bam", "-pos", "-ref", "-out") %in% res$args))
  expect_equal(res$args[match("-pos", res$args) + 1L], "barcode09:1-12")
  expect_equal(res$args[match("-imagetype", res$args) + 1L], "png")
  expect_true("-save_image_only" %in% res$args)
  expect_true("-silence" %in% res$args)
})

test_that("bamsnap_snapshot can generate JPG snapshot names", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- bamsnap_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    format = "jpg",
    dry_run = TRUE
  )

  expect_equal(res$format, "jpg")
  expect_match(res$snapshot, "contig1_1-4\\.jpg$")
})

test_that("bamsnap_snapshot can show a whole contig when start and end are omitted", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- bamsnap_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    dry_run = TRUE
  )

  expect_equal(res$region, "contig1")
  expect_equal(res$args[match("-pos", res$args) + 1L], "contig1")
  expect_match(res$snapshot, "contig1\\.png$")
})

test_that("bamsnap_snapshot passes display-related options", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- bamsnap_snapshot(
    genome_fasta = ref,
    bam = bam,
    chr = "contig1",
    start = 1,
    end = 4,
    title = "barcode09",
    draw = c("coordinates", "bamplot"),
    bamplot = c("coverage", "read"),
    width = 1200,
    height = 500,
    margin = 20,
    no_title = TRUE,
    no_target_line = TRUE,
    process = 2,
    extra_args = c("-read_group", "sample"),
    dry_run = TRUE,
    quiet = FALSE
  )

  expect_equal(res$args[match("-title", res$args) + 1L], "barcode09")
  expect_equal(res$args[match("-width", res$args) + 1L], "1200")
  expect_equal(res$args[match("-height", res$args) + 1L], "500")
  expect_equal(res$args[match("-margin", res$args) + 1L], "20")
  expect_equal(res$args[match("-process", res$args) + 1L], "2")
  expect_true("-no_title" %in% res$args)
  expect_true("-no_target_line" %in% res$args)
  expect_true("-read_group" %in% res$args)
  expect_false("-silence" %in% res$args)
  expect_equal(res$stdout, "")
  expect_equal(res$stderr, "")
})

test_that("bamsnap_snapshot creates a missing FASTA index with Rsamtools", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), ref)
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  res <- bamsnap_snapshot(
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

test_that("bamsnap_snapshot validates inputs", {
  ref <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGT"), ref)
  file.create(paste0(ref, ".fai"))
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_error(
    bamsnap_snapshot(ref, bam, chr = "contig1", start = 1, dry_run = TRUE),
    "both be supplied"
  )
  expect_error(
    bamsnap_snapshot(
      ref,
      bam,
      chr = "contig1",
      start = 1,
      end = 4,
      format = "jpg",
      snapshot_name = "contig1.png",
      dry_run = TRUE
    ),
    "extension must match"
  )
  expect_error(
    bamsnap_snapshot(ref, bam, chr = "contig1", process = 0, dry_run = TRUE),
    "positive integer"
  )
})
