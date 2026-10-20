read_abif_test_tags <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- rawToChar(readBin(con, "raw", n = 4L))
  expect_equal(magic, "ABIF")
  readBin(con, "integer", n = 1L, size = 2L, endian = "big")
  root <- readBin(con, "raw", n = 28L)
  count <- read_uint32_test(root[13:16])
  directory_offset <- read_uint32_test(root[21:24])
  seek(con, directory_offset, origin = "start")
  vapply(
    seq_len(count),
    function(i) {
      entry <- readBin(con, "raw", n = 28L)
      paste0(rawToChar(entry[1:4]), ":", read_uint32_test(entry[5:8]))
    },
    character(1)
  )
}

read_uint32_test <- function(bytes) {
  sum(as.integer(bytes) * c(16777216, 65536, 256, 1))
}

test_that("synthetic_ab1_from_bam validates arguments", {
  consensus <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  ab1 <- tempfile(fileext = ".ab1")
  writeLines(c(">contig1", "ACGTACGT"), consensus)
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  expect_error(
    synthetic_ab1_from_bam(consensus, bam, ab1, spacing = 0, echo = FALSE),
    "spacing"
  )
  expect_error(
    synthetic_ab1_from_bam(consensus, bam, ab1, noise_fraction = 2, echo = FALSE),
    "noise_fraction"
  )
  expect_error(
    synthetic_ab1_from_bam(
      consensus, bam, ab1, peak_value_mode = "unknown", echo = FALSE
    ),
    "arg"
  )
  expect_error(
    synthetic_ab1_from_bam("missing.fasta", bam, ab1, echo = FALSE),
    "consensus"
  )
})

test_that("pileup parser counts reference bases, mismatches, and indels", {
  counts <- parse_pileup_bases(".,Aa+2tt-1c^]G$", "A")

  expect_equal(unname(counts[c("A", "G")]), c(4L, 1L))
  expect_equal(sum(counts), 5L)
})

test_that("synthetic ABIF writer creates expected core tags", {
  ab1 <- tempfile(fileext = ".ab1")
  sequence <- "ACGTACGT"
  counts <- lapply(seq_len(nchar(sequence)), function(i) {
    stats::setNames(c(3L, 0L, 0L, 0L), c("A", "C", "G", "T"))
  })
  names(counts) <- as.character(seq_along(counts))
  trace <- synthetic_ab1_trace(sequence, counts, spacing = 12, sigma = 2,
                               peak_height = 900, baseline = 15,
                               noise_fraction = 0.03)

  write_abif_file(ab1, trace$sequence, trace$channels, trace$peak_locations,
                  trace$qualities, "sample1")

  tags <- read_abif_test_tags(ab1)
  expect_true(all(
    c("PBAS:2", "PLOC:2", "PCON:2", "DATA:9", "DATA:10", "DATA:11",
      "DATA:12", "FWO_:1") %in% tags
  ))
  expect_gt(file.info(ab1)$size, 500)
})

test_that("pileup support controls peak ratios and base calls", {
  counts <- list(
    `1` = stats::setNames(c(80L, 0L, 20L, 0L), c("A", "C", "G", "T")),
    `2` = stats::setNames(c(5L, 0L, 5L, 0L), c("A", "C", "G", "T"))
  )
  fractions <- ab1_base_values(counts[["1"]], "G", 0.03, "relative")
  trace <- synthetic_ab1_trace(
    sequence = "GGN",
    pileup_counts = counts,
    spacing = 12,
    sigma = 2,
    peak_height = 900,
    baseline = 15,
    noise_fraction = 0.03
  )

  expect_equal(unname(fractions), c(1, 0.03, 0.25, 0.03))
  expect_equal(trace$sequence, "AGN")
  expect_equal(trace$qualities, c(7L, 3L, 20L))
  expect_gt(
    trace$channels$A[trace$peak_locations[[1L]] + 1L],
    trace$channels$G[trace$peak_locations[[1L]] + 1L]
  )
})

test_that("FASTA bases are used only when pileup coverage is absent", {
  fractions <- ab1_base_values(NULL, "C", 0.03, "relative")
  called <- pileup_called_bases(
    c("C", "T"),
    list(`1` = NULL, `2` = stats::setNames(
      c(4L, 0L, 0L, 6L),
      c("A", "C", "G", "T")
    ))
  )

  expect_equal(unname(fractions), c(0.03, 1, 0.03, 0.03))
  expect_equal(called, c("C", "T"))
})

test_that("peak value modes have distinct documented meanings", {
  counts <- stats::setNames(
    c(100L, 24L, 0L, 1000L),
    c("A", "C", "G", "T")
  )

  relative <- ab1_base_values(counts, "T", 0.03, "relative")
  count <- ab1_base_values(counts, "T", 0.03, "count")
  fraction <- ab1_base_values(counts, "T", 0.03, "fraction")

  expect_equal(unname(relative), c(0.1, 0.03, 0.03, 1))
  expect_equal(unname(count), c(100, 24, 0, 1000))
  expect_equal(unname(fraction), c(100, 24, 0, 1000) / 1124)
  expect_equal(
    unname(ab1_base_values(NULL, "T", 0.03, "count")),
    rep(0, 4)
  )
  expect_equal(
    unname(ab1_base_values(NULL, "T", 0.03, "fraction")),
    rep(0, 4)
  )
})

test_that("count mode stores read support as amplitude above baseline", {
  counts <- list(`1` = stats::setNames(
    c(100L, 24L, 0L, 1000L),
    c("A", "C", "G", "T")
  ))
  trace <- synthetic_ab1_trace(
    sequence = "T",
    pileup_counts = counts,
    spacing = 12,
    sigma = 2,
    peak_height = 900,
    baseline = 1,
    noise_fraction = 0.03,
    peak_value_mode = "count"
  )
  center <- trace$peak_locations[[1L]] + 1L

  expect_equal(trace$channels$A[[center]] - 1L, 100L)
  expect_equal(trace$channels$C[[center]] - 1L, 24L)
  expect_equal(trace$channels$G[[center]] - 1L, 0L)
  expect_equal(trace$channels$T[[center]] - 1L, 1000L)
})

test_that("synthetic_ab1_from_bam runs on a tiny BAM", {
  skip_if(Sys.which("samtools") == "", "samtools is not installed")

  work <- tempfile("synthetic-ab1-bam-")
  dir.create(work)
  consensus <- file.path(work, "consensus.fasta")
  sam <- file.path(work, "reads.sam")
  bam <- file.path(work, "reads.bam")
  ab1 <- file.path(work, "synthetic.ab1")
  sequence <- "ACGTACGTACGT"
  writeLines(c(">contig1", sequence), consensus)
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      paste0("@SQ\tSN:contig1\tLN:", nchar(sequence)),
      paste0("read1\t0\tcontig1\t1\t60\t12M\t*\t0\t0\t", sequence, "\tIIIIIIIIIIII"),
      "read2\t0\tcontig1\t1\t60\t12M\t*\t0\t0\tACGTGCGTACGT\tIIIIIIIIIIII"
    ),
    sam
  )
  system2("samtools", c("view", "-bS", sam), stdout = bam)
  system2("samtools", c("index", bam))

  res <- synthetic_ab1_from_bam(consensus, bam, ab1, sample = "tiny",
                                echo = FALSE, stderr = FALSE)

  expect_equal(res$status, 0L)
  expect_equal(res$reference_name, "contig1")
  expect_equal(res$called_sequence, sequence)
  expect_equal(res$sequence_length, nchar(sequence))
  expect_equal(res$peak_value_mode, "relative")
  expect_true(file.exists(ab1))
  expect_false(file.exists(paste0(consensus, ".fai")))
  expect_true("PBAS:2" %in% read_abif_test_tags(ab1))
})
