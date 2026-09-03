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

  write_abif_file(ab1, sequence, trace$channels, trace$peak_locations,
                  trace$qualities, "sample1")

  tags <- read_abif_test_tags(ab1)
  expect_true(all(
    c("PBAS:2", "PLOC:2", "PCON:2", "DATA:9", "DATA:10", "DATA:11",
      "DATA:12", "FWO_:1") %in% tags
  ))
  expect_gt(file.info(ab1)$size, 500)
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
  expect_equal(res$sequence_length, nchar(sequence))
  expect_true(file.exists(ab1))
  expect_false(file.exists(paste0(consensus, ".fai")))
  expect_true("PBAS:2" %in% read_abif_test_tags(ab1))
})

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
