test_that("normalize_linked_variants splits multi-allelic candidates", {
  variants <- data.frame(
    CHROM = "chr1",
    POS = 10,
    REF = "A",
    ALT = "G,T"
  )

  out <- normalize_linked_variants(variants)

  expect_equal(nrow(out), 2)
  expect_equal(out$variant_id, c("chr1:10:A>G", "chr1:10:A>T"))
})

test_that("classify_read_variants captures linked SNP insertion and deletion calls", {
  read <- list(
    qname = "read1",
    rname = "chr1",
    pos = 1L,
    cigar = "4M2I3M2D3M",
    seq = "ACGTGGAAAGGG",
    qual = rep(30L, 12),
    mapq = 60L,
    strand = "+"
  )
  variants <- normalize_linked_variants(data.frame(
    CHROM = "chr1",
    POS = c(5, 4, 7),
    REF = c("C", "T", "ACC"),
    ALT = c("A", "TGG", "A")
  ))

  calls <- classify_read_variants(read, variants, min_baseq = 10, no_call_label = "NO_CALL")

  expect_equal(calls$call, c("ALT", "ALT", "ALT"))
})

test_that("summarize_linked_haplotype_table reports frequencies", {
  read_table <- data.frame(
    haplotype = c("WT", "WT", "chr1:10:A>G", "chr1:10:A>G;chr1:20:C>T"),
    stringsAsFactors = FALSE
  )

  out <- summarize_linked_haplotype_table(read_table)

  expect_equal(out$read_count[match("WT", out$haplotype)], 2)
  expect_equal(out$frequency_percent[match("WT", out$haplotype)], 50)
})

test_that("parse_linked_region parses reference and interval", {
  region <- parse_linked_region("chr1:100-250")

  expect_equal(region$reference, "chr1")
  expect_equal(region$start, 100L)
  expect_equal(region$end, 250L)
})
