test_that("read_vcf parses single-sample VCF files", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode001",
    "chr1\t10\t.\tA\tG\t60\tPASS\tDP=30;AF=0.4;SOMATIC\tGT:GQ:DP:AD:AF\t0/1:99:30:18,12:0.4",
    "chr1\t20\trs1\tTT\tT\t.\tq10\tDP=12;NOTE=low\tGT:DP\t1/1:12"
  ), vcf)

  res <- read_vcf(vcf)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$CHROM, c("chr1", "chr1"))
  expect_equal(res$POS, c(10L, 20L))
  expect_equal(res$QUAL, c(60, NA))
  expect_equal(res$Type, c("SNP", "INDEL"))
  expect_equal(res$DP, c(30, 12))
  expect_equal(res$AF, c(0.4, NA))
  expect_equal(res$SOMATIC, c(TRUE, NA))
  expect_equal(res$sampleID, c("barcode001", "barcode001"))
  expect_equal(res$GT, c("0/1", "1/1"))
  expect_equal(res$GQ, c(99, NA))
  expect_equal(res$FORMAT_DP, c(30, 12))
  expect_equal(res$FORMAT_AF, c(0.4, NA))
  expect_equal(res$AD, c("18,12", NA))
})

test_that("vcf_variant_type infers wf-amplicon style variant types", {
  expect_equal(
    vcf_variant_type(
      ref = c("A", "AT", "A", "AC", "A", "A", "CC"),
      alt = c("G", "A", "AT", "GT", "A", ".", "G")
    ),
    c("SNP", "INDEL", "INDEL", "MNP", "REF", NA, "OTHER")
  )

  expect_equal(vcf_variant_type("A", "G,AT"), "MIXED")
  expect_equal(vcf_variant_type("CC", "G"), "OTHER")
  expect_equal(vcf_variant_type("A", "G,T"), "SNP")
  expect_equal(vcf_variant_type("A", "G,AT", collapse_multiallelic = FALSE)[[1]],
               c("SNP", "INDEL"))
})

test_that("read_vcf can omit inferred variant Type", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tA\tG\t60\tPASS\tDP=30"
  ), vcf)

  res <- read_vcf(vcf, add_variant_type = FALSE)

  expect_false("Type" %in% names(res))
})

test_that("read_vcf derives medaka allele depth columns from SR and AR", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode279",
    paste0(
      "PLA3_A_\t76\t.\tCG\tC\t1.858\tPASS\t",
      "DP=37;DPS=13,24;DPSP=9;SR=0,1,4,3;AR=0,1;",
      "SC=681,782,693,790\tGT:GQ\t1:2"
    )
  ), vcf)

  res <- read_vcf(vcf)

  expect_equal(res$ref_fwd_depth, 0)
  expect_equal(res$ref_rev_depth, 1)
  expect_equal(res$alt_fwd_depth, 4)
  expect_equal(res$alt_rev_depth, 3)
  expect_equal(res$ref_depth, 1)
  expect_equal(res$alt_depth, 7)
  expect_equal(res$variant_percent, 87.5)
  expect_equal(res$ambiguous_fwd_depth, 0)
  expect_equal(res$ambiguous_rev_depth, 1)
  expect_equal(res$ambiguous_depth, 1)
})

test_that("read_vcf can derive medaka allele depth columns without parsing INFO output", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tA\tG\t60\tPASS\tSR=2,3,4,5;AR=1,2"
  ), vcf)

  res <- read_vcf(vcf, parse_info = FALSE)

  expect_false("SR" %in% names(res))
  expect_equal(res$ref_depth, 5)
  expect_equal(res$alt_depth, 9)
  expect_equal(res$variant_percent, 9 / 14 * 100)
  expect_equal(res$ambiguous_depth, 3)
})

test_that("read_vcf sums SR depth across multiple ALT alleles", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tA\tG,AT\t60\tPASS\tSR=1,2,3,4,5,6"
  ), vcf)

  res <- read_vcf(vcf)

  expect_equal(res$ref_depth, 3)
  expect_equal(res$alt_fwd_depth, 8)
  expect_equal(res$alt_rev_depth, 10)
  expect_equal(res$alt_depth, 18)
  expect_equal(res$variant_percent, 18 / 21 * 100)
})

test_that("read_vcf adds reference flanking sequences", {
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(">chr1", "ACGTACGTACGT"), ref)
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t5\t.\tA\tG\t60\tPASS\tDP=30",
    "chr1\t5\t.\tAC\tA\t60\tPASS\tDP=30"
  ), vcf)

  res <- read_vcf(vcf, reference_fasta = ref, flank_width = 4)

  expect_equal(res$ref_upstream_4bp, c("ACGT", "ACGT"))
  expect_equal(res$ref_downstream_4bp, c("CGTA", "GTAC"))
})

test_that("read_vcf reference flanks are clipped at sequence boundaries", {
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(">chr1", "ACGTACGT"), ref)
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1\t.\tA\tG\t60\tPASS\tDP=30",
    "chr1\t8\t.\tT\tC\t60\tPASS\tDP=30"
  ), vcf)

  res <- read_vcf(vcf, reference_fasta = ref, flank_width = 4)

  expect_equal(res$ref_upstream_4bp, c("", "TACG"))
  expect_equal(res$ref_downstream_4bp, c("CGTA", ""))
})

test_that("read_vcf validates CHROM values when adding reference flanks", {
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(">chr2", "ACGTACGT"), ref)
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1\t.\tA\tG\t60\tPASS\tDP=30"
  ), vcf)

  expect_error(
    read_vcf(vcf, reference_fasta = ref),
    "not found in `reference_fasta`"
  )
})

test_that("read_vcf parses gzipped VCF files", {
  vcf <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    "chr2\t5\t.\tC\tT\t50\tPASS\tDP=9\tGT:DP\t0/1:9"
  ), con)
  close(con)
  on.exit(NULL)

  res <- read_vcf(vcf)

  expect_equal(nrow(res), 1)
  expect_equal(res$CHROM, "chr2")
  expect_equal(res$DP, 9)
  expect_equal(res$FORMAT_DP, 9)
  expect_equal(res$GT, "0/1")
})

test_that("read_vcf supports multi-sample wide and long output", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t10\t.\tA\tG\t60\tPASS\tDP=30\tGT:DP\t0/1:12\t1/1:18"
  ), vcf)

  wide <- read_vcf(vcf)
  long <- read_vcf(vcf, sample_format = "long")
  s2 <- read_vcf(vcf, sample_format = "long", samples = "s2")

  expect_true(all(c("s1_GT", "s1_DP", "s2_GT", "s2_DP") %in% names(wide)))
  expect_equal(wide$s1_DP, 12)
  expect_equal(wide$s2_GT, "1/1")
  expect_equal(nrow(long), 2)
  expect_equal(long$sampleID, c("s1", "s2"))
  expect_equal(long$DP, c(30, 30))
  expect_equal(long$FORMAT_DP, c(12, 18))
  expect_equal(nrow(s2), 1)
  expect_equal(s2$sampleID, "s2")
})

test_that("read_vcf can keep raw INFO and FORMAT fields", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    "chr1\t10\t.\tA\tG\t60\tPASS\tDP=30\tGT:DP\t0/1:30"
  ), vcf)

  res <- read_vcf(vcf, keep_info = TRUE, keep_format = TRUE)

  expect_true(all(c("INFO", "FORMAT", "sample1") %in% names(res)))
  expect_equal(res$INFO, "DP=30")
  expect_equal(res$sample1, "0/1:30")
})

test_that("read_vcf handles empty VCF files", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ), vcf)

  res <- read_vcf(vcf)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
  expect_true(all(c("CHROM", "POS", "INFO") %in% names(res)))
})

test_that("parse_vcf_info supports flags, missing fields, and prefixes", {
  res <- parse_vcf_info(c("DP=10;AF=0.2;SOMATIC", ".", "NOTE=low"), prefix = "INFO_")

  expect_equal(res$INFO_DP, c(10, NA, NA))
  expect_equal(res$INFO_AF, c(0.2, NA, NA))
  expect_equal(res$INFO_SOMATIC, c(TRUE, NA, NA))
  expect_equal(res$INFO_NOTE, c(NA, NA, "low"))
})

test_that("read_vcf validates inputs", {
  expect_error(read_vcf("missing.vcf"), "vcf_file")
  expect_error(read_vcf(tempfile()), "vcf_file")
})
