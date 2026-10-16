test_that("extract_aligned_region applies CIGAR operations consistently", {
  out <- extract_aligned_region(
    reference = "ACGTACGT",
    read_seq = "ACGTTACGT",
    cigar = "4M1I4M",
    ref_start = 1,
    x = 2,
    y = 7
  )

  expect_equal(out$reference_aligned, "CGT-ACG")
  expect_equal(out$read_aligned, "CGTTACG")
  expect_error(
    extract_aligned_region("ACGTACGT", "ACGTACGT", "4Mbad4M", 1, 1, 8),
    "complete CIGAR"
  )
})

test_that("call_read_edit_events calls substitutions and normalizes case", {
  wt <- call_read_edit_events(
    "wt", "acgtacgt", "8M", 1, "ACGTACGT", 1, 8
  )
  snp <- call_read_edit_events(
    "snp", "ACATACGT", "8M", 1, "acgtacgt", 1, 8
  )

  expect_equal(wt$events$type, "WT")
  expect_equal(snp$events$type, "SNP")
  expect_equal(snp$events$label, "SNP_3_3:G>A")
})

test_that("call_read_edit_events handles insertion and deletion boundaries", {
  insertion <- call_read_edit_events(
    "ins", "ACGTTACGT", "4M1I4M", 1, "ACGTACGT", 1, 8
  )
  deletion <- call_read_edit_events(
    "del", "ACGTCGT", "4M1D3M", 1, "ACGTACGT", 1, 8
  )
  outside <- call_read_edit_events(
    "outside", "ACGTTACGT", "4M1I4M", 1, "ACGTACGT", 1, 4
  )

  expect_equal(insertion$events$label, "IN_after_4:T")
  expect_equal(deletion$events$label, "DEL_5_5:A")
  expect_equal(outside$events$type, "WT")
})

test_that("only related indels are marked as complex events", {
  out <- call_read_edit_events(
    "complex",
    "ATGTCCGT",
    "1M1I1D6M",
    1,
    "ACGTACGT",
    1,
    8
  )$events

  expect_equal(out$type, c("OTHER", "OTHER", "SNP"))
  expect_match(out$label[[1L]], "IN_after_1")
  expect_match(out$label[[2L]], "DEL_2_2")
  expect_equal(out$label[[3L]], "SNP_5_5:A>C")
})

test_that("skipped reference regions are non-callable", {
  skipped <- call_read_edit_events(
    "skip", "ACGCGT", "3M2N3M", 1, "ACGTACGT", 1, 8
  )

  expect_equal(skipped$status, "NO_CALL")
  expect_equal(skipped$events$type, "NO_CALL")
  expect_equal(skipped$events$label, "NO_CALL_N_skip")
})

test_that("filter_reads_by_region respects anchors", {
  alignments <- data.frame(
    id = c("full", "left_short", "right_short"),
    pos = c(1L, 3L, 1L),
    cigar = c("12M", "10M", "9M")
  )

  out <- filter_reads_by_region(alignments, 4, 9, min_anchor = 2)
  expect_equal(out$id, "full")
})

test_that("make_variant_events returns typed events and typed empty output", {
  alignments <- data.frame(
    qname = c("wt", "snp"),
    seq = c("ACGTACGT", "ACATACGT"),
    cigar = c("8M", "8M"),
    pos = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  events <- make_variant_events(
    alignments,
    reference = c(amplicon = "ACGTACGT"),
    region_start = 1,
    region_end = 8,
    min_anchor = 0,
    mc.cores = 1
  )
  empty <- make_variant_events(
    transform(alignments, pos = 4L),
    reference = c(amplicon = "ACGTACGTACGT"),
    region_start = 1,
    region_end = 8,
    min_anchor = 0,
    mc.cores = 1
  )

  expect_equal(events$type, c("WT", "SNP"))
  expect_equal(unique(events$Ref_name), "amplicon")
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("read_id", "sample_id", "Ref_name") %in% names(empty)))
})

test_that("summarise_variant_events excludes every NO_CALL subtype", {
  alignments <- data.frame(
    qname = c("wt", "snp"),
    seq = c("ACGTACGT", "ACATACGT"),
    cigar = c("8M", "8M"),
    pos = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  events <- make_variant_events(
    alignments,
    reference = c(amplicon = "ACGTACGT"),
    region_start = 1,
    region_end = 8,
    min_anchor = 0,
    mc.cores = 1
  )
  no_call <- events[1L, , drop = FALSE]
  no_call$read_id <- "skip"
  no_call$type <- "NO_CALL"
  no_call$label <- "NO_CALL_N_skip"

  summary <- summarise_variant_events(rbind(events, no_call), min_ratio = 1)

  expect_equal(unique(summary$event_frequency$Callable_reads), 2L)
  expect_false("NO_CALL_N_skip" %in% summary$event_frequency$label)
  expect_equal(sum(summary$var_info$Read_count), 2L)
  expect_error(summarise_variant_events(events, min_ratio = 101),
               "between 0 and 100")
})

test_that("editing caches are invalidated by signature changes", {
  bam <- tempfile(fileext = ".bam")
  writeLines("placeholder", bam)
  cache <- tempfile(fileext = ".rds")
  signature <- editing_cache_signature(
    bam,
    reference = "ACGT",
    reference_name = "amplicon",
    region_start = 1,
    region_end = 4,
    min_anchor = 0,
    min_mapq = 20,
    reference_name_map = NULL
  )
  events <- data.frame(read_id = "read1")
  write_editing_cache(events, cache, signature)

  expect_equal(read_editing_cache(cache, signature), events)
  changed <- signature
  changed$reference_sequence <- "TGCA"
  expect_null(read_editing_cache(cache, changed))

  saveRDS(events, cache)
  expect_null(read_editing_cache(cache, signature))
})

test_that("sample sheet normalization accepts legacy names", {
  legacy <- data.frame(
    barcode = "barcode01",
    Reference_name = "amplicon",
    Reference_Region_Start = 2,
    Reference_Region_End = 7
  )

  out <- normalize_editing_sample_info(legacy)
  expect_equal(out$ref_name, "amplicon")
  expect_equal(out$region_start, 2L)
  expect_equal(out$region_end, 7L)
})
