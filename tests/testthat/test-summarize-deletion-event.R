test_that("parse_deletion_event_cigar extracts aligned and deleted intervals", {
  parsed <- parse_deletion_event_cigar(1L, "99M61D100M")

  expect_equal(parsed$aligned$start, c(1L, 161L))
  expect_equal(parsed$aligned$end, c(99L, 260L))
  expect_equal(parsed$deletions$start, 100L)
  expect_equal(parsed$deletions$end, 160L)
  expect_equal(parsed$read_end, 260L)
})

test_that("classify_deletion_event_read identifies deletion support", {
  read <- list(
    qname = "read-del",
    rname = "PLA3_B_",
    pos = 1L,
    cigar = "99M61D100M",
    mapq = 60L,
    strand = "+"
  )

  out <- classify_deletion_event_read(
    read,
    deletion_start = 100,
    deletion_end = 160,
    breakpoint_tolerance = 5,
    min_flank_coverage = 20,
    max_wt_deletion_bases = 0
  )

  expect_equal(out$event_class, "DEL_SUPPORT")
  expect_equal(out$reason, "matched_deletion")
  expect_equal(out$observed_deletion_start, 100L)
  expect_equal(out$observed_deletion_end, 160L)
})

test_that("classify_deletion_event_read identifies WT support", {
  read <- list(
    qname = "read-wt",
    rname = "PLA3_B_",
    pos = 1L,
    cigar = "260M",
    mapq = 60L,
    strand = "+"
  )

  out <- classify_deletion_event_read(
    read,
    deletion_start = 100,
    deletion_end = 160,
    breakpoint_tolerance = 5,
    min_flank_coverage = 20,
    max_wt_deletion_bases = 0
  )

  expect_equal(out$event_class, "WT_SUPPORT")
  expect_equal(out$reason, "covers_region_no_target_deletion")
})

test_that("classify_deletion_event_read marks missing flanks as ambiguous", {
  read <- list(
    qname = "read-partial",
    rname = "PLA3_B_",
    pos = 120L,
    cigar = "120M",
    mapq = 60L,
    strand = "+"
  )

  out <- classify_deletion_event_read(
    read,
    deletion_start = 100,
    deletion_end = 160,
    breakpoint_tolerance = 5,
    min_flank_coverage = 20,
    max_wt_deletion_bases = 0
  )

  expect_equal(out$event_class, "AMBIGUOUS")
  expect_match(out$reason, "missing_left_flank")
})

test_that("classify_deletion_event_read marks boundary mismatch as ambiguous", {
  read <- list(
    qname = "read-shifted-del",
    rname = "PLA3_B_",
    pos = 1L,
    cigar = "99M51D120M",
    mapq = 60L,
    strand = "+"
  )

  out <- classify_deletion_event_read(
    read,
    deletion_start = 100,
    deletion_end = 160,
    breakpoint_tolerance = 5,
    min_flank_coverage = 20,
    max_wt_deletion_bases = 0
  )

  expect_equal(out$event_class, "AMBIGUOUS")
  expect_equal(out$reason, "deletion_boundary_mismatch")
  expect_equal(out$observed_deletion_start, 100L)
  expect_equal(out$observed_deletion_end, 150L)
})

test_that("summarize_deletion_event_table reports class frequencies", {
  read_table <- data.frame(
    event_class = c("DEL_SUPPORT", "DEL_SUPPORT", "WT_SUPPORT", "AMBIGUOUS"),
    stringsAsFactors = FALSE
  )

  out <- summarize_deletion_event_table(read_table)

  expect_equal(out$read_count[out$event_class == "DEL_SUPPORT"], 2L)
  expect_equal(out$frequency_percent[out$event_class == "DEL_SUPPORT"], 50)
})
