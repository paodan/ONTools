test_that("get_primers extracts a single primer pair by project and size", {
  sample_info <- data.frame(
    Project = c("P1", "P1", "P2"),
    Expected_Amplicon_Size_bp = c("1600bp", "1600bp", "3500bp"),
    Primer_Sequence_5 = c("agagtttgatcmtggctcag", "agagtttgatcmtggctcag", "AAA"),
    Primer_Sequence_3 = c("tacggytaccttgttacgactt", "tacggytaccttgttacgactt", "TTT"),
    stringsAsFactors = FALSE
  )

  primers <- get_primers(
    sample_info,
    size = "1600bp",
    project = "P1",
    project_col = "Project"
  )

  expect_equal(primers$f_primer, "AGAGTTTGATCMTGGCTCAG")
  expect_equal(primers$r_primer, "TACGGYTACCTTGTTACGACTT")
  expect_equal(primers$Primer_Sequence_5, primers$f_primer)
  expect_equal(primers$Primer_Sequence_3, primers$r_primer)
  expect_equal(nrow(primers$primer_table), 1L)
  expect_equal(primers$n_rows, 2L)
  expect_equal(primers$filters$size, "1600bp")
  expect_equal(primers$filters$project, "P1")
})

test_that("get_primers reads sample information from CSV", {
  sample_info <- data.frame(
    Expected_Amplicon_Size_bp = c(1600, 3500),
    Primer_Sequence_5 = c("AAAA", "CCCC"),
    Primer_Sequence_3 = c("TTTT", "GGGG"),
    stringsAsFactors = FALSE
  )
  sample_csv <- tempfile(fileext = ".csv")
  utils::write.csv(sample_info, sample_csv, row.names = FALSE)

  primers <- get_primers(sample_csv, size = 3500)

  expect_equal(primers$f_primer, "CCCC")
  expect_equal(primers$r_primer, "GGGG")
})

test_that("get_primers can return multiple primer pairs when requested", {
  sample_info <- data.frame(
    Expected_Amplicon_Size_bp = c("1600bp", "1600bp"),
    Primer_Sequence_5 = c("AAAA", "CCCC"),
    Primer_Sequence_3 = c("TTTT", "GGGG"),
    stringsAsFactors = FALSE
  )

  expect_error(
    get_primers(sample_info, size = "1600bp"),
    "More than one unique primer pair"
  )

  primers <- get_primers(sample_info, size = "1600bp", allow_multiple = TRUE)

  expect_equal(primers$f_primer, c("AAAA", "CCCC"))
  expect_equal(primers$r_primer, c("TTTT", "GGGG"))
  expect_equal(nrow(primers$primer_table), 2L)
})

test_that("get_primers preserves primer case when requested", {
  sample_info <- data.frame(
    Primer_Sequence_5 = "aaaa",
    Primer_Sequence_3 = "tttt",
    stringsAsFactors = FALSE
  )

  primers <- get_primers(sample_info, uppercase = FALSE)

  expect_equal(primers$f_primer, "aaaa")
  expect_equal(primers$r_primer, "tttt")
})

test_that("get_primers validates inputs", {
  sample_info <- data.frame(
    Expected_Amplicon_Size_bp = "1600bp",
    Primer_Sequence_5 = "AAAA",
    stringsAsFactors = FALSE
  )

  expect_error(get_primers(sample_info), "missing required column")
  expect_error(get_primers(sample_info, project = "P1"), "project_col")
  expect_error(get_primers(sample_info, size = NA), "size")
  expect_error(get_primers(sample_info, allow_multiple = "yes"), "allow_multiple")

  sample_info2 <- data.frame(
    Expected_Amplicon_Size_bp = "1600bp",
    Primer_Sequence_5 = "AAAA",
    Primer_Sequence_3 = "",
    stringsAsFactors = FALSE
  )

  expect_warning(
    expect_error(get_primers(sample_info2), "No complete primer pairs"),
    "missing primer sequences"
  )
})
