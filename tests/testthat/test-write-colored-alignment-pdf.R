test_that("write_colored_alignment_pdf writes a PDF", {
  alignment <- Biostrings::DNAStringSet(c(
    Sample_1 = "ACTGAA--CTGA",
    Sample_2 = "ACTGAA--CTGA",
    Sample_3 = "ACTGAATTCT-A"
  ))
  pdf_file <- tempfile(fileext = ".pdf")

  expect_message(
    result <- write_colored_alignment_pdf(
      alignment,
      pdf_file,
      bases_per_block = 5,
      blocks_per_page = 2,
      character_spacing = 1.2,
      line_spacing = 1.1
    ),
    "PDF file:"
  )

  expect_true(file.exists(pdf_file))
  expect_gt(file.info(pdf_file)$size, 0)
  expect_equal(normalizePath(result), normalizePath(pdf_file))
})

test_that("write_colored_alignment_pdf supports coordinate offsets", {
  alignment <- Biostrings::DNAStringSet(c(
    Sample_1 = "AAAAACCCCC",
    Sample_2 = "AAAATCCCCC"
  ))
  pdf_file <- tempfile(fileext = ".pdf")

  expect_silent(
    result <- write_colored_alignment_pdf(
      alignment,
      pdf_file,
      coordinate_start = 1750,
      show_conservation = FALSE,
      quiet = TRUE
    )
  )

  expect_true(file.exists(result))
})

test_that("write_colored_alignment_pdf validates its inputs", {
  alignment <- Biostrings::DNAStringSet(c(
    Sample_1 = "ACTG",
    Sample_2 = "ACTG"
  ))

  expect_error(
    write_colored_alignment_pdf("ACTG", tempfile(fileext = ".pdf")),
    "DNAStringSet"
  )
  expect_error(
    write_colored_alignment_pdf(
      Biostrings::DNAStringSet(c("ACTG", "ACT")),
      tempfile(fileext = ".pdf")
    ),
    "same width"
  )
  expect_error(
    write_colored_alignment_pdf(
      alignment,
      tempfile(fileext = ".pdf"),
      character_spacing = 0
    ),
    "character_spacing"
  )
  expect_error(
    write_colored_alignment_pdf(
      alignment,
      tempfile(fileext = ".pdf"),
      base_colors = c(A = "not-a-color")
    ),
    "valid R colors"
  )

  pdf_file <- tempfile(fileext = ".pdf")
  write_colored_alignment_pdf(alignment, pdf_file, quiet = TRUE)
  expect_error(
    write_colored_alignment_pdf(alignment, pdf_file, quiet = TRUE),
    "already exists"
  )
  expect_silent(
    write_colored_alignment_pdf(
      alignment,
      pdf_file,
      overwrite = TRUE,
      quiet = TRUE
    )
  )
})
