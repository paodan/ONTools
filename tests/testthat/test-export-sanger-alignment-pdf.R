test_that("export_sanger_alignment_pdf writes an HTML alignment", {
  alignment <- Biostrings::DNAStringSet(c(
    reference = "AAACCCGGGTTT",
    read1 = "AAACCCGGGTTT",
    read2 = "AAACCCGGTTTT"
  ))
  html_file <- tempfile(fileext = ".html")
  pdf_file <- tempfile(fileext = ".pdf")

  res <- export_sanger_alignment_pdf(
    alignment = alignment,
    pdf_file = pdf_file,
    html_file = html_file,
    render_pdf = FALSE,
    quiet = TRUE
  )

  expect_true(file.exists(html_file))
  expect_equal(normalizePath(res$html), normalizePath(html_file))
  expect_true(is.na(res$pdf))
  expect_false(file.exists(pdf_file))
  expect_equal(res$alignment_names, names(alignment))
})

test_that("export_sanger_alignment_pdf accepts SangerMultipleAlignment results", {
  reference <- Biostrings::DNAStringSet(c(
    reference = "ATGCGTACCAAAGGGTTTCCCA"
  ))
  sanger <- Biostrings::DNAStringSet(c(
    read_forward = "GTACCAAAGG"
  ))
  alignment_result <- align_sanger_with_full_reference(reference, sanger)
  html_file <- tempfile(fileext = ".html")

  res <- export_sanger_alignment_pdf(
    alignment = alignment_result,
    html_file = html_file,
    render_pdf = FALSE,
    quiet = TRUE
  )

  expect_true(file.exists(html_file))
  expect_equal(res$alignment_names, names(alignment_result$alignment))
})

test_that("export_sanger_alignment_pdf validates arguments", {
  alignment <- Biostrings::DNAStringSet(c(reference = "AAACCCGGGTTT"))
  html_file <- tempfile(fileext = ".html")
  writeLines("old", html_file)

  expect_error(
    export_sanger_alignment_pdf(
      alignment = alignment,
      html_file = html_file,
      render_pdf = FALSE,
      quiet = TRUE
    ),
    "already exists"
  )
  expect_error(
    export_sanger_alignment_pdf(
      alignment = alignment,
      html_file = tempfile(fileext = ".html"),
      block_width = 75,
      render_pdf = FALSE,
      quiet = TRUE
    ),
    "multiple of 20"
  )
  expect_error(
    export_sanger_alignment_pdf(
      alignment = "not alignment",
      html_file = tempfile(fileext = ".html"),
      render_pdf = FALSE,
      quiet = TRUE
    ),
    "XStringSet"
  )
  expect_error(
    export_sanger_alignment_pdf(
      alignment = alignment,
      html_file = tempfile(fileext = ".html"),
      render_pdf = "no",
      quiet = TRUE
    ),
    "render_pdf"
  )

  res <- export_sanger_alignment_pdf(
    alignment = alignment,
    html_file = html_file,
    render_pdf = FALSE,
    overwrite = TRUE,
    quiet = TRUE
  )
  expect_true(file.exists(res$html))
})
