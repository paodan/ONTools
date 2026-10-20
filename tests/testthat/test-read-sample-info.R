make_sectioned_workbook <- function(path) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "2_metadata")
  metadata <- matrix(NA_character_, nrow = 10L, ncol = 4L)
  metadata[1L, 1L] <- "[header]"
  metadata[2L, 1:3] <- c("Project_ID", "Project identifier", "P001")
  metadata[4L, 1L] <- "[columns]"
  metadata[5L, 1:2] <- c("Sample_ID", "Sample identifier")
  metadata[7L, 1L] <- "[data]"
  metadata[8L, ] <- c("Sample_ID", "Barcode_ID", "Barcode_ID", NA)
  metadata[9L, ] <- c("S1", "NB01", "duplicate", "extra value")
  metadata[10L, ] <- c("S2", "NB02", "duplicate 2", NA)
  openxlsx::writeData(
    wb, "2_metadata", metadata,
    colNames = FALSE, rowNames = FALSE, keepNA = FALSE
  )

  openxlsx::addWorksheet(wb, "3_gel")
  gel <- matrix(NA_character_, nrow = 2L, ncol = 5L)
  gel[1L, ] <- c("[Example]", NA, "Gel_lane", "Gel_Quality", "Clear_Band")
  openxlsx::writeData(
    wb, "3_gel", gel,
    colNames = FALSE, rowNames = FALSE, keepNA = FALSE
  )

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

test_that("read_sample_info reads sheets and parses known block types", {
  skip_if_not_installed("openxlsx")
  workbook <- tempfile(fileext = ".xlsx")
  make_sectioned_workbook(workbook)

  info <- suppressWarnings(read_sample_info(workbook))

  expect_named(info, c("2_metadata", "3_gel"))
  expect_named(info[["2_metadata"]], c("header", "columns", "data"))
  expect_named(
    info[["2_metadata"]][["header"]],
    c("field", "description", "value")
  )
  expect_named(
    info[["2_metadata"]][["columns"]],
    c("column", "description")
  )
  expect_named(
    info[["2_metadata"]][["data"]],
    c("Sample_ID", "Barcode_ID", "Barcode_ID_1", ".extra_4")
  )
  expect_equal(info[["2_metadata"]][["data"]]$.extra_4[[1L]], "extra value")
  expect_named(info[["3_gel"]], "Example")
  expect_equal(
    unname(as.character(info[["3_gel"]][["Example"]][1L, ])),
    c("Gel_lane", "Gel_Quality", "Clear_Band")
  )
  expect_equal(attr(info, "source"), normalizePath(workbook))
})

test_that("read_sample_info supports sheet names and positions", {
  skip_if_not_installed("openxlsx")
  workbook <- tempfile(fileext = ".xlsx")
  make_sectioned_workbook(workbook)

  by_name <- suppressWarnings(read_sample_info(workbook, sheets = "3_gel"))
  by_position <- suppressWarnings(read_sample_info(workbook, sheets = 2L))

  expect_named(by_name, "3_gel")
  expect_identical(by_name, by_position)
  expect_error(
    read_sample_info(workbook, sheets = "missing"),
    "does not contain"
  )
})

test_that("read_sample_info strict mode validates workbook structure", {
  skip_if_not_installed("openxlsx")

  unmarked <- tempfile(fileext = ".xlsx")
  openxlsx::write.xlsx(data.frame(value = "plain content"), unmarked,
                       colNames = FALSE)
  expect_error(read_sample_info(unmarked), "no \\[name\\] block markers")
  expect_warning(
    relaxed <- read_sample_info(unmarked, strict = FALSE),
    "returning one `unmarked` block"
  )
  expect_named(relaxed[[1L]], "unmarked")

  preamble_cells <- matrix(
    c("intro", NA, "[data]", NA, "id", "value", "S1", "one"),
    ncol = 2L,
    byrow = TRUE
  )
  preamble <- tempfile(fileext = ".xlsx")
  openxlsx::write.xlsx(preamble_cells, preamble, colNames = FALSE)
  expect_error(read_sample_info(preamble), "before its first block marker")
  expect_warning(
    relaxed <- read_sample_info(preamble, strict = FALSE),
    "returning it as `preamble`"
  )
  expect_named(relaxed[[1L]], c("preamble", "data"))
})

test_that("read_sample_info can extract embedded workbook media", {
  skip_if_not_installed("openxlsx")
  skip_if_not(capabilities("png"))

  image <- tempfile(fileext = ".png")
  grDevices::png(image, width = 400L, height = 400L)
  graphics::plot.new()
  grDevices::dev.off()

  workbook <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "image_sheet")
  openxlsx::writeData(wb, "image_sheet", "[Example]", colNames = FALSE)
  openxlsx::insertImage(wb, "image_sheet", image, startCol = 2L, startRow = 1L)
  openxlsx::saveWorkbook(wb, workbook, overwrite = TRUE)

  expect_warning(read_sample_info(workbook), "embedded image file")
  info <- read_sample_info(workbook, include_images = TRUE)
  expect_true(length(attr(info, "images")) >= 1L)
  expect_true(all(file.exists(attr(info, "images"))))
  expect_true(dir.exists(attr(info, "image_dir")))
})
