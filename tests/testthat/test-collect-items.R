test_that("collect_items vectorizes mixed files and directories", {
  source_root <- tempfile("collect-items-source-")
  destination_root <- tempfile("collect-items-destination-")
  dir.create(file.path(source_root, "barcode01"), recursive = TRUE)
  writeLines("file", file.path(source_root, "report.txt"))
  writeLines(
    "nested",
    file.path(source_root, "barcode01", "alignment.txt")
  )

  result <- collect_items(
    source = c(
      file.path(source_root, "report.txt"),
      file.path(source_root, "barcode01")
    ),
    destination_dir = destination_root
  )

  expect_equal(result$type, c("file", "directory"))
  expect_true(all(result$exists))
  expect_true(all(result$copied))
  expect_equal(
    result$destination,
    file.path(destination_root, c("report.txt", "barcode01"))
  )
  expect_true(file.exists(file.path(destination_root, "report.txt")))
  expect_true(file.exists(file.path(
    destination_root,
    "barcode01",
    "alignment.txt"
  )))
})

test_that("collect_items supports one destination parent per source", {
  source_root <- tempfile("collect-items-vector-parent-")
  destination_root <- tempfile("collect-items-vector-output-")
  dir.create(source_root)
  sources <- file.path(source_root, c("a.txt", "b.txt"))
  writeLines("a", sources[[1L]])
  writeLines("b", sources[[2L]])

  result <- collect_items(
    source = sources,
    destination_dir = file.path(destination_root, c("one", "two")),
    label = c("first", "second")
  )

  expect_equal(result$label, c("first", "second"))
  expect_true(file.exists(file.path(destination_root, "one", "a.txt")))
  expect_true(file.exists(file.path(destination_root, "two", "b.txt")))
})

test_that("collect_items reports missing sources without stopping", {
  source_root <- tempfile("collect-items-missing-source-")
  destination_root <- tempfile("collect-items-missing-output-")
  dir.create(source_root)
  existing <- file.path(source_root, "present.txt")
  missing <- file.path(source_root, "missing.txt")
  writeLines("present", existing)

  result <- collect_items(
    source = c(existing, missing),
    destination_dir = destination_root
  )

  expect_equal(result$type[[1L]], "file")
  expect_true(is.na(result$type[[2L]]))
  expect_equal(result$exists, c(TRUE, FALSE))
  expect_equal(result$copied, c(TRUE, FALSE))
  expect_false(file.exists(file.path(destination_root, "missing.txt")))
})

test_that("collect_items honors overwrite", {
  source_root <- tempfile("collect-items-overwrite-source-")
  destination_root <- tempfile("collect-items-overwrite-output-")
  dir.create(source_root)
  dir.create(destination_root)
  source <- file.path(source_root, "report.txt")
  destination <- file.path(destination_root, "report.txt")
  writeLines("new", source)
  writeLines("old", destination)

  skipped <- collect_items(source, destination_root, overwrite = FALSE)
  expect_false(skipped$copied)
  expect_equal(readLines(destination), "old")

  replaced <- collect_items(source, destination_root, overwrite = TRUE)
  expect_true(replaced$copied)
  expect_equal(readLines(destination), "new")
})

test_that("collect_items preserves file modification time", {
  source_root <- tempfile("collect-items-time-source-")
  destination_root <- tempfile("collect-items-time-output-")
  dir.create(source_root)
  source <- file.path(source_root, "dated.txt")
  writeLines("dated", source)
  expected_time <- as.POSIXct("2020-01-02 03:04:05", tz = "UTC")
  Sys.setFileTime(source, expected_time)

  result <- collect_items(source, destination_root)
  actual_time <- file.info(result$destination)$mtime

  expect_true(result$copied)
  expect_equal(as.numeric(actual_time), as.numeric(file.info(source)$mtime))
})

test_that("collect_items rejects ambiguous or unsafe copy plans", {
  root <- tempfile("collect-items-plan-")
  dir.create(file.path(root, "one"), recursive = TRUE)
  dir.create(file.path(root, "two"), recursive = TRUE)
  first <- file.path(root, "one", "same.txt")
  second <- file.path(root, "two", "same.txt")
  writeLines("one", first)
  writeLines("two", second)

  expect_error(
    collect_items(c(first, second), file.path(root, "output")),
    "same destination"
  )
  expect_error(
    collect_items(first, dirname(first)),
    "must not be the same path"
  )
  expect_error(
    collect_items(root, file.path(root, "inside")),
    "cannot be copied into itself"
  )
})

test_that("collect_items validates vector lengths and arguments", {
  expect_error(collect_items(character(), "output"), "source")
  expect_error(collect_items(c("a", "b"), c("x", "y", "z")), "length")
  expect_error(collect_items("a", "output", overwrite = NA), "overwrite")
  expect_error(collect_items("a", "output", label = character()), "label")
})
