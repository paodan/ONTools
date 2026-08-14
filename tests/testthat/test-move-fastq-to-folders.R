test_that("move_fastq_to_folders plans moves with project and amplicon size", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(file.path(fastq_dir, "barcode01"), recursive = TRUE)
  dir.create(file.path(fastq_dir, "barcode02"), recursive = TRUE)
  sample_info <- data.frame(
    barcode = c("barcode01", "barcode02"),
    Project = c("Project A", "Project B"),
    Expected_Amplicon_Size_bp = c("3500bp", "1600 bp")
  )

  plan <- move_fastq_to_folders(
    fastq_dir = fastq_dir,
    sample_info = sample_info,
    project_col = "Project",
    dry_run = TRUE
  )

  expect_equal(plan$group, c("Project_A_3500bp", "Project_B_1600_bp"))
  expect_equal(plan$status, c("dry_run", "dry_run"))
  expect_true(all(plan$source_exists))
  expect_false(dir.exists(file.path(fastq_dir, "Project_A_3500bp", "barcode01")))
})

test_that("move_fastq_to_folders moves barcode folders", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(file.path(fastq_dir, "barcode01"), recursive = TRUE)
  writeLines("x", file.path(fastq_dir, "barcode01", "reads.fastq.gz"))
  sample_info <- data.frame(
    barcode = "barcode01",
    Expected_Amplicon_Size_bp = "3500bp"
  )

  plan <- move_fastq_to_folders(
    fastq_dir = fastq_dir,
    sample_info = sample_info
  )

  expect_true(plan$moved[[1]])
  expect_equal(plan$status[[1]], "moved")
  expect_false(dir.exists(file.path(fastq_dir, "barcode01")))
  expect_true(dir.exists(file.path(fastq_dir, "3500bp", "barcode01")))
  expect_true(file.exists(file.path(fastq_dir, "3500bp", "barcode01", "reads.fastq.gz")))
})

test_that("move_fastq_to_folders reads sample info from CSV", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(file.path(fastq_dir, "barcode01"), recursive = TRUE)
  sample_csv <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      barcode = "barcode01",
      Expected_Amplicon_Size_bp = "3500bp"
    ),
    sample_csv,
    row.names = FALSE
  )

  plan <- move_fastq_to_folders(
    fastq_dir = fastq_dir,
    sample_info = sample_csv,
    dry_run = TRUE
  )

  expect_equal(plan$barcode, "barcode01")
  expect_equal(plan$group, "3500bp")
})

test_that("move_fastq_to_folders reports missing barcode folders", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(fastq_dir)
  sample_info <- data.frame(
    barcode = "barcode01",
    Expected_Amplicon_Size_bp = "3500bp"
  )

  expect_warning(
    plan <- move_fastq_to_folders(
      fastq_dir = fastq_dir,
      sample_info = sample_info,
      dry_run = TRUE
    ),
    "not found"
  )
  expect_false(plan$source_exists[[1]])
  expect_equal(plan$status[[1]], "missing_source")
})

test_that("move_fastq_to_folders protects existing destinations", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(file.path(fastq_dir, "barcode01"), recursive = TRUE)
  dir.create(file.path(fastq_dir, "3500bp", "barcode01"), recursive = TRUE)
  sample_info <- data.frame(
    barcode = "barcode01",
    Expected_Amplicon_Size_bp = "3500bp"
  )

  expect_error(
    move_fastq_to_folders(
      fastq_dir = fastq_dir,
      sample_info = sample_info,
      dry_run = TRUE
    ),
    "already exist"
  )
})

test_that("move_fastq_to_folders validates required columns", {
  fastq_dir <- tempfile("fastq-pass-trim-")
  dir.create(fastq_dir)
  sample_info <- data.frame(barcode = "barcode01")

  expect_error(
    move_fastq_to_folders(
      fastq_dir = fastq_dir,
      sample_info = sample_info,
      dry_run = TRUE
    ),
    "missing required column"
  )
})
