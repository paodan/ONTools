test_that("standard_paths returns the documented standard layout", {
  paths <- standard_paths(
    input_path_root = "analysis/run01",
    delivery_path_root = "delivery",
    delivery_folder = "project01"
  )

  expect_equal(paths$input_path_root, "analysis/run01")
  expect_equal(paths$work_path, file.path("analysis/run01", "work"))
  expect_equal(paths$work_bam_path, file.path("analysis/run01", "work", "Bam"))
  expect_equal(paths$work_fastq_path, file.path("analysis/run01", "work", "0_fastq"))
  expect_equal(paths$work_qc_path, file.path("analysis/run01", "work", "QC"))
  expect_equal(paths$work_ab1_path, file.path("analysis/run01", "work", "Ab1"))
  expect_equal(paths$delivery_path, file.path("delivery", "project01"))
  expect_equal(paths$delivery_bam_path, file.path("delivery", "project01", "Bam"))
  expect_equal(paths$delivery_fastq_path, file.path("delivery", "project01", "0_fastq"))
  expect_equal(paths$delivery_qc_path, file.path("delivery", "project01", "QC"))
  expect_equal(paths$delivery_ab1_path, file.path("delivery", "project01", "Ab1"))
})

test_that("standard_paths accepts an explicit work root and named extensions", {
  paths <- standard_paths(
    input_path_root = "analysis/run01",
    delivery_path_root = "delivery",
    delivery_folder = "project01",
    work_path = "scratch/run01",
    report_file = "scratch/run01/report.xlsx"
  )

  expect_equal(paths$work_path, "scratch/run01")
  expect_equal(paths$work_qc_path, file.path("scratch/run01", "QC"))
  expect_equal(paths$report_file, "scratch/run01/report.xlsx")
})

test_that("standard_paths validates standard and additional paths", {
  expect_error(
    standard_paths(NA_character_, "delivery", "project01"),
    "input_path_root"
  )
  expect_error(
    standard_paths("input", "delivery", "project01", "work", "x"),
    "must be named"
  )
  expect_error(
    standard_paths(
      "input", "delivery", "project01",
      work_bam_path = "another/Bam"
    ),
    "duplicate standard path"
  )
  expect_error(
    standard_paths("input", "delivery", "project01", report = character()),
    "report"
  )
})

test_that("standard_paths_gene_editing extends the base layout", {
  paths <- standard_paths_gene_editing(
    input_path_root = "variant_results",
    delivery_path_root = "delivery"
  )

  expect_equal(paths$delivery_folder, "gene_editing_delivery")
  expect_equal(
    paths$ref_sanitized,
    file.path("variant_results", "reference_sanitized_seqID.fasta")
  )
  expect_equal(
    paths$xlsx_file_full,
    file.path("variant_results", "work", "variant_table_full.xlsx")
  )
  expect_equal(
    paths$xlsx_file_filtered,
    file.path("variant_results", "work", "variant_table_filtered.xlsx")
  )
})

test_that("standard_paths_gene_editing supports custom locations and filenames", {
  paths <- standard_paths_gene_editing(
    input_path_root = "variant_results",
    delivery_path_root = "delivery",
    delivery_folder = "project01",
    work_path = "scratch/project01",
    ref_sanitized_file = "reference.fa",
    xlsx_file_full = "all.xlsx",
    xlsx_file_filtered = "filtered.xlsx"
  )

  expect_equal(paths$ref_sanitized, file.path("variant_results", "reference.fa"))
  expect_equal(paths$xlsx_file_full, file.path("scratch/project01", "all.xlsx"))
  expect_equal(
    paths$xlsx_file_filtered,
    file.path("scratch/project01", "filtered.xlsx")
  )
})
