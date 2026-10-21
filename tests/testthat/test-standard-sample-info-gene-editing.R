test_that("standard_sample_info_gene_editing creates canonical columns and paths", {
  sample_info <- data.frame(
    Sample_HiTOM = c(" sample_A ", "sample_A", "sample_B"),
    barcode = c("barcode001", "barcode001", "barcode002"),
    Reference_name = c("amp1", "amp2", "amp1"),
    Reference_Region_Start = c("10", "20", "30"),
    Reference_Region_End = c(15, 25, 35),
    ignored = 1:3,
    stringsAsFactors = FALSE
  )

  result <- standard_sample_info_gene_editing(
    sample_info = sample_info,
    input_path_root = "variant_results",
    delivery_path_root = "delivery"
  )

  expect_equal(
    names(result),
    c(
      "sample_id", "barcode", "ref_name", "region_start", "region_end",
      "bam_file_origin", "bam_file", "output_fastq", "reference_file",
      "ab1_file"
    )
  )
  expect_equal(result$sample_id, c("sample_A", "sample_A", "sample_B"))
  expect_type(result$region_start, "integer")
  expect_type(result$region_end, "integer")
  expect_equal(
    result$bam_file_origin[[1L]],
    file.path(
      "variant_results", "barcode001", "alignments",
      "barcode001.aligned.sorted.bam"
    )
  )
  expect_equal(
    result$bam_file[[1L]],
    file.path("variant_results", "work", "Bam", "sample_A.aligned.sorted.bam")
  )
  expect_equal(
    result$output_fastq[[1L]],
    file.path("variant_results", "work", "0_fastq", "sample_A_filtered.fastq.gz")
  )
  expect_equal(
    unique(result$reference_file),
    file.path("variant_results", "reference_sanitized_seqIDs.fasta")
  )
  expect_equal(
    result$ab1_file[[1L]],
    file.path("variant_results", "work", "Ab1", "sample_A_filtered.ab1")
  )
})

test_that("standard_sample_info_gene_editing accepts custom source columns", {
  sample_info <- data.frame(
    specimen = "S01",
    barcode_id = "barcode010",
    target = "amp10",
    edit_start = 100,
    edit_end = 150
  )

  result <- standard_sample_info_gene_editing(
    sample_info = sample_info,
    input_path_root = "input",
    delivery_path_root = "delivery",
    work_path = "scratch/project",
    sample_id_column = "specimen",
    barcode_column = "barcode_id",
    ref_name_column = "target",
    ref_region_start_column = "edit_start",
    ref_region_end_column = "edit_end",
    ref_sanitized_file = "custom_reference.fasta"
  )

  expect_equal(result$sample_id, "S01")
  expect_equal(result$barcode, "barcode010")
  expect_equal(result$ref_name, "amp10")
  expect_equal(
    result$bam_file,
    file.path("scratch/project", "Bam", "S01.aligned.sorted.bam")
  )
  expect_equal(
    result$reference_file,
    file.path("input", "custom_reference.fasta")
  )
})

test_that("standard_sample_info_gene_editing reports missing and repeated columns", {
  incomplete <- data.frame(Sample_HiTOM = "S01", barcode = "barcode001")
  expect_error(
    standard_sample_info_gene_editing(incomplete, "input", "delivery"),
    "Reference_name.*Reference_Region_Start.*Reference_Region_End"
  )

  complete <- data.frame(
    Sample_HiTOM = "S01",
    barcode = "barcode001",
    Reference_name = "amp1",
    Reference_Region_Start = 1,
    Reference_Region_End = 2
  )
  expect_error(
    standard_sample_info_gene_editing(
      complete, "input", "delivery",
      sample_id_column = "barcode"
    ),
    "distinct columns"
  )
})

test_that("standard_sample_info_gene_editing validates values and intervals", {
  make_info <- function(sample_id = "S01", barcode = "barcode001",
                        start = 10, end = 20) {
    data.frame(
      Sample_HiTOM = sample_id,
      barcode = barcode,
      Reference_name = "amp1",
      Reference_Region_Start = start,
      Reference_Region_End = end
    )
  }

  expect_error(
    standard_sample_info_gene_editing(make_info(start = 10.5), "input", "delivery"),
    "positive whole numbers"
  )
  expect_error(
    standard_sample_info_gene_editing(make_info(start = 21, end = 20), "input", "delivery"),
    "start must be less"
  )
  expect_error(
    standard_sample_info_gene_editing(make_info(sample_id = "../S01"), "input", "delivery"),
    "sample_id.*unsafe path"
  )
  expect_error(
    standard_sample_info_gene_editing(make_info(barcode = "../barcode001"), "input", "delivery"),
    "barcode.*unsafe path"
  )
  expect_error(
    standard_sample_info_gene_editing(make_info(barcode = " "), "input", "delivery"),
    "barcode.*missing or empty"
  )
})

test_that("standard_sample_info_gene_editing prevents sample output collisions", {
  sample_info <- data.frame(
    Sample_HiTOM = c("S01", "S01"),
    barcode = c("barcode001", "barcode002"),
    Reference_name = c("amp1", "amp2"),
    Reference_Region_Start = c(10, 20),
    Reference_Region_End = c(15, 25)
  )

  expect_error(
    standard_sample_info_gene_editing(sample_info, "input", "delivery"),
    "exactly one barcode.*S01"
  )

  sample_info$Sample_HiTOM <- c("S01", "S02")
  sample_info$barcode <- c("barcode001", "barcode001")
  expect_error(
    standard_sample_info_gene_editing(sample_info, "input", "delivery"),
    "exactly one `sample_id`.*barcode001"
  )
})

test_that("standard_sample_info_gene_editing validates the input object", {
  expect_error(
    standard_sample_info_gene_editing(list(), "input", "delivery"),
    "must be a data frame"
  )
  expect_error(
    standard_sample_info_gene_editing(data.frame(), "input", "delivery"),
    "at least one row"
  )
})
