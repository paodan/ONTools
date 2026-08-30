test_that("make_consensus_delivery builds dry-run plans and commands", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("amplicon-delivery-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = c("PBC001-001", "PBC001-002"),
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 1200,
      Max_Read_Length = 1800,
      Primer_F = "ACGT",
      Primer_R = "TGCA"
    ),
    sample_info,
    row.names = FALSE
  )

  res <- make_consensus_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    path_delivery = delivery_dir,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$stat$status, NA_integer_)
  expect_true("PROJECT001_1600" %in% names(res$workflow))
  expect_match(res$workflow$PROJECT001_1600$command_string, "wf-amplicon")
  expect_match(res$workflow$PROJECT001_1600$command_string, "--min_read_length")
  expect_match(res$workflow$PROJECT001_1600$command_string, "1200")
  expect_equal(res$move_plans$PROJECT001_1600$status, c("dry_run", "dry_run"))
  expect_equal(res$move_plans$PROJECT001_1600$barcode, c("barcode001", "barcode002"))
  expect_equal(res$trim$PROJECT001_1600$status, "dry_run")
  expect_equal(res$delivery$PROJECT001_1600$destination, file.path(normalizePath(delivery_dir), "PROJECT001_1600"))
})

test_that("make_consensus_delivery can reconstruct skipped demux paths", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("amplicon-delivery-")
  fastq_dir <- file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01",
    "fastq_pass_trim"
  )
  dir.create(file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01",
    "bam_pass",
    "barcode001"
  ), recursive = TRUE)
  dir.create(file.path(fastq_dir, "barcode001"), recursive = TRUE)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = "PBC001-001",
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 1200,
      Max_Read_Length = 1800
    ),
    sample_info,
    row.names = FALSE
  )

  res <- make_consensus_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    collect_results_step = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(
    res$workflow$PROJECT001_1600$paths$fastq,
    file.path(normalizePath(fastq_dir), "PROJECT001_1600")
  )
  expect_match(res$workflow$PROJECT001_1600$paths$out_dir, "wf_amplicon_denovo")
})

test_that("make_consensus_delivery auto mode reuses grouped FASTQ folders", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("amplicon-delivery-")
  fastq_dir <- file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01",
    "fastq_pass_trim"
  )
  dir.create(file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01",
    "bam_pass",
    "barcode001"
  ), recursive = TRUE)
  dir.create(file.path(fastq_dir, "PROJECT001_1600", "barcode001"), recursive = TRUE)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = "PBC001-001",
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 1200,
      Max_Read_Length = 1800
    ),
    sample_info,
    row.names = FALSE
  )

  res <- make_consensus_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = TRUE,
    move_fastq_mode = "auto",
    run_amplicon_step = FALSE,
    trim_consensus_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    collect_results_step = FALSE,
    echo = FALSE
  )

  expect_equal(res$move_plans$PROJECT001_1600$status, "reused")
  expect_equal(
    res$workflow$PROJECT001_1600$paths$fastq,
    file.path(normalizePath(fastq_dir), "PROJECT001_1600")
  )
})

test_that("make_consensus_delivery validates sample information inputs", {
  proj <- tempfile("ont-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = "PBC001-001",
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 1200
    ),
    sample_info,
    row.names = FALSE
  )

  expect_error(
    make_consensus_delivery(
      path_proj = proj,
      path_sampleInfo_file_list = sample_info,
      path_delivery = tempfile("amplicon-delivery-"),
      dry_run = TRUE,
      echo = FALSE
    ),
    "must be named"
  )

  expect_error(
    make_consensus_delivery(
      path_proj = proj,
      path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
      path_delivery = tempfile("amplicon-delivery-"),
      dry_run = TRUE,
      echo = FALSE
    ),
    "missing required column"
  )
})
