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
  expect_equal(
    res$delivery$PROJECT001_1600$destination,
    file.path(normalizePath(delivery_dir), "PROJECT001_1600", "consensus_results")
  )
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
    make_ab1 = FALSE,
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
    make_ab1 = FALSE,
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
      make_ab1 = "yes",
      dry_run = TRUE,
      echo = FALSE
    ),
    "make_ab1"
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

test_that("make_consensus_delivery writes AB1 files to delivered barcode folders", {
  skip_if(Sys.which("samtools") == "", "samtools is not installed")

  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("amplicon-delivery-")
  run_root <- file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01"
  )
  fastq_root <- file.path(run_root, "fastq_pass_trim")
  result_dir <- file.path(run_root, "results", "wf_amplicon_denovo", "PROJECT001_1600")
  dir.create(file.path(fastq_root, "PROJECT001_1600", "barcode001"), recursive = TRUE)
  dir.create(file.path(run_root, "bam_pass", "barcode001"), recursive = TRUE)
  dir.create(file.path(result_dir, "barcode001", "alignments"), recursive = TRUE)
  dir.create(file.path(result_dir, "barcode001", "consensus"), recursive = TRUE)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = "PBC001-001",
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 8,
      Max_Read_Length = 20
    ),
    sample_info,
    row.names = FALSE
  )

  sequence <- "ACGTACGTACGT"
  consensus <- file.path(result_dir, "combined-consensus.fasta")
  sam <- file.path(result_dir, "reads.sam")
  bam <- file.path(result_dir, "barcode001", "alignments", "barcode001.aligned.sorted.bam")
  writeLines(c(">barcode001", sequence), consensus)
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      paste0("@SQ\tSN:barcode001\tLN:", nchar(sequence)),
      paste0("read1\t0\tbarcode001\t1\t60\t12M\t*\t0\t0\t", sequence, "\tIIIIIIIIIIII"),
      "read2\t0\tbarcode001\t1\t60\t12M\t*\t0\t0\tACGTGCGTACGT\tIIIIIIIIIIII"
    ),
    sam
  )
  system2("samtools", c("faidx", consensus))
  system2("samtools", c("view", "-bS", sam), stdout = bam)
  system2("samtools", c("index", bam))
  writeLines("fastq", file.path(result_dir, "barcode001", "consensus", "consensus.fastq"))

  res <- make_consensus_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = FALSE,
    run_amplicon_step = FALSE,
    trim_consensus_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    collect_results_step = TRUE,
    consensus_file = "combined-consensus.fasta",
    consensus_index_file = "combined-consensus.fasta.fai",
    ab1_name_template = "{barcode}.trace.ab1",
    readme_name = "README.custom.txt",
    chinese_readme_name = NULL,
    echo = FALSE,
    stderr = FALSE
  )

  delivered_ab1 <- file.path(
    delivery_dir,
    "PROJECT001_1600",
    "consensus_results",
    "barcode001",
    "barcode001.trace.ab1"
  )
  source_ab1 <- file.path(result_dir, "barcode001", "barcode001.trace.ab1")

  expect_true(file.exists(delivered_ab1))
  expect_true(file.exists(source_ab1))
  expect_equal(res$ab1$PROJECT001_1600$status, "generated")
  expect_equal(res$ab1$PROJECT001_1600$barcode, "barcode001")
  expect_equal(res$ab1$PROJECT001_1600$ab1, normalizePath(source_ab1))
  expect_true(any(grepl(
    "barcode*/<barcode>.trace.ab1",
    readLines(file.path(
      delivery_dir,
      "PROJECT001_1600",
      "consensus_results",
      "README.custom.txt"
    )),
    fixed = TRUE
  )))
  expect_false(file.exists(file.path(
    delivery_dir,
    "PROJECT001_1600",
    "consensus_results",
    "README.zh-CN.txt"
  )))
})
