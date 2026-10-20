test_that("make_variant_delivery builds dry-run variant workflow plans", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("variant-delivery-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGTACGTACGTACGT"), reference)

  sample_info <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Barcode_ID = c("PBC001-001", "PBC001-002"),
      Project_ID = "PROJECT001",
      Expected_Size_bp = "1600",
      Min_Read_Length = 1200,
      Max_Read_Length = 1800
    ),
    sample_info,
    row.names = FALSE
  )

  res <- make_variant_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    reference = reference,
    path_delivery = delivery_dir,
    amplicon_extra_args = "--threads 20",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$stat$status, NA_integer_)
  expect_true("PROJECT001_1600" %in% names(res$workflow))
  expect_match(res$workflow$PROJECT001_1600$command_string, "wf-amplicon")
  expect_match(res$workflow$PROJECT001_1600$command_string, "--reference")
  expect_match(
    res$workflow$PROJECT001_1600$command_string,
    normalizePath(reference),
    fixed = TRUE
  )
  expect_match(res$workflow$PROJECT001_1600$command_string, "--threads 20", fixed = TRUE)
  expect_equal(res$variant_tables$PROJECT001_1600$status, "dry_run")
  expect_equal(
    res$delivery$PROJECT001_1600$destination,
    file.path(normalizePath(delivery_dir), "PROJECT001_1600", "variant_results")
  )
})

test_that("make_variant_delivery skips filtered-read QC in dry-run mode", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("variant-delivery-")
  run_root <- file.path(proj, "demux_out_YS-NB576", "run01")
  fastq_root <- file.path(run_root, "fastq_pass_trim")
  result_dir <- file.path(run_root, "results", "wf_amplicon_variant", "PROJECT001_1600")

  dir.create(file.path(run_root, "bam_pass", "barcode001"), recursive = TRUE)
  dir.create(file.path(fastq_root, "PROJECT001_1600", "barcode001"), recursive = TRUE)
  dir.create(result_dir, recursive = TRUE)
  writeLines("placeholder sequencing summary", file.path(proj, "demux_out_YS-NB576", "sequencing_summary.txt"))

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGTACGTACGTACGT"), reference)

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

  res <- make_variant_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    reference = reference,
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = TRUE,
    move_fastq_mode = "reuse",
    run_amplicon_step = FALSE,
    make_variant_table_step = FALSE,
    run_filtered_QC_step = TRUE,
    run_igv_step = FALSE,
    make_ab1 = FALSE,
    collect_results_step = FALSE,
    dry_run = TRUE,
    echo = FALSE,
    stderr = FALSE
  )

  expect_null(res$g2$PROJECT001_1600)
  expect_false(file.exists(file.path(result_dir, "numRead_perSample.png")))
  expect_false(file.exists(file.path(result_dir, "lenReadDistribution.png")))
  expect_false(file.exists(file.path(result_dir, "numRead_perSample.csv")))
})

test_that("make_variant_delivery writes per-barcode and merged variant tables", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("variant-delivery-")
  run_root <- file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01"
  )
  fastq_root <- file.path(run_root, "fastq_pass_trim")
  result_dir <- file.path(run_root, "results", "wf_amplicon_variant", "PROJECT001_1600")
  barcode_dir <- file.path(result_dir, "barcode001")
  vcf_dir <- file.path(barcode_dir, "variants")

  dir.create(file.path(run_root, "bam_pass", "barcode001"), recursive = TRUE)
  dir.create(file.path(fastq_root, "PROJECT001_1600", "barcode001"), recursive = TRUE)
  dir.create(vcf_dir, recursive = TRUE)

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp 1", "ACGTACGTACGTACGTACGTACGTACGT"), reference)
  sanitized_reference <- file.path(result_dir, "reference_sanitized_seqID.fasta")
  writeLines(c(">amp_1", "ACGTACGTACGTACGTACGTACGTACGT"), sanitized_reference)

  vcf <- file.path(vcf_dir, "medaka.annotated.vcf.gz")
  con <- gzfile(vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode001",
    "amp_1\t5\t.\tA\tG\t60\tPASS\tDP=30;SR=1,2,3,4;AR=0,1\tGT:GQ\t1:20"
  ), con)
  close(con)
  on.exit(NULL)

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

  res <- make_variant_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    reference = reference,
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = TRUE,
    move_fastq_mode = "reuse",
    run_amplicon_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    make_ab1 = FALSE,
    collect_results_step = TRUE,
    echo = FALSE,
    stderr = FALSE
  )

  selected_cols <- c(
    "CHROM", "reference_name_original", "reference_name_sanitized",
    "POS", "REF", "ALT", "QUAL", "FILTER", "Type",
    "ref_depth_downsampled", "alt_depth_downsampled",
    "variant_percent", "sampleID", "ref_upstream_20bp", "ref_downstream_20bp"
  )
  barcode_tsv <- file.path(barcode_dir, "variant.tsv")
  merged_tsv <- file.path(result_dir, "variant_all.tsv")
  delivered_tsv <- file.path(
    delivery_dir,
    "PROJECT001_1600",
    "variant_results",
    "variant_all.tsv"
  )
  delivered_readme <- file.path(
    delivery_dir,
    "PROJECT001_1600",
    "variant_results",
    "README.txt"
  )
  delivered_readme_zh <- file.path(
    delivery_dir,
    "PROJECT001_1600",
    "variant_results",
    "README.zh-CN.txt"
  )

  expect_true(file.exists(barcode_tsv))
  expect_true(file.exists(merged_tsv))
  expect_true(file.exists(delivered_tsv))
  expect_true(file.exists(delivered_readme))
  expect_true(file.exists(delivered_readme_zh))
  expect_equal(res$variant_tables$PROJECT001_1600$status, "written")
  expect_equal(res$variant_tables$PROJECT001_1600$files$status, "written")
  expect_equal(
    res$variant_references$PROJECT001_1600,
    normalizePath(sanitized_reference)
  )

  variant <- utils::read.delim(barcode_tsv, check.names = FALSE)
  merged <- utils::read.delim(merged_tsv, check.names = FALSE)
  expect_equal(names(variant), selected_cols)
  expect_equal(names(merged), selected_cols)
  expect_equal(variant$CHROM, "amp_1")
  expect_equal(variant$reference_name_original, "amp 1")
  expect_equal(variant$reference_name_sanitized, "amp_1")
  expect_equal(variant$Type, "SNP")
  expect_equal(variant$ref_depth_downsampled, 3)
  expect_equal(variant$alt_depth_downsampled, 7)
  expect_equal(variant$variant_percent, 70)
  expect_equal(variant$sampleID, "barcode001")
  expect_equal(merged, variant)
  expect_true(any(grepl("Variant table columns", readLines(delivered_readme), fixed = TRUE)))
  expect_true(any(grepl("variant_percent", readLines(delivered_readme), fixed = TRUE)))
  expect_true(any(grepl("变异表字段说明", readLines(delivered_readme_zh), fixed = TRUE)))
  expect_true(any(grepl("variant_percent", readLines(delivered_readme_zh), fixed = TRUE)))
})

test_that("make_variant_delivery uses workflow-generated reference matching VCF CHROM values", {
  proj <- tempfile("ont-project-")
  delivery_dir <- tempfile("variant-delivery-")
  run_root <- file.path(
    proj,
    "demux_out_YS-NB576",
    "run01",
    "sample01",
    "readset01"
  )
  fastq_root <- file.path(run_root, "fastq_pass_trim")
  result_dir <- file.path(run_root, "results", "wf_amplicon_variant", "PROJECT001_1600")
  barcode_dir <- file.path(result_dir, "barcode001")
  vcf_dir <- file.path(barcode_dir, "variants")

  dir.create(file.path(run_root, "bam_pass", "barcode001"), recursive = TRUE)
  dir.create(file.path(fastq_root, "PROJECT001_1600", "barcode001"), recursive = TRUE)
  dir.create(vcf_dir, recursive = TRUE)

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">PLA3 A",
    "ACGTACGTACGTACGTACGTACGTACGT",
    ">PLA3 B",
    "TGCATGCATGCATGCATGCATGCATGCA"
  ), reference)
  workflow_reference <- file.path(result_dir, "workflow_final_reference.fasta")
  writeLines(c(
    ">PLA3_A_",
    "ACGTACGTACGTACGTACGTACGTACGT",
    ">PLA3_B_",
    "TGCATGCATGCATGCATGCATGCATGCA"
  ), workflow_reference)

  vcf <- file.path(vcf_dir, "medaka.annotated.vcf.gz")
  con <- gzfile(vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode001",
    "PLA3_A_\t5\t.\tA\tG\t60\tPASS\tDP=30;SR=1,2,3,4;AR=0,1\tGT:GQ\t1:20",
    "PLA3_B_\t6\t.\tG\tT\t60\tPASS\tDP=40;SR=2,3,4,5;AR=0,1\tGT:GQ\t1:20"
  ), con)
  close(con)
  on.exit(NULL)

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

  res <- make_variant_delivery(
    path_proj = proj,
    path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
    reference = reference,
    path_delivery = delivery_dir,
    run_basecalling_demux_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = TRUE,
    move_fastq_mode = "reuse",
    run_amplicon_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    make_ab1 = FALSE,
    collect_results_step = FALSE,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(
    res$variant_references$PROJECT001_1600,
    normalizePath(workflow_reference)
  )
  variant <- utils::read.delim(
    file.path(barcode_dir, "variant.tsv"),
    check.names = FALSE
  )
  expect_equal(variant$CHROM, c("PLA3_A_", "PLA3_B_"))
  expect_equal(variant$reference_name_original, c("PLA3 A", "PLA3 B"))
})

test_that("make_variant_delivery writes empty variant tables for empty VCFs", {
  result_dir <- tempfile("wf-amplicon-result-")
  vcf_dir <- file.path(result_dir, "barcode001", "variants")
  dir.create(vcf_dir, recursive = TRUE)

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGTACGT"), reference)

  vcf <- file.path(vcf_dir, "medaka.annotated.vcf.gz")
  con <- gzfile(vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode001"
  ), con)
  close(con)
  on.exit(NULL)

  selected_cols <- c(
    "CHROM", "reference_name_original", "reference_name_sanitized",
    "POS", "REF", "ALT", "QUAL", "FILTER", "Type",
    "ref_depth_downsampled", "alt_depth_downsampled",
    "variant_percent", "sampleID", "ref_upstream_20bp", "ref_downstream_20bp"
  )

  res <- make_wf_amplicon_variant_tables(
    result_dir = result_dir,
    reference = reference,
    variant_vcf_name = "medaka.annotated.vcf.gz",
    variant_tsv_name = "variant.tsv",
    variant_all_name = "variant_all.tsv",
    variant_columns = selected_cols,
    flank_width = 20,
    dry_run = FALSE
  )

  variant <- utils::read.delim(
    file.path(result_dir, "barcode001", "variant.tsv"),
    check.names = FALSE
  )
  merged <- utils::read.delim(file.path(result_dir, "variant_all.tsv"), check.names = FALSE)

  expect_equal(res$files$status, "no_variants")
  expect_equal(names(variant), selected_cols)
  expect_equal(names(merged), selected_cols)
  expect_equal(nrow(variant), 0)
  expect_equal(nrow(merged), 0)
})

test_that("make_variant_delivery plans IGV snapshots and AB1 files in dry-run mode", {
  result_dir <- tempfile("wf-amplicon-result-")
  align_dir <- file.path(result_dir, "barcode001", "alignments")
  dir.create(align_dir, recursive = TRUE)
  bam <- file.path(align_dir, "barcode001.aligned.sorted.bam")
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGT", ">amp2", "TGCATGCATGCA"), reference)
  file.create(paste0(reference, ".fai"))

  igv <- generate_variant_igv_snapshots(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    format = "png",
    igv = "igv.sh",
    sort = "base",
    display = "collapse",
    max_panel_height = NULL,
    use_xvfb = FALSE,
    igv_dir_name = "IGV",
    dry_run = TRUE
  )
  ab1 <- generate_variant_ab1_files(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    ab1_name_template = "{barcode}.{reference}.synthetic.ab1",
    ab1_dir_name = "AB1",
    samtools = "samtools",
    overwrite = TRUE,
    dry_run = TRUE,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(nrow(igv$files), 2)
  expect_equal(nrow(ab1), 2)
  expect_equal(igv$files$reference_name, c("amp1", "amp2"))
  expect_equal(ab1$reference_name, c("amp1", "amp2"))
  expect_true(all(grepl("barcode001/IGV/barcode001[.]amp[12][.]png$", igv$files$snapshot)))
  expect_true(all(grepl("barcode001/AB1/barcode001[.]amp[12][.]synthetic[.]ab1$", ab1$ab1)))
  expect_equal(igv$files$status, c("dry_run", "dry_run"))
  expect_equal(igv$workers, 1L)
  expect_equal(ab1$status, c("dry_run", "dry_run"))
  expect_equal(attr(ab1, "workers"), 1L)
})

test_that("variant IGV worker count is capped by available barcodes", {
  detected <- parallel::detectCores(logical = TRUE)
  available <- if (is.na(detected) || detected < 2L) 1L else detected - 1L
  expected_auto <- min(as.integer(available), 3L)
  if (.Platform$OS.type == "windows") expected_auto <- 1L

  expect_equal(resolve_variant_igv_workers(NULL, 3L), expected_auto)
  expect_equal(resolve_variant_igv_workers(8L, 3L),
               if (.Platform$OS.type == "windows") 1L else 3L)
  expect_equal(resolve_variant_igv_workers(1L, 10L), 1L)
  expect_equal(resolve_variant_igv_workers(NULL, 0L), 1L)
  expect_error(resolve_variant_igv_workers(0L, 3L), "positive integer")
})

test_that("generate_variant_igv_snapshots plans barcodes with multiple workers", {
  result_dir <- tempfile("wf-amplicon-result-")
  for (barcode in c("barcode001", "barcode002", "barcode003")) {
    align_dir <- file.path(result_dir, barcode, "alignments")
    dir.create(align_dir, recursive = TRUE)
    bam <- file.path(align_dir, paste0(barcode, ".aligned.sorted.bam"))
    file.create(bam)
    file.create(paste0(bam, ".bai"))
  }

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGT"), reference)

  result <- generate_variant_igv_snapshots(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    format = "png",
    igv = "igv.sh",
    sort = "base",
    display = "collapse",
    max_panel_height = NULL,
    use_xvfb = FALSE,
    igv_dir_name = "IGV",
    dry_run = TRUE,
    igv_workers = 2L
  )

  expect_equal(result$workers,
               if (.Platform$OS.type == "windows") 1L else 2L)
  expect_equal(result$files$barcode,
               c("barcode001", "barcode002", "barcode003"))
  expect_true(all(result$files$status == "dry_run"))
})

test_that("generate_variant_igv_snapshots defaults to primary-only filtering", {
  result_dir <- tempfile("wf-amplicon-result-")
  align_dir <- file.path(result_dir, "barcode001", "alignments")
  dir.create(align_dir, recursive = TRUE)
  bam <- file.path(align_dir, "barcode001.aligned.sorted.bam")
  file.create(bam)
  file.create(paste0(bam, ".bai"))

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGT"), reference)
  file.create(paste0(reference, ".fai"))

  seen <- new.env(parent = emptyenv())
  seen$igv_primary_only <- NULL
  testthat::local_mocked_bindings(
    igv_snapshot = function(out_dir,
                            snapshot_name,
                            igv_primary_only,
                            ...) {
      seen$igv_primary_only <- igv_primary_only
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      file.create(file.path(out_dir, snapshot_name))
      list(status = 0L)
    }
  )

  result <- generate_variant_igv_snapshots(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    format = "png",
    igv = "igv.sh",
    sort = "base",
    display = "collapse",
    max_panel_height = NULL,
    use_xvfb = FALSE,
    igv_dir_name = "IGV",
    dry_run = FALSE,
    igv_workers = 1L
  )

  expect_true(seen$igv_primary_only)
  expect_equal(result$files$status, "generated")
})

test_that("variant AB1 worker count is capped by available barcodes", {
  detected <- parallel::detectCores(logical = TRUE)
  available <- if (is.na(detected) || detected < 2L) 1L else detected - 1L
  expected_auto <- min(as.integer(available), 3L)
  if (.Platform$OS.type == "windows") expected_auto <- 1L

  expect_equal(resolve_variant_ab1_workers(NULL, 3L), expected_auto)
  expect_equal(resolve_variant_ab1_workers(8L, 3L),
               if (.Platform$OS.type == "windows") 1L else 3L)
  expect_equal(resolve_variant_ab1_workers(1L, 10L), 1L)
  expect_equal(resolve_variant_ab1_workers(NULL, 0L), 1L)
  expect_error(resolve_variant_ab1_workers(0L, 3L), "positive integer")
})

test_that("generate_variant_ab1_files plans barcodes with multiple workers", {
  result_dir <- tempfile("wf-amplicon-result-")
  for (barcode in c("barcode001", "barcode002", "barcode003")) {
    align_dir <- file.path(result_dir, barcode, "alignments")
    dir.create(align_dir, recursive = TRUE)
    bam <- file.path(align_dir, paste0(barcode, ".aligned.sorted.bam"))
    file.create(bam)
    file.create(paste0(bam, ".bai"))
  }

  reference <- tempfile(fileext = ".fasta")
  writeLines(c(">amp1", "ACGTACGTACGT"), reference)

  result <- generate_variant_ab1_files(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    ab1_name_template = "{barcode}.{reference}.synthetic.ab1",
    ab1_dir_name = "AB1",
    samtools = "samtools",
    overwrite = TRUE,
    dry_run = TRUE,
    echo = FALSE,
    stderr = FALSE,
    ab1_workers = 2L
  )

  expect_equal(attr(result, "workers"),
               if (.Platform$OS.type == "windows") 1L else 2L)
  expect_equal(result$barcode, c("barcode001", "barcode002", "barcode003"))
  expect_true(all(result$status == "dry_run"))
})

test_that("generate_variant_ab1_files writes AB1 files in barcode workers", {
  skip_if(Sys.which("samtools") == "", "samtools is not installed")
  skip_if(.Platform$OS.type == "windows", "forked workers are unavailable")

  result_dir <- tempfile("wf-amplicon-result-")
  reference <- tempfile(fileext = ".fasta")
  sequence <- "ACGTACGTACGT"
  writeLines(c(">amp1", sequence), reference)

  for (barcode in c("barcode001", "barcode002")) {
    align_dir <- file.path(result_dir, barcode, "alignments")
    dir.create(align_dir, recursive = TRUE)
    sam <- file.path(align_dir, paste0(barcode, ".sam"))
    bam <- file.path(align_dir, paste0(barcode, ".aligned.sorted.bam"))
    writeLines(
      c(
        "@HD\tVN:1.6\tSO:coordinate",
        paste0("@SQ\tSN:amp1\tLN:", nchar(sequence)),
        paste0(
          barcode, "_read\t0\tamp1\t1\t60\t12M\t*\t0\t0\t",
          sequence, "\tIIIIIIIIIIII"
        )
      ),
      sam
    )
    system2("samtools", c("view", "-bS", sam), stdout = bam)
    system2("samtools", c("index", bam))
  }

  result <- generate_variant_ab1_files(
    result_dir = result_dir,
    reference = reference,
    barcode_pattern = "^barcode[0-9]+$",
    ab1_name_template = "{barcode}.{reference}.synthetic.ab1",
    ab1_dir_name = "AB1",
    samtools = "samtools",
    overwrite = TRUE,
    dry_run = FALSE,
    echo = FALSE,
    stderr = FALSE,
    ab1_workers = 2L
  )

  expect_equal(attr(result, "workers"), 2L)
  expect_equal(result$status, c("generated", "generated"))
  expect_true(all(file.exists(result$ab1)))
  expect_false(file.exists(paste0(reference, ".fai")))
})

test_that("collect_variant_results documents IGV snapshots and AB1 files", {
  result_dir <- tempfile("wf-amplicon-result-")
  output_dir <- tempfile("variant-delivery-")
  barcode_dir <- file.path(result_dir, "barcode001")
  dir.create(file.path(barcode_dir, "alignments"), recursive = TRUE)
  dir.create(file.path(barcode_dir, "IGV"), recursive = TRUE)
  dir.create(file.path(barcode_dir, "AB1"), recursive = TRUE)
  writeLines(
    "CHROM\treference_name_original\treference_name_sanitized\tPOS\tREF\tALT\tQUAL\tFILTER\tType\tref_depth_downsampled\talt_depth_downsampled\tvariant_percent\tsampleID\tref_upstream_20bp\tref_downstream_20bp",
    file.path(result_dir, "variant_all.tsv")
  )
  writeLines("x", file.path(barcode_dir, "IGV", "barcode001.amp1.png"))
  writeLines("x", file.path(barcode_dir, "AB1", "barcode001.amp1.synthetic.ab1"))

  collect_variant_results(
    result_dir = result_dir,
    output_dir = output_dir,
    include_igv = TRUE,
    include_ab1 = TRUE,
    ab1_name_template = "{barcode}.{reference}.synthetic.ab1",
    igv_dir_name = "IGV",
    ab1_dir_name = "AB1"
  )

  readme <- readLines(file.path(output_dir, "README.txt"))
  readme_zh <- readLines(file.path(output_dir, "README.zh-CN.txt"))

  expect_true(any(grepl("IGV alignment snapshot", readme, fixed = TRUE)))
  expect_true(any(grepl("synthetic Sanger-style AB1", readme, fixed = TRUE)))
  expect_true(any(grepl("barcode*/IGV/", readme, fixed = TRUE)))
  expect_true(any(grepl("barcode*/AB1/", readme, fixed = TRUE)))
  expect_true(any(grepl("IGV 比对截图", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("AB1 峰图文件", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("barcode*/IGV/", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("barcode*/AB1/", readme_zh, fixed = TRUE)))
})
