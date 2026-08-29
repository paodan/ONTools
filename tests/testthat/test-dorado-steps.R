test_that("dorado_basecall builds expected dry-run command", {
  proj <- tempfile("dorado-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  res <- dorado_basecall(
    proj = proj,
    model = "sup",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$command, "dorado")
  expect_equal(res$args[1:2], c("basecaller", "sup"))
  expect_equal(res$paths$output_bam, file.path(normalizePath(proj), "bam", "calls_sup.bam"))
  expect_match(res$command_string, "basecaller", fixed = TRUE)
})

test_that("dorado_basecall supports conda_env", {
  proj <- tempfile("dorado-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  res <- dorado_basecall(
    proj = proj,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$command, "conda")
  expect_equal(res$args[1:4], c("run", "-n", "ont-tools", "dorado"))
})

test_that("dorado_demux_bam builds kit and custom barcode commands", {
  calls_bam <- tempfile(fileext = ".bam")
  demux_dir <- tempfile("demux-")
  config_dir <- tempfile("barcode-config-")
  dir.create(config_dir)
  writeLines("bam", calls_bam)
  writeLines("[arrangement]\nname = \"YS-NB576\"", file.path(config_dir, "yisheng_576.toml"))
  writeLines(c(">barcode001", "ACGT"), file.path(config_dir, "yisheng_576.fasta"))

  standard <- dorado_demux_bam(
    calls_bam = calls_bam,
    demux_dir = demux_dir,
    kit_name = "EXP-NBD196",
    barcode_both_ends = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_true(all(c("demux", "--output-dir", "--emit-summary", "--kit-name", "EXP-NBD196") %in% standard$args))
  expect_false("--barcode-both-ends" %in% standard$args)

  custom <- dorado_demux_bam(
    calls_bam = calls_bam,
    demux_dir = demux_dir,
    kit_name = "YS-NB576",
    barcode_config_dir = config_dir,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_true("--barcode-arrangement" %in% custom$args)
  expect_true("--barcode-sequences" %in% custom$args)
  expect_true("--barcode-both-ends" %in% custom$args)
  expect_equal(custom$paths$barcode_arrangement, normalizePath(file.path(config_dir, "yisheng_576.toml")))
})

test_that("dorado_bam_to_fastq plans barcode conversions", {
  demux_dir <- tempfile("demux-")
  bam_pass <- file.path(demux_dir, "run", "sample", "readset", "bam_pass")
  dir.create(file.path(bam_pass, "barcode001"), recursive = TRUE)
  dir.create(file.path(bam_pass, "barcode002"), recursive = TRUE)
  dir.create(file.path(bam_pass, "unclassified"), recursive = TRUE)
  writeLines("bam1", file.path(bam_pass, "barcode001", "a.bam"))
  writeLines("bam2", file.path(bam_pass, "barcode001", "b.bam"))
  writeLines("bam3", file.path(bam_pass, "barcode002", "c.bam"))
  writeLines("bam4", file.path(bam_pass, "unclassified", "u.bam"))

  res <- dorado_bam_to_fastq(
    demux_dir = demux_dir,
    include_unclassified = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(nrow(res$conversions), 2L)
  expect_equal(res$conversions$barcode, c("barcode001", "barcode002"))
  expect_equal(res$conversions$n_bam, c(2L, 1L))
  expect_true(all(grepl("fastq_pass_trim", res$conversions$output_fastq, fixed = TRUE)))
  expect_true(all(grepl("[.]fastq[.]gz$", res$conversions$output_fastq)))
  expect_equal(length(res$commands), 2L)
})

test_that("dorado_bam_to_fastq supports conda_env command prefix", {
  demux_dir <- tempfile("demux-")
  bam_pass <- file.path(demux_dir, "run", "sample", "readset", "bam_pass")
  dir.create(file.path(bam_pass, "barcode001"), recursive = TRUE)
  writeLines("bam1", file.path(bam_pass, "barcode001", "a.bam"))

  res <- dorado_bam_to_fastq(
    demux_dir = demux_dir,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_match(res$commands[[1]], "'conda' run -n 'ont-tools' 'samtools' fastq", fixed = TRUE)
})

test_that("dorado_bam_to_fastq runs conversion commands", {
  demux_dir <- tempfile("demux-")
  bam_pass <- file.path(demux_dir, "run", "sample", "readset", "bam_pass")
  fake_bin <- tempfile("dorado-step-bin-")
  dir.create(file.path(bam_pass, "barcode001"), recursive = TRUE)
  dir.create(fake_bin)
  writeLines("bam1", file.path(bam_pass, "barcode001", "a.bam"))

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "if [[ \"${1:-}\" == 'fastq' ]]; then",
      "  printf '@read1\\nACGT\\n+\\n!!!!\\n'",
      "else",
      "  exit 1",
      "fi"
    ),
    file.path(fake_bin, "samtools")
  )
  Sys.chmod(file.path(fake_bin, "samtools"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- dorado_bam_to_fastq(
    demux_dir = demux_dir,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(file.exists(res$conversions$output_fastq[[1]]))
  expect_gt(file.info(res$conversions$output_fastq[[1]])$size, 0)
})

test_that("dorado step functions validate inputs", {
  expect_error(dorado_basecall("missing-project", dry_run = TRUE, echo = FALSE), "proj")

  calls_bam <- tempfile(fileext = ".bam")
  writeLines("bam", calls_bam)
  expect_error(
    dorado_demux_bam(calls_bam, tempfile("demux-"), conda_env = "",
                     dry_run = TRUE, echo = FALSE),
    "conda_env"
  )

  demux_dir <- tempfile("demux-")
  dir.create(demux_dir)
  expect_error(
    dorado_bam_to_fastq(demux_dir, dry_run = TRUE, echo = FALSE),
    "No barcode BAM folders"
  )
})
