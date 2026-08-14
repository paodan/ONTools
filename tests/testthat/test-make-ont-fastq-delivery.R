test_that("make_ont_fastq_delivery builds the default command", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  sample_sheet <- tempfile(fileext = ".csv")
  dir.create(fastq_dir)
  dir.create(output_dir)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq_dir, "sample.fastq"))
  writeLines("barcode,sample\nbarcode01,sampleA", sample_sheet)

  res <- make_ont_fastq_delivery(
    input = fastq_dir,
    output = output_dir,
    project = "PROJECT001",
    sample_sheet = sample_sheet,
    threads = 4,
    run_nanoplot = FALSE,
    run_multiqc = FALSE,
    overwrite = TRUE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$execution_command, "bash")
  expect_equal(res$status, NA_integer_)
  expect_true(any(res$args == "--input"))
  expect_true(any(res$args == normalizePath(fastq_dir)))
  expect_true(any(res$args == "--output"))
  expect_true(any(res$args == normalizePath(output_dir)))
  expect_true(any(res$args == "--project"))
  expect_true(any(res$args == "PROJECT001"))
  expect_true(any(res$args == "--threads"))
  expect_true(any(res$args == "4"))
  expect_true(any(res$args == "--sample-sheet"))
  expect_true(any(res$args == normalizePath(sample_sheet)))
  expect_true(any(res$args == "--overwrite"))
  expect_false(any(res$args == "--no-reuse-nanoplot"))
  expect_true(any(res$args == "--skip-nanoplot"))
  expect_true(any(res$args == "--skip-multiqc"))
  expect_match(res$command_string, "make_ont_fastq_delivery.sh", fixed = TRUE)
  expect_equal(res$paths$delivery_dir, file.path(normalizePath(output_dir), "PROJECT001_delivery"))
  expect_equal(res$paths$archive, file.path(normalizePath(output_dir), "PROJECT001_delivery.tar.gz"))
})

test_that("make_ont_fastq_delivery validates arguments", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  dir.create(fastq_dir)

  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      threads = 0,
      dry_run = TRUE,
      echo = FALSE
    ),
    "threads"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      run_nanoplot = "no",
      dry_run = TRUE,
      echo = FALSE
    ),
    "run_nanoplot"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      reuse_nanoplot = "yes",
      dry_run = TRUE,
      echo = FALSE
    ),
    "reuse_nanoplot"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      overwrite = "yes",
      dry_run = TRUE,
      echo = FALSE
    ),
    "overwrite"
  )
  expect_error(
    make_ont_fastq_delivery(
      input = fastq_dir,
      output = output_dir,
      project = "PROJECT001",
      sample_sheet = file.path(output_dir, "missing.csv"),
      dry_run = TRUE,
      echo = FALSE
    ),
    "sample_sheet"
  )
})

test_that("make_ont_fastq_delivery gives MultiQC per-sample NanoStats names", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  fake_bin <- tempfile("ont-delivery-bin-")
  dir.create(fastq_dir)
  dir.create(output_dir)
  dir.create(fake_bin)
  dir.create(file.path(output_dir, "PROJECT001_delivery"))
  dir.create(
    file.path(output_dir, "PROJECT001_delivery", "02_qc_report", "nanoplot", "barcode01"),
    recursive = TRUE
  )
  writeLines("stale", file.path(output_dir, "PROJECT001_delivery", "stale.txt"))
  writeLines("stale", file.path(output_dir, "PROJECT001_delivery.tar.gz"))
  writeLines(
    c("metric\tvalue", "reused\t1"),
    file.path(
      output_dir,
      "PROJECT001_delivery",
      "02_qc_report",
      "nanoplot",
      "barcode01",
      "barcode01_NanoStats.txt"
    )
  )

  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq_dir, "barcode01.fastq"))
  writeLines(c("@read1", "TGCA", "+", "!!!!"), file.path(fastq_dir, "barcode02.fastq"))

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\n'",
      "for f in \"$@\"; do",
      "  case \"$f\" in",
      "    --*) ;;",
      "    *) [[ -f \"$f\" ]] && printf '%s\\tFASTQ\\tDNA\\t1\\t4\\n' \"$f\" ;;",
      "  esac",
      "done"
    ),
    file.path(fake_bin, "seqkit")
  )
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "outdir=''",
      "prefix=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --outdir) outdir=\"$2\"; shift 2 ;;",
      "    --prefix) prefix=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "mkdir -p \"$outdir\"",
      "printf 'metric\\tvalue\\nreads\\t1\\n' > \"$outdir/${prefix}NanoStats.txt\""
    ),
    file.path(fake_bin, "NanoPlot")
  )
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "input=\"$1\"",
      "outdir=''",
      "filename='multiqc_report.html'",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --outdir|-o) outdir=\"$2\"; shift 2 ;;",
      "    --filename) filename=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "mkdir -p \"$outdir/multiqc_data\"",
      "find \"$input\" -type f -name '*_NanoStats.txt' -exec basename {} \\; | sort > \"$outdir/multiqc_input_files.txt\"",
      "printf '<html></html>\\n' > \"$outdir/$filename\""
    ),
    file.path(fake_bin, "multiqc")
  )
  Sys.chmod(file.path(fake_bin, c("seqkit", "NanoPlot", "multiqc")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  make_ont_fastq_delivery(
    input = fastq_dir,
    output = output_dir,
    project = "PROJECT001",
    script = system.file("scripts", "make_ont_fastq_delivery.sh", package = "ONTools"),
    overwrite = TRUE,
    echo = FALSE,
    stdout = FALSE,
    stderr = FALSE
  )

  multiqc_inputs <- readLines(file.path(
    output_dir,
    "PROJECT001_delivery",
    "02_qc_report",
    "multiqc_input_files.txt"
  ))
  fastq_stats <- utils::read.delim(
    file.path(
      output_dir,
      "PROJECT001_delivery",
      "02_qc_report",
      "fastq_stats.tsv"
    ),
    check.names = FALSE
  )
  readme_zh <- readLines(file.path(
    output_dir,
    "PROJECT001_delivery",
    "README.zh-CN.txt"
  ))

  expect_equal(multiqc_inputs, c("barcode01_NanoStats.txt", "barcode02_NanoStats.txt"))
  expect_equal(fastq_stats$file, c("01_fastq/barcode01.fastq", "01_fastq/barcode02.fastq"))
  expect_true(any(grepl("项目：PROJECT001", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("01_fastq/", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("MD5 校验值", readme_zh, fixed = TRUE)))
  expect_equal(
    readLines(file.path(
      output_dir,
      "PROJECT001_delivery",
      "02_qc_report",
      "nanoplot",
      "barcode01",
      "barcode01_NanoStats.txt"
    )),
    c("metric\tvalue", "reused\t1")
  )
  expect_false(file.exists(file.path(output_dir, "PROJECT001_delivery", "stale.txt")))
  expect_false(file.exists(file.path(
    output_dir,
    "PROJECT001_delivery",
    "02_qc_report",
    "fastq_files.txt"
  )))
  expect_false(dir.exists(file.path(
    output_dir,
    "PROJECT001_delivery",
    "02_qc_report",
    "nanoplot_for_multiqc"
  )))
})

test_that("make_ont_fastq_delivery can disable NanoPlot reuse", {
  fastq_dir <- tempfile("ont-fastq-")
  output_dir <- tempfile("ont-delivery-")
  dir.create(fastq_dir)
  dir.create(output_dir)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq_dir, "barcode01.fastq"))

  res <- make_ont_fastq_delivery(
    input = fastq_dir,
    output = output_dir,
    project = "PROJECT001",
    reuse_nanoplot = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true(any(res$args == "--no-reuse-nanoplot"))
})
