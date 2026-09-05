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

test_that("dorado_basecall protects existing outputs", {
  proj <- tempfile("dorado-project-")
  output_bam <- tempfile(fileext = ".bam")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)
  writeLines("old-bam", output_bam)

  expect_error(
    dorado_basecall(
      proj = proj,
      output_bam = output_bam,
      overwrite = FALSE,
      dry_run = TRUE,
      echo = FALSE
    ),
    "overwrite"
  )
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

test_that("dorado_demux_bam supports threads", {
  calls_bam <- tempfile(fileext = ".bam")
  demux_dir <- tempfile("demux-")
  writeLines("bam", calls_bam)

  res <- dorado_demux_bam(
    calls_bam = calls_bam,
    demux_dir = demux_dir,
    threads = 8,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true(all(c("--threads", "8") %in% res$args))
})

test_that("dorado_basecall runs dorado and rejects empty output", {
  proj <- tempfile("dorado-project-")
  fake_bin <- tempfile("dorado-basecall-bin-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "if [[ \"${1:-}\" == 'basecaller' ]]; then",
      "  printf 'fake-bam\\n'",
      "else",
      "  exit 1",
      "fi"
    ),
    file.path(fake_bin, "dorado")
  )
  Sys.chmod(file.path(fake_bin, "dorado"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- dorado_basecall(proj = proj, echo = FALSE, stderr = FALSE)
  expect_equal(res$status, 0L)
  expect_true(file.exists(res$paths$output_bam))
  expect_gt(file.info(res$paths$output_bam)$size, 0)
})

test_that("dorado_basecall removes failed partial output", {
  proj <- tempfile("dorado-project-")
  fake_bin <- tempfile("dorado-basecall-bin-")
  output_bam <- tempfile(fileext = ".bam")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf 'partial-bam\\n'",
      "exit 17"
    ),
    file.path(fake_bin, "dorado")
  )
  Sys.chmod(file.path(fake_bin, "dorado"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  expect_error(
    dorado_basecall(
      proj = proj,
      output_bam = output_bam,
      echo = FALSE,
      stderr = FALSE
    ),
    "basecaller failed"
  )
  expect_false(file.exists(output_bam))
})

test_that("dorado_demux_bam runs dorado and checks bam_pass output", {
  calls_bam <- tempfile(fileext = ".bam")
  demux_dir <- tempfile("demux-")
  fake_bin <- tempfile("dorado-demux-bin-")
  dir.create(fake_bin)
  writeLines("bam", calls_bam)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "if [[ \"${1:-}\" != 'demux' ]]; then exit 1; fi",
      "outdir=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --output-dir) outdir=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "mkdir -p \"$outdir/run/sample/readset/bam_pass/barcode001\"",
      "printf 'bam\\n' > \"$outdir/run/sample/readset/bam_pass/barcode001/a.bam\""
    ),
    file.path(fake_bin, "dorado")
  )
  Sys.chmod(file.path(fake_bin, "dorado"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- dorado_demux_bam(
    calls_bam = calls_bam,
    demux_dir = demux_dir,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(dorado_demux_has_bam_pass(res$paths$demux_dir))
})

test_that("dorado_demux_bam rejects successful runs without BAM output", {
  calls_bam <- tempfile(fileext = ".bam")
  demux_dir <- tempfile("demux-")
  fake_bin <- tempfile("dorado-demux-bin-")
  dir.create(fake_bin)
  writeLines("bam", calls_bam)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "exit 0"
    ),
    file.path(fake_bin, "dorado")
  )
  Sys.chmod(file.path(fake_bin, "dorado"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  expect_error(
    dorado_demux_bam(
      calls_bam = calls_bam,
      demux_dir = demux_dir,
      echo = FALSE,
      stderr = FALSE
    ),
    "no BAM files"
  )
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
  expect_true(all(grepl("[.]fastq[.]gz[.]md5$", res$conversions$md5_file)))
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

  expect_match(res$commands[[1]], "set -o pipefail", fixed = TRUE)
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
  expect_true(file.exists(res$conversions$md5_file[[1]]))
  expect_gt(file.info(res$conversions$output_fastq[[1]])$size, 0)
  expect_equal(readLines(gzfile(res$conversions$output_fastq[[1]])), c("@read1", "ACGT", "+", "!!!!"))
  expect_equal(
    readLines(res$conversions$md5_file[[1]]),
    paste(unname(tools::md5sum(res$conversions$output_fastq[[1]])), "barcode001.fastq.gz")
  )
})

test_that("dorado_bam_to_fastq surfaces samtools failures in the pipe", {
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
      "exit 23"
    ),
    file.path(fake_bin, "samtools")
  )
  Sys.chmod(file.path(fake_bin, "samtools"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  expect_error(
    dorado_bam_to_fastq(
      demux_dir = demux_dir,
      echo = FALSE,
      stderr = FALSE
    ),
    "BAM to FASTQ conversion failed"
  )
})

test_that("dorado_bam_to_fastq rejects empty converted FASTQ by default", {
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
      "  exit 0",
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

  expect_error(
    dorado_bam_to_fastq(
      demux_dir = demux_dir,
      echo = FALSE,
      stderr = FALSE
    ),
    "empty or invalid"
  )
})

test_that("dorado_bam_to_fastq can allow empty converted FASTQ", {
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
      "  exit 0",
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
    allow_empty = TRUE,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(file.exists(res$conversions$output_fastq[[1]]))
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
