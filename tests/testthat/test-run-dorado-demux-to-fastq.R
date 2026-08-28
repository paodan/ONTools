test_that("run_dorado_demux_to_fastq exposes barcode-both-ends control", {
  proj <- tempfile("dorado-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  default_cmd <- run_dorado_demux_to_fastq(
    proj = proj,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_true("--barcode-both-ends" %in% default_cmd$args)

  relaxed_cmd <- run_dorado_demux_to_fastq(
    proj = proj,
    barcode_both_ends = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_false("--barcode-both-ends" %in% relaxed_cmd$args)
})

test_that("run_dorado_demux_to_fastq accepts the YS-NB576 kit alias", {
  proj <- tempfile("dorado-project-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)

  cmd <- run_dorado_demux_to_fastq(
    proj = proj,
    kit_name = "YS-NB576",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true("--kit-name" %in% cmd$args)
  expect_true("YS-NB576" %in% cmd$args)
  expect_equal(cmd$paths$demux_dir, file.path(normalizePath(proj), "demux_out_YS-NB576"))
})

test_that("run_dorado_demux_to_fastq shell maps YS-NB576 to custom barcode files", {
  script <- system.file("scripts", "run_dorado_demux_to_fastq", package = "ONTools")
  skip_if(script == "", "run_dorado_demux_to_fastq shell script is not installed")

  proj <- tempfile("dorado-project-")
  fake_bin <- tempfile("dorado-bin-")
  demux_args_file <- tempfile("dorado-demux-args-")
  dir.create(file.path(proj, "pod5"), recursive = TRUE)
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cmd=\"$1\"",
      "shift",
      "case \"$cmd\" in",
      "  basecaller)",
      "    printf 'fake-bam\\n'",
      "    ;;",
      "  demux)",
      "    printf '%s\\n' \"$@\" > \"$DORADO_DEMUX_ARGS_FILE\"",
      "    outdir=''",
      "    while [[ $# -gt 0 ]]; do",
      "      case \"$1\" in",
      "        --output-dir) outdir=\"$2\"; shift 2 ;;",
      "        *) shift ;;",
      "      esac",
      "    done",
      "    mkdir -p \"$outdir/run/sample/readset/bam_pass/barcode01\"",
      "    printf 'bam\\n' > \"$outdir/run/sample/readset/bam_pass/barcode01/barcode01.bam\"",
      "    ;;",
      "  *)",
      "    exit 1",
      "    ;;",
      "esac"
    ),
    file.path(fake_bin, "dorado")
  )
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
  Sys.chmod(file.path(fake_bin, c("dorado", "samtools")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  old_args_file <- Sys.getenv("DORADO_DEMUX_ARGS_FILE", unset = NA)
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  on.exit({
    if (is.na(old_args_file)) {
      Sys.unsetenv("DORADO_DEMUX_ARGS_FILE")
    } else {
      Sys.setenv(DORADO_DEMUX_ARGS_FILE = old_args_file)
    }
  }, add = TRUE)
  Sys.setenv(
    PATH = paste(fake_bin, old_path, sep = .Platform$path.sep),
    DORADO_DEMUX_ARGS_FILE = demux_args_file
  )

  status <- system2(
    "bash",
    c(script, "--proj", proj, "--kit-name", "YS-NB576"),
    stdout = FALSE,
    stderr = FALSE
  )

  demux_args <- readLines(demux_args_file)
  expect_identical(status, 0L)
  expect_true("--barcode-arrangement" %in% demux_args)
  expect_true("--barcode-sequences" %in% demux_args)
  expect_true(any(grepl("yisheng_576[.]toml$", demux_args)))
  expect_true(any(grepl("yisheng_576[.]fasta$", demux_args)))
  expect_false("--kit-name" %in% demux_args)
  expect_false("YS-NB576" %in% demux_args)
})
