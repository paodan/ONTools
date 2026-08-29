test_that("run_flye_assembly builds expected dry-run command", {
  reads <- tempfile(fileext = ".fastq")
  out_dir <- tempfile("flye-")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- run_flye_assembly(
    reads = reads,
    genome_size = "180k",
    min_overlap = 3000,
    threads = 22,
    out_dir = out_dir,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$command, "flye")
  expect_equal(res$read_type, "nano_raw")
  expect_true(all(c(
    "--nano-raw", normalizePath(reads),
    "--genome-size", "180k",
    "--threads", "22",
    "--out-dir", normalizePath(out_dir, mustWork = FALSE),
    "--min-overlap", "3000"
  ) %in% res$args))
  expect_match(res$command_string, "--nano-raw", fixed = TRUE)
})

test_that("run_flye_assembly supports optional flags and conda_env", {
  reads <- tempfile(fileext = ".fastq")
  out_dir <- tempfile("flye-")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- run_flye_assembly(
    reads = reads,
    genome_size = "180k",
    read_type = "nano_hq",
    iterations = 2,
    meta = TRUE,
    plasmids = TRUE,
    trestle = TRUE,
    keep_haplotypes = TRUE,
    no_alt_contigs = TRUE,
    polisher = "medaka",
    read_error = 0.03,
    extra_args = c("--debug"),
    out_dir = out_dir,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$command, "conda")
  expect_equal(res$args[1:4], c("run", "-n", "ont-tools", "flye"))
  expect_true(all(c(
    "--nano-hq",
    "--iterations", "2",
    "--meta",
    "--plasmids",
    "--trestle",
    "--keep-haplotypes",
    "--no-alt-contigs",
    "--polisher", "medaka",
    "--read-error", "0.03",
    "--debug"
  ) %in% res$args))
})

test_that("run_flye_assembly runs flye command", {
  reads <- tempfile(fileext = ".fastq")
  out_dir <- tempfile("flye-")
  fake_bin <- tempfile("flye-bin-")
  dir.create(fake_bin)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "out=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --out-dir) out=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "mkdir -p \"$out\"",
      "printf '>contig_1\\nACGT\\n' > \"$out/assembly.fasta\"",
      "printf 'seq_name\\tlength\\tcov.\\ncontig_1\\t4\\t1\\n' > \"$out/assembly_info.txt\"",
      "printf 'flye done\\n' > \"$out/flye.log\""
    ),
    file.path(fake_bin, "flye")
  )
  Sys.chmod(file.path(fake_bin, "flye"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- run_flye_assembly(
    reads = reads,
    genome_size = "180k",
    out_dir = out_dir,
    echo = FALSE,
    stdout = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(file.exists(res$paths$assembly))
  expect_true(file.exists(res$paths$assembly_info))
  expect_true(file.exists(res$paths$log))
})

test_that("run_flye_assembly validates arguments", {
  reads <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  expect_error(
    run_flye_assembly(reads, tempfile("flye-"), "180k", threads = 0,
                      dry_run = TRUE, echo = FALSE),
    "threads"
  )
  expect_error(
    run_flye_assembly(reads, tempfile("flye-"), "180k", min_overlap = 0,
                      dry_run = TRUE, echo = FALSE),
    "min_overlap"
  )
  expect_error(
    run_flye_assembly(reads, tempfile("flye-"), "180k", read_error = 2,
                      dry_run = TRUE, echo = FALSE),
    "read_error"
  )
  expect_error(
    run_flye_assembly(reads, tempfile("flye-"), "180k", conda_env = "",
                      dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
  expect_error(
    run_flye_assembly("missing.fastq", tempfile("flye-"), "180k",
                      dry_run = TRUE, echo = FALSE),
    "reads"
  )
})
