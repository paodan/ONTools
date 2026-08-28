test_that("dehost_fastq builds expected dry-run commands and paths", {
  reference <- tempfile(fileext = ".fasta")
  fastq <- tempfile(fileext = ".fastq.gz")
  out_dir <- tempfile("dehost-")
  writeLines(c(">ref", "ACGT"), reference)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- dehost_fastq(
    reference = reference,
    fastq = fastq,
    out_dir = out_dir,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_match(res$commands$minimap2, "minimap2", fixed = TRUE)
  expect_match(res$commands$minimap2, "--secondary=no", fixed = TRUE)
  expect_match(res$commands$seqkit_grep, "grep", fixed = TRUE)
  expect_match(res$commands$seqkit_grep, "-v", fixed = TRUE)
  expect_equal(res$paths$output_fastq, file.path(normalizePath(out_dir), paste0(basename(sub("[.]fastq[.]gz$", "", fastq)), ".dehost.fastq.gz")))
})

test_that("dehost_fastq dry-run supports conda_env", {
  reference <- tempfile(fileext = ".fasta")
  fastq <- tempfile(fileext = ".fastq.gz")
  writeLines(c(">ref", "ACGT"), reference)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- dehost_fastq(
    reference = reference,
    fastq = fastq,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$conda_env, "ont-tools")
  expect_match(res$commands$minimap2, "'conda' 'run' '-n' 'ont-tools' 'minimap2'", fixed = TRUE)
  expect_match(res$commands$seqkit_grep, "'conda' 'run' '-n' 'ont-tools' 'seqkit'", fixed = TRUE)
})

test_that("dehost_fastq_read_ids filters PAF by MAPQ and aligned fraction", {
  paf <- tempfile(fileext = ".paf")
  writeLines(
    c(
      paste("read1", 10000, 0, 9000, "+", "ref", 5000000, 1, 9001, 9000, 9000, 60, sep = "\t"),
      paste("read2", 10000, 0, 7000, "+", "ref", 5000000, 1, 7001, 7000, 7000, 60, sep = "\t"),
      paste("read3", 10000, 0, 9000, "+", "ref", 5000000, 1, 9001, 9000, 9000, 10, sep = "\t")
    ),
    paf
  )

  ids <- dehost_fastq_read_ids(paf, min_mapq = 20, min_aln_frac = 0.8)

  expect_equal(ids, "read1")
})

test_that("dehost_fastq runs minimap2 and seqkit commands", {
  reference <- tempfile(fileext = ".fasta")
  fastq <- tempfile(fileext = ".fastq.gz")
  out_dir <- tempfile("dehost-")
  fake_bin <- tempfile("dehost-bin-")
  dir.create(fake_bin)
  writeLines(c(">ref", "ACGT"), reference)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf 'read1\\t10000\\t0\\t9000\\t+\\tref\\t5000\\t0\\t9000\\t9000\\t9000\\t60\\n'",
      "printf 'read2\\t10000\\t0\\t7000\\t+\\tref\\t5000\\t0\\t7000\\t7000\\t7000\\t60\\n'"
    ),
    file.path(fake_bin, "minimap2")
  )
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cmd=\"$1\"",
      "shift",
      "case \"$cmd\" in",
      "  grep)",
      "    out=''",
      "    while [[ $# -gt 0 ]]; do",
      "      case \"$1\" in",
      "        -o) out=\"$2\"; shift 2 ;;",
      "        *) shift ;;",
      "      esac",
      "    done",
      "    printf '@read2\\nACGT\\n+\\n!!!!\\n' > \"$out\"",
      "    ;;",
      "  stats)",
      "    printf 'file\\tformat\\ttype\\tnum_seqs\\n%s\\tFASTQ\\tDNA\\t1\\n' \"$2\"",
      "    ;;",
      "  *)",
      "    exit 1",
      "    ;;",
      "esac"
    ),
    file.path(fake_bin, "seqkit")
  )
  Sys.chmod(file.path(fake_bin, c("minimap2", "seqkit")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- dehost_fastq(
    reference = reference,
    fastq = fastq,
    out_dir = out_dir,
    echo = FALSE
  )

  expect_equal(res$status, 0L)
  expect_equal(res$n_reference_like_reads, 1L)
  expect_equal(readLines(res$paths$reference_like_read_ids), "read1")
  expect_true(file.exists(res$paths$output_fastq))
  expect_true(file.exists(res$paths$stats))
})

test_that("dehost_fastq runs commands through conda_env", {
  reference <- tempfile(fileext = ".fasta")
  fastq <- tempfile(fileext = ".fastq.gz")
  out_dir <- tempfile("dehost-")
  fake_bin <- tempfile("dehost-bin-")
  conda_log <- tempfile("conda-log-")
  dir.create(fake_bin)
  writeLines(c(">ref", "ACGT"), reference)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf 'read1\\t10000\\t0\\t9000\\t+\\tref\\t5000\\t0\\t9000\\t9000\\t9000\\t60\\n'"
    ),
    file.path(fake_bin, "minimap2")
  )
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cmd=\"$1\"",
      "shift",
      "case \"$cmd\" in",
      "  grep)",
      "    out=''",
      "    while [[ $# -gt 0 ]]; do",
      "      case \"$1\" in",
      "        -o) out=\"$2\"; shift 2 ;;",
      "        *) shift ;;",
      "      esac",
      "    done",
      "    printf '@read2\\nACGT\\n+\\n!!!!\\n' > \"$out\"",
      "    ;;",
      "  stats)",
      "    printf 'file\\tformat\\ttype\\tnum_seqs\\n%s\\tFASTQ\\tDNA\\t1\\n' \"$2\"",
      "    ;;",
      "  *)",
      "    exit 1",
      "    ;;",
      "esac"
    ),
    file.path(fake_bin, "seqkit")
  )
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf '%s\\n' \"$*\" >> \"$CONDA_LOG\"",
      "[[ \"$1\" == 'run' ]]",
      "[[ \"$2\" == '-n' ]]",
      "[[ \"$3\" == 'ont-tools' ]]",
      "shift 3",
      "exec \"$@\""
    ),
    file.path(fake_bin, "conda")
  )
  Sys.chmod(file.path(fake_bin, c("minimap2", "seqkit", "conda")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  old_conda_log <- Sys.getenv("CONDA_LOG", unset = NA)
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  on.exit({
    if (is.na(old_conda_log)) {
      Sys.unsetenv("CONDA_LOG")
    } else {
      Sys.setenv(CONDA_LOG = old_conda_log)
    }
  }, add = TRUE)
  Sys.setenv(
    PATH = paste(fake_bin, old_path, sep = .Platform$path.sep),
    CONDA_LOG = conda_log
  )

  res <- dehost_fastq(
    reference = reference,
    fastq = fastq,
    out_dir = out_dir,
    conda_env = "ont-tools",
    echo = FALSE
  )

  expect_equal(res$status, 0L)
  expect_equal(res$conda_env, "ont-tools")
  expect_true(file.exists(res$paths$output_fastq))
  expect_true(all(grepl("^run -n ont-tools ", readLines(conda_log))))
})

test_that("dehost_fastq validates arguments", {
  reference <- tempfile(fileext = ".fasta")
  fastq <- tempfile(fileext = ".fastq.gz")
  writeLines(c(">ref", "ACGT"), reference)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  expect_error(
    dehost_fastq(reference, fastq, threads = 0, dry_run = TRUE, echo = FALSE),
    "threads"
  )
  expect_error(
    dehost_fastq(reference, fastq, min_aln_frac = 1.2, dry_run = TRUE, echo = FALSE),
    "min_aln_frac"
  )
  expect_error(
    dehost_fastq(reference, fastq, conda_env = "", dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
  expect_error(
    dehost_fastq("missing.fasta", fastq, dry_run = TRUE, echo = FALSE),
    "reference"
  )
})
