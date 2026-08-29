test_that("fastq_stats builds direct and length-filtered commands", {
  fastq1 <- tempfile(fileext = ".fastq.gz")
  fastq2 <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq1)
  writeLines(c("@read2", "ACGT", "+", "!!!!"), fastq2)

  res <- fastq_stats(
    fastq = c(fastq1, fastq2),
    min_lengths = c(1000, 5000),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(length(res$commands), 5L)
  expect_match(res$commands[[1L]], "'seqkit' 'stats' '-a' '-T'", fixed = TRUE)
  expect_true(any(grepl("'seqkit' 'seq' '-m' '1000'", res$commands, fixed = TRUE)))
  expect_true(any(grepl("| 'seqkit' 'stats' '-a' '-T'", res$commands, fixed = TRUE)))
  expect_equal(res$command_meta$stat_type, c("original", rep("min_length", 4)))
})

test_that("fastq_stats dry-run supports conda_env", {
  fastq <- tempfile(fileext = ".fastq.gz")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  res <- fastq_stats(
    fastq = fastq,
    min_lengths = 1000,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$conda_env, "ont-tools")
  expect_match(res$commands[[1L]], "'conda' 'run' '-n' 'ont-tools' 'seqkit'", fixed = TRUE)
  expect_match(res$commands[[2L]], "| 'conda' 'run' '-n' 'ont-tools' 'seqkit' 'stats' '-a' '-T'", fixed = TRUE)
})

test_that("fastq_stats parses seqkit stats output", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_tsv <- tempfile(fileext = ".tsv")
  fake_bin <- tempfile("seqkit-bin-")
  dir.create(fake_bin)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cmd=\"${1:-}\"",
      "[[ -n \"$cmd\" ]] || exit 0",
      "shift",
      "case \"$cmd\" in",
      "  stats)",
      "    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tavg_len\\tQ20(%%)\\n'",
      "    if [[ \"$#\" -eq 2 && \"$1\" == '-a' && \"$2\" == '-T' ]]; then",
      "      cat >/dev/null",
      "      printf 'stdin\\tFASTQ\\tDNA\\t1\\t4\\t4.0\\t100.00\\n'",
      "    else",
      "      for arg in \"$@\"; do",
      "        [[ \"$arg\" == '-a' ]] && continue",
      "        [[ \"$arg\" == '-T' ]] && continue",
      "        printf '%s\\tFASTQ\\tDNA\\t1,234\\t4,936\\t4.0\\t100.00\\n' \"$arg\"",
      "      done",
      "    fi",
      "    ;;",
      "  seq)",
      "    cat \"${@: -1}\"",
      "    ;;",
      "  *)",
      "    exit 1",
      "    ;;",
      "esac"
    ),
    file.path(fake_bin, "seqkit")
  )
  Sys.chmod(file.path(fake_bin, "seqkit"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- fastq_stats(
    fastq = fastq,
    min_lengths = c(1000, 3000),
    output_tsv = output_tsv,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_s3_class(res$stats, "data.frame")
  expect_equal(nrow(res$stats), 3L)
  expect_equal(res$stats$stat_type, c("original", "min_length", "min_length"))
  expect_equal(res$stats$min_length, c(NA_integer_, 1000L, 3000L))
  expect_equal(res$stats$num_seqs[[1L]], 1234L)
  expect_equal(res$stats$sum_len[[1L]], 4936L)
  expect_equal(res$stats$avg_len[[1L]], 4)
  expect_equal(res$stats$q20pct[[1L]], 100)
  expect_false("source_fastq" %in% names(res$stats))
  expect_equal(res$stats$file[[2L]], normalizePath(fastq))
  expect_true(file.exists(output_tsv))
})

test_that("fastq_stats validates arguments", {
  fastq <- tempfile(fileext = ".fastq.gz")
  writeLines(c("@read1", "ACGT", "+", "!!!!"), fastq)

  expect_error(
    fastq_stats(character(), dry_run = TRUE, echo = FALSE),
    "fastq"
  )
  expect_error(
    fastq_stats(fastq, min_lengths = 0, dry_run = TRUE, echo = FALSE),
    "min_lengths"
  )
  expect_error(
    fastq_stats(fastq, include_original = FALSE, dry_run = TRUE, echo = FALSE),
    "At least one"
  )
  expect_error(
    fastq_stats(fastq, conda_env = "", dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
})
