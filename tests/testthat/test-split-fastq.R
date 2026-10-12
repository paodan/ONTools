write_test_fastq <- function(path, lengths, compressed = FALSE) {
  records <- unlist(lapply(seq_along(lengths), function(i) {
    sequence <- paste(rep("A", lengths[[i]]), collapse = "")
    c(paste0("@read", i), sequence, "+", chartr("A", "!", sequence))
  }), use.names = FALSE)

  if (isTRUE(compressed)) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(records, con)
  } else {
    writeLines(records, path)
  }
  invisible(records)
}

test_that("split_fastq plans the expected range and counts both boundaries", {
  fastq <- tempfile(fileext = ".fastq")
  output_dir <- tempfile("split-fastq-")
  write_test_fastq(fastq, c(4L, 5L, 6L, 8L))

  result <- split_fastq(
    fastq,
    out_fastq_path = output_dir,
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(result$read_count, 4L)
  expect_equal(result$peak_counts$source, "expected")
  expect_equal(result$peak_counts$read_count, 3L)
  expect_equal(result$peak_counts$read_percent, 75)
  expect_match(result$output_fastq, "_peak_1_5bp[.]fastq$")
  expect_equal(result$filter_results[[1L]]$status, NA_integer_)
  expect_true(all(c("-m", "4", "-M", "6", "-o") %in%
                    result$filter_results[[1L]]$args))
  expect_null(result$plot)
})

test_that("split_fastq reads gzip input and writes genuine gzip output", {
  fastq <- tempfile(fileext = ".fastq.gz")
  output_dir <- tempfile("split-fastq-")
  fake_bin <- tempfile("seqkit-bin-")
  dir.create(fake_bin)
  records <- write_test_fastq(fastq, c(4L, 5L, 6L), compressed = TRUE)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "out=''",
      "input=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    seq|--remove-gaps) shift ;;",
      "    -m|-M) shift 2 ;;",
      "    -o) out=\"$2\"; shift 2 ;;",
      "    *) input=\"$1\"; shift ;;",
      "  esac",
      "done",
      "gzip -cd \"$input\" | gzip -c > \"$out\""
    ),
    file.path(fake_bin, "seqkit")
  )
  Sys.chmod(file.path(fake_bin, "seqkit"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  result <- split_fastq(
    fastq,
    out_fastq_path = output_dir,
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    echo = FALSE
  )

  output <- unname(result$output_fastq[[1L]])
  con <- file(output, open = "rb")
  magic <- readBin(con, what = "raw", n = 2L)
  close(con)
  expect_equal(magic, as.raw(c(0x1f, 0x8b)))
  expect_equal(readLines(gzfile(output)), records)
  expect_equal(result$filter_results[[1L]]$status, 0L)
})

test_that("split_fastq adds detected peak ranges", {
  fastq <- tempfile(fileext = ".fastq")
  lengths <- c(
    95:105, 97:103, 96:104,
    295:305, 297:303, 296:304
  )
  write_test_fastq(fastq, lengths)

  result <- split_fastq(
    fastq,
    out_fastq_path = tempfile("split-fastq-"),
    expected_size = 200,
    min_read_length = 180,
    max_read_length = 220,
    by_peak = TRUE,
    min_distance = 100,
    density_n = 1024,
    x_breaks = NULL,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_s3_class(result$plot, "ggplot")
  expect_gt(nrow(result$peak_counts), 1L)
  expect_equal(result$peak_counts$source[[1L]], "expected")
  expect_true(all(result$peak_counts$source[-1L] == "detected"))
  expect_equal(length(result$output_fastq), nrow(result$peak_counts))
  expect_equal(length(result$filter_results), nrow(result$peak_counts))
})

test_that("split_fastq validates file format, records, and ranges", {
  wrong_suffix <- tempfile(fileext = ".fq")
  write_test_fastq(wrong_suffix, 4L)
  expect_error(split_fastq(wrong_suffix, dry_run = TRUE), "must end")

  malformed <- tempfile(fileext = ".fastq")
  writeLines(c("@read1", "ACGT", "+", "!!!"), malformed)
  expect_error(split_fastq(malformed, dry_run = TRUE), "lengths differ")

  valid <- tempfile(fileext = ".fastq")
  write_test_fastq(valid, 4L)
  expect_error(
    split_fastq(
      valid,
      expected_size = 5,
      min_read_length = 6,
      max_read_length = 8,
      dry_run = TRUE
    ),
    "expected_size"
  )
})

test_that("split_fastq warns and skips an empty FASTQ", {
  fastq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fastq, open = "wt")
  close(con)

  expect_warning(
    result <- split_fastq(fastq, dry_run = TRUE, echo = FALSE),
    "contains no FASTQ records; skipping"
  )
  expect_null(result)
})
