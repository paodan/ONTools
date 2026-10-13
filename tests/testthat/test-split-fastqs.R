write_split_fastqs_test_file <- function(path, lengths, compressed = FALSE) {
  records <- unlist(lapply(seq_along(lengths), function(i) {
    sequence <- paste(rep("A", lengths[[i]]), collapse = "")
    quality <- paste(rep("!", lengths[[i]]), collapse = "")
    c(paste0("@read", i), sequence, "+", quality)
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

test_that("split_fastqs plans multiple files in a shared directory", {
  input_dir <- tempfile("split-fastqs-input-")
  output_dir <- tempfile("split-fastqs-output-")
  dir.create(input_dir)
  fastqs <- file.path(input_dir, c("a.fastq", "b.fastq.gz"))
  write_split_fastqs_test_file(fastqs[[1L]], c(4L, 5L, 6L, 8L))
  write_split_fastqs_test_file(
    fastqs[[2L]],
    c(4L, 5L, 5L, 7L),
    compressed = TRUE
  )

  results <- split_fastqs(
    fastqs,
    out_fastq_path = output_dir,
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_named(results, c("a.fastq", "b.fastq.gz"))
  expect_length(results, 2L)
  expect_true(all(vapply(results, is.list, logical(1))))
  expect_equal(results[[1L]]$peak_counts$read_count, 3L)
  expect_equal(results[[2L]]$peak_counts$read_count, 3L)
  expect_true(all(vapply(
    results,
    function(x) identical(x$output_dir, normalizePath(output_dir)),
    logical(1)
  )))
  expect_match(results[[1L]]$output_fastq, "[.]fastq$")
  expect_match(results[[2L]]$output_fastq, "[.]fastq[.]gz$")
})

test_that("split_fastqs supports default and per-file output directories", {
  root <- tempfile("split-fastqs-root-")
  dir.create(root)
  input_dirs <- file.path(root, c("one", "two"))
  dir.create(input_dirs[[1L]])
  dir.create(input_dirs[[2L]])
  fastqs <- file.path(input_dirs, c("a.fastq", "b.fastq"))
  write_split_fastqs_test_file(fastqs[[1L]], c(4L, 5L, 6L))
  write_split_fastqs_test_file(fastqs[[2L]], c(4L, 5L, 6L))

  default_results <- split_fastqs(
    fastqs,
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_equal(
    unname(vapply(default_results, `[[`, character(1), "output_dir")),
    normalizePath(file.path(input_dirs, "peaks"))
  )

  output_dirs <- file.path(root, c("output-one", "output-two"))
  separate_results <- split_fastqs(
    fastqs,
    out_fastq_path = output_dirs,
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    dry_run = TRUE,
    echo = FALSE
  )
  expect_equal(
    unname(vapply(separate_results, `[[`, character(1), "output_dir")),
    normalizePath(output_dirs)
  )
})

test_that("split_fastqs preserves NULL for an empty input and continues", {
  input_dir <- tempfile("split-fastqs-input-")
  dir.create(input_dir)
  empty <- file.path(input_dir, "empty.fastq.gz")
  con <- gzfile(empty, open = "wt")
  close(con)
  populated <- file.path(input_dir, "populated.fastq")
  write_split_fastqs_test_file(populated, c(4L, 5L, 6L))

  expect_warning(
    results <- split_fastqs(
      c(empty, populated),
      expected_size = 5,
      min_read_length = 4,
      max_read_length = 6,
      by_peak = FALSE,
      dry_run = TRUE,
      echo = FALSE
    ),
    "contains no FASTQ records; skipping"
  )

  expect_length(results, 2L)
  expect_null(results[[1L]])
  expect_true(is.list(results[[2L]]))
})

test_that("split_fastqs detects duplicated inputs and output collisions", {
  first_dir <- tempfile("split-fastqs-first-")
  second_dir <- tempfile("split-fastqs-second-")
  output_dir <- tempfile("split-fastqs-output-")
  dir.create(first_dir)
  dir.create(second_dir)
  first <- file.path(first_dir, "sample.fastq")
  second <- file.path(second_dir, "sample.fastq")
  write_split_fastqs_test_file(first, c(4L, 5L, 6L))
  write_split_fastqs_test_file(second, c(4L, 5L, 6L))

  expect_error(
    split_fastqs(c(first, first), dry_run = TRUE, echo = FALSE),
    "duplicated input"
  )
  expect_error(
    split_fastqs(
      c(first, second),
      out_fastq_path = output_dir,
      dry_run = TRUE,
      echo = FALSE
    ),
    "colliding outputs"
  )

  expect_no_error(split_fastqs(
    c(first, second),
    out_fastq_path = c(file.path(output_dir, "one"),
                       file.path(output_dir, "two")),
    expected_size = 5,
    min_read_length = 4,
    max_read_length = 6,
    by_peak = FALSE,
    dry_run = TRUE,
    echo = FALSE
  ))
})

test_that("split_fastqs validates vector arguments before processing", {
  fastq <- tempfile(fileext = ".fastq")
  write_split_fastqs_test_file(fastq, c(4L, 5L, 6L))

  expect_error(split_fastqs(character()), "fastqs")
  expect_error(split_fastqs("missing.fastq"), "does not exist")
  expect_error(
    split_fastqs(
      fastq,
      out_fastq_path = c("one", "two"),
      dry_run = TRUE,
      echo = FALSE
    ),
    "out_fastq_path"
  )
})
