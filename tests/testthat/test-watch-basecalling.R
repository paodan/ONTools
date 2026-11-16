test_that("watch_file lists matching regular files", {
  path <- tempfile("watch-file-")
  dir.create(path)
  writeBin(as.raw(1:4), file.path(path, "batch01.pod5"))
  writeLines("bam", file.path(path, "batch01.bam"))
  dir.create(file.path(path, "not-a-file.pod5"))

  result <- watch_file(path, ".pod5")

  expect_equal(result$name, "batch01")
  expect_equal(result$size, 4)
  expect_equal(result$path, normalizePath(file.path(path, "batch01.pod5")))
})

test_that("watch_file validates its extension", {
  path <- tempfile("watch-file-")
  dir.create(path)

  expect_error(watch_file(path, "pod5"), "must start")
  expect_error(watch_file(path, "."), "include a suffix")
})

test_that("watch_file uses natural numeric ordering", {
  path <- tempfile("watch-natural-order-")
  dir.create(path)
  stems <- c("batch11", "batch2", "batch0", "batch10", "batch1", "batch9")
  for (stem in stems) {
    writeLines("pod5", file.path(path, paste0(stem, ".pod5")))
  }

  result <- watch_file(path, ".pod5")

  expect_equal(
    result$name,
    c("batch0", "batch1", "batch2", "batch9", "batch10", "batch11")
  )
})

test_that("basecall_one builds direct and conda dry-run commands", {
  pod5 <- tempfile(fileext = ".pod5")
  writeLines("pod5", pod5)

  direct <- basecall_one(
    pod5,
    tempfile(fileext = ".bam"),
    model = "hac",
    dry_run = TRUE,
    echo = FALSE
  )
  expect_equal(direct$status, NA_integer_)
  expect_equal(direct$command, "dorado")
  expect_equal(direct$args[1:2], c("basecaller", "hac"))

  conda <- basecall_one(
    pod5,
    tempfile(fileext = ".bam"),
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )
  expect_equal(conda$command, "conda")
  expect_equal(
    conda$args[1:6],
    c("run", "--no-capture-output", "-n", "ont-tools", "--", "dorado")
  )
})

test_that("basecall_one adds custom barcode files and extra parameters", {
  pod5 <- tempfile(fileext = ".pod5")
  arrangement <- tempfile(fileext = ".toml")
  sequences <- tempfile(fileext = ".fasta")
  writeLines("pod5", pod5)
  writeLines(c("[arrangement]", "name = \"MY_CUSTOM_KIT\""), arrangement)
  writeLines(c(">BC01", "ACGT"), sequences)

  result <- basecall_one(
    pod5_file = pod5,
    bam_file = tempfile(fileext = ".bam"),
    kit_name = "MY_CUSTOM_KIT",
    barcode_arrangement = arrangement,
    barcode_sequences = sequences,
    extra_basecaller_params = c("--device", "cuda:all", "--batchsize", "0"),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(
    tail(result$args, 10L),
    c(
      "--kit-name", "MY_CUSTOM_KIT",
      "--barcode-arrangement", normalizePath(arrangement),
      "--barcode-sequences", normalizePath(sequences),
      "--device", "cuda:all", "--batchsize", "0"
    )
  )
  expect_equal(result$paths$barcode_arrangement, normalizePath(arrangement))
  expect_equal(result$paths$barcode_sequences, normalizePath(sequences))
})

test_that("basecall_one validates custom barcode and extra parameters", {
  pod5 <- tempfile(fileext = ".pod5")
  arrangement <- tempfile(fileext = ".toml")
  sequences <- tempfile(fileext = ".fasta")
  writeLines("pod5", pod5)
  writeLines("[arrangement]", arrangement)
  writeLines(c(">BC01", "ACGT"), sequences)

  expect_error(
    basecall_one(
      pod5, tempfile(fileext = ".bam"),
      barcode_arrangement = arrangement, dry_run = TRUE, echo = FALSE
    ),
    "requires `kit_name`"
  )
  expect_error(
    basecall_one(
      pod5, tempfile(fileext = ".bam"), kit_name = "CUSTOM",
      barcode_sequences = sequences, dry_run = TRUE, echo = FALSE
    ),
    "requires `barcode_arrangement`"
  )
  expect_error(
    basecall_one(
      pod5, tempfile(fileext = ".bam"),
      extra_basecaller_params = c("--kit-name", "OTHER"),
      dry_run = TRUE, echo = FALSE
    ),
    "cannot contain"
  )
  expect_error(
    basecall_one(
      pod5, tempfile(fileext = ".bam"),
      extra_basecaller_params = "--emit-fastq",
      dry_run = TRUE, echo = FALSE
    ),
    "output-changing"
  )
  expect_error(
    basecall_one(
      pod5, tempfile(fileext = ".bam"),
      extra_basecaller_params = c("--device", NA_character_),
      dry_run = TRUE, echo = FALSE
    ),
    "without missing"
  )
})

test_that("basecall_one runs Dorado and publishes a completed BAM", {
  pod5 <- tempfile(fileext = ".pod5")
  bam <- tempfile(fileext = ".bam")
  dorado <- tempfile("fake-dorado-")
  writeLines("pod5", pod5)
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "[[ \"${1:-}\" == 'basecaller' ]]",
      "printf 'fake-bam\\n'"
    ),
    dorado
  )
  Sys.chmod(dorado, "0755")

  result <- basecall_one(
    pod5,
    bam,
    dorado = dorado,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(result$status, 0L)
  expect_true(file.exists(bam))
  expect_equal(readLines(bam), "fake-bam")
  expect_length(list.files(dirname(bam), paste0(basename(bam), ".*partial")), 0L)
})

test_that("basecall_one removes partial output after failure", {
  pod5 <- tempfile(fileext = ".pod5")
  bam <- tempfile(fileext = ".bam")
  dorado <- tempfile("fake-dorado-")
  writeLines("pod5", pod5)
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "printf 'partial\\n'",
      "exit 9"
    ),
    dorado
  )
  Sys.chmod(dorado, "0755")

  expect_error(
    basecall_one(pod5, bam, dorado = dorado, echo = FALSE, stderr = FALSE),
    "exit status: 9"
  )
  expect_false(file.exists(bam))
})

test_that("watch_basecalling plans only POD5 files without BAM files", {
  pod5_dir <- tempfile("pod5-")
  bam_dir <- tempfile("bam-")
  dir.create(pod5_dir)
  dir.create(bam_dir)
  writeLines("pod5", file.path(pod5_dir, "new.pod5"))
  writeLines("pod5", file.path(pod5_dir, "done.pod5"))
  writeLines("bam", file.path(bam_dir, "done.bam"))

  result <- watch_basecalling(
    pod5_dir,
    bam_dir,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_named(result$basecalls, "new")
  expect_equal(result$basecalls$new$status, NA_integer_)
  expect_equal(result$basecalls$new$paths$bam_file, file.path(normalizePath(bam_dir), "new.bam"))
})

test_that("watch_basecalling passes custom basecaller parameters through", {
  pod5_dir <- tempfile("pod5-")
  bam_dir <- tempfile("bam-")
  arrangement <- tempfile(fileext = ".toml")
  sequences <- tempfile(fileext = ".fasta")
  dir.create(pod5_dir)
  writeLines("pod5", file.path(pod5_dir, "batch01.pod5"))
  writeLines("[arrangement]", arrangement)
  writeLines(c(">BC01", "ACGT"), sequences)

  result <- watch_basecalling(
    pod5_dir,
    bam_dir,
    kit_name = "CUSTOM",
    barcode_arrangement = arrangement,
    barcode_sequences = sequences,
    extra_basecaller_params = c("--device", "cpu"),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true("--barcode-arrangement" %in% result$basecalls$batch01$args)
  expect_equal(tail(result$basecalls$batch01$args, 2L), c("--device", "cpu"))
})

test_that("watch_basecalling invokes the real basecaller", {
  pod5_dir <- tempfile("pod5-")
  bam_dir <- tempfile("bam-")
  dorado <- tempfile("fake-dorado-")
  dir.create(pod5_dir)
  writeLines("pod5", file.path(pod5_dir, "batch01.pod5"))
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf 'fake-bam\\n'"
    ),
    dorado
  )
  Sys.chmod(dorado, "0755")

  result <- watch_basecalling(
    pod5_dir,
    bam_dir,
    max_run_time = 0.5,
    check_interval = 0.01,
    idle_interval = 0.01,
    dorado = dorado,
    echo = FALSE,
    stderr = FALSE
  )

  expect_true(file.exists(file.path(bam_dir, "batch01.bam")))
  expect_named(result$basecalls, "batch01")
  expect_equal(result$basecalls$batch01$status, 0L)
})
