test_that("run_wf_16s builds expected dry-run command", {
  res <- run_wf_16s(
    fastq = "reads",
    out_dir = "results",
    work_dir = "work",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$command, "nextflow")
  expect_true(all(c(
    "run", "epi2me-labs/wf-16s",
    "--fastq", "reads",
    "--out_dir", "results",
    "-work-dir", "work",
    "-profile", "standard",
    "-resume"
  ) %in% res$args))
  expect_false("--work-dir" %in% res$args)
  expect_true("NXF_SYNTAX_PARSER=v1" %in% res$env)
  expect_true("NXF_ANSI_LOG=false" %in% res$env)
  expect_equal(res$paths, list(fastq = "reads", out_dir = "results", work_dir = "work"))
})

test_that("run_wf_16s supports quiet, env, and extra args", {
  res <- run_wf_16s(
    fastq = "reads",
    out_dir = "results",
    work_dir = "work",
    quiet = TRUE,
    resume = FALSE,
    extra_args = "--minimap2_by_reference --abundance_threshold 1",
    syntax_parser = "v1",
    ansi_log = FALSE,
    nextflow_env = c("NXF_OPTS=-Xms1g"),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$args[[1]], "-q")
  expect_false("-resume" %in% res$args)
  expect_true(res$uses_shell)
  expect_equal(res$execution_command, "sh")
  expect_equal(res$execution_args, "<temporary shell script>")
  expect_true("NXF_SYNTAX_PARSER=v1" %in% res$env)
  expect_true("NXF_ANSI_LOG=false" %in% res$env)
  expect_true("NXF_OPTS=-Xms1g" %in% res$env)
  expect_match(res$command_string, "--minimap2_by_reference", fixed = TRUE)
})

test_that("run_wf_16s validates arguments", {
  expect_error(run_wf_16s(fastq = "", dry_run = TRUE, echo = FALSE), "fastq")
  expect_error(run_wf_16s(work_dir = "", dry_run = TRUE, echo = FALSE), "work_dir")
  expect_error(run_wf_16s(resume = NA, dry_run = TRUE, echo = FALSE), "resume")
  expect_error(run_wf_16s(extra_args = character(), dry_run = TRUE, echo = FALSE), "extra_args")
  expect_error(run_wf_16s(nextflow_env = "BAD", dry_run = TRUE, echo = FALSE), "nextflow_env")
})

test_that("run_wf_16s lets nextflow_env override syntax parser", {
  res <- run_wf_16s(
    fastq = "reads",
    out_dir = "results",
    work_dir = "work",
    syntax_parser = NULL,
    nextflow_env = "NXF_SYNTAX_PARSER=v2",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$env, c("NXF_ANSI_LOG=false", "NXF_SYNTAX_PARSER=v2"))
})

test_that("run_wf_16s runs nextflow command", {
  fake_bin <- tempfile("wf16s-bin-")
  log_file <- tempfile("wf16s-log-")
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf '%s\\n' \"$*\" > \"$WF16S_LOG\""
    ),
    file.path(fake_bin, "nextflow")
  )
  Sys.chmod(file.path(fake_bin, "nextflow"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  old_log <- Sys.getenv("WF16S_LOG", unset = NA)
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  on.exit({
    if (is.na(old_log)) {
      Sys.unsetenv("WF16S_LOG")
    } else {
      Sys.setenv(WF16S_LOG = old_log)
    }
  }, add = TRUE)
  Sys.setenv(
    PATH = paste(fake_bin, old_path, sep = .Platform$path.sep),
    WF16S_LOG = log_file
  )

  res <- run_wf_16s(
    fastq = "reads",
    out_dir = "results",
    work_dir = "work",
    echo = FALSE,
    stdout = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_match(readLines(log_file), "-work-dir work", fixed = TRUE)
})
