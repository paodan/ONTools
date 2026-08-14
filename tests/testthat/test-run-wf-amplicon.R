test_that("run_wf_amplicon builds the default Nextflow command", {
  res <- run_wf_amplicon(dry_run = TRUE, echo = FALSE)

  expect_equal(res$command, "nextflow")
  expect_equal(res$execution_command, "nextflow")
  expect_equal(res$execution_args, res$args)
  expect_false(res$uses_shell)
  expect_equal(res$env, "NXF_SYNTAX_PARSER=v1")
  expect_equal(res$status, NA_integer_)
  expect_true(any(res$args == "julibeg/wf-amplicon"))
  expect_true(any(res$args == "--fastq"))
  expect_true(any(res$args == "./fastq_pass_trim"))
  expect_true(any(res$args == "--out_dir"))
  expect_true(any(res$args == "./results/wf_amplicon_denovo"))
  expect_true(any(res$args == "--min_read_length"))
  expect_true(any(res$args == "2000"))
  expect_true(any(res$args == "--max_read_length"))
  expect_true(any(res$args == "3300"))
  expect_true(any(res$args == "-resume"))
  expect_match(res$command_string, "nextflow")
  expect_match(res$command_string, "julibeg/wf-amplicon", fixed = TRUE)
})

test_that("run_wf_amplicon appends raw extra arguments", {
  res <- run_wf_amplicon(
    fastq = "fastq",
    out_dir = "results",
    profile = "standard",
    resume = FALSE,
    extra_args = "--threads 16 --custom_param 'raw value'",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_false(any(res$args == "-resume"))
  expect_equal(res$extra_args, "--threads 16 --custom_param 'raw value'")
  expect_equal(res$execution_command, "sh")
  expect_equal(res$execution_args, "<temporary shell script>")
  expect_true(res$uses_shell)
  expect_match(res$command_string, "--threads 16 --custom_param 'raw value'", fixed = TRUE)
})

test_that("run_wf_amplicon can override Nextflow environment", {
  res <- run_wf_amplicon(
    syntax_parser = NULL,
    nextflow_env = c("NXF_OFFLINE=true", "NXF_ANSI_LOG=false"),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$env, c("NXF_OFFLINE=true", "NXF_ANSI_LOG=false"))

  res2 <- run_wf_amplicon(
    syntax_parser = "v1",
    nextflow_env = "NXF_SYNTAX_PARSER=v2",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res2$env, "NXF_SYNTAX_PARSER=v2")
})

test_that("run_wf_amplicon validates numeric arguments", {
  expect_error(
    run_wf_amplicon(
      min_read_length = 3300,
      max_read_length = 2000,
      dry_run = TRUE,
      echo = FALSE
    ),
    "less than or equal"
  )
  expect_error(
    run_wf_amplicon(min_n_reads = 0, dry_run = TRUE, echo = FALSE),
    "min_n_reads"
  )
  expect_error(
    run_wf_amplicon(min_read_qual = -1, dry_run = TRUE, echo = FALSE),
    "min_read_qual"
  )
  expect_error(
    run_wf_amplicon(nextflow_env = "bad-env", dry_run = TRUE, echo = FALSE),
    "nextflow_env"
  )
})
