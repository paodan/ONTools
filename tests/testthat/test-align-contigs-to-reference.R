test_that("align_contigs_to_reference builds expected dry-run commands", {
  reference <- tempfile(fileext = ".fasta")
  assembly <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">ref", "ACGTACGT"), reference)
  writeLines(c(">contig1", "ACGT"), assembly)

  res <- align_contigs_to_reference(
    reference = reference,
    assembly = assembly,
    align_bam = bam,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$preset, "asm5")
  expect_match(res$commands$minimap2, "'minimap2' '-x' 'asm5' '-a' '-t' '16'", fixed = TRUE)
  expect_match(res$commands$minimap2, "--secondary=no", fixed = TRUE)
  expect_match(res$commands$samtools_sort, "'samtools' 'sort' '-@' '16' '-o'", fixed = TRUE)
  expect_match(res$commands$samtools_index, "'samtools' 'index'", fixed = TRUE)
  expect_equal(res$paths$bai, paste0(normalizePath(bam, mustWork = FALSE), ".bai"))
})

test_that("align_contigs_to_reference supports conda_env and secondary alignments", {
  reference <- tempfile(fileext = ".fasta")
  assembly <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">ref", "ACGTACGT"), reference)
  writeLines(c(">contig1", "ACGT"), assembly)

  res <- align_contigs_to_reference(
    reference = reference,
    assembly = assembly,
    align_bam = bam,
    preset = "asm10",
    threads = 8,
    secondary = TRUE,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_match(res$commands$minimap2, "'conda' 'run' '-n' 'ont-tools' 'minimap2'", fixed = TRUE)
  expect_match(res$commands$samtools_sort, "'conda' 'run' '-n' 'ont-tools' 'samtools'", fixed = TRUE)
  expect_match(res$commands$samtools_index, "'conda' 'run' '-n' 'ont-tools' 'samtools'", fixed = TRUE)
  expect_match(res$commands$minimap2, "'-x' 'asm10'", fixed = TRUE)
  expect_match(res$commands$minimap2, "'-t' '8'", fixed = TRUE)
  expect_false(grepl("--secondary=no", res$commands$minimap2, fixed = TRUE))
})

test_that("align_contigs_to_reference runs minimap2 and samtools commands", {
  reference <- tempfile(fileext = ".fasta")
  assembly <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  fake_bin <- tempfile("align-bin-")
  dir.create(fake_bin)
  writeLines(c(">ref", "ACGTACGT"), reference)
  writeLines(c(">contig1", "ACGT"), assembly)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf '@HD\\tVN:1.6\\tSO:unsorted\\n'",
      "printf '@SQ\\tSN:ref\\tLN:8\\n'",
      "printf 'contig1\\t0\\tref\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\t*\\n'"
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
      "  sort)",
      "    out=''",
      "    input=''",
      "    while [[ $# -gt 0 ]]; do",
      "      case \"$1\" in",
      "        -o) out=\"$2\"; shift 2 ;;",
      "        -@) shift 2 ;;",
      "        *) input=\"$1\"; shift ;;",
      "      esac",
      "    done",
      "    cp \"$input\" \"$out\"",
      "    ;;",
      "  index)",
      "    touch \"$1.bai\"",
      "    ;;",
      "  *)",
      "    exit 1",
      "    ;;",
      "esac"
    ),
    file.path(fake_bin, "samtools")
  )
  Sys.chmod(file.path(fake_bin, c("minimap2", "samtools")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- align_contigs_to_reference(
    reference = reference,
    assembly = assembly,
    align_bam = bam,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(file.exists(res$paths$bam))
  expect_true(file.exists(res$paths$bai))
})

test_that("align_contigs_to_reference validates inputs", {
  reference <- tempfile(fileext = ".fasta")
  assembly <- tempfile(fileext = ".fasta")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">ref", "ACGTACGT"), reference)
  writeLines(c(">contig1", "ACGT"), assembly)

  expect_error(
    align_contigs_to_reference(reference, assembly, bam, threads = 0,
                               dry_run = TRUE, echo = FALSE),
    "threads"
  )
  expect_error(
    align_contigs_to_reference(reference, assembly, bam, secondary = NA,
                               dry_run = TRUE, echo = FALSE),
    "secondary"
  )
  expect_error(
    align_contigs_to_reference(reference, assembly, bam, conda_env = "",
                               dry_run = TRUE, echo = FALSE),
    "conda_env"
  )
  expect_error(
    align_contigs_to_reference("missing.fasta", assembly, bam,
                               dry_run = TRUE, echo = FALSE),
    "reference"
  )
})
