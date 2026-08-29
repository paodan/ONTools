test_that("map_reads_to_assembly builds expected dry-run commands", {
  assembly <- tempfile(fileext = ".fasta")
  reads <- tempfile(fileext = ".fastq")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), assembly)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- map_reads_to_assembly(
    assembly = assembly,
    reads = reads,
    align_bam = bam,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_equal(res$preset, "map-ont")
  expect_match(res$commands$minimap2, "'minimap2' '-a' '-x' 'map-ont' '-t' '16'", fixed = TRUE)
  expect_match(res$commands$samtools_sort, "'samtools' 'sort' '-@' '16' '-o'", fixed = TRUE)
  expect_match(res$commands$samtools_index, "'samtools' 'index'", fixed = TRUE)
  expect_match(res$commands$samtools_depth, "'samtools' 'depth'", fixed = TRUE)
  expect_match(res$commands$plot_depth, "Plot depth in R", fixed = TRUE)
  expect_equal(res$paths$depth, paste0(sub("[.]bam$", "", normalizePath(bam, mustWork = FALSE)), ".depth.txt"))
  expect_equal(res$paths$depth_plot, paste0(res$paths$depth, ".png"))
})

test_that("map_reads_to_assembly supports conda_env and depth options", {
  assembly <- tempfile(fileext = ".fasta")
  reads <- tempfile(fileext = ".fastq")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), assembly)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- map_reads_to_assembly(
    assembly = assembly,
    reads = reads,
    align_bam = bam,
    threads = 8,
    secondary = FALSE,
    depth_all_positions = TRUE,
    min_mapping_quality = 20,
    min_base_quality = 10,
    plot_depth = FALSE,
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_match(res$commands$minimap2, "'conda' 'run' '-n' 'ont-tools' 'minimap2'", fixed = TRUE)
  expect_match(res$commands$samtools_depth, "'conda' 'run' '-n' 'ont-tools' 'samtools'", fixed = TRUE)
  expect_match(res$commands$minimap2, "--secondary=no", fixed = TRUE)
  expect_match(res$commands$samtools_depth, "'-a'", fixed = TRUE)
  expect_match(res$commands$samtools_depth, "'-Q' '20'", fixed = TRUE)
  expect_match(res$commands$samtools_depth, "'-q' '10'", fixed = TRUE)
  expect_null(res$paths$depth_plot)
  expect_null(res$commands$plot_depth)
})

test_that("map_reads_to_assembly prefers -aa over -a for depth", {
  assembly <- tempfile(fileext = ".fasta")
  reads <- tempfile(fileext = ".fastq")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), assembly)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  res <- map_reads_to_assembly(
    assembly = assembly,
    reads = reads,
    align_bam = bam,
    depth_all_positions = TRUE,
    depth_all_references = TRUE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_match(res$commands$samtools_depth, "'-aa'", fixed = TRUE)
  expect_false(grepl("'-a'", res$commands$samtools_depth, fixed = TRUE))
})

test_that("map_reads_to_assembly runs commands and plots depth", {
  assembly <- tempfile(fileext = ".fasta")
  reads <- tempfile(fileext = ".fastq")
  bam <- tempfile(fileext = ".bam")
  depth <- tempfile(fileext = ".depth.txt")
  depth_png <- tempfile(fileext = ".png")
  fake_bin <- tempfile("map-reads-bin-")
  dir.create(fake_bin)
  writeLines(c(">contig1", "ACGTACGT"), assembly)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "printf '@HD\\tVN:1.6\\tSO:unsorted\\n'",
      "printf '@SQ\\tSN:contig1\\tLN:8\\n'",
      "printf 'read1\\t0\\tcontig1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\t!!!!\\n'"
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
      "  depth)",
      "    printf 'contig1\\t1\\t7\\ncontig1\\t2\\t8\\ncontig1\\t3\\t9\\n'",
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

  res <- map_reads_to_assembly(
    assembly = assembly,
    reads = reads,
    align_bam = bam,
    depth_file = depth,
    depth_plot = depth_png,
    echo = FALSE,
    stderr = FALSE
  )

  expect_equal(res$status, 0L)
  expect_true(file.exists(res$paths$bam))
  expect_true(file.exists(res$paths$bai))
  expect_true(file.exists(res$paths$depth))
  expect_true(file.exists(res$paths$depth_plot))
  expect_s3_class(res$plot, "ggplot")
})

test_that("map_reads_to_assembly validates arguments", {
  assembly <- tempfile(fileext = ".fasta")
  reads <- tempfile(fileext = ".fastq")
  bam <- tempfile(fileext = ".bam")
  writeLines(c(">contig1", "ACGTACGT"), assembly)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)

  expect_error(
    map_reads_to_assembly(assembly, reads, bam, threads = 0,
                          dry_run = TRUE, echo = FALSE),
    "threads"
  )
  expect_error(
    map_reads_to_assembly(assembly, reads, bam, min_mapping_quality = -1,
                          dry_run = TRUE, echo = FALSE),
    "min_mapping_quality"
  )
  expect_error(
    map_reads_to_assembly(assembly, reads, bam, plot_width = 0,
                          dry_run = TRUE, echo = FALSE),
    "plot_width"
  )
  expect_error(
    map_reads_to_assembly("missing.fasta", reads, bam,
                          dry_run = TRUE, echo = FALSE),
    "assembly"
  )
})
