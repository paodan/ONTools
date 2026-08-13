test_that("collect_amplicon_results copies deliverables and writes README", {
  result_dir <- tempfile("amplicon-result-")
  output_dir <- tempfile("amplicon-final-")
  dir.create(result_dir)
  dir.create(file.path(result_dir, "barcode57", "alignments"), recursive = TRUE)
  dir.create(file.path(result_dir, "barcode57", "consensus"), recursive = TRUE)
  dir.create(file.path(result_dir, "barcode58", "alignments"), recursive = TRUE)
  dir.create(file.path(result_dir, "execution"), recursive = TRUE)
  fastq_pass_trim_dir <- tempfile("fastq-pass-trim-")
  dir.create(file.path(fastq_pass_trim_dir, "barcode57"), recursive = TRUE)
  dir.create(file.path(fastq_pass_trim_dir, "barcode58"), recursive = TRUE)

  writeLines(">barcode57\nACGT", file.path(result_dir, "all-consensus-seqs.fasta"))
  writeLines(
    "barcode57\t4\t0\t4\t5",
    file.path(result_dir, "all-consensus-seqs.fasta.fai")
  )
  writeLines(
    ">barcode57\nACGT",
    file.path(result_dir, "all-consensus-seqs_trimmed.fasta")
  )
  writeLines("bam", file.path(result_dir, "barcode57", "alignments", "barcode57.aligned.sorted.bam"))
  writeLines("bai", file.path(result_dir, "barcode57", "alignments", "barcode57.aligned.sorted.bam.bai"))
  writeLines("png", file.path(result_dir, "barcode57", "alignments", "barcode57.png"))
  writeLines("fastq", file.path(result_dir, "barcode57", "consensus", "consensus.fastq"))
  writeLines("report", file.path(result_dir, "wf-amplicon-report.html"))
  writeLines("{}", file.path(result_dir, "params.json"))
  writeLines("versions", file.path(result_dir, "versions.txt"))
  writeLines("trace", file.path(result_dir, "execution", "trace.txt"))
  writeLines(
    "read length png",
    file.path(
      fastq_pass_trim_dir,
      "barcode57",
      "Distribution_seqLength__sup__barcode57.png"
    )
  )
  writeLines(
    "read length png",
    file.path(
      fastq_pass_trim_dir,
      "barcode58",
      "Distribution_seqLength__sup__barcode58.png"
    )
  )
  sample_map <- file.path(result_dir, "sample_map.tsv")
  writeLines("barcode\tsample\nbarcode57\tsampleA", sample_map)

  copied <- collect_amplicon_results(
    result_dir = result_dir,
    output_dir = output_dir,
    sample_map = sample_map,
    fastq_pass_trim_dir = fastq_pass_trim_dir,
    include_execution = TRUE
  )

  expect_true(file.exists(file.path(output_dir, "all-consensus-seqs.fasta")))
  expect_true(file.exists(file.path(output_dir, "all-consensus-seqs.fasta.fai")))
  expect_true(file.exists(file.path(output_dir, "all-consensus-seqs_trimmed.fasta")))
  expect_true(file.exists(file.path(output_dir, "sample_map.tsv")))
  expect_true(dir.exists(file.path(output_dir, "barcode57")))
  expect_true(dir.exists(file.path(output_dir, "barcode58")))
  expect_true(file.exists(file.path(
    output_dir,
    "barcode57",
    "Distribution_seqLength__sup__barcode57.png"
  )))
  expect_true(file.exists(file.path(
    output_dir,
    "barcode58",
    "Distribution_seqLength__sup__barcode58.png"
  )))
  expect_true(dir.exists(file.path(output_dir, "execution")))
  expect_false(file.exists(file.path(output_dir, "wf-amplicon-report.html")))
  expect_false(file.exists(file.path(output_dir, "params.json")))
  expect_false(file.exists(file.path(output_dir, "versions.txt")))
  expect_true(file.exists(file.path(output_dir, "README.txt")))
  expect_true(any(copied$label == "README"))

  readme <- readLines(file.path(output_dir, "README.txt"))
  expect_true(any(grepl("Recommended file", readme)))
  expect_true(any(grepl("all-consensus-seqs_trimmed.fasta", readme)))
  expect_true(any(grepl("sample_map.tsv", readme)))
  expect_true(any(grepl("barcode57", readme)))
  expect_true(any(grepl("read length distribution", readme)))
  expect_true(any(grepl("Distribution_seqLength", readme)))
  expect_false(any(grepl("Workflow reports", readme)))
  expect_false(any(grepl("params.json", readme)))
})

test_that("collect_amplicon_results reports missing optional files", {
  result_dir <- tempfile("amplicon-result-")
  output_dir <- tempfile("amplicon-final-")
  dir.create(result_dir)
  writeLines(
    ">barcode01\nACGT",
    file.path(result_dir, "all-consensus-seqs_trimmed.fasta")
  )

  copied <- collect_amplicon_results(
    result_dir = result_dir,
    output_dir = output_dir
  )

  expect_false(copied$exists[copied$label == "untrimmed_consensus"])
  expect_true(copied$exists[copied$label == "trimmed_consensus"])
  expect_true(file.exists(file.path(output_dir, "README.txt")))
})

test_that("collect_amplicon_results validates inputs", {
  expect_error(
    collect_amplicon_results("missing-dir", tempfile()),
    "result_dir"
  )
})
