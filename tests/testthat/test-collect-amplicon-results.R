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
  expect_true(file.exists(file.path(output_dir, "README.zh-CN.txt")))
  expect_true(any(copied$label == "README"))
  expect_true(any(copied$label == "README_zh_CN"))

  readme <- readLines(file.path(output_dir, "README.txt"))
  expect_true(any(grepl("Recommended file", readme)))
  expect_true(any(grepl("all-consensus-seqs_trimmed.fasta", readme)))
  expect_true(any(grepl("sample_map.tsv", readme)))
  expect_true(any(grepl("barcode57", readme)))
  expect_true(any(grepl("read length distribution", readme)))
  expect_true(any(grepl("Distribution_seqLength", readme)))
  expect_false(any(grepl("synthetic.ab1", readme)))
  expect_false(any(grepl("Workflow reports", readme)))
  expect_false(any(grepl("params.json", readme)))

  readme_zh <- readLines(file.path(output_dir, "README.zh-CN.txt"))
  expect_true(any(grepl("推荐使用的文件", readme_zh)))
  expect_true(any(grepl("all-consensus-seqs_trimmed.fasta", readme_zh)))
  expect_true(any(grepl("sample_map.tsv", readme_zh)))
  expect_true(any(grepl("barcode57", readme_zh)))
  expect_true(any(grepl("读长分布图", readme_zh)))
  expect_true(any(grepl("Distribution_seqLength", readme_zh)))
  expect_false(any(grepl("synthetic.ab1", readme_zh)))
})

test_that("collect_amplicon_results can mention synthetic AB1 files in README", {
  skip_if(Sys.which("samtools") == "", "samtools is not installed")

  result_dir <- tempfile("amplicon-result-")
  output_dir <- tempfile("amplicon-final-")
  dir.create(file.path(result_dir, "barcode57", "alignments"), recursive = TRUE)
  sequence <- "ACGTACGTACGT"
  consensus <- file.path(result_dir, "all-consensus-seqs.fasta")
  sam <- file.path(result_dir, "reads.sam")
  bam <- file.path(result_dir, "barcode57", "alignments", "barcode57.aligned.sorted.bam")
  writeLines(c(">barcode57", sequence), consensus)
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      paste0("@SQ\tSN:barcode57\tLN:", nchar(sequence)),
      paste0("read1\t0\tbarcode57\t1\t60\t12M\t*\t0\t0\t", sequence, "\tIIIIIIIIIIII")
    ),
    sam
  )
  system2("samtools", c("faidx", consensus))
  system2("samtools", c("view", "-bS", sam), stdout = bam)
  system2("samtools", c("index", bam))

  copied <- collect_amplicon_results(
    result_dir = result_dir,
    output_dir = output_dir,
    include_ab1 = TRUE,
    ab1_name_template = "{barcode}.synthetic.ab1",
    ab1_echo = FALSE,
    ab1_stderr = FALSE
  )

  readme <- readLines(file.path(output_dir, "README.txt"))
  readme_zh <- readLines(file.path(output_dir, "README.zh-CN.txt"))

  expect_true(any(grepl("barcode*/<barcode>.synthetic.ab1", readme, fixed = TRUE)))
  expect_true(any(grepl("synthetic Sanger-style AB1 chromatogram", readme, fixed = TRUE)))
  expect_true(any(grepl("barcode*/<barcode>.synthetic.ab1", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("模拟 Sanger AB1 峰图文件", readme_zh, fixed = TRUE)))
  expect_true(file.exists(file.path(output_dir, "barcode57", "barcode57.synthetic.ab1")))
  expect_true(any(copied$label == "barcode57_synthetic_ab1"))
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
  expect_true(file.exists(file.path(output_dir, "README.zh-CN.txt")))
})

test_that("collect_amplicon_results can skip the Chinese README", {
  result_dir <- tempfile("amplicon-result-")
  output_dir <- tempfile("amplicon-final-")
  dir.create(result_dir)

  copied <- collect_amplicon_results(
    result_dir = result_dir,
    output_dir = output_dir,
    chinese_readme_name = NULL
  )

  expect_true(file.exists(file.path(output_dir, "README.txt")))
  expect_false(file.exists(file.path(output_dir, "README.zh-CN.txt")))
  expect_true(any(copied$label == "README"))
  expect_false(any(copied$label == "README_zh_CN"))
})

test_that("collect_amplicon_results validates inputs", {
  expect_error(
    collect_amplicon_results("missing-dir", tempfile()),
    "result_dir"
  )
})
