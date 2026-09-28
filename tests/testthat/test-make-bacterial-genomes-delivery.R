test_that("make_bacterial_genomes_delivery copies files without renaming", {
  input <- tempfile("bg-input-")
  output <- tempfile("bg-output-")
  dir.create(file.path(input, "Bam"), recursive = TRUE)
  dir.create(file.path(input, "QC"), recursive = TRUE)
  dir.create(file.path(input, "Sequence"), recursive = TRUE)
  dir.create(file.path(input, "Var"), recursive = TRUE)
  dir.create(output)

  writeLines("bam", file.path(input, "Bam", "G22510280341-G418.bam"))
  writeLines("bai", file.path(input, "Bam", "G22510280341-G418.bam.bai"))
  writeLines("png", file.path(input, "QC", "G22510280341-G418.coverage.png"))
  writeLines(c(">consensus", "ACGT"), file.path(input, "Sequence", "G22510280341-G418.consensus.fasta"))
  writeLines("ab1", file.path(input, "Sequence", "G22510280341-G418.ab1"))
  writeLines(
    "CHROM\tPOS\tREF\tALT\tvariant_percent",
    file.path(input, "Var", "G22510280341-G418.var.xls")
  )
  writeLines(
    "CHROM\tPOS\tREF\tALT\tvariant_percent",
    file.path(input, "Var", "G22510280341-G418.filt.var.xls")
  )

  res <- make_bacterial_genomes_delivery(
    input_dir = input,
    output = output,
    project = "PROJECT001",
    overwrite = TRUE,
    make_archive = FALSE,
    echo = FALSE
  )

  delivery <- res$paths$delivery_dir
  expect_true(file.exists(file.path(delivery, "Bam", "G22510280341-G418.bam")))
  expect_true(file.exists(file.path(delivery, "Bam", "G22510280341-G418.bam.bai")))
  expect_true(file.exists(file.path(delivery, "QC", "G22510280341-G418.coverage.png")))
  expect_true(file.exists(file.path(delivery, "Sequence", "G22510280341-G418.consensus.fasta")))
  expect_true(file.exists(file.path(delivery, "Sequence", "G22510280341-G418.ab1")))
  expect_true(file.exists(file.path(delivery, "Var", "G22510280341-G418.var.xls")))
  expect_true(file.exists(file.path(delivery, "Var", "G22510280341-G418.filt.var.xls")))
  expect_true(file.exists(file.path(delivery, "merged_data.var.xls")))
  expect_true(file.exists(file.path(delivery, "merged_data.filt.var.xls")))
  expect_true(file.exists(file.path(delivery, "manifest.tsv")))
  expect_true(file.exists(file.path(delivery, "04_md5", "md5.txt")))

  manifest <- utils::read.delim(file.path(delivery, "manifest.tsv"), check.names = FALSE)
  expect_true("variant_table" %in% manifest$label)
  expect_true("merged_variant_table" %in% manifest$label)
})

test_that("make_bacterial_genomes_delivery generates variant tables from VCF", {
  input <- tempfile("bg-input-")
  output <- tempfile("bg-output-")
  dir.create(input)
  dir.create(output)

  vcf <- file.path(input, "G22510280342-PKNTG.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tG22510280342-PKNTG",
    "contig1\t10\t.\tA\tG\t60\tPASS\tSR=8,2,5,5;AR=0,0\tGT:DP\t0/1:20",
    "contig1\t20\t.\tC\tT\t60\tPASS\tSR=99,1,0,1;AR=0,0\tGT:DP\t0/1:101"
  ), vcf)

  res <- make_bacterial_genomes_delivery(
    input_dir = input,
    output = output,
    project = "PROJECT002",
    overwrite = TRUE,
    make_archive = FALSE,
    echo = FALSE
  )

  delivery <- res$paths$delivery_dir
  var_file <- file.path(delivery, "Var", "G22510280342-PKNTG.var.xls")
  filt_file <- file.path(delivery, "Var", "G22510280342-PKNTG.filt.var.xls")
  expect_true(file.exists(var_file))
  expect_true(file.exists(filt_file))

  var <- utils::read.delim(var_file, check.names = FALSE)
  filt <- utils::read.delim(filt_file, check.names = FALSE)
  expect_equal(nrow(var), 2)
  expect_equal(nrow(filt), 1)
  expect_equal(
    names(var),
    c("Chr", "Pos", "Ref", "Alt", "DP", "Ref_dp", "Alt_dp", "Freq", "DP4", "Seq", "Supplement")
  )
  expect_equal(var$Chr, c("contig1", "contig1"))
  expect_equal(var$Ref, c("A", "C"))
  expect_equal(var$Alt, c("G", "T"))
  expect_equal(var$Ref_dp, c(10L, 100L))
  expect_equal(var$Alt_dp, c(10L, 1L))
  expect_equal(var$Freq, c(0.5, 0.0099), tolerance = 0.00001)
  expect_equal(var$Seq, c("[A/G]", "[C/T]"))
  expect_true(file.exists(file.path(delivery, "merged_data.var.xls")))
  expect_true(file.exists(file.path(delivery, "merged_data.filt.var.xls")))
})

test_that("make_bacterial_genomes_delivery supports dry-run and archive", {
  input <- tempfile("bg-input-")
  output <- tempfile("bg-output-")
  dir.create(input)
  dir.create(output)
  writeLines(c(">consensus", "ACGT"), file.path(input, "sample.consensus.fasta"))

  dry <- make_bacterial_genomes_delivery(
    input_dir = input,
    output = output,
    project = "PROJECT003",
    dry_run = TRUE,
    echo = FALSE
  )
  expect_false(dir.exists(dry$paths$delivery_dir))
  expect_equal(dry$status, NA_integer_)
  expect_true(any(dry$plan$destination == file.path(dry$paths$delivery_dir, "Sequence", "sample.consensus.fasta")))

  res <- make_bacterial_genomes_delivery(
    input_dir = input,
    output = output,
    project = "PROJECT003",
    overwrite = TRUE,
    make_archive = TRUE,
    echo = FALSE
  )
  expect_true(file.exists(res$paths$archive))
})

test_that("make_bacterial_genomes_delivery validates duplicate delivery names", {
  input <- tempfile("bg-input-")
  output <- tempfile("bg-output-")
  dir.create(file.path(input, "a"), recursive = TRUE)
  dir.create(file.path(input, "b"), recursive = TRUE)
  dir.create(output)
  writeLines("bam-a", file.path(input, "a", "sample.bam"))
  writeLines("bam-b", file.path(input, "b", "sample.bam"))

  expect_error(
    make_bacterial_genomes_delivery(
      input_dir = input,
      output = output,
      project = "PROJECT004",
      dry_run = TRUE,
      echo = FALSE
    ),
    "same delivery filename"
  )
})

test_that("make_bacterial_genomes_delivery can collect workflow outputs after FASTQ analysis", {
  fastq <- tempfile("bg-fastq-")
  wf_out <- tempfile("bg-wf-")
  output <- tempfile("bg-output-")
  dir.create(file.path(fastq, "barcode001"), recursive = TRUE)
  dir.create(file.path(wf_out, "sample1"), recursive = TRUE)
  dir.create(output)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq, "barcode001", "barcode001.fastq.gz"))

  fasta_gz <- file.path(wf_out, "sample1", "sample1.medaka.fasta.gz")
  con <- gzfile(fasta_gz, open = "wt")
  writeLines(c(">sample1", "ACGTACGT"), con)
  close(con)
  writeLines("bam", file.path(wf_out, "sample1", "sample1.bam"))
  writeLines("bai", file.path(wf_out, "sample1", "sample1.bam.bai"))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    "sample1\t2\t.\tC\tT\t60\tPASS\tSR=8,2,5,5;AR=0,0\tGT:DP\t0/1:20"
  ), file.path(wf_out, "sample1", "sample1.vcf"))

  res <- make_bacterial_genomes_delivery(
    fastq = fastq,
    output = output,
    project = "PROJECT005",
    run_wf = FALSE,
    wf_out_dir = wf_out,
    overwrite = TRUE,
    make_archive = FALSE,
    echo = FALSE
  )

  delivery <- res$paths$delivery_dir
  expect_true(file.exists(file.path(delivery, "Sequence", "sample1.consensus.fasta")))
  expect_true(file.exists(file.path(delivery, "Bam", "sample1.bam")))
  expect_true(file.exists(file.path(delivery, "Bam", "sample1.bam.bai")))
  expect_true(file.exists(file.path(delivery, "Var", "sample1.var.xls")))
  expect_true(file.exists(file.path(delivery, "Var", "sample1.filt.var.xls")))
  expect_true(file.exists(file.path(delivery, "merged_data.var.xls")))
  expect_equal(res$collection$status, 0L)
})

test_that("make_bacterial_genomes_delivery can call variants from consensus and BAM", {
  fastq <- tempfile("bg-fastq-")
  wf_out <- tempfile("bg-wf-")
  output <- tempfile("bg-output-")
  bin <- tempfile("bg-bin-")
  dir.create(file.path(fastq, "barcode001"), recursive = TRUE)
  dir.create(file.path(wf_out, "sample1"), recursive = TRUE)
  dir.create(output)
  dir.create(bin)
  writeLines(c("@read1", "ACGT", "+", "!!!!"), file.path(fastq, "barcode001", "barcode001.fastq.gz"))

  fasta_gz <- file.path(wf_out, "sample1", "sample1.medaka.fasta.gz")
  con <- gzfile(fasta_gz, open = "wt")
  writeLines(c(">sample1", "ACGTACGT"), con)
  close(con)
  writeLines("bam", file.path(wf_out, "sample1", "sample1.bam"))
  writeLines("bai", file.path(wf_out, "sample1", "sample1.bam.bai"))

  fake_samtools <- file.path(bin, "samtools")
  writeLines(c(
    "#!/bin/sh",
    "if [ \"$1\" = \"faidx\" ]; then",
    "  printf 'sample1\\t8\\t9\\t8\\t9\\n' > \"$2.fai\"",
    "  exit 0",
    "fi",
    "exit 0"
  ), fake_samtools)
  Sys.chmod(fake_samtools, "0755")

  fake_bcftools <- file.path(bin, "bcftools")
  writeLines(c(
    "#!/bin/sh",
    "cmd=\"$1\"",
    "shift",
    "if [ \"$cmd\" = \"mpileup\" ]; then",
    "  out=''",
    "  while [ \"$#\" -gt 0 ]; do",
    "    if [ \"$1\" = \"-o\" ]; then out=\"$2\"; shift 2; else shift; fi",
    "  done",
    "  printf 'bcf\\n' > \"$out\"",
    "  exit 0",
    "fi",
    "if [ \"$cmd\" = \"call\" ]; then",
    "  out=''",
    "  while [ \"$#\" -gt 0 ]; do",
    "    if [ \"$1\" = \"-o\" ]; then out=\"$2\"; shift 2; else shift; fi",
    "  done",
    "  {",
    "    printf '##fileformat=VCFv4.2\\n'",
    "    printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tsample1\\n'",
    "    printf 'sample1\\t2\\t.\\tC\\tT\\t60\\tPASS\\tDP4=6,4,5,5\\tGT:DP\\t0/1:20\\n'",
    "  } > \"$out\"",
    "  exit 0",
    "fi",
    "if [ \"$cmd\" = \"index\" ]; then",
    "  printf 'index\\n' > \"$1.csi\"",
    "  exit 0",
    "fi",
    "exit 1"
  ), fake_bcftools)
  Sys.chmod(fake_bcftools, "0755")

  res <- make_bacterial_genomes_delivery(
    fastq = fastq,
    output = output,
    project = "PROJECT006",
    run_wf = FALSE,
    wf_out_dir = wf_out,
    run_variant_calling = TRUE,
    samtools = fake_samtools,
    bcftools = fake_bcftools,
    overwrite = TRUE,
    make_archive = FALSE,
    echo = FALSE
  )

  delivery <- res$paths$delivery_dir
  expect_true(file.exists(file.path(delivery, "Var", "sample1.var.xls")))
  expect_true(file.exists(file.path(delivery, "Var", "sample1.filt.var.xls")))
  var <- utils::read.delim(file.path(delivery, "Var", "sample1.var.xls"), check.names = FALSE)
  expect_equal(var$DP4, "DP4=6,4,5,5")
  expect_equal(var$Ref_dp, 10L)
  expect_equal(var$Alt_dp, 10L)
  expect_equal(var$Freq, 0.5)
  expect_equal(res$collection$variants$sample1$status, 0L)
})
