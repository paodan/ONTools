test_that("make_circular_consensus_delivery builds the default command", {
  consensus <- tempfile(fileext = ".fasta")
  output_dir <- tempfile("consensus-delivery-")
  sample_sheet <- tempfile(fileext = ".csv")
  bam <- tempfile(fileext = ".bam")
  bai <- tempfile(fileext = ".bam.bai")
  depth <- tempfile(fileext = ".depth.txt")
  variants_vcf <- tempfile(fileext = ".vcf.gz")
  variants_vcf_index <- tempfile(fileext = ".vcf.gz.csi")
  variants_af <- tempfile(fileext = ".af.tsv")
  variants_af_filtered <- tempfile(fileext = ".gt0.05.af.tsv")
  assembly_info <- tempfile("assembly_info-")
  flye_log <- tempfile("flye-", fileext = ".log")
  notes <- tempfile("notes-", fileext = ".txt")

  dir.create(output_dir)
  writeLines(c(">consensus1", "ACGTACGT"), consensus)
  writeLines("sample,barcode\nsampleA,barcode01", sample_sheet)
  writeLines("bam", bam)
  writeLines("bai", bai)
  writeLines("depth", depth)
  writeLines("vcf", variants_vcf)
  writeLines("vcf index", variants_vcf_index)
  writeLines("CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD\tALT_DEPTH\tALT_AF", variants_af)
  writeLines("CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD\tALT_DEPTH\tALT_AF", variants_af_filtered)
  writeLines("assembly", assembly_info)
  writeLines("flye log", flye_log)
  writeLines("notes", notes)

  res <- make_circular_consensus_delivery(
    consensus = consensus,
    output = output_dir,
    project = "PROJECT001",
    sequence_type = "plasmid",
    sample_sheet = sample_sheet,
    bam = bam,
    bai = bai,
    depth = depth,
    variants_vcf = variants_vcf,
    variants_vcf_index = variants_vcf_index,
    variants_af = variants_af,
    variants_af_filtered = variants_af_filtered,
    assembly_info = assembly_info,
    flye_log = flye_log,
    notes = notes,
    overwrite = TRUE,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$execution_command, "bash")
  expect_equal(res$status, NA_integer_)
  expect_true(any(res$args == "--consensus"))
  expect_true(any(res$args == normalizePath(consensus)))
  expect_true(any(res$args == "--output"))
  expect_true(any(res$args == normalizePath(output_dir)))
  expect_true(any(res$args == "--project"))
  expect_true(any(res$args == "PROJECT001"))
  expect_true(any(res$args == "--sequence-type"))
  expect_true(any(res$args == "plasmid"))
  expect_true(any(res$args == "--sample-sheet"))
  expect_true(any(res$args == normalizePath(sample_sheet)))
  expect_true(any(res$args == "--bam"))
  expect_true(any(res$args == normalizePath(bam)))
  expect_true(any(res$args == "--bai"))
  expect_true(any(res$args == normalizePath(bai)))
  expect_true(any(res$args == "--depth"))
  expect_true(any(res$args == normalizePath(depth)))
  expect_true(any(res$args == "--variants-vcf"))
  expect_true(any(res$args == normalizePath(variants_vcf)))
  expect_true(any(res$args == "--variants-vcf-index"))
  expect_true(any(res$args == normalizePath(variants_vcf_index)))
  expect_true(any(res$args == "--variants-af"))
  expect_true(any(res$args == normalizePath(variants_af)))
  expect_true(any(res$args == "--variants-af-filtered"))
  expect_true(any(res$args == normalizePath(variants_af_filtered)))
  expect_true(any(res$args == "--assembly-info"))
  expect_true(any(res$args == normalizePath(assembly_info)))
  expect_true(any(res$args == "--flye-log"))
  expect_true(any(res$args == normalizePath(flye_log)))
  expect_true(any(res$args == "--notes"))
  expect_true(any(res$args == normalizePath(notes)))
  expect_true(any(res$args == "--overwrite"))
  expect_match(res$command_string, "make_circular_consensus_delivery.sh", fixed = TRUE)
  expect_equal(
    res$paths$delivery_dir,
    file.path(normalizePath(output_dir), "PROJECT001_consensus_delivery")
  )
  expect_equal(
    res$paths$archive,
    file.path(normalizePath(output_dir), "PROJECT001_consensus_delivery.tar.gz")
  )
})

test_that("make_circular_consensus_delivery validates arguments", {
  consensus <- tempfile(fileext = ".fasta")
  output_dir <- tempfile("consensus-delivery-")
  dir.create(output_dir)
  writeLines(c(">consensus1", "ACGTACGT"), consensus)

  expect_error(
    make_circular_consensus_delivery(
      consensus = file.path(output_dir, "missing.fasta"),
      output = output_dir,
      project = "PROJECT001",
      dry_run = TRUE,
      echo = FALSE
    ),
    "consensus"
  )
  expect_error(
    make_circular_consensus_delivery(
      consensus = consensus,
      output = output_dir,
      project = "PROJECT001",
      sequence_type = "linear",
      dry_run = TRUE,
      echo = FALSE
    ),
    "arg"
  )
  expect_error(
    make_circular_consensus_delivery(
      consensus = consensus,
      output = output_dir,
      project = "PROJECT001",
      overwrite = "yes",
      dry_run = TRUE,
      echo = FALSE
    ),
    "overwrite"
  )
})

test_that("make_circular_consensus_delivery builds a delivery package", {
  consensus <- tempfile(fileext = ".fasta")
  output_dir <- tempfile("consensus-delivery-")
  fake_bin <- tempfile("consensus-delivery-bin-")
  sample_sheet <- tempfile(fileext = ".csv")
  notes <- tempfile(fileext = ".txt")
  variants_af <- tempfile(fileext = ".af.tsv")
  variants_af_filtered <- tempfile(fileext = ".gt0.05.af.tsv")

  dir.create(output_dir)
  dir.create(fake_bin)
  writeLines(c(">consensus1", "ACGTACGT"), consensus)
  writeLines("sample,barcode\nsampleA,barcode01", sample_sheet)
  writeLines("A short project note.", notes)
  writeLines(
    "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD\tALT_DEPTH\tALT_AF",
    variants_af
  )
  writeLines(
    "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD\tALT_DEPTH\tALT_AF",
    variants_af_filtered
  )

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\n'",
      "for f in \"$@\"; do",
      "  case \"$f\" in",
      "    --*) ;;",
      "    *) [[ -f \"$f\" ]] && printf '%s\\tFASTA\\tDNA\\t1\\t8\\n' \"$f\" ;;",
      "  esac",
      "done"
    ),
    file.path(fake_bin, "seqkit")
  )
  Sys.chmod(file.path(fake_bin, "seqkit"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  make_circular_consensus_delivery(
    consensus = consensus,
    output = output_dir,
    project = "PROJECT001",
    sequence_type = "virus",
    sample_sheet = sample_sheet,
    variants_af = variants_af,
    variants_af_filtered = variants_af_filtered,
    notes = notes,
    script = system.file(
      "scripts",
      "make_circular_consensus_delivery.sh",
      package = "ONTools"
    ),
    overwrite = TRUE,
    echo = FALSE,
    stdout = FALSE,
    stderr = FALSE
  )

  delivery_dir <- file.path(output_dir, "PROJECT001_consensus_delivery")
  expect_true(dir.exists(delivery_dir))
  expect_true(file.exists(file.path(output_dir, "PROJECT001_consensus_delivery.tar.gz")))
  expect_equal(
    readLines(file.path(delivery_dir, "01_consensus", "consensus.fasta")),
    c(">consensus1", "ACGTACGT")
  )
  expect_true(file.exists(file.path(delivery_dir, "00_metadata", basename(sample_sheet))))
  expect_true(file.exists(file.path(delivery_dir, "00_metadata", basename(notes))))
  expect_true(file.exists(file.path(delivery_dir, "03_qc", "consensus_stats.tsv")))
  expect_true(file.exists(file.path(delivery_dir, "03_qc", "variants.af.tsv")))
  expect_true(file.exists(file.path(delivery_dir, "03_qc", "variants_gt0.05.af.tsv")))
  expect_true(file.exists(file.path(delivery_dir, "04_md5", "md5.txt")))

  manifest <- utils::read.delim(file.path(delivery_dir, "manifest.tsv"))
  readme_zh <- readLines(file.path(delivery_dir, "README.zh-CN.txt"))

  expect_true("consensus_fasta" %in% manifest$label)
  expect_true("consensus_stats" %in% manifest$label)
  expect_true("variants_af" %in% manifest$label)
  expect_true("variants_af_filtered" %in% manifest$label)
  expect_true(any(grepl("项目：PROJECT001", readme_zh, fixed = TRUE)))
  expect_true(any(grepl("序列类型：virus", readme_zh, fixed = TRUE)))
})
