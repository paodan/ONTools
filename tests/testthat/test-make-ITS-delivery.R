make_fake_ITS_inputs <- function() {
  consensus_delivery <- tempfile("consensus-delivery-")
  its_result <- tempfile("its-result-")

  dir.create(file.path(consensus_delivery, "project", "barcode303"), recursive = TRUE)
  writeLines(
    c(">barcode303_original_consensus", "ACGTACGTACGT"),
    file.path(consensus_delivery, "project", "all-consensus-seqs.fasta")
  )
  writeLines(
    "barcode303_original_consensus\t12\t0\t12\t12",
    file.path(consensus_delivery, "project", "all-consensus-seqs.fasta.fai")
  )
  writeLines(
    c(">barcode303_trimmed_consensus", "ACGT"),
    file.path(consensus_delivery, "project", "all-consensus-seqs_trimmed.fasta")
  )
  writeLines(
    "barcode303_trimmed_consensus\t4\t0\t4\t4",
    file.path(consensus_delivery, "project", "all-consensus-seqs_trimmed.fasta.fai")
  )
  writeLines(
    "consensus summary",
    file.path(consensus_delivery, "project", "barcode303", "summary.txt")
  )

  dir.create(file.path(its_result, "alignment_tables"), recursive = TRUE)
  writeLines(c(
    "tax\tbarcode303\ttotal",
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t10\t10",
    "Eukaryota;Fungi;Mucoromycota;Mucoromycetes;Mucorales;Choanephoraceae;Poitrasia\t90\t90"
  ), file.path(its_result, "abundance_table_genus.tsv"))
  writeLines(
    paste(
      "reference", "startpos", "ref length", "number of reads", "covbases",
      "% coverage", "meandepth", "meanbaseq", "meanmapq", "taxid",
      "superkingdom", "kingdom", "phylum", "class", "order", "family",
      "genus", "species", "mean", "sd", "Coefficient of Variance", "pcreads",
      sep = "\t"
    ),
    file.path(its_result, "alignment_tables", "barcode303-alignment-stats.tsv")
  )

  list(consensus_delivery = consensus_delivery, its_result = its_result)
}

make_fake_grouped_ITS_inputs <- function() {
  consensus_delivery <- tempfile("consensus-delivery-")
  its_result <- tempfile("its-result-")
  sample_info_dir <- tempfile("sample-info-")
  dir.create(sample_info_dir)
  sample_info_303 <- file.path(sample_info_dir, "PROJECT001_ITS.csv")
  sample_info_304 <- file.path(sample_info_dir, "PROJECT002_ITS.csv")

  for (group in c("PROJECT001_ITS", "PROJECT002_ITS")) {
    barcode <- if (identical(group, "PROJECT001_ITS")) "barcode303" else "barcode304"
    dir.create(
      file.path(consensus_delivery, group, "consensus_results", barcode),
      recursive = TRUE
    )
    writeLines(
      paste("summary", barcode),
      file.path(consensus_delivery, group, "consensus_results", barcode, "summary.txt")
    )
    writeLines(
      c(paste0(">", barcode, "_original_consensus"), "ACGTACGTACGT"),
      file.path(
        consensus_delivery,
        group,
        "consensus_results",
        "all-consensus-seqs.fasta"
      )
    )
    writeLines(
      paste0(barcode, "_original_consensus\t12\t0\t12\t12"),
      file.path(
        consensus_delivery,
        group,
        "consensus_results",
        "all-consensus-seqs.fasta.fai"
      )
    )
    writeLines(
      c(paste0(">", barcode, "_trimmed_consensus"), "ACGT"),
      file.path(
        consensus_delivery,
        group,
        "consensus_results",
        "all-consensus-seqs_trimmed.fasta"
      )
    )
    writeLines(
      paste0(barcode, "_trimmed_consensus\t4\t0\t4\t4"),
      file.path(
        consensus_delivery,
        group,
        "consensus_results",
        "all-consensus-seqs_trimmed.fasta.fai"
      )
    )
  }

  dir.create(file.path(its_result, "alignment_tables"), recursive = TRUE)
  writeLines(c(
    "tax\tbarcode303\tbarcode304\ttotal",
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t10\t20\t30",
    "Eukaryota;Fungi;Mucoromycota;Mucoromycetes;Mucorales;Choanephoraceae;Poitrasia\t90\t80\t170"
  ), file.path(its_result, "abundance_table_genus.tsv"))
  for (barcode in c("barcode303", "barcode304")) {
    writeLines(
      paste(
        "reference", "startpos", "ref length", "number of reads", "covbases",
        "% coverage", "meandepth", "meanbaseq", "meanmapq", "taxid",
        "superkingdom", "kingdom", "phylum", "class", "order", "family",
        "genus", "species", "mean", "sd", "Coefficient of Variance", "pcreads",
        sep = "\t"
      ),
      file.path(its_result, "alignment_tables", paste0(barcode, "-alignment-stats.tsv"))
    )
  }

  utils::write.csv(
    data.frame(Barcode_ID = "PBC001-303", Sample_ID = "sample303"),
    sample_info_303,
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(Barcode_ID = "PBC001-304", Sample_ID = "sample304"),
    sample_info_304,
    row.names = FALSE
  )

  list(
    consensus_delivery = consensus_delivery,
    its_result = its_result,
    sample_info = c(
      PROJECT001_ITS = sample_info_303,
      PROJECT002_ITS = sample_info_304
    )
  )
}

test_that("make_ITS_delivery reports a dry-run plan", {
  plan <- make_ITS_delivery(
    path_ITS_result = "wf_its",
    path_delivery = "delivery_its",
    consensus_delivery_path = NULL,
    run_unite_annotation = TRUE,
    dry_run = TRUE
  )

  expect_equal(plan$status, "dry_run")
  expect_true(plan$need_ITS)
  expect_true(plan$need_consensus)
  expect_match(plan$ITS_plan$command_string, "--database_set' 'ncbi_16s_18s_28s_ITS", fixed = TRUE)
  expect_equal(plan$paths$effective_ITS_result, file.path(".", "delivery_its_wf_ITS"))
  expect_equal(plan$paths$samples, file.path("delivery_its", "samples"))
  expect_equal(plan$paths$readme$readme, file.path("delivery_its", "README.txt"))
  expect_equal(plan$paths$readme$readme_zh, file.path("delivery_its", "README.zh-CN.txt"))
  expect_true(plan$consensus_steps$run_basecalling_demux_step)
  expect_true(plan$consensus_steps$run_dorado_basecall_step)
  expect_true(plan$consensus_steps$run_dorado_demux_step)
  expect_true(plan$consensus_steps$run_dorado_fastq_step)
  expect_true(plan$consensus_steps$run_QC_step)
  expect_true(plan$ITS_steps$run_ITS_step)
  expect_true(plan$ITS_steps$move_ITS_step)
})

test_that("make_ITS_delivery passes unified ITS workflow arguments", {
  plan <- make_ITS_delivery(
    path_ITS_result = NULL,
    path_delivery = "delivery_its",
    consensus_delivery_path = "missing-consensus",
    fastq_out = "fastq_pass_trim",
    path_work = "work_root",
    out_dir = "custom_its_out",
    work_dir = "custom_its_work",
    profile = "docker",
    resume = FALSE,
    database_set = "ncbi_16s_18s_28s_ITS",
    min_len = 300,
    max_len = 900,
    extra_args = "--minimap2_by_reference",
    run_basecalling_demux_step = FALSE,
    run_dorado_basecall_step = FALSE,
    run_dorado_demux_step = FALSE,
    run_dorado_fastq_step = FALSE,
    run_QC_step = FALSE,
    move_fastq_step = FALSE,
    move_fastq_mode = "auto",
    run_amplicon_step = FALSE,
    trim_consensus_step = FALSE,
    run_filtered_QC_step = FALSE,
    run_igv_step = FALSE,
    collect_results_step = FALSE,
    make_ab1 = FALSE,
    dry_run = TRUE
  )

  expect_equal(plan$ITS_plan$paths$fastq, "fastq_pass_trim")
  expect_equal(plan$ITS_plan$paths$out_dir, "custom_its_out")
  expect_equal(plan$ITS_plan$paths$work_dir, "custom_its_work")
  expect_match(plan$ITS_plan$command_string, "-profile' 'docker", fixed = TRUE)
  expect_match(plan$ITS_plan$command_string, "--max_len' '900", fixed = TRUE)
  expect_false(any(plan$ITS_plan$args == "-resume"))
  expect_false(plan$consensus_steps$run_basecalling_demux_step)
  expect_false(plan$consensus_steps$run_dorado_basecall_step)
  expect_false(plan$consensus_steps$run_dorado_demux_step)
  expect_false(plan$consensus_steps$run_dorado_fastq_step)
  expect_false(plan$consensus_steps$run_QC_step)
  expect_false(plan$consensus_steps$move_fastq_step)
  expect_equal(plan$consensus_steps$move_fastq_mode, "auto")
  expect_false(plan$consensus_steps$run_amplicon_step)
  expect_false(plan$consensus_steps$trim_consensus_step)
  expect_false(plan$consensus_steps$run_filtered_QC_step)
  expect_false(plan$consensus_steps$run_igv_step)
  expect_false(plan$consensus_steps$collect_results_step)
  expect_false(plan$consensus_steps$make_ab1)
  expect_true(plan$ITS_steps$run_ITS_step)
  expect_true(plan$ITS_steps$move_ITS_step)
})

test_that("make_ITS_delivery rejects removed output_dir argument", {
  expect_error(
    make_ITS_delivery(output_dir = "old", dry_run = TRUE),
    "unused argument"
  )
})

test_that("make_ITS_delivery rejects removed ITS-prefixed workflow arguments", {
  expect_error(
    make_ITS_delivery(ITS_fastq = "fastq_pass_trim", dry_run = TRUE),
    "unused argument"
  )
  expect_error(
    make_ITS_delivery(ITS_out_dir = "wf_its", dry_run = TRUE),
    "unused argument"
  )
})

test_that("make_ITS_delivery can run ITS when path_ITS_result is missing", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_ITS_inputs()
  path_delivery <- tempfile("its-delivery-")
  its_out_dir <- tempfile("wf-its-")
  its_work_dir <- tempfile("wf-its-work-")
  fake_bin <- tempfile("nextflow-bin-")
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "out=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --out_dir) out=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "mkdir -p \"$out/alignment_tables\"",
      "printf 'tax\\tbarcode303\\ttotal\\nUnclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\\t1\\t1\\nEukaryota;Fungi;Mucoromycota;Mucoromycetes;Mucorales;Choanephoraceae;Poitrasia\\t9\\t9\\n' > \"$out/abundance_table_genus.tsv\"",
      "printf 'reference\\tstartpos\\tref length\\tnumber of reads\\tcovbases\\t% coverage\\tmeandepth\\tmeanbaseq\\tmeanmapq\\ttaxid\\tsuperkingdom\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tmean\\tsd\\tCoefficient of Variance\\tpcreads\\n' > \"$out/alignment_tables/barcode303-alignment-stats.tsv\""
    ),
    file.path(fake_bin, "nextflow")
  )
  Sys.chmod(file.path(fake_bin, "nextflow"), mode = "0755")

  res <- NULL
  expect_warning(
    res <- make_ITS_delivery(
      path_ITS_result = NULL,
      path_delivery = path_delivery,
      consensus_delivery_path = inputs$consensus_delivery,
      fastq_out = "fastq_pass_trim",
      out_dir = its_out_dir,
      work_dir = its_work_dir,
      nextflow = file.path(fake_bin, "nextflow"),
      tax_levels = "Genus",
      width = 4,
      height = 3,
      echo = FALSE
    ),
    "Skipping UNITE annotation"
  )

  expect_true(res$ITS_generated)
  expect_equal(res$ITS_result$status, 0L)
  expect_true(file.exists(file.path(its_out_dir, "abundance_table_genus.tsv")))
  expect_true(file.exists(file.path(path_delivery, "abundance_table_genus.tsv")))
  expect_true(file.exists(file.path(
    path_delivery,
    "samples",
    "barcode303",
    "ITS_results",
    "identification_tables",
    "barcode303-alignment-stats.tsv"
  )))
})

test_that("make_ITS_delivery requires wait when ITS must be generated", {
  inputs <- make_fake_ITS_inputs()

  expect_error(
    make_ITS_delivery(
      path_ITS_result = NULL,
      path_delivery = tempfile("its-delivery-"),
      consensus_delivery_path = inputs$consensus_delivery,
      wait = FALSE
    ),
    "Use `wait = TRUE`"
  )
})

test_that("make_ITS_delivery creates grouped delivery folders from sample info", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_grouped_ITS_inputs()
  path_delivery <- tempfile("its-delivery-")
  dir.create(path_delivery)

  res <- NULL
  expect_warning(
    res <- make_ITS_delivery(
      path_ITS_result = inputs$its_result,
      path_delivery = path_delivery,
      consensus_delivery_path = inputs$consensus_delivery,
      path_sampleInfo_file_list = inputs$sample_info,
      tax_levels = "Genus",
      width = 4,
      height = 3,
      run_unite_annotation = FALSE,
      echo = FALSE
    ),
    NA
  )

  expect_true(res$grouped_delivery)
  expect_true(all(c("PROJECT001_ITS", "PROJECT002_ITS") %in% names(res$delivery)))
  group1 <- file.path(path_delivery, "PROJECT001_ITS")
  group2 <- file.path(path_delivery, "PROJECT002_ITS")
  expect_true(file.exists(file.path(group1, "abundance_table_genus.tsv")))
  expect_true(file.exists(file.path(group2, "abundance_table_genus.tsv")))
  expect_true(file.exists(file.path(group1, "PROJECT001_ITS.csv")))
  expect_true(file.exists(file.path(group1, "all-consensus-seqs.fasta")))
  expect_true(file.exists(file.path(group1, "all-consensus-seqs.fasta.fai")))
  expect_true(file.exists(file.path(group1, "all-consensus-seqs_trimmed.fasta")))
  expect_true(file.exists(file.path(group1, "all-consensus-seqs_trimmed.fasta.fai")))
  expect_match(
    paste(readLines(file.path(group1, "all-consensus-seqs.fasta")), collapse = "\n"),
    "barcode303_trimmed_consensus"
  )
  expect_match(
    paste(readLines(file.path(group1, "all-consensus-seqs_trimmed.fasta")), collapse = "\n"),
    "barcode303_trimmed_consensus"
  )
  expect_true(file.exists(file.path(
    group1,
    "samples",
    "barcode303",
    "consensus_results",
    "summary.txt"
  )))
  expect_false(dir.exists(file.path(group1, "samples", "barcode304")))
  expect_true(file.exists(file.path(
    group1,
    "samples",
    "barcode303",
    "ITS_results",
    "identification_tables",
    "barcode303-alignment-stats.tsv"
  )))
  expect_false(file.exists(file.path(
    group1,
    "samples",
    "barcode304",
    "ITS_results",
    "identification_tables",
    "barcode304-alignment-stats.tsv"
  )))

  abundance1 <- utils::read.delim(
    file.path(group1, "abundance_table_genus.tsv"),
    sep = "\t",
    check.names = FALSE
  )
  expect_equal(names(abundance1), c("tax", "barcode303", "total"))
  expect_equal(abundance1$total, abundance1$barcode303)
})

test_that("make_ITS_delivery organizes consensus and ITS results under samples", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_ITS_inputs()
  path_delivery <- tempfile("its-delivery-")

  res <- NULL
  expect_warning(
    res <- make_ITS_delivery(
      path_ITS_result = inputs$its_result,
      path_delivery = path_delivery,
      consensus_delivery_path = inputs$consensus_delivery,
      tax_levels = "Genus",
      width = 4,
      height = 3,
      echo = FALSE
    ),
    "Skipping UNITE annotation"
  )

  expect_true(file.exists(file.path(path_delivery, "abundance_table_genus.tsv")))
  expect_true(dir.exists(file.path(path_delivery, "figures")))
  expect_true(file.exists(file.path(
    path_delivery,
    "samples",
    "barcode303",
    "consensus_results",
    "summary.txt"
  )))
  expect_true(file.exists(file.path(
    path_delivery,
    "samples",
    "barcode303",
    "ITS_results",
    "identification_tables",
    "barcode303-alignment-stats.tsv"
  )))
  expect_true(file.exists(file.path(path_delivery, "all-consensus-seqs.fasta")))
  expect_true(file.exists(file.path(path_delivery, "all-consensus-seqs.fasta.fai")))
  expect_true(file.exists(file.path(path_delivery, "all-consensus-seqs_trimmed.fasta")))
  expect_true(file.exists(file.path(path_delivery, "all-consensus-seqs_trimmed.fasta.fai")))
  expect_match(
    paste(readLines(file.path(path_delivery, "all-consensus-seqs.fasta")), collapse = "\n"),
    "barcode303_trimmed_consensus"
  )
  expect_match(
    paste(readLines(file.path(path_delivery, "all-consensus-seqs_trimmed.fasta")), collapse = "\n"),
    "barcode303_trimmed_consensus"
  )
  expect_true(file.exists(file.path(path_delivery, "README.txt")))
  expect_true(file.exists(file.path(path_delivery, "README.zh-CN.txt")))
  expect_null(res$unite_annotation)
})

test_that("make_ITS_delivery writes UNITE detailed and root top-hit tables", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_ITS_inputs()
  path_delivery <- tempfile("its-delivery-")
  fake_bin <- tempfile("blast-bin-")
  dir.create(fake_bin)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "out=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    -out) out=\"$2\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "printf 'barcode303_consensus\\t1000\\t1\\t990\\tSH001\\t1000\\t1\\t990\\t990\\t99.2\\t99\\t1\\t0\\t1e-100\\t500\\tSH001|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae|reps\\t4932\\n' > \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- NULL
  expect_warning(
    res <- make_ITS_delivery(
      path_ITS_result = inputs$its_result,
      path_delivery = path_delivery,
      consensus_delivery_path = inputs$consensus_delivery,
      tax_levels = "Genus",
      unite_db = "unite_eukaryotes",
      width = 4,
      height = 3,
      echo = FALSE
    ),
    NA
  )

  expect_true(file.exists(file.path(
    path_delivery,
    "unite_consensus_annotation",
    "consensus.blast.tsv"
  )))
  expect_true(file.exists(file.path(
    path_delivery,
    "unite_consensus_annotation",
    "consensus.top_hits.tsv"
  )))
  expect_true(file.exists(file.path(path_delivery, "unite_consensus_top_hits.tsv")))
  expect_s3_class(res$unite_annotation$top_hits, "data.frame")
  expect_equal(res$unite_annotation$top_hits$genus, "Saccharomyces")
})

test_that("write_ITS_delivery_readme can skip Chinese README", {
  path_delivery <- tempfile("its-delivery-")
  dir.create(path_delivery)

  files <- write_ITS_delivery_readme(path_delivery, chinese_readme_name = NULL)

  expect_true(file.exists(file.path(path_delivery, "README.txt")))
  expect_false(file.exists(file.path(path_delivery, "README.zh-CN.txt")))
  expect_equal(names(files), "README")
})
