make_fake_ITS_inputs <- function() {
  consensus_delivery <- tempfile("consensus-delivery-")
  its_result <- tempfile("its-result-")

  dir.create(file.path(consensus_delivery, "project", "barcode303"), recursive = TRUE)
  writeLines(
    c(">barcode303_consensus", "ACGTACGTACGT"),
    file.path(consensus_delivery, "project", "all-consensus-seqs_trimmed.fasta")
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

test_that("make_ITS_delivery reports a dry-run plan", {
  plan <- make_ITS_delivery(
    path_ITS_result = "wf_its",
    output_dir = "delivery_its",
    consensus_delivery_path = NULL,
    run_unite_annotation = TRUE,
    dry_run = TRUE
  )

  expect_equal(plan$status, "dry_run")
  expect_true(plan$need_consensus)
  expect_equal(plan$paths$samples, file.path("delivery_its", "samples"))
  expect_equal(plan$paths$readme$readme, file.path("delivery_its", "README.txt"))
  expect_equal(plan$paths$readme$readme_zh, file.path("delivery_its", "README.zh-CN.txt"))
})

test_that("make_ITS_delivery organizes consensus and ITS results under samples", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_ITS_inputs()
  output_dir <- tempfile("its-delivery-")

  res <- expect_warning(
    make_ITS_delivery(
      path_ITS_result = inputs$its_result,
      output_dir = output_dir,
      consensus_delivery_path = inputs$consensus_delivery,
      tax_levels = "Genus",
      width = 4,
      height = 3,
      echo = FALSE
    ),
    "Skipping UNITE annotation"
  )

  expect_true(file.exists(file.path(output_dir, "abundance_table_genus.tsv")))
  expect_true(dir.exists(file.path(output_dir, "figures")))
  expect_true(file.exists(file.path(
    output_dir,
    "samples",
    "barcode303",
    "consensus_results",
    "summary.txt"
  )))
  expect_true(file.exists(file.path(
    output_dir,
    "samples",
    "barcode303",
    "identification_tables",
    "barcode303-alignment-stats.tsv"
  )))
  expect_true(file.exists(file.path(output_dir, "consensus_for_unite.fasta")))
  expect_true(file.exists(file.path(output_dir, "README.txt")))
  expect_true(file.exists(file.path(output_dir, "README.zh-CN.txt")))
  expect_null(res$unite_annotation)
})

test_that("make_ITS_delivery writes UNITE detailed and root top-hit tables", {
  skip_if_not(capabilities("png"))

  inputs <- make_fake_ITS_inputs()
  output_dir <- tempfile("its-delivery-")
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

  res <- expect_warning(
    make_ITS_delivery(
      path_ITS_result = inputs$its_result,
      output_dir = output_dir,
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
    output_dir,
    "unite_consensus_annotation",
    "consensus.blast.tsv"
  )))
  expect_true(file.exists(file.path(
    output_dir,
    "unite_consensus_annotation",
    "consensus.top_hits.tsv"
  )))
  expect_true(file.exists(file.path(output_dir, "unite_consensus_top_hits.tsv")))
  expect_s3_class(res$unite_annotation$top_hits, "data.frame")
  expect_equal(res$unite_annotation$top_hits$genus, "Saccharomyces")
})

test_that("write_ITS_delivery_readme can skip Chinese README", {
  output_dir <- tempfile("its-delivery-")
  dir.create(output_dir)

  files <- write_ITS_delivery_readme(output_dir, chinese_readme_name = NULL)

  expect_true(file.exists(file.path(output_dir, "README.txt")))
  expect_false(file.exists(file.path(output_dir, "README.zh-CN.txt")))
  expect_equal(names(files), "README")
})
