test_that("annotate_consensus_blast builds expected dry-run command", {
  query <- tempfile(fileext = ".fasta")
  writeLines(c(">consensus1", "ACGTACGT"), query)

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "unite_eukaryotes",
    out_dir = tempfile("blast-"),
    threads = 4,
    max_target_seqs = 10,
    task = "blastn",
    perc_identity = 80,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$status, NA_integer_)
  expect_match(res$commands$blastn, "'blastn'", fixed = TRUE)
  expect_match(res$commands$blastn, "'-db' 'unite_eukaryotes'", fixed = TRUE)
  expect_match(res$commands$blastn, "'-num_threads' '4'", fixed = TRUE)
  expect_match(res$commands$blastn, "'-max_target_seqs' '10'", fixed = TRUE)
  expect_match(res$commands$blastn, "'-perc_identity' '80'", fixed = TRUE)
  expect_null(res$blast)
})

test_that("annotate_consensus_blast dry-run supports db_fasta and conda_env", {
  query <- tempfile(fileext = ".fasta")
  db_fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">consensus1", "ACGTACGT"), query)
  writeLines(c(">ref1", "ACGTACGT"), db_fasta)

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db_fasta = db_fasta,
    out_dir = tempfile("blast-"),
    conda_env = "ont-tools",
    dry_run = TRUE,
    echo = FALSE
  )

  expect_equal(res$conda_env, "ont-tools")
  expect_match(res$commands$makeblastdb, "'conda' 'run' '-n' 'ont-tools' 'makeblastdb'", fixed = TRUE)
  expect_match(res$commands$blastn, "'conda' 'run' '-n' 'ont-tools' 'blastn'", fixed = TRUE)
  expect_match(res$commands$makeblastdb, "'-dbtype' 'nucl'", fixed = TRUE)
})

test_that("annotate_consensus_blast parses BLAST output and taxonomy", {
  query <- tempfile(fileext = ".fasta")
  out_dir <- tempfile("blast-")
  fake_bin <- tempfile("blast-bin-")
  dir.create(fake_bin)
  writeLines(c(">consensus1", "ACGTACGT"), query)

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
      "printf 'consensus1\\t1000\\t1\\t980\\tSH1\\t1000\\t1\\t980\\t980\\t99.1\\t98\\t2\\t0\\t1e-100\\t500\\tSH1|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae|reps\\t4932\\n' > \"$out\"",
      "printf 'consensus1\\t1000\\t1\\t900\\tSH2\\t1000\\t1\\t900\\t900\\t96.0\\t90\\t10\\t1\\t1e-80\\t400\\tSH2|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Kazachstania;s__Kazachstania_unispora|reps\\t564\\n' >> \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "unite_eukaryotes",
    out_dir = out_dir,
    echo = FALSE
  )

  expect_equal(res$status, 0L)
  expect_s3_class(res$blast, "data.frame")
  expect_equal(nrow(res$blast), 2L)
  expect_equal(res$blast$rank, c(1L, 2L))
  expect_equal(res$top_hits$sseqid, "SH1")
  expect_equal(res$top_hits$genus, "Saccharomyces")
  expect_equal(res$top_hits$species, "Saccharomyces cerevisiae")
  expect_false("annotation_level" %in% names(res$top_hits))
  expect_false("novel_candidate" %in% names(res$top_hits))
  expect_true(file.exists(res$paths$output_tsv))
  expect_match(readLines(res$paths$output_tsv, n = 1), "qseqid\tqlen\tqstart", fixed = TRUE)
})

test_that("annotate_consensus_blast fills taxonomy from NCBI taxdump staxids", {
  query <- tempfile(fileext = ".fasta")
  out_dir <- tempfile("blast-")
  fake_bin <- tempfile("blast-bin-")
  taxdump <- tempfile("taxdump-")
  dir.create(fake_bin)
  dir.create(taxdump)
  writeLines(c(">consensus1", "ACGTACGT"), query)
  writeLines(c(
    "1\t|\t1\t|\tno rank\t|",
    "2\t|\t1\t|\tsuperkingdom\t|",
    "1239\t|\t2\t|\tphylum\t|",
    "91061\t|\t1239\t|\tclass\t|",
    "1385\t|\t91061\t|\torder\t|",
    "186817\t|\t1385\t|\tfamily\t|",
    "1386\t|\t186817\t|\tgenus\t|",
    "1423\t|\t1386\t|\tspecies\t|"
  ), file.path(taxdump, "nodes.dmp"))
  writeLines(c(
    "1\t|\troot\t|\t\t|\tscientific name\t|",
    "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
    "1239\t|\tBacillota\t|\t\t|\tscientific name\t|",
    "91061\t|\tBacilli\t|\t\t|\tscientific name\t|",
    "1385\t|\tBacillales\t|\t\t|\tscientific name\t|",
    "186817\t|\tBacillaceae\t|\t\t|\tscientific name\t|",
    "1386\t|\tBacillus\t|\t\t|\tscientific name\t|",
    "1423\t|\tBacillus subtilis\t|\t\t|\tscientific name\t|"
  ), file.path(taxdump, "names.dmp"))

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
      "printf 'consensus1\\t1500\\t1\\t1490\\tNR_001\\t1500\\t1\\t1490\\t1490\\t99.5\\t99\\t1\\t0\\t1e-120\\t600\\tBacillus subtilis 16S ribosomal RNA\\t1423\\n' > \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "ncbi_16s_rRNA",
    out_dir = out_dir,
    taxdump_dir = taxdump,
    echo = FALSE
  )

  expect_equal(res$top_hits$kingdom, "Bacteria")
  expect_equal(res$top_hits$phylum, "Bacillota")
  expect_equal(res$top_hits$class, "Bacilli")
  expect_equal(res$top_hits$order, "Bacillales")
  expect_equal(res$top_hits$family, "Bacillaceae")
  expect_equal(res$top_hits$genus, "Bacillus")
  expect_equal(res$top_hits$species, "Bacillus subtilis")
  expect_match(res$top_hits$taxonomy_path, "Bacteria;Bacillota;Bacilli", fixed = TRUE)
})

test_that("annotate_consensus_blast passes BLAST outfmt as one argument", {
  query <- tempfile(fileext = ".fasta")
  out_dir <- tempfile("blast-")
  fake_bin <- tempfile("blast-bin-")
  captured <- file.path(fake_bin, "captured_args.txt")
  dir.create(fake_bin)
  writeLines(c(">consensus1", "ACGTACGT"), query)

  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "out=''",
      "outfmt_seen=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    -out) out=\"$2\"; shift 2 ;;",
      "    -outfmt) outfmt_seen=\"$2\"; printf '%s\\n' \"$2\" > \"$CAPTURED_ARGS\"; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "if [[ \"$outfmt_seen\" != *'qseqid qlen qstart'* ]]; then",
      "  echo \"bad outfmt: $outfmt_seen\" >&2",
      "  exit 9",
      "fi",
      ": > \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  old_capture <- Sys.getenv("CAPTURED_ARGS")
  on.exit({
    Sys.setenv(PATH = old_path)
    Sys.setenv(CAPTURED_ARGS = old_capture)
  }, add = TRUE)
  Sys.setenv(
    PATH = paste(fake_bin, old_path, sep = .Platform$path.sep),
    CAPTURED_ARGS = captured
  )

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "unite_eukaryotes",
    out_dir = out_dir,
    echo = FALSE
  )

  expect_equal(res$status, 0L)
  expect_match(readLines(captured), "qseqid qlen qstart", fixed = TRUE)
})

test_that("annotate_consensus_blast keeps low-identity top-hit metrics without automatic labels", {
  query <- tempfile(fileext = ".fasta")
  out_dir <- tempfile("blast-")
  fake_bin <- tempfile("blast-bin-")
  dir.create(fake_bin)
  writeLines(c(">consensus1", "ACGTACGT"), query)

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
      "printf 'consensus1\\t1000\\t1\\t950\\tSH3\\t1000\\t1\\t950\\t950\\t94.5\\t95\\t20\\t2\\t1e-60\\t300\\tSH3|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_sp|reps\\t0\\n' > \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "unite_eukaryotes",
    out_dir = out_dir,
    echo = FALSE
  )

  expect_equal(res$top_hits$sseqid, "SH3")
  expect_equal(res$top_hits$pident, 94.5)
  expect_equal(res$top_hits$query_coverage, 95)
  expect_false("annotation_level" %in% names(res$top_hits))
  expect_false("novel_candidate" %in% names(res$top_hits))
})

test_that("read_consensus_blast_table returns empty data frame for empty BLAST output", {
  query <- tempfile(fileext = ".fasta")
  out_dir <- tempfile("blast-")
  fake_bin <- tempfile("blast-bin-")
  dir.create(fake_bin)
  writeLines(c(">consensus1", "ACGTACGT"), query)

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
      ": > \"$out\""
    ),
    file.path(fake_bin, "blastn")
  )
  Sys.chmod(file.path(fake_bin, "blastn"), mode = "0755")

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  res <- annotate_consensus_blast(
    consensus_fasta = query,
    db = "unite_eukaryotes",
    out_dir = out_dir,
    echo = FALSE
  )

  expect_s3_class(res$blast, "data.frame")
  expect_equal(nrow(res$blast), 0L)
  expect_equal(nrow(res$top_hits), 0L)
})

test_that("annotate_consensus_blast validates arguments", {
  query <- tempfile(fileext = ".fasta")
  writeLines(c(">consensus1", "ACGTACGT"), query)

  expect_error(
    annotate_consensus_blast(query, dry_run = TRUE, echo = FALSE),
    "Supply either `db` or `db_fasta`"
  )
  expect_error(
    annotate_consensus_blast(query, db = "db", threads = 0,
                             dry_run = TRUE, echo = FALSE),
    "threads"
  )
  expect_error(
    annotate_consensus_blast(query, db = "db", extra_args = c("--x", NA),
                             dry_run = TRUE, echo = FALSE),
    "extra_args"
  )
})
