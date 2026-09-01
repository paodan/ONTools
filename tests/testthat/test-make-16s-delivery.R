test_that("parse_abundance_table returns long and wide data frames", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "tax\tsample1\tsample2",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus\t10\t0",
    "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia\t30\t20"
  ), tmp)

  long <- parse_abundance_table(tmp)
  wide <- parse_abundance_table(tmp, format = "wide")

  expect_s3_class(long, "data.frame")
  expect_s3_class(wide, "data.frame")
  expect_equal(nrow(long), 4)
  expect_true(all(c("tax", "samples", "count", "pct", "Kingdom", "Genus") %in% names(long)))
  expect_equal(long$pct[long$samples == "sample1" & long$Genus == "Bacillus"], 0.25)
  expect_equal(long$pct[long$samples == "sample1" & long$Genus == "Escherichia"], 0.75)
  expect_true(all(c("sample1", "sample2", "Phylum", "Genus") %in% names(wide)))
})

test_that("plot_abundance_bar writes an image", {
  skip_if_not(capabilities("png"))

  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "tax\tsample1\tsample2",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus\t10\t0",
    "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia\t30\t20"
  ), tmp)
  output <- tempfile(fileext = ".png")

  plot <- plot_abundance_bar(
    abundance_table_genus = tmp,
    output = output,
    fill = "Genus",
    width = 4,
    height = 3
  )

  expect_s3_class(plot, "ggplot")
  expect_true(file.exists(output))
})

test_that("move_16s creates delivery files without file.copy directory warnings", {
  skip_if_not(capabilities("png"))

  path_result <- tempfile()
  path_delivery <- tempfile()
  dir.create(file.path(path_result, "alignment_tables"), recursive = TRUE)
  writeLines(c(
    "tax\tsample1\tsample2",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus\t10\t0",
    "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia\t30\t20"
  ), file.path(path_result, "abundance_table_genus.tsv"))
  writeLines(
    "reference\tstartpos\tref length\tnumber of reads",
    file.path(path_result, "alignment_tables", "barcode001-alignment-stats.tsv")
  )

  res <- expect_warning(
    move_16s(
      path_result = path_result,
      path_delivery = path_delivery,
      overwrite = TRUE,
      tax_levels = c("Phylum", "Genus"),
      width = 4,
      height = 3
    ),
    NA
  )

  path_16s <- file.path(normalizePath(path_delivery), "16s")
  expect_true(file.exists(file.path(path_16s, "abundance_table_genus.tsv")))
  expect_true(file.exists(file.path(
    path_16s,
    "identification_tables",
    "barcode001-alignment-stats.tsv"
  )))
  expect_equal(length(res$plot_files), 4)
  expect_true(all(file.exists(res$plot_files)))
})

test_that("move_16s protects existing output and supports overwrite", {
  skip_if_not(capabilities("png"))

  path_result <- tempfile()
  path_delivery <- tempfile()
  dir.create(file.path(path_result, "alignment_tables"), recursive = TRUE)
  writeLines(c(
    "tax\tsample1",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus\t10"
  ), file.path(path_result, "abundance_table_genus.tsv"))
  writeLines(
    "reference\tstartpos",
    file.path(path_result, "alignment_tables", "barcode001-alignment-stats.tsv")
  )
  dir.create(file.path(path_delivery, "16s"), recursive = TRUE)

  expect_error(
    move_16s(path_result = path_result, path_delivery = path_delivery),
    "exists"
  )
  expect_warning(
    move_16s(
      path_result = path_result,
      path_delivery = path_delivery,
      overwrite = TRUE,
      tax_levels = "Genus",
      width = 4,
      height = 3
    ),
    NA
  )
})

test_that("make_16s_delivery points users to move_16s", {
  expect_error(make_16s_delivery(), "move_16s")
})
