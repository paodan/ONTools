test_that("plot_species_pie aggregates species and groups rare entries", {
  alignment_table <- data.frame(
    `number of reads` = c(80, 15, 3, 2),
    pcreads = c(80, 15, 3, 2),
    genuspecies = c("Species A", "Species B", "Species C", "Species D"),
    check.names = FALSE
  )

  plot <- plot_species_pie(alignment_table, otherPct = 0.05)
  plot_data <- attr(plot, "plot_data")

  expect_s3_class(plot, "ggplot")
  expect_equal(as.character(plot_data$Species), c("Species A", "Species B", "Others"))
  expect_equal(plot_data$value, c(80, 15, 5))
  expect_equal(plot_data$reads, c(80, 15, 5))
})

test_that("plot_species_pie supports percent mode and explicit columns", {
  alignment_table <- data.frame(
    reads = c(20, 10),
    pct = c(66.7, 33.3),
    taxon = c("Species A", "Species B"),
    check.names = FALSE
  )

  plot <- plot_species_pie(
    alignment_table,
    value = "percent",
    species_col = "taxon",
    reads_col = "reads",
    percent_col = "pct",
    top_n = 1
  )
  plot_data <- attr(plot, "plot_data")

  expect_equal(as.character(plot_data$Species), c("Species A", "Others"))
  expect_equal(plot_data$value, c(66.7, 33.3))
})

test_that("plot_species_pie reads TSV files and writes output images", {
  skip_if_not(capabilities("png"))

  alignment_file <- tempfile(fileext = ".tsv")
  output <- tempfile(fileext = ".png")
  writeLines(c(
    "reference\tnumber of reads\tpcreads\tgenuspecies",
    "NR_1\t10\t90\tSpecies A",
    "NR_2\t1\t10\tSpecies B"
  ), alignment_file)

  plot <- plot_species_pie(
    alignment_file,
    showPercent = TRUE,
    title = "barcode001",
    label_slices = TRUE,
    output = output,
    width = 4,
    height = 3
  )

  expect_s3_class(plot, "ggplot")
  expect_true(file.exists(output))
})

test_that("plot_species_pie validates inputs", {
  expect_error(plot_species_pie(data.frame(species = "A")), "reads_col")
  expect_error(
    plot_species_pie(
      data.frame(species = "A", `number of reads` = 1, check.names = FALSE),
      value = "percent"
    ),
    "percent_col"
  )
  expect_error(plot_species_pie(data.frame(species = "A", reads = 0)), "No positive")
})
