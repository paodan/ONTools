#' Build a 16S delivery folder from wf-16s results
#'
#' `move_16s()` copies key wf-16s output files into a delivery folder and
#' generates abundance bar plots.
#'
#' @param path_result Path to a completed wf-16s result directory.
#' @param path_delivery Delivery root directory. The function creates a `16s/`
#'   subdirectory inside this path.
#' @param overwrite Logical. If `TRUE`, replace an existing `16s/` delivery
#'   directory.
#' @param tax_levels Taxonomic levels to plot.
#' @param abundance_table Filename of the wf-16s genus abundance table under
#'   `path_result`.
#' @param alignment_tables_dir Directory name of wf-16s per-barcode alignment
#'   tables under `path_result`.
#' @param figure_dir Name of the figures directory under the `16s/` delivery
#'   directory.
#' @param identification_dir Name of the copied alignment-table directory under
#'   the `16s/` delivery directory.
#' @param cutoff Minimum relative abundance kept in abundance plots.
#' @param width,height Plot width and height in inches.
#' @param readme_name English README filename written under the `16s/`
#'   delivery directory.
#' @param chinese_readme_name Chinese README filename written under the `16s/`
#'   delivery directory. Set to `NULL` to skip writing it.
#'
#' @return Invisibly returns a list with input paths, output paths, and generated
#'   ggplot objects.
#'
#' @export
move_16s <- function(path_result,
                     path_delivery = "/data/project_delivery",
                     overwrite = FALSE,
                     tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                     abundance_table = "abundance_table_genus.tsv",
                     alignment_tables_dir = "alignment_tables",
                     figure_dir = "figures",
                     identification_dir = "identification_tables",
                     cutoff = 0.01,
                     width = 12,
                     height = 6,
                     readme_name = "README.txt",
                     chinese_readme_name = "README.zh-CN.txt") {
  check_dir_arg(path_result, "path_result")
  check_scalar_character(path_delivery, "path_delivery")
  check_logical_scalar(overwrite, "overwrite")
  if (!is.character(tax_levels) || length(tax_levels) == 0L || anyNA(tax_levels)) {
    stop("`tax_levels` must be a non-empty character vector.", call. = FALSE)
  }
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(alignment_tables_dir, "alignment_tables_dir")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  cutoff <- validate_fraction(cutoff, "cutoff")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }

  path_result <- normalizePath(path_result, mustWork = TRUE)
  dir.create(path_delivery, recursive = TRUE, showWarnings = FALSE)
  path_delivery <- normalizePath(path_delivery, mustWork = TRUE)

  path_16s <- file.path(path_delivery, "16s")
  if (dir.exists(path_16s)) {
    if (!isTRUE(overwrite)) {
      stop(path_16s, " exists. Use `overwrite = TRUE` to replace it.",
           call. = FALSE)
    }
    unlink(path_16s, recursive = TRUE)
  }

  path_fig <- file.path(path_16s, figure_dir)
  path_identification <- file.path(path_16s, identification_dir)
  dir.create(path_fig, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_identification, showWarnings = FALSE, recursive = TRUE)

  abun <- file.path(path_result, abundance_table)
  tbls <- file.path(path_result, alignment_tables_dir)
  check_file_arg(abun, "abundance_table")
  check_dir_arg(tbls, "alignment_tables_dir")

  # abun_name <- tools::file_path_sans_ext(basename(abun))
  plots <- list()
  plot_files <- character()
  for (level in tax_levels) {
    percentage_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_percentage.png")
    )
    count_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_count.png")
    )

    plots[[paste0(level, "_percentage")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = percentage_file,
      fill = level,
      cutoff = cutoff,
      position = "fill",
      width = width,
      height = height
    )
    plots[[paste0(level, "_count")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = count_file,
      fill = level,
      cutoff = cutoff,
      position = "stack",
      width = width,
      height = height
    )
    plot_files <- c(
      plot_files,
      stats::setNames(percentage_file, paste0(level, "_percentage")),
      stats::setNames(count_file, paste0(level, "_count"))
    )
  }

  copied_tables <- copy_directory_contents(tbls, path_identification, overwrite = TRUE)
  copied_abundance <- file.copy(
    abun,
    file.path(path_16s, basename(abun)),
    overwrite = TRUE
  )
  if (!isTRUE(copied_abundance)) {
    stop("Failed to copy abundance table: ", abun, call. = FALSE)
  }

  readme_files <- write_16s_readme(
    path_16s = path_16s,
    abundance_table = basename(abun),
    figure_dir = figure_dir,
    identification_dir = identification_dir,
    readme_name = readme_name,
    chinese_readme_name = chinese_readme_name
  )

  invisible(list(
    path_result = path_result,
    path_delivery = path_delivery,
    path_16s = path_16s,
    path_fig = path_fig,
    path_identification = path_identification,
    abundance_table = file.path(path_16s, basename(abun)),
    copied_tables = copied_tables,
    plot_files = plot_files,
    readme_files = readme_files,
    plot = plots
  ))
}


#' Write README files for a 16S delivery folder
#'
#' @param path_16s Path to the 16S delivery directory.
#' @param abundance_table Name of the copied abundance table.
#' @param figure_dir Name of the figure directory under `path_16s`.
#' @param identification_dir Name of the per-barcode identification table
#'   directory under `path_16s`.
#' @param readme_name English README filename written under `path_16s`.
#' @param chinese_readme_name Chinese README filename written under `path_16s`.
#'   Set to `NULL` to skip writing it.
#'
#' @return Invisibly returns the generated README file paths.
#'
#' @export
write_16s_readme <- function(path_16s,
                             abundance_table = "abundance_table_genus.tsv",
                             figure_dir = "figures",
                             identification_dir = "identification_tables",
                             readme_name = "README.txt",
                             chinese_readme_name = "README.zh-CN.txt") {
  check_dir_arg(path_16s, "path_16s")
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }

  path_16s <- normalizePath(path_16s, mustWork = TRUE)
  output <- file.path(path_16s, readme_name)
  writeLines(
    s16_results_readme(
      abundance_table = abundance_table,
      figure_dir = figure_dir,
      identification_dir = identification_dir
    ),
    output,
    useBytes = TRUE
  )

  files <- c(README = output)
  if (!is.null(chinese_readme_name)) {
    output_zh <- file.path(path_16s, chinese_readme_name)
    writeLines(
      s16_results_readme_zh(
        abundance_table = abundance_table,
        figure_dir = figure_dir,
        identification_dir = identification_dir
      ),
      output_zh,
      useBytes = TRUE
    )
    files <- c(files, README_zh_CN = output_zh)
  }

  invisible(files)
}


#' Build a ITS delivery folder from wf-16s results
#'
#' `move_ITS()` copies key wf-16s output files (using specific params for ITS sequences)
#' into a delivery folder and generates abundance bar plots.
#'
#' @param path_result Path to a completed wf-16s result directory.
#' @param path_delivery Delivery root directory. The function creates a `ITS/`
#'   subdirectory inside this path.
#' @param overwrite Logical. If `TRUE`, replace an existing `ITS/` delivery
#'   directory.
#' @param tax_levels Taxonomic levels to plot.
#' @param abundance_table Filename of the wf-16s genus abundance table under
#'   `path_result`.
#' @param alignment_tables_dir Directory name of wf-16s per-barcode alignment
#'   tables under `path_result`.
#' @param figure_dir Name of the figures directory under the `ITS/` delivery
#'   directory.
#' @param identification_dir Name of the copied alignment-table directory under
#'   the `ITS/` delivery directory.
#' @param cutoff Minimum relative abundance kept in abundance plots.
#' @param width,height Plot width and height in inches.
#' @param readme_name English README filename written under the `ITS/`
#'   delivery directory.
#' @param chinese_readme_name Chinese README filename written under the `ITS/`
#'   delivery directory. Set to `NULL` to skip writing it.
#'
#' @return Invisibly returns a list with input paths, output paths, and generated
#'   ggplot objects.
#'
#' @export
move_ITS <- function(path_result,
                     path_delivery = "/data/project_delivery",
                     overwrite = FALSE,
                     tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                     abundance_table = "abundance_table_genus.tsv",
                     alignment_tables_dir = "alignment_tables",
                     figure_dir = "figures",
                     identification_dir = "identification_tables",
                     cutoff = 0.01,
                     width = 12,
                     height = 6,
                     readme_name = "README.txt",
                     chinese_readme_name = "README.zh-CN.txt") {
  check_dir_arg(path_result, "path_result")
  check_scalar_character(path_delivery, "path_delivery")
  check_logical_scalar(overwrite, "overwrite")
  if (!is.character(tax_levels) || length(tax_levels) == 0L || anyNA(tax_levels)) {
    stop("`tax_levels` must be a non-empty character vector.", call. = FALSE)
  }
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(alignment_tables_dir, "alignment_tables_dir")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  cutoff <- validate_fraction(cutoff, "cutoff")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }

  path_result <- normalizePath(path_result, mustWork = TRUE)
  dir.create(path_delivery, recursive = TRUE, showWarnings = FALSE)
  path_delivery <- normalizePath(path_delivery, mustWork = TRUE)

  path_ITS <- file.path(path_delivery, "ITS")
  if (dir.exists(path_ITS)) {
    if (!isTRUE(overwrite)) {
      stop(path_ITS, " exists. Use `overwrite = TRUE` to replace it.",
           call. = FALSE)
    }
    unlink(path_ITS, recursive = TRUE)
  }

  path_fig <- file.path(path_ITS, figure_dir)
  path_identification <- file.path(path_ITS, identification_dir)
  dir.create(path_fig, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_identification, showWarnings = FALSE, recursive = TRUE)

  abun <- file.path(path_result, abundance_table)
  tbls <- file.path(path_result, alignment_tables_dir)
  check_file_arg(abun, "abundance_table")
  check_dir_arg(tbls, "alignment_tables_dir")

  # abun_name <- tools::file_path_sans_ext(basename(abun))
  plots <- list()
  plot_files <- character()
  for (level in tax_levels) {
    percentage_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_percentage.png")
    )
    count_file <- file.path(
      path_fig,
      paste0("abundance_", level, "_count.png")
    )

    plots[[paste0(level, "_percentage")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = percentage_file,
      fill = level,
      cutoff = cutoff,
      position = "fill",
      width = width,
      height = height
    )
    plots[[paste0(level, "_count")]] <- plot_abundance_bar(
      abundance_table_genus = abun,
      output = count_file,
      fill = level,
      cutoff = cutoff,
      position = "stack",
      width = width,
      height = height
    )
    plot_files <- c(
      plot_files,
      stats::setNames(percentage_file, paste0(level, "_percentage")),
      stats::setNames(count_file, paste0(level, "_count"))
    )
  }

  copied_tables <- copy_directory_contents(tbls, path_identification, overwrite = TRUE)
  copied_abundance <- file.copy(
    abun,
    file.path(path_ITS, basename(abun)),
    overwrite = TRUE
  )
  if (!isTRUE(copied_abundance)) {
    stop("Failed to copy abundance table: ", abun, call. = FALSE)
  }

  readme_files <- write_ITS_readme(
    path_ITS = path_ITS,
    abundance_table = basename(abun),
    figure_dir = figure_dir,
    identification_dir = identification_dir,
    readme_name = readme_name,
    chinese_readme_name = chinese_readme_name
  )

  invisible(list(
    path_result = path_result,
    path_delivery = path_delivery,
    path_ITS = path_ITS,
    path_fig = path_fig,
    path_identification = path_identification,
    abundance_table = file.path(path_ITS, basename(abun)),
    copied_tables = copied_tables,
    plot_files = plot_files,
    readme_files = readme_files,
    plot = plots
  ))
}

#' Write README files for an ITS delivery folder
#'
#' @param path_ITS Path to the ITS delivery directory.
#' @param abundance_table Name of the copied abundance table.
#' @param figure_dir Name of the figure directory under `path_ITS`.
#' @param identification_dir Name of the per-barcode identification table
#'   directory under `path_ITS`.
#' @param readme_name English README filename written under `path_ITS`.
#' @param chinese_readme_name Chinese README filename written under `path_ITS`.
#'   Set to `NULL` to skip writing it.
#'
#' @return Invisibly returns the generated README file paths.
#'
#' @export
write_ITS_readme <- function(path_ITS,
                             abundance_table = "abundance_table_genus.tsv",
                             figure_dir = "figures",
                             identification_dir = "identification_tables",
                             readme_name = "README.txt",
                             chinese_readme_name = "README.zh-CN.txt") {
  check_dir_arg(path_ITS, "path_ITS")
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }

  path_ITS <- normalizePath(path_ITS, mustWork = TRUE)
  output <- file.path(path_ITS, readme_name)
  writeLines(
    its_results_readme(
      abundance_table = abundance_table,
      figure_dir = figure_dir,
      identification_dir = identification_dir
    ),
    output,
    useBytes = TRUE
  )

  files <- c(README = output)
  if (!is.null(chinese_readme_name)) {
    output_zh <- file.path(path_ITS, chinese_readme_name)
    writeLines(
      its_results_readme_zh(
        abundance_table = abundance_table,
        figure_dir = figure_dir,
        identification_dir = identification_dir
      ),
      output_zh,
      useBytes = TRUE
    )
    files <- c(files, README_zh_CN = output_zh)
  }

  invisible(files)
}


#' Parse a wf-16s genus abundance table
#'
#' @param abundance_table_genus Path to `abundance_table_genus.tsv`.
#' @param format Output format: `"long"` or `"wide"`.
#' @param include_pct Logical. If `TRUE` and `format = "long"`, add per-sample
#'   relative abundance in `pct`.
#' @param levels Taxonomic level names used to split the `tax` column.
#'
#' @return A data frame.
#'
#' @export
parse_abundance_table <- function(abundance_table_genus,
                                  format = c("long", "wide"),
                                  include_pct = TRUE,
                                  levels = c(
                                    "Superkingdom", "Kingdom", "Phylum",
                                    "Class", "Order", "Family", "Genus",
                                    "Species"
                                  )[1:7]) {
  check_file_arg(abundance_table_genus, "abundance_table_genus")
  format <- match.arg(format)
  check_logical_scalar(include_pct, "include_pct")
  if (!is.character(levels) || length(levels) == 0L || anyNA(levels)) {
    stop("`levels` must be a non-empty character vector.", call. = FALSE)
  }

  abun_data <- utils::read.delim(
    abundance_table_genus,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"tax" %in% names(abun_data)) {
    stop("`abundance_table_genus` must contain a `tax` column.", call. = FALSE)
  }

  sample_cols <- setdiff(names(abun_data), "tax")
  if (length(sample_cols) == 0L) {
    stop("`abundance_table_genus` must contain at least one sample column.",
         call. = FALSE)
  }

  annotations <- split_taxonomy(abun_data$tax, levels)
  if (identical(format, "wide")) {
    return(cbind(abun_data, annotations))
  }

  count_data <- stack_abundance_table(abun_data, sample_cols)
  if (isTRUE(include_pct)) {
    totals <- stats::setNames(
      colSums(abun_data[, sample_cols, drop = FALSE], na.rm = TRUE),
      sample_cols
    )
    count_data$pct <- ifelse(
      totals[count_data$samples] > 0,
      count_data$count / totals[count_data$samples],
      NA_real_
    )
  }

  annotation_index <- match(count_data$tax, abun_data$tax)
  cbind(count_data, annotations[annotation_index, , drop = FALSE])
}

#' Plot wf-16s abundance bars
#'
#' @param abundance_table_genus Path to `abundance_table_genus.tsv`.
#' @param output Output image path.
#' @param fill Taxonomic level used for bar fill.
#' @param cutoff Minimum relative abundance kept in the plot.
#' @param position Bar position. `"fill"` shows relative abundance and
#'   `"stack"` shows read counts.
#' @param levels Taxonomic level names used to split the `tax` column.
#' @param x_angle,x_hjust,x_vjust X-axis text angle and justification.
#' @param width,height Plot width and height in inches.
#'
#' @return Invisibly returns a ggplot object.
#'
#' @export
plot_abundance_bar <- function(abundance_table_genus,
                               output = paste0(
                                 tools::file_path_sans_ext(abundance_table_genus),
                                 ".png"
                               ),
                               fill = "Order",
                               cutoff = 0.01,
                               position = c("fill", "stack"),
                               levels = c(
                                 "Superkingdom", "Kingdom", "Phylum",
                                 "Class", "Order", "Family", "Genus",
                                 "Species"
                               )[1:7],
                               x_angle = 60,
                               x_hjust = 1,
                               x_vjust = 1,
                               width = 12,
                               height = 6) {
  check_scalar_character(output, "output")
  check_scalar_character(fill, "fill")
  position <- match.arg(position)
  cutoff <- validate_fraction(cutoff, "cutoff")
  width <- validate_positive_number(width, "width")
  height <- validate_positive_number(height, "height")

  abun_data <- parse_abundance_table(
    abundance_table_genus = abundance_table_genus,
    format = "long",
    include_pct = TRUE,
    levels = levels
  )
  if (!fill %in% names(abun_data)) {
    stop("`fill` must be one of: ", paste(names(abun_data), collapse = ", "),
         call. = FALSE)
  }

  abun_data <- abun_data[!is.na(abun_data$pct) & abun_data$pct >= cutoff, , drop = FALSE]
  plot <- ggplot2::ggplot(
    abun_data,
    ggplot2::aes(.data$samples, .data$count, fill = .data[[fill]])
  ) +
    ggplot2::geom_col(position = position) +
    seqqc_theme(x_angle = x_angle, x_hjust = x_hjust, x_vjust = x_vjust) +
    ggplot2::xlab(NULL)

  if (identical(position, "fill")) {
    plot <- plot +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
      ggplot2::ylab("Relative abundance (%)")
  } else {
    plot <- plot + ggplot2::ylab("Reads")
  }

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = output, plot = plot, width = width, height = height)
  invisible(plot)
}

make_16s_delivery <- function() {
  stop("`make_16s_delivery()` is not implemented. Use `move_16s()` instead.",
       call. = FALSE)
}

copy_directory_contents <- function(from, to, overwrite) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(from, all.files = FALSE, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0L) {
    return(character())
  }

  copied <- file.copy(files, to, recursive = TRUE, overwrite = overwrite)
  if (any(!copied)) {
    stop("Failed to copy: ", files[which(!copied)[1L]], call. = FALSE)
  }

  file.path(to, basename(files))
}

split_taxonomy <- function(tax, levels) {
  parts <- strsplit(as.character(tax), ";", fixed = TRUE)
  max_depth <- max(lengths(parts))
  if (max_depth > length(levels)) {
    stop(
      "There are ",
      max_depth,
      " levels in `tax`, but `levels` has only ",
      length(levels),
      " values.",
      call. = FALSE
    )
  }

  mat <- matrix(NA_character_, nrow = length(parts), ncol = length(levels))
  for (i in seq_along(parts)) {
    if (length(parts[[i]]) > 0L) {
      mat[i, seq_along(parts[[i]])] <- trimws(parts[[i]])
    }
  }
  colnames(mat) <- levels
  as.data.frame(mat, stringsAsFactors = FALSE)
}

stack_abundance_table <- function(abun_data, sample_cols) {
  rows <- lapply(sample_cols, function(sample) {
    data.frame(
      tax = abun_data$tax,
      samples = sample,
      count = suppressWarnings(as.numeric(abun_data[[sample]])),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

s16_results_readme <- function(abundance_table,
                              figure_dir,
                              identification_dir) {
  c(
    "16S Taxonomic Profiling Delivery",
    "",
    "Overview",
    "This directory contains the delivery files for 16S rRNA amplicon taxonomic profiling generated from wf-16s results. The workflow filters reads, assigns them to a reference taxonomy database, aggregates genus-level abundance, and reports per-barcode reference alignment summaries.",
    "",
    "Directory Contents",
    paste0("- ", abundance_table, ": genus-level abundance table. Rows are taxonomic paths and columns are sample/barcode read counts plus a total column when present."),
    paste0("- ", figure_dir, "/: abundance bar plots generated from the abundance table. Files ending in `_percentage.png` show relative abundance, and files ending in `_count.png` show read counts."),
    paste0("- ", identification_dir, "/: per-barcode alignment summary tables copied from wf-16s `alignment_tables/`. Each `barcode*-alignment-stats.tsv` file summarizes reference-level alignment evidence."),
    "",
    "Abundance Table",
    "The `tax` column is a semicolon-separated taxonomy path, usually in this order: superkingdom; kingdom; phylum; class; order; family; genus.",
    "Sample/barcode columns contain read counts assigned to each taxonomic path. Relative abundance can be calculated within each sample as: taxon reads / total reads in that sample.",
    "`Unclassified;Unknown;Unknown;...` represents reads that were not assigned to a genus in the abundance aggregation step, commonly because they did not pass the configured identity, reference-coverage, read-length, or read-quality filters.",
    "",
    "Per-Barcode Identification Tables",
    "The alignment tables are useful for reviewing which database references were hit by each barcode. Important columns include:",
    "- reference: database reference sequence identifier.",
    "- ref length: length of the reference sequence.",
    "- number of reads: reads aligned to that reference.",
    "- covbases and % coverage: covered bases and percentage of the reference covered by aligned reads.",
    "- meandepth: average depth across the reference length.",
    "- meanbaseq and meanmapq: average base quality and mapping quality.",
    "- taxid and taxonomy columns: taxonomy assigned to the database reference.",
    "- pcreads: percentage of sample reads represented by that reference hit.",
    "",
    "Reference Database",
    "When this delivery is generated from the default epi2me-labs/wf-16s settings, the workflow commonly uses the `ncbi_16s_18s` database. This database is suitable for routine 16S/18S composition profiling, initial taxonomic screening, and results that need to remain compatible with wf-16s outputs.",
    "`ncbi_16s_18s` should not be treated as the only evidence for strict species confirmation or novel-species assessment. For species-level confirmation, review the consensus/representative sequence independently against curated 16S resources, type-strain/type-material records, or other validated databases, and inspect the top 10-20 hits rather than only the single best hit.",
    "",
    "Notes for 16S Interpretation",
    "The abundance table is the main file for sample-level composition summaries. The alignment tables help with manual review, but they should not be interpreted as abundance tables because they summarize reference hits rather than final per-read taxonomic assignments.",
    "Species-level calls from 16S should be interpreted cautiously, especially when several closely related species have similar 16S sequences, coverage is low, or mapping quality is poor.",
    "For potential novel-species assessment, useful evidence includes a high-quality consensus sequence, high query and reference coverage, percent identity to the closest known references, whether multiple near-identical species are tied among top hits, and whether the closest hits are from type strains or curated reference material.",
    "",
    "Recommended Use",
    "Use the abundance table and figures for routine reporting. Use the per-barcode identification tables to review candidate taxa, low-abundance hits, high `Unknown` samples, and references with low coverage or low mapping quality."
  )
}

s16_results_readme_zh <- function(abundance_table,
                                 figure_dir,
                                 identification_dir) {
  c(
    "16S 物种注释结果说明",
    "",
    "概述",
    "本目录为 16S rRNA 扩增子物种注释结果交付目录。结果来自 wf-16s 流程：对 reads 进行质控过滤，基于参考分类数据库进行注释，汇总属水平丰度，并输出每个 barcode 的参考序列比对统计表。",
    "",
    "目录内容",
    paste0("- ", abundance_table, "：属水平丰度表。每一行为一个分类路径，每个样本或 barcode 对应一列 read 数；如果存在 total 列，则表示所有样本的合计 read 数。"),
    paste0("- ", figure_dir, "/：由丰度表生成的丰度柱状图。`_percentage.png` 表示相对丰度，`_count.png` 表示 read 数。"),
    paste0("- ", identification_dir, "/：每个 barcode 的参考序列比对统计表，来自 wf-16s 的 `alignment_tables/`。每个 `barcode*-alignment-stats.tsv` 文件汇总该 barcode 比对到各参考序列的证据。"),
    "",
    "丰度表说明",
    "`tax` 列为分号分隔的分类路径，通常顺序为：superkingdom; kingdom; phylum; class; order; family; genus。",
    "样本或 barcode 列为每个分类路径对应的 read 数。每个样本内的相对丰度可按以下方式计算：该分类 read 数 / 该样本总 read 数。",
    "`Unclassified;Unknown;Unknown;...` 表示这些 reads 在丰度汇总阶段没有被注释到属水平。常见原因包括未通过流程设置的 identity、reference coverage、读长或 read 质量阈值。",
    "",
    "每个 barcode 的注释表说明",
    "`identification_tables/` 中的表可用于查看每个 barcode 命中了哪些数据库参考序列。重要字段包括：",
    "- reference：数据库参考序列 ID。",
    "- ref length：参考序列长度。",
    "- number of reads：比对到该参考序列的 reads 数。",
    "- covbases 和 % coverage：参考序列上被覆盖的碱基数及覆盖比例。",
    "- meandepth：参考序列全长范围内的平均深度。",
    "- meanbaseq 和 meanmapq：平均碱基质量和平均比对质量。",
    "- taxid 及 taxonomy 相关列：该数据库参考序列对应的分类注释。",
    "- pcreads：该参考命中的 reads 在样本中的占比。",
    "",
    "参考数据库说明",
    "如果该结果来自 epi2me-labs/wf-16s 的默认设置，流程通常使用 `ncbi_16s_18s` 数据库。该数据库适合常规 16S/18S 组成分析、分类初筛，以及保持结果与 wf-16s 输出格式兼容。",
    "`ncbi_16s_18s` 不建议作为严格种水平确认或新物种判断的唯一证据。如果需要确认到种，建议把共识序列或代表序列单独与 curated 16S 数据库、type strain/type material 相关记录，或其他经过验证的数据库进行复核，并查看 top 10-20 hits，而不是只看单个 best hit。",
    "",
    "16S 结果解读注意事项",
    "丰度表是样本整体组成分析的主要结果。注释表适合人工复核候选分类，但不能直接当作丰度表使用，因为它汇总的是参考序列命中情况，而不是最终逐条 read 分类后的丰度。",
    "16S 的种水平注释需要谨慎解释，尤其是在近缘物种 16S 序列非常相似、覆盖度较低或 mapping quality 较低的情况下。",
    "如果需要评估潜在新物种，建议重点查看：共识序列质量、query coverage、reference coverage、与最接近已知参考序列的 percent identity、top hits 中是否有多个近缘种并列，以及最近命中是否来自 type strain 或高质量 curated reference。",
    "",
    "推荐使用方式",
    "常规报告建议使用丰度表和图片；当样本 Unknown 比例较高、存在低丰度命中，或某些参考序列覆盖度和比对质量较低时，再结合每个 barcode 的注释表进行人工复核。"
  )
}

its_results_readme <- function(abundance_table,
                               figure_dir,
                               identification_dir) {
  c(
    "ITS Taxonomic Profiling Delivery",
    "",
    "Overview",
    "This directory contains the delivery files for ITS amplicon taxonomic profiling generated from wf-16s-style results. The workflow aligns reads to a reference database, assigns taxonomy, aggregates genus-level abundance, and reports per-barcode reference alignment summaries.",
    "",
    "Directory Contents",
    paste0("- ", abundance_table, ": genus-level abundance table. Rows are taxonomic paths and columns are sample/barcode read counts plus a total column when present."),
    paste0("- ", figure_dir, "/: abundance bar plots generated from the abundance table. Files ending in `_percentage.png` show relative abundance, and files ending in `_count.png` show read counts."),
    paste0("- ", identification_dir, "/: per-barcode alignment summary tables copied from wf-16s `alignment_tables/`. Each `barcode*-alignment-stats.tsv` file summarizes reference-level alignment evidence."),
    "",
    "Abundance Table",
    "The `tax` column is a semicolon-separated taxonomy path, usually in this order: superkingdom; kingdom; phylum; class; order; family; genus.",
    "Sample/barcode columns contain read counts assigned to each taxonomic path. Relative abundance can be calculated within each sample as: taxon reads / total reads in that sample.",
    "`Unclassified;Unknown;Unknown;...` represents reads that were not assigned to a genus in the abundance aggregation step, often because no sufficiently confident database match passed the configured identity and reference-coverage thresholds.",
    "",
    "Per-Barcode Identification Tables",
    "The alignment tables are useful for reviewing which database references were hit by each barcode. Important columns include:",
    "- reference: database reference sequence identifier.",
    "- ref length: length of the reference sequence.",
    "- number of reads: reads aligned to that reference.",
    "- covbases and % coverage: covered bases and percentage of the reference covered by aligned reads.",
    "- meandepth: average depth across the reference length.",
    "- meanbaseq and meanmapq: average base quality and mapping quality.",
    "- taxid and taxonomy columns: taxonomy assigned to the database reference.",
    "- pcreads: percentage of sample reads represented by that reference hit.",
    "",
    "Reference Database",
    "For ITS analysis, UNITE or other curated ITS databases are generally more appropriate than generic 16S-focused databases. A broad eukaryote ITS database can help detect non-fungal eukaryotic sequences and possible contamination, while a fungi-focused ITS database is smaller, more targeted, and often easier to interpret for fungal samples.",
    "Database-based ITS annotation is suitable for routine composition summaries and candidate taxon screening, but it should not be used as the only evidence for strict species confirmation or novel-species assessment. For species confirmation, review consensus/representative sequences against curated ITS resources such as UNITE, inspect the top 10-20 hits, and consider whether the closest hits belong to the same Species Hypothesis or a set of closely related taxa.",
    "",
    "Notes for ITS Interpretation",
    "ITS lengths vary substantially across fungi and database records, so a read may align to a named reference but still fail the abundance-table filters if identity or reference coverage is below the selected thresholds.",
    "The alignment tables do not contain per-read identity. For novel-species assessment, inspect read-level statistics, representative sequences or contigs, alignment coverage, percent identity, top-hit versus second-hit separation, and database completeness rather than relying on the abundance table alone.",
    "For potential novel-species assessment, useful evidence includes a high-quality consensus sequence, high query coverage, adequate reference coverage for the matched database sequence, percent identity to the closest known references, top-hit versus second-hit separation, and whether the best matches are database representatives, type material, or well-curated records.",
    "",
    "Recommended Use",
    "Use the abundance table and figures for sample-level composition summaries. Use the per-barcode identification tables for manual review of candidate taxa, especially when `Unknown` is high or when evaluating potentially novel organisms."
  )
}

its_results_readme_zh <- function(abundance_table,
                                  figure_dir,
                                  identification_dir) {
  c(
    "ITS 物种注释结果说明",
    "",
    "概述",
    "本目录为 ITS 扩增子物种注释结果交付目录。结果来自 wf-16s 风格流程：将 reads 比对到参考数据库，进行分类注释，汇总属水平丰度，并输出每个 barcode 的参考序列比对统计表。",
    "",
    "目录内容",
    paste0("- ", abundance_table, "：属水平丰度表。每一行为一个分类路径，每个样本或 barcode 对应一列 read 数；如果存在 total 列，则表示所有样本的合计 read 数。"),
    paste0("- ", figure_dir, "/：由丰度表生成的丰度柱状图。`_percentage.png` 表示相对丰度，`_count.png` 表示 read 数。"),
    paste0("- ", identification_dir, "/：每个 barcode 的参考序列比对统计表，来自 wf-16s 的 `alignment_tables/`。每个 `barcode*-alignment-stats.tsv` 文件汇总该 barcode 比对到各参考序列的证据。"),
    "",
    "丰度表说明",
    "`tax` 列为分号分隔的分类路径，通常顺序为：superkingdom; kingdom; phylum; class; order; family; genus。",
    "样本或 barcode 列为每个分类路径对应的 read 数。每个样本内的相对丰度可按以下方式计算：该分类 read 数 / 该样本总 read 数。",
    "`Unclassified;Unknown;Unknown;...` 表示这些 reads 在丰度汇总阶段没有被注释到属水平。常见原因是没有足够可信的数据库匹配，或比对结果没有通过流程设置的 identity 与 reference coverage 阈值。",
    "",
    "每个 barcode 的注释表说明",
    "`identification_tables/` 中的表可用于查看每个 barcode 命中了哪些数据库参考序列。重要字段包括：",
    "- reference：数据库参考序列 ID。",
    "- ref length：参考序列长度。",
    "- number of reads：比对到该参考序列的 reads 数。",
    "- covbases 和 % coverage：参考序列上被覆盖的碱基数及覆盖比例。",
    "- meandepth：参考序列全长范围内的平均深度。",
    "- meanbaseq 和 meanmapq：平均碱基质量和平均比对质量。",
    "- taxid 及 taxonomy 相关列：该数据库参考序列对应的分类注释。",
    "- pcreads：该参考命中的 reads 在样本中的占比。",
    "",
    "参考数据库说明",
    "对于 ITS 分析，UNITE 或其他 curated ITS 数据库通常比偏 16S 的通用数据库更合适。较宽泛的真核 ITS 数据库有助于发现非真菌真核序列或潜在污染；而真菌专用 ITS 数据库更小、更聚焦，在真菌样本中通常更容易解释。",
    "基于数据库的 ITS 注释适合常规组成分析和候选分类筛查，但不建议作为严格种水平确认或新物种判断的唯一证据。如果需要确认到种，建议将共识序列或代表序列单独与 UNITE 等 curated ITS 资源进行复核，查看 top 10-20 hits，并关注最近命中是否属于同一个 Species Hypothesis 或同一组近缘分类。",
    "",
    "ITS 结果解读注意事项",
    "ITS 区域在不同真菌及数据库记录之间长度差异较大。因此，某条 read 即使能比对到带有属种名称的参考序列，也可能因为 identity 或 reference coverage 未达到阈值，而在丰度表中被归为 Unknown。",
    "`alignment-stats.tsv` 表不包含逐条 read 的 identity 信息。如果需要判断是否可能为新物种，建议结合逐 read 统计、代表序列或组装 contig、比对覆盖度、percent identity、最佳命中与次佳命中的差距，以及参考数据库完整性进行综合判断。",
    "如果需要评估潜在新物种，建议重点查看：共识序列质量、query coverage、对命中参考序列的 coverage、与最接近已知参考序列的 percent identity、最佳命中与次佳命中的差距，以及最近命中是否来自数据库代表序列、type material 或高质量 curated record。",
    "",
    "推荐使用方式",
    "丰度表和图片适合用于样本整体组成展示；每个 barcode 的注释表适合用于人工复核候选分类，尤其适用于 Unknown 比例较高或需要评估潜在新物种的样本。"
  )
}
