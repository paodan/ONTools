#' Collect amplicon workflow deliverables
#'
#' `collect_amplicon_results()` collects the main deliverables from a dorado /
#' wf_amplicon_denovo result directory into a clean final-results directory. It
#' copies consensus FASTA files, barcode result directories, a barcode-to-sample
#' map, and writes a README explaining the files.
#'
#' @param result_dir Path to a wf_amplicon_denovo result directory.
#' @param output_dir Path to the final-results directory to create.
#' @param sample_map Optional path to a barcode-to-sample mapping file. The file
#'   is copied into `output_dir`.
#' @param consensus_file Name of the untrimmed consensus FASTA in `result_dir`.
#' @param consensus_index_file Name of the FASTA index for `consensus_file`.
#' @param trimmed_consensus_file Name of the trimmed consensus FASTA in
#'   `result_dir`. The generated README recommends this file for downstream use.
#' @param barcode_pattern Regular expression used to identify barcode
#'   directories.
#' @param include_execution Logical. If `TRUE`, copy the `execution` directory
#'   when present.
#' @param include_reports Logical. If `TRUE`, copy common workflow report files
#'   such as `wf-amplicon-report.html`, `params.json`, and `versions.txt` when
#'   present.
#' @param readme_name README filename written into `output_dir`.
#' @param overwrite Logical. If `TRUE`, replace existing files in `output_dir`.
#'
#' @return Invisibly returns a data frame listing attempted copy operations,
#'   source paths, destination paths, item type, and whether each source existed.
#'
#' @examples
#' result_dir <- tempfile("amplicon-result-")
#' output_dir <- tempfile("amplicon-final-")
#' dir.create(result_dir)
#' dir.create(file.path(result_dir, "barcode01", "alignments"), recursive = TRUE)
#' writeLines(">barcode01\nACGT", file.path(result_dir, "all-consensus-seqs.fasta"))
#' writeLines("barcode01\t4\t0\t4\t5", file.path(result_dir, "all-consensus-seqs.fasta.fai"))
#' writeLines(">barcode01\nACGT", file.path(result_dir, "all-consensus-seqs_trimmed.fasta"))
#' writeLines("bam", file.path(result_dir, "barcode01", "alignments", "barcode01.bam"))
#' sample_map <- file.path(result_dir, "sample_map.tsv")
#' writeLines("barcode\tsample\nbarcode01\tsampleA", sample_map)
#'
#' collect_amplicon_results(result_dir, output_dir, sample_map = sample_map)
#'
#' @export
collect_amplicon_results <- function(result_dir,
                                     output_dir,
                                     sample_map = NULL,
                                     consensus_file = "all-consensus-seqs.fasta",
                                     consensus_index_file = paste0(consensus_file, ".fai"),
                                     trimmed_consensus_file = "all-consensus-seqs_trimmed.fasta",
                                     barcode_pattern = "^barcode[0-9]+$",
                                     include_execution = FALSE,
                                     include_reports = TRUE,
                                     readme_name = "README.txt",
                                     overwrite = TRUE) {
  check_dir_arg(result_dir, "result_dir")
  check_scalar_character(output_dir, "output_dir")
  check_scalar_character(consensus_file, "consensus_file")
  check_scalar_character(consensus_index_file, "consensus_index_file")
  check_scalar_character(trimmed_consensus_file, "trimmed_consensus_file")
  check_scalar_character(barcode_pattern, "barcode_pattern")
  check_scalar_character(readme_name, "readme_name")
  check_logical_scalar(include_execution, "include_execution")
  check_logical_scalar(include_reports, "include_reports")
  check_logical_scalar(overwrite, "overwrite")

  result_dir <- normalizePath(result_dir, mustWork = TRUE)
  if (!is.null(sample_map)) {
    check_file_arg(sample_map, "sample_map")
    sample_map <- normalizePath(sample_map, mustWork = TRUE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  barcode_dirs <- list.dirs(result_dir, recursive = FALSE, full.names = TRUE)
  barcode_dirs <- barcode_dirs[grepl(barcode_pattern, basename(barcode_dirs))]
  barcode_dirs <- sort(barcode_dirs)

  items <- data.frame(
    label = character(),
    source = character(),
    destination = character(),
    type = character(),
    exists = logical(),
    copied = logical(),
    stringsAsFactors = FALSE
  )

  items <- rbind(
    items,
    collect_item(
      label = "untrimmed_consensus",
      source = file.path(result_dir, consensus_file),
      destination = file.path(output_dir, consensus_file),
      type = "file",
      overwrite = overwrite
    ),
    collect_item(
      label = "consensus_fasta_index",
      source = file.path(result_dir, consensus_index_file),
      destination = file.path(output_dir, consensus_index_file),
      type = "file",
      overwrite = overwrite
    ),
    collect_item(
      label = "trimmed_consensus",
      source = file.path(result_dir, trimmed_consensus_file),
      destination = file.path(output_dir, trimmed_consensus_file),
      type = "file",
      overwrite = overwrite
    )
  )

  if (!is.null(sample_map)) {
    items <- rbind(
      items,
      collect_item(
        label = "sample_map",
        source = sample_map,
        destination = file.path(output_dir, basename(sample_map)),
        type = "file",
        overwrite = overwrite
      )
    )
  }

  for (barcode_dir in barcode_dirs) {
    items <- rbind(
      items,
      collect_item(
        label = basename(barcode_dir),
        source = barcode_dir,
        destination = file.path(output_dir, basename(barcode_dir)),
        type = "directory",
        overwrite = overwrite
      )
    )
  }

  if (isTRUE(include_execution)) {
    items <- rbind(
      items,
      collect_item(
        label = "execution",
        source = file.path(result_dir, "execution"),
        destination = file.path(output_dir, "execution"),
        type = "directory",
        overwrite = overwrite
      )
    )
  }

  if (isTRUE(include_reports)) {
    report_files <- c(
      "wf-amplicon-report.html",
      "params.json",
      "versions.txt"
    )
    for (report_file in report_files) {
      items <- rbind(
        items,
        collect_item(
          label = report_file,
          source = file.path(result_dir, report_file),
          destination = file.path(output_dir, report_file),
          type = "file",
          overwrite = overwrite
        )
      )
    }
  }

  readme_path <- file.path(output_dir, readme_name)
  writeLines(
    amplicon_results_readme(
      trimmed_consensus_file = trimmed_consensus_file,
      consensus_file = consensus_file,
      consensus_index_file = consensus_index_file,
      sample_map_name = if (is.null(sample_map)) NULL else basename(sample_map),
      barcode_names = basename(barcode_dirs),
      include_execution = include_execution,
      include_reports = include_reports
    ),
    readme_path
  )

  items <- rbind(
    items,
    data.frame(
      label = "README",
      source = NA_character_,
      destination = readme_path,
      type = "file",
      exists = TRUE,
      copied = TRUE,
      stringsAsFactors = FALSE
    )
  )

  invisible(items)
}

collect_item <- function(label, source, destination, type, overwrite) {
  exists <- if (identical(type, "directory")) {
    dir.exists(source)
  } else {
    file.exists(source)
  }

  copied <- FALSE
  if (exists) {
    if (identical(type, "directory")) {
      if (dir.exists(destination) && isTRUE(overwrite)) {
        unlink(destination, recursive = TRUE, force = TRUE)
      }
      copied <- file.copy(source, dirname(destination), recursive = TRUE)
    } else {
      copied <- file.copy(source, destination, overwrite = overwrite)
    }
  }

  data.frame(
    label = label,
    source = source,
    destination = destination,
    type = type,
    exists = exists,
    copied = copied,
    stringsAsFactors = FALSE
  )
}

amplicon_results_readme <- function(trimmed_consensus_file,
                                    consensus_file,
                                    consensus_index_file,
                                    sample_map_name,
                                    barcode_names,
                                    include_execution,
                                    include_reports) {
  lines <- c(
    "Amplicon sequencing final results",
    "=================================",
    "",
    "Recommended file",
    "----------------",
    paste0(
      "- ",
      trimmed_consensus_file,
      ": trimmed consensus FASTA. Use this file for downstream analyses; ",
      "extra bases outside the primer-bounded amplicon region have been removed while primer sequences are kept."
    ),
    "",
    "Other files",
    "-----------",
    paste0(
      "- ",
      consensus_file,
      ": original untrimmed consensus FASTA generated by the amplicon workflow."
    ),
    paste0(
      "- ",
      consensus_index_file,
      ": FASTA index for ",
      consensus_file,
      "."
    )
  )

  if (!is.null(sample_map_name)) {
    lines <- c(
      lines,
      paste0(
        "- ",
        sample_map_name,
        ": barcode-to-sample mapping table. Use this file to match barcode folders to sample names."
      )
    )
  }

  lines <- c(
    lines,
    "",
    "Barcode folders",
    "---------------"
  )

  if (length(barcode_names) > 0L) {
    lines <- c(
      lines,
      "- barcode*/: per-barcode results copied from the workflow output.",
      "- barcode*/consensus/consensus.fastq: per-barcode consensus FASTQ.",
      "- barcode*/alignments/*.bam: read alignments for the barcode.",
      "- barcode*/alignments/*.bam.bai: BAM index files.",
      "- barcode*/alignments/*.png: alignment snapshot images when generated.",
      "",
      "Included barcode folders:",
      paste0("- ", barcode_names)
    )
  } else {
    lines <- c(lines, "- No barcode folders were found in the source result directory.")
  }

  if (isTRUE(include_reports)) {
    lines <- c(
      lines,
      "",
      "Workflow reports",
      "----------------",
      "- wf-amplicon-report.html: workflow HTML report.",
      "- params.json: workflow parameter record.",
      "- versions.txt: software version record."
    )
  }

  if (isTRUE(include_execution)) {
    lines <- c(
      lines,
      "",
      "Execution metadata",
      "------------------",
      "- execution/: Nextflow execution report, timeline, and trace files when present."
    )
  }

  c(lines, "")
}

check_dir_arg <- function(x, name) {
  check_scalar_character(x, name)
  if (!dir.exists(x)) {
    stop("`", name, "` does not exist or is not a directory: ", x,
         call. = FALSE)
  }
}

check_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
}
