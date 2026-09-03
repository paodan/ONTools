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
#' @param fastq_pass_trim_dir Optional directory containing per-barcode
#'   `fastq_pass_trim` folders. When supplied, per-barcode read-length
#'   distribution images are copied into the corresponding final barcode
#'   folders.
#' @param sample_length_plot_pattern Regular expression used to identify
#'   per-barcode read-length distribution image files inside each
#'   `fastq_pass_trim/<barcode>` folder.
#' @param include_execution Logical. If `TRUE`, copy the `execution` directory
#'   when present.
#' @param include_ab1 Logical. If `TRUE`, mention synthetic AB1 chromatogram
#'   files in the generated README files. The files are expected to already be
#'   present in source barcode directories and are copied as part of those
#'   directories.
#' @param ab1_name_template AB1 filename template used in README files. The
#'   placeholder `{barcode}` is replaced with each barcode directory name.
#' @param readme_name README filename written into `output_dir`.
#' @param chinese_readme_name Chinese README filename written into `output_dir`.
#'   Set to `NULL` to skip writing the Chinese README.
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
                                     fastq_pass_trim_dir = NULL,
                                     sample_length_plot_pattern = "^Distribution_seqLength__.*\\.png$",
                                     include_execution = FALSE,
                                     include_ab1 = FALSE,
                                     ab1_name_template = "{barcode}.synthetic.ab1",
                                     readme_name = "README.txt",
                                     chinese_readme_name = "README.zh-CN.txt",
                                     overwrite = TRUE) {
  check_dir_arg(result_dir, "result_dir")
  check_scalar_character(output_dir, "output_dir")
  check_scalar_character(consensus_file, "consensus_file")
  check_scalar_character(consensus_index_file, "consensus_index_file")
  check_scalar_character(trimmed_consensus_file, "trimmed_consensus_file")
  check_scalar_character(barcode_pattern, "barcode_pattern")
  check_scalar_character(sample_length_plot_pattern, "sample_length_plot_pattern")
  check_scalar_character(ab1_name_template, "ab1_name_template")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }
  check_logical_scalar(include_execution, "include_execution")
  check_logical_scalar(include_ab1, "include_ab1")
  check_logical_scalar(overwrite, "overwrite")
  if (isTRUE(include_ab1) && !grepl("[{]barcode[}]", ab1_name_template)) {
    stop("`ab1_name_template` must contain `{barcode}`.", call. = FALSE)
  }

  result_dir <- normalizePath(result_dir, mustWork = TRUE)
  if (!is.null(sample_map)) {
    check_file_arg(sample_map, "sample_map")
    sample_map <- normalizePath(sample_map, mustWork = TRUE)
  }
  if (!is.null(fastq_pass_trim_dir)) {
    check_dir_arg(fastq_pass_trim_dir, "fastq_pass_trim_dir")
    fastq_pass_trim_dir <- normalizePath(fastq_pass_trim_dir, mustWork = TRUE)
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
    barcode_name <- basename(barcode_dir)
    items <- rbind(
      items,
      collect_item(
        label = barcode_name,
        source = barcode_dir,
        destination = file.path(output_dir, barcode_name),
        type = "directory",
        overwrite = overwrite
      )
    )

    if (!is.null(fastq_pass_trim_dir)) {
      plot_files <- find_sample_length_plots(
        fastq_pass_trim_dir = fastq_pass_trim_dir,
        barcode_name = barcode_name,
        pattern = sample_length_plot_pattern
      )
      for (plot_file in plot_files) {
        items <- rbind(
          items,
          collect_item(
            label = paste0(barcode_name, "_read_length_plot"),
            source = plot_file,
            destination = file.path(output_dir, barcode_name, basename(plot_file)),
            type = "file",
            overwrite = overwrite
          )
        )
      }
    }
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

  readme_path <- file.path(output_dir, readme_name)
  writeLines(
    amplicon_results_readme(
      trimmed_consensus_file = trimmed_consensus_file,
      consensus_file = consensus_file,
      consensus_index_file = consensus_index_file,
      sample_map_name = if (is.null(sample_map)) NULL else basename(sample_map),
      barcode_names = basename(barcode_dirs),
      include_sample_length_plots = !is.null(fastq_pass_trim_dir),
      include_execution = include_execution,
      include_ab1 = include_ab1,
      ab1_name_template = ab1_name_template
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

  if (!is.null(chinese_readme_name)) {
    chinese_readme_path <- file.path(output_dir, chinese_readme_name)
    writeLines(
      amplicon_results_readme_zh(
        trimmed_consensus_file = trimmed_consensus_file,
        consensus_file = consensus_file,
        consensus_index_file = consensus_index_file,
        sample_map_name = if (is.null(sample_map)) NULL else basename(sample_map),
        barcode_names = basename(barcode_dirs),
        include_sample_length_plots = !is.null(fastq_pass_trim_dir),
        include_execution = include_execution,
        include_ab1 = include_ab1,
        ab1_name_template = ab1_name_template
      ),
      chinese_readme_path
    )

    items <- rbind(
      items,
      data.frame(
        label = "README_zh_CN",
        source = NA_character_,
        destination = chinese_readme_path,
        type = "file",
        exists = TRUE,
        copied = TRUE,
        stringsAsFactors = FALSE
      )
    )
  }

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
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
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

find_sample_length_plots <- function(fastq_pass_trim_dir, barcode_name, pattern) {
  barcode_dir <- file.path(fastq_pass_trim_dir, barcode_name)
  if (!dir.exists(barcode_dir)) {
    return(character())
  }

  list.files(
    barcode_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE
  )
}

amplicon_results_readme <- function(trimmed_consensus_file,
                                    consensus_file,
                                    consensus_index_file,
                                    sample_map_name,
                                    barcode_names,
                                    include_sample_length_plots,
                                    include_execution,
                                    include_ab1,
                                    ab1_name_template) {
  lines <- c(
    "Amplicon sequencing final results",
    "=================================",
    "",
    "Recommended file",
    "----------------",
    paste0(
      "- ",
      trimmed_consensus_file,
      ": trimmed consensus FASTA. ",
      "This file is generated only when the primer sequences are provided. ",
      "Use this file (if exists) or ", consensus_file, " for downstream analyses; ",
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
      if (isTRUE(include_ab1)) {
        paste0(
          "- barcode*/",
          sub("[{]barcode[}]", "<barcode>", ab1_name_template),
          ": synthetic Sanger-style AB1 chromatogram generated from ",
          consensus_file,
          " and the per-barcode BAM alignment."
        )
      },
      if (isTRUE(include_sample_length_plots)) {
        "- barcode*/Distribution_seqLength__*.png: per-barcode read length distribution plot generated from filtered reads."
      },
      "",
      "Included barcode folders:",
      paste0("- ", barcode_names)
    )
  } else {
    lines <- c(lines, "- No barcode folders were found in the source result directory.")
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

amplicon_results_readme_zh <- function(trimmed_consensus_file,
                                       consensus_file,
                                       consensus_index_file,
                                       sample_map_name,
                                       barcode_names,
                                       include_sample_length_plots,
                                       include_execution,
                                       include_ab1,
                                       ab1_name_template) {
  lines <- c(
    "扩增子测序最终结果",
    "==================",
    "",
    "推荐使用的文件",
    "--------------",
    paste0(
      "- ",
      trimmed_consensus_file,
      "：去除引物外侧序列后的共识序列 FASTA。只有在提供引物序列时才会生成此文件。下游分析优先使用该文件；如果该文件不存在，可使用 ",
      consensus_file,
      "。修剪过程中会去除扩增子引物边界外的额外碱基，并保留引物序列。"
    ),
    "",
    "其他文件",
    "--------",
    paste0(
      "- ",
      consensus_file,
      "：扩增子分析流程生成的原始未修剪共识序列 FASTA。"
    ),
    paste0(
      "- ",
      consensus_index_file,
      "：",
      consensus_file,
      " 的 FASTA 索引文件。"
    )
  )

  if (!is.null(sample_map_name)) {
    lines <- c(
      lines,
      paste0(
        "- ",
        sample_map_name,
        "：barcode 与样本名称的对应表，可用于把 barcode 文件夹对应到具体样本。"
      )
    )
  }

  lines <- c(
    lines,
    "",
    "Barcode 文件夹",
    "--------------"
  )

  if (length(barcode_names) > 0L) {
    lines <- c(
      lines,
      "- barcode*/：从流程结果目录复制得到的各 barcode 结果文件夹。",
      "- barcode*/consensus/consensus.fastq：每个 barcode 的共识序列 FASTQ。",
      "- barcode*/alignments/*.bam：该 barcode 的 reads 比对结果。",
      "- barcode*/alignments/*.bam.bai：BAM 索引文件。",
      "- barcode*/alignments/*.png：已生成的比对结果截图。",
      if (isTRUE(include_ab1)) {
        paste0(
          "- barcode*/",
          sub("[{]barcode[}]", "<barcode>", ab1_name_template),
          "：由 ",
          consensus_file,
          " 和该 barcode 的 BAM 比对结果生成的模拟 Sanger AB1 峰图文件。"
        )
      },
      if (isTRUE(include_sample_length_plots)) {
        "- barcode*/Distribution_seqLength__*.png：基于过滤后 reads 生成的每个 barcode 的读长分布图。"
      },
      "",
      "包含的 barcode 文件夹：",
      paste0("- ", barcode_names)
    )
  } else {
    lines <- c(lines, "- 在源结果目录中没有找到 barcode 文件夹。")
  }

  if (isTRUE(include_execution)) {
    lines <- c(
      lines,
      "",
      "运行信息",
      "--------",
      "- execution/：如果源结果中存在，则包含 Nextflow 的运行报告、时间线和 trace 文件。"
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
