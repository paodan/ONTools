#' Build amplicon consensus delivery packages
#'
#' `make_consensus_delivery()` runs the common ONT amplicon consensus delivery
#' workflow: optional Dorado basecalling and demultiplexing, run-level QC,
#' grouping demultiplexed barcode FASTQ folders by sample information, running
#' `wf-amplicon`, optional primer-boundary trimming, filtered-read QC, optional
#' IGV snapshots, and final result collection.
#'
#' @param path_proj MinKNOW project directory.
#' @param path_sampleInfo_file_list Named character vector or list of sample
#'   information CSV files. Names are used as grouped FASTQ folder names.
#' @param path_delivery Root directory where per-group final results are copied.
#' @param kit_name Dorado barcode kit name.
#' @param model Dorado basecalling model or model alias.
#' @param demux_out Demux output directory name under `path_proj`. If `NULL`,
#'   defaults to `demux_out_<kit_name>`.
#' @param fastq_out FASTQ output directory name created by Dorado conversion.
#' @param barcode_both_ends Logical. Passed to [run_dorado_demux_to_fastq()].
#' @param run_basecalling_demux_step Logical. If `FALSE`, reuse
#'   `work/stat.Rdata` when present or reconstruct output paths from
#'   `path_proj`, `demux_out`, and `fastq_out`.
#' @param run_QC_step Logical. Generate run-level QC plots and source FASTQ
#'   statistics.
#' @param move_fastq_step Logical. Move barcode FASTQ folders into sample-info
#'   groups before running amplicon analysis.
#' @param run_amplicon_step Logical. Run [run_wf_amplicon()] for each group.
#' @param trim_consensus_step Logical. Trim consensus FASTA files when primer
#'   columns are present and non-empty.
#' @param run_filtered_QC_step Logical. Generate QC plots for filtered reads for
#'   each group.
#' @param run_igv_step Logical. Generate IGV snapshots for per-barcode BAM files.
#' @param collect_results_step Logical. Copy final deliverables with
#'   [collect_amplicon_results()].
#' @param project_col,amplicon_size_col,barcode_col Column names in each sample
#'   information file.
#' @param min_read_length_col,max_read_length_col Primer-independent read-length
#'   filter column names in each sample information file.
#' @param f_primer_col,r_primer_col Forward and reverse primer column names.
#' @param min_read_qual,min_n_reads,force_spoa_length_threshold,override_basecaller_cfg,profile,resume
#'   Parameters passed to [run_wf_amplicon()].
#' @param barcode_digits Number of barcode digits used by [plot_seqQC()] for
#'   filtered-read QC.
#' @param igv IGV executable path passed to [igv_snapshot()].
#' @param overwrite_fastq Logical. Passed to [move_fastq_to_folders()].
#' @param overwrite_delivery Logical. Passed to [collect_amplicon_results()].
#' @param include_execution Logical. Passed to [collect_amplicon_results()].
#' @param dry_run Logical. If `TRUE`, build plans and commands without running
#'   external analysis steps or copying final deliverables.
#' @param echo Logical. If `TRUE`, print commands from wrapped command runners.
#' @param wait Logical. Passed to wrapped command runners.
#' @param stdout,stderr Passed to wrapped command runners.
#'
#' @return Invisibly returns a list with run-level QC (`g1`), per-group QC
#'   (`g2`), Dorado/basecalling status (`stat`), FASTQ move plans, amplicon
#'   workflow results, trim summaries, IGV results, delivery copy summaries, and
#'   important output paths.
#'
#' @examples
#' proj <- tempfile("ont-project-")
#' dir.create(file.path(proj, "pod5"), recursive = TRUE)
#' sample_info <- tempfile(fileext = ".csv")
#' utils::write.csv(
#'   data.frame(
#'     Barcode_ID = "PBC001-001",
#'     Project_ID = "PROJECT001",
#'     Expected_Size_bp = "1600",
#'     Min_Read_Length = 1200,
#'     Max_Read_Length = 1800
#'   ),
#'   sample_info,
#'   row.names = FALSE
#' )
#'
#' res <- make_consensus_delivery(
#'   path_proj = proj,
#'   path_sampleInfo_file_list = c(PROJECT001_1600 = sample_info),
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' res$workflow$PROJECT001_1600$command_string
#'
#' @export
make_consensus_delivery <- function(path_proj,
                                    path_sampleInfo_file_list,
                                    path_delivery = "/data/project_delivery",
                                    kit_name = "YS-NB576",
                                    model = "sup",
                                    demux_out = NULL,
                                    fastq_out = "fastq_pass_trim",
                                    barcode_both_ends = FALSE,
                                    run_basecalling_demux_step = TRUE,
                                    run_QC_step = TRUE,
                                    move_fastq_step = TRUE,
                                    run_amplicon_step = TRUE,
                                    trim_consensus_step = TRUE,
                                    run_filtered_QC_step = TRUE,
                                    run_igv_step = TRUE,
                                    collect_results_step = TRUE,
                                    project_col = "Project_ID",
                                    amplicon_size_col = "Expected_Size_bp",
                                    barcode_col = "Barcode_ID",
                                    min_read_length_col = "Min_Read_Length",
                                    max_read_length_col = "Max_Read_Length",
                                    f_primer_col = "Primer_F",
                                    r_primer_col = "Primer_R",
                                    min_read_qual = 10,
                                    min_n_reads = 40,
                                    force_spoa_length_threshold = 2000,
                                    override_basecaller_cfg = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
                                    profile = "standard",
                                    resume = TRUE,
                                    barcode_digits = 3,
                                    igv = "/usr/local/bin/IGV_Linux_2.19.8/igv.sh",
                                    overwrite_fastq = FALSE,
                                    overwrite_delivery = TRUE,
                                    include_execution = FALSE,
                                    dry_run = FALSE,
                                    echo = TRUE,
                                    wait = TRUE,
                                    stdout = "",
                                    stderr = "") {
  check_dir_arg(path_proj, "path_proj")
  check_scalar_character(path_delivery, "path_delivery")
  check_scalar_character(kit_name, "kit_name")
  check_scalar_character(model, "model")
  check_scalar_character(fastq_out, "fastq_out")
  check_scalar_character(project_col, "project_col")
  check_scalar_character(amplicon_size_col, "amplicon_size_col")
  check_scalar_character(barcode_col, "barcode_col")
  check_scalar_character(min_read_length_col, "min_read_length_col")
  check_scalar_character(max_read_length_col, "max_read_length_col")
  check_scalar_character(f_primer_col, "f_primer_col")
  check_scalar_character(r_primer_col, "r_primer_col")
  check_scalar_character(override_basecaller_cfg, "override_basecaller_cfg")
  check_scalar_character(profile, "profile")
  check_scalar_character(igv, "igv")
  check_logical_scalar(barcode_both_ends, "barcode_both_ends")
  check_logical_scalar(run_basecalling_demux_step, "run_basecalling_demux_step")
  check_logical_scalar(run_QC_step, "run_QC_step")
  check_logical_scalar(move_fastq_step, "move_fastq_step")
  check_logical_scalar(run_amplicon_step, "run_amplicon_step")
  check_logical_scalar(trim_consensus_step, "trim_consensus_step")
  check_logical_scalar(run_filtered_QC_step, "run_filtered_QC_step")
  check_logical_scalar(run_igv_step, "run_igv_step")
  check_logical_scalar(collect_results_step, "collect_results_step")
  check_logical_scalar(resume, "resume")
  check_logical_scalar(overwrite_fastq, "overwrite_fastq")
  check_logical_scalar(overwrite_delivery, "overwrite_delivery")
  check_logical_scalar(include_execution, "include_execution")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")

  min_n_reads <- validate_positive_integer(min_n_reads, "min_n_reads")
  force_spoa_length_threshold <- validate_positive_integer(
    force_spoa_length_threshold,
    "force_spoa_length_threshold"
  )
  barcode_digits <- validate_positive_integer(barcode_digits, "barcode_digits")
  min_read_qual <- validate_nonnegative_number(min_read_qual, "min_read_qual")

  if (!is.null(demux_out)) {
    check_scalar_character(demux_out, "demux_out")
  }

  path_proj <- normalizePath(path_proj, mustWork = TRUE)
  path_sampleInfo_file_list <- validate_sample_info_file_list(
    path_sampleInfo_file_list
  )
  demux_out_name <- if (is.null(demux_out)) paste0("demux_out_", kit_name) else demux_out

  path_work <- file.path(path_proj, "work")
  dir.create(path_work, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_delivery, showWarnings = FALSE, recursive = TRUE)
  if (dir.exists(path_delivery)) {
    path_delivery <- normalizePath(path_delivery, mustWork = TRUE)
  } else if (isTRUE(dry_run)) {
    path_delivery <- normalizePath(path_delivery, mustWork = FALSE)
  } else {
    stop("Could not create `path_delivery`: ", path_delivery, call. = FALSE)
  }

  if (isTRUE(run_basecalling_demux_step)) {
    message("Step 1: basecalling and demultiplexing")
    stat <- run_dorado_demux_to_fastq(
      proj = path_proj,
      kit_name = kit_name,
      model = model,
      demux_out = demux_out_name,
      fastq_out = fastq_out,
      barcode_both_ends = barcode_both_ends,
      dry_run = dry_run,
      echo = echo,
      wait = wait,
      stdout = stdout,
      stderr = stderr
    )
    if (!isTRUE(dry_run) && isTRUE(wait)) {
      save(stat, file = file.path(path_work, "stat.Rdata"))
    }
  } else {
    message("Step 1: basecalling and demultiplexing (skipped)")
    stat_file <- file.path(path_work, "stat.Rdata")
    if (file.exists(stat_file)) {
      load(stat_file)
    } else {
      stat <- list(
        command = NA_character_,
        args = character(),
        status = NA_integer_,
        paths = dorado_demux_to_fastq_paths(
          proj = path_proj,
          model = model,
          demux_out = demux_out_name,
          fastq_out = fastq_out,
          scan_dynamic = TRUE
        )
      )
    }
  }

  fastq_root <- select_single_path(
    stat$paths$fastq_dirs,
    "FASTQ output directory",
    fallback = if (isTRUE(dry_run)) file.path(stat$paths$demux_dir, "<run_dir>", fastq_out) else NA_character_
  )
  run_root <- select_single_path(
    stat$paths$run_dirs,
    "Dorado run directory",
    fallback = if (isTRUE(dry_run)) file.path(stat$paths$demux_dir, "<run_dir>") else NA_character_
  )
  path_seq_summary <- file.path(stat$paths$demux_dir, "sequencing_summary.txt")
  message("Sequencing summary is saved here: ", path_seq_summary)

  fig_path <- file.path(stat$paths$demux_dir, "figs", "png")
  stat_path <- file.path(stat$paths$demux_dir, "stats")
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(stat_path, showWarnings = FALSE, recursive = TRUE)

  g1 <- NULL
  source_fastq_stats <- NULL
  if (isTRUE(run_QC_step)) {
    if (!file.exists(path_seq_summary)) {
      if (isTRUE(dry_run)) {
        message("Step 2: Making QC plots for all samples (dry-run; sequencing summary not found)")
      } else {
        stop("Sequencing summary file not found: ", path_seq_summary, call. = FALSE)
      }
    } else if (is.na(fastq_root)) {
      stop("No FASTQ output directory was found for run-level QC.", call. = FALSE)
    } else {
    message("Step 2: Making QC plots for all samples")
    g1 <- plot_seqQC(
      filePath = path_seq_summary,
      runName = model,
      device = "png",
      barcodes = NULL,
      sample_len_dir_mode = "fastq_pass_trim",
      fastq_pass_trim_dir = fastq_root
    )
    message("Demultiplexing rate is ", g1$recovery)

    ggplot2::ggsave(
      filename = file.path(fig_path, "numRead_perSample.png"),
      plot = g1$numRead,
      width = 16,
      height = 8
    )
    ggplot2::ggsave(
      filename = file.path(fig_path, "lenReadDistribution.png"),
      plot = g1$lenRead,
      width = 22,
      height = 22
    )
    utils::write.csv(g1$read_counts, file = file.path(stat_path, "numRead_perSample.csv"))

    source_fastqs <- list.files(
      file.path(fastq_root, list.files(fastq_root)),
      pattern = "[.]fastq[.]gz$",
      full.names = TRUE
    )
    source_fastq_stats <- fastq_stats(
      fastq = source_fastqs,
      stderr = FALSE,
      echo = FALSE,
      conda_env = NULL
    )$stats
    utils::write.csv(source_fastq_stats, file = file.path(stat_path, "stat_source_fastq.csv"))
    }
  } else {
    message("Step 2: Making QC plots for all samples (skipped)")
  }

  move_plans <- list()
  workflow <- list()
  trim <- list()
  g2 <- list()
  igv_results <- list()
  delivery <- list()

  for (folder in names(path_sampleInfo_file_list)) {
    message(folder)
    sample_info_file <- path_sampleInfo_file_list[[folder]]
    sample_info <- read_sample_info_table(sample_info_file)
    sample_barcode_col <- resolve_barcode_col(sample_info, barcode_col)
    validate_amplicon_delivery_sample_info(
      sample_info = sample_info,
      required_cols = c(sample_barcode_col, min_read_length_col, max_read_length_col)
    )
    sample_info_sub <- sample_info[1, , drop = FALSE]
    min_read_length <- validate_positive_integer(
      sample_info_sub[[min_read_length_col]],
      min_read_length_col
    )
    max_read_length <- validate_positive_integer(
      sample_info_sub[[max_read_length_col]],
      max_read_length_col
    )
    if (min_read_length > max_read_length) {
      stop("`", min_read_length_col, "` must be less than or equal to `",
           max_read_length_col, "` for group `", folder, "`.", call. = FALSE)
    }

    if (isTRUE(move_fastq_step)) {
      if (is.na(fastq_root)) {
        stop("No FASTQ output directory was found for moving barcode folders.",
             call. = FALSE)
      }
      message("Step 3: Move fastq files to a folder")
      if (isTRUE(dry_run) && !dir.exists(fastq_root)) {
        move_plans[[folder]] <- make_fastq_move_plan_without_source_dir(
          fastq_dir = fastq_root,
          sample_info = sample_info,
          project_col = project_col,
          amplicon_size_col = amplicon_size_col,
          barcode_col = sample_barcode_col
        )
      } else {
        move_plans[[folder]] <- move_fastq_to_folders(
          fastq_dir = fastq_root,
          sample_info = sample_info_file,
          project_col = project_col,
          amplicon_size_col = amplicon_size_col,
          barcode_col = barcode_col,
          overwrite = overwrite_fastq,
          dry_run = dry_run
        )
      }
    } else {
      message("Step 3: Move fastq files to a folder (skipped)")
      move_plans[[folder]] <- NULL
    }

    group_fastq <- if (is.na(fastq_root)) NA_character_ else file.path(fastq_root, folder)
    group_out_dir <- if (is.na(run_root)) {
      file.path(stat$paths$demux_dir, "results", "wf_amplicon_denovo", folder)
    } else {
      file.path(run_root, "results", "wf_amplicon_denovo", folder)
    }

    if (isTRUE(run_amplicon_step)) {
      message("Step 4: Run amplicon workflow")
      workflow[[folder]] <- run_wf_amplicon(
        fastq = group_fastq,
        out_dir = group_out_dir,
        min_read_length = min_read_length,
        max_read_length = max_read_length,
        min_read_qual = min_read_qual,
        min_n_reads = min_n_reads,
        force_spoa_length_threshold = force_spoa_length_threshold,
        override_basecaller_cfg = override_basecaller_cfg,
        profile = profile,
        resume = resume,
        dry_run = dry_run,
        echo = echo,
        wait = wait,
        stdout = stdout,
        stderr = stderr
      )
    } else {
      message("Step 4: Run amplicon workflow (skipped)")
      workflow[[folder]] <- list(
        status = NA_integer_,
        paths = list(fastq = group_fastq, out_dir = group_out_dir)
      )
    }

    if (isTRUE(trim_consensus_step)) {
      message("Step 5: Trim extra bases")
      if (has_primer_pair(sample_info_sub, f_primer_col, r_primer_col)) {
        input_fasta <- file.path(workflow[[folder]]$paths$out_dir, "all-consensus-seqs.fasta")
        output_fasta <- file.path(workflow[[folder]]$paths$out_dir, "all-consensus-seqs_trimmed.fasta")
        if (isTRUE(dry_run)) {
          trim[[folder]] <- list(
            input_fasta = input_fasta,
            output_fasta = output_fasta,
            f_primer = sample_info_sub[[f_primer_col]],
            r_primer = sample_info_sub[[r_primer_col]],
            status = "dry_run"
          )
        } else {
          trim[[folder]] <- trim_fasta_keep_primers(
            input_fasta = input_fasta,
            f_primer = sample_info_sub[[f_primer_col]],
            r_primer = sample_info_sub[[r_primer_col]],
            output_fasta = output_fasta
          )
        }
      } else {
        trim[[folder]] <- NULL
      }
    } else {
      message("Step 5: Trim extra bases (skipped)")
      trim[[folder]] <- NULL
    }

    if (isTRUE(run_filtered_QC_step)) {
      if (!file.exists(path_seq_summary)) {
        if (isTRUE(dry_run)) {
          message("Step 6: Run QC plots for filtered reads (dry-run; sequencing summary not found)")
          g2[[folder]] <- NULL
        } else {
          stop("Sequencing summary file not found: ", path_seq_summary, call. = FALSE)
        }
      } else {
      message("Step 6: Run QC plots for filtered reads")
      barcode_numbers <- barcode_ids_to_numbers(sample_info[[sample_barcode_col]])
      g2[[folder]] <- plot_seqQC(
        filePath = path_seq_summary,
        runName = model,
        device = "png",
        barcodes = barcode_numbers,
        barcode_digits = barcode_digits,
        min_read_length = min_read_length,
        max_read_length = max_read_length,
        filename_suffix = "_filtered",
        sample_len_dir_mode = "fastq_pass_trim",
        fastq_pass_trim_dir = workflow[[folder]]$paths$fastq,
        plot_unclassified = FALSE
      )

      ggplot2::ggsave(
        filename = file.path(workflow[[folder]]$paths$out_dir, "numRead_perSample.png"),
        plot = g2[[folder]]$numRead,
        width = 12,
        height = 8
      )
      ggplot2::ggsave(
        filename = file.path(workflow[[folder]]$paths$out_dir, "lenReadDistribution.png"),
        plot = g2[[folder]]$lenRead,
        width = 15,
        height = 15
      )
      utils::write.csv(g2[[folder]]$read_counts, file = file.path(workflow[[folder]]$paths$out_dir, "numRead_perSample.csv"))
      }
    } else {
      message("Step 6: Run QC plots for filtered reads (skipped)")
      g2[[folder]] <- NULL
    }

    if (isTRUE(run_igv_step)) {
      message("Step 7: Generate IGV plots")
      if (isTRUE(dry_run)) {
        igv_results[[folder]] <- list(status = "dry_run")
      } else {
        barcode_dirs <- list.files(
          workflow[[folder]]$paths$out_dir,
          pattern = "^barcode",
          full.names = FALSE
        )
        igv_results[[folder]] <- lapply(barcode_dirs, function(barcode_dir) {
          igv_snapshot(
            genome_fasta = file.path(workflow[[folder]]$paths$out_dir, "all-consensus-seqs.fasta"),
            bam = file.path(
              workflow[[folder]]$paths$out_dir,
              barcode_dir,
              "alignments",
              paste0(barcode_dir, ".aligned.sorted.bam")
            ),
            out_dir = file.path(workflow[[folder]]$paths$out_dir, barcode_dir, "alignments"),
            chr = barcode_dir,
            igv = igv
          )
        })
        names(igv_results[[folder]]) <- barcode_dirs
      }
    } else {
      message("Step 7: Generate IGV plots (skipped)")
      igv_results[[folder]] <- NULL
    }

    if (isTRUE(collect_results_step)) {
      message("Step 8: Deliver the results")
      if (isTRUE(dry_run)) {
        delivery[[folder]] <- data.frame(
          label = "delivery",
          source = workflow[[folder]]$paths$out_dir,
          destination = file.path(path_delivery, folder),
          type = "directory",
          exists = dir.exists(workflow[[folder]]$paths$out_dir),
          copied = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        delivery[[folder]] <- collect_amplicon_results(
          result_dir = workflow[[folder]]$paths$out_dir,
          output_dir = file.path(path_delivery, folder),
          sample_map = sample_info_file,
          fastq_pass_trim_dir = workflow[[folder]]$paths$fastq,
          include_execution = include_execution,
          overwrite = overwrite_delivery
        )
        utils::write.csv(
          delivery[[folder]],
          file = file.path(path_work, paste0("delivery_stats_", folder, ".csv")),
          row.names = FALSE
        )
      }
    } else {
      message("Step 8: Deliver the results (skipped)")
      delivery[[folder]] <- NULL
    }
  }

  invisible(list(
    g1 = g1,
    g2 = g2,
    stat = stat,
    move_plans = move_plans,
    workflow = workflow,
    trim = trim,
    igv = igv_results,
    delivery = delivery,
    source_fastq_stats = source_fastq_stats,
    path_work = path_work,
    path_delivery = path_delivery,
    path_seq_summary = path_seq_summary,
    fig_path = fig_path,
    stat_path = stat_path
  ))
}

validate_sample_info_file_list <- function(path_sampleInfo_file_list) {
  if (is.character(path_sampleInfo_file_list)) {
    files <- path_sampleInfo_file_list
  } else if (is.list(path_sampleInfo_file_list) &&
             all(vapply(path_sampleInfo_file_list, function(x) {
               is.character(x) && length(x) == 1L
             }, logical(1)))) {
    files <- unlist(path_sampleInfo_file_list, use.names = TRUE)
  } else {
    stop("`path_sampleInfo_file_list` must be a character vector or list of file paths.",
         call. = FALSE)
  }

  if (length(files) == 0L || anyNA(files) || any(!nzchar(files))) {
    stop("`path_sampleInfo_file_list` must contain at least one non-empty file path.",
         call. = FALSE)
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0L) {
    stop("Sample information file(s) not found: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  if (is.null(names(files)) || anyNA(names(files)) || any(!nzchar(names(files)))) {
    stop("`path_sampleInfo_file_list` must be named; names are used as FASTQ group folders.",
         call. = FALSE)
  }

  stats::setNames(normalizePath(files, mustWork = TRUE), names(files))
}

validate_amplicon_delivery_sample_info <- function(sample_info, required_cols) {
  if (nrow(sample_info) == 0L) {
    stop("`sample_info` must contain at least one row.", call. = FALSE)
  }
  missing_cols <- setdiff(required_cols, names(sample_info))
  if (length(missing_cols) > 0L) {
    stop(
      "`sample_info` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

select_single_path <- function(paths, label, fallback = NA_character_) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0L) {
    return(fallback)
  }
  if (length(paths) > 1L) {
    warning(
      "Multiple ", label, " values were found; using the first: ", paths[[1]],
      call. = FALSE
    )
  }

  paths[[1]]
}

has_primer_pair <- function(sample_info_row, f_primer_col, r_primer_col) {
  if (!all(c(f_primer_col, r_primer_col) %in% names(sample_info_row))) {
    return(FALSE)
  }

  f_primer <- sample_info_row[[f_primer_col]]
  r_primer <- sample_info_row[[r_primer_col]]
  !is.na(f_primer) && !is.na(r_primer) && nzchar(f_primer) && nzchar(r_primer)
}

barcode_ids_to_numbers <- function(barcode_ids) {
  barcode_ids <- as.character(barcode_ids)
  out <- sub("^.*?([0-9]+)$", "\\1", barcode_ids)
  out[!grepl("[0-9]+$", barcode_ids)] <- NA_character_
  out <- as.integer(out)
  if (anyNA(out)) {
    stop(
      "`Barcode_ID` values must end with barcode numbers for filtered QC: ",
      paste(barcode_ids[is.na(out)], collapse = ", "),
      call. = FALSE
    )
  }

  out
}

make_fastq_move_plan_without_source_dir <- function(fastq_dir,
                                                    sample_info,
                                                    project_col,
                                                    amplicon_size_col,
                                                    barcode_col) {
  barcode_col <- resolve_barcode_col(sample_info, barcode_col)
  required_cols <- c(barcode_col, amplicon_size_col)
  if (!is.null(project_col)) {
    required_cols <- c(required_cols, project_col)
  }
  validate_amplicon_delivery_sample_info(sample_info, required_cols)

  if (barcode_col == "Barcode_ID") {
    barcodes <- paste0("barcode", sub(".+-([0-9]+$)", "\\1", sample_info[[barcode_col]]))
  } else {
    barcodes <- as.character(sample_info[[barcode_col]])
  }
  amplicon_sizes <- as.character(sample_info[[amplicon_size_col]])
  projects <- if (is.null(project_col)) {
    rep(NA_character_, nrow(sample_info))
  } else {
    as.character(sample_info[[project_col]])
  }
  groups <- if (is.null(project_col)) {
    safe_path_component(amplicon_sizes)
  } else {
    paste(safe_path_component(projects), safe_path_component(amplicon_sizes), sep = "_")
  }

  data.frame(
    barcode = barcodes,
    project = projects,
    amplicon_size = amplicon_sizes,
    group = groups,
    source = file.path(fastq_dir, barcodes),
    destination = file.path(fastq_dir, groups, barcodes),
    source_exists = FALSE,
    destination_exists = FALSE,
    moved = FALSE,
    status = "dry_run",
    stringsAsFactors = FALSE
  )
}
