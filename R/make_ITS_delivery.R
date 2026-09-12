#' Build an integrated ITS delivery folder
#'
#' `make_ITS_delivery()` combines amplicon consensus delivery results, wf-16s
#' ITS taxonomic profiling results, and optional consensus-vs-UNITE BLAST
#' review results into one customer-facing delivery folder.
#'
#' The default layout keeps per-sample files under `samples/barcode*/` and
#' keeps project-level summaries, figures, and UNITE top-hit results directly
#' under `path_delivery`. This makes the root directory useful for quick review
#' while keeping each barcode's detailed evidence in one place.
#'
#' @param path_ITS_result Path to a completed wf-16s result directory for ITS
#'   data. If `NULL` or if the directory does not exist, [run_ITS()] is run
#'   first. When an existing directory is supplied, the function expects it to
#'   contain `abundance_table_genus.tsv` and `alignment_tables/`.
#' @param path_delivery Final ITS delivery root. In grouped mode, per-group ITS
#'   deliveries are created below this directory.
#' @param consensus_delivery_path Existing consensus delivery directory. If
#'   `NULL` or if the directory does not exist, [make_consensus_delivery()] is
#'   run first using `path_proj`, `path_sampleInfo_file_list`, and `...`.
#' @param path_proj,path_sampleInfo_file_list Arguments passed to
#'   [make_consensus_delivery()] when `consensus_delivery_path` is missing.
#'   When `path_sampleInfo_file_list` is supplied, its names are also used as
#'   delivery group names, matching [make_consensus_delivery()].
#' @param consensus_delivery_output Directory used to store newly generated
#'   consensus delivery results. The default `NULL` creates a sibling work
#'   directory named `basename(path_delivery)_consensus_work`.
#' @param ITS_fastq,ITS_out_dir,ITS_work_dir,ITS_profile,ITS_resume,ITS_database_set,ITS_min_len,ITS_max_len,ITS_workflow,ITS_nextflow,ITS_quiet,ITS_extra_args,ITS_syntax_parser,ITS_ansi_log,ITS_nextflow_env
#'   Parameters passed to [run_ITS()] when `path_ITS_result` is `NULL` or
#'   missing. Defaults mirror [run_ITS()] except that `ITS_out_dir` and
#'   `ITS_work_dir` default to sibling work directories derived from
#'   `path_delivery`.
#' @param group_delivery_dir Directory name created under each sample-info
#'   group when `path_sampleInfo_file_list` is supplied. Default is
#'   `"ITS_results"`, producing `path_delivery/<group>/ITS_results/`.
#' @param barcode_col Column name in each sample information file used to infer
#'   barcode directory names. Default is `"Barcode_ID"`, matching
#'   [make_consensus_delivery()]. Values like `PBC001-097` are converted to
#'   `barcode097`; otherwise values are used as-is after resolving the barcode
#'   column.
#' @param samples_dir Directory under `path_delivery` that stores per-barcode
#'   sample results. The default `"samples"` creates
#'   `path_delivery/samples/barcode*/`. Set to `NULL` to put `barcode*/`
#'   directly under `path_delivery`.
#' @param abundance_table,alignment_tables_dir,figure_dir,identification_dir
#'   File/directory names used by [move_ITS()] and the final delivery layout.
#'   Defaults are `"abundance_table_genus.tsv"`, `"alignment_tables"`,
#'   `"figures"`, and `"identification_tables"`. In the final delivery,
#'   `abundance_table` and `figure_dir` are written at the root, while each
#'   barcode alignment table is copied to
#'   `samples/barcode*/identification_tables/`.
#' @param tax_levels,cutoff,width,height Plotting parameters passed to
#'   [move_ITS()]. By default, plots are generated for Kingdom, Phylum, Class,
#'   Order, Family, and Genus; taxa below `cutoff = 0.01` relative abundance
#'   are grouped by the plotting helper; figures are saved at 12 x 6 inches.
#' @param consensus_file,trimmed_consensus_file Consensus FASTA names searched
#'   under the consensus delivery. The default prefers
#'   `"all-consensus-seqs_trimmed.fasta"` for UNITE review and falls back to
#'   `"all-consensus-seqs.fasta"` when trimmed consensus files are absent.
#' @param barcode_pattern Regular expression used to identify barcode
#'   directories and alignment-stat files. The default `"^barcode[0-9]+$"`
#'   matches names such as `barcode097` and `barcode303`.
#' @param run_unite_annotation Logical. If `TRUE`, run
#'   [annotate_consensus_blast()] on the collected consensus FASTA. Default is
#'   `TRUE`; if neither `unite_db` nor `unite_db_fasta` is supplied, this step
#'   is skipped with a warning.
#' @param unite_db,unite_db_fasta BLAST database prefix or FASTA passed to
#'   [annotate_consensus_blast()]. When both are `NULL`, UNITE annotation is
#'   skipped with a warning. Use `unite_db` for an already-built BLAST database
#'   prefix, or `unite_db_fasta` to build a database from a FASTA file using
#'   `makeblastdb`.
#' @param unite_dir Directory under `path_delivery` for detailed UNITE BLAST
#'   results. Default is `"unite_consensus_annotation"`.
#' @param unite_top_hits_name Filename copied to `path_delivery` for the top-hit
#'   summary. Default is `"unite_consensus_top_hits.tsv"`, placed beside
#'   `abundance_table_genus.tsv` for quick review.
#' @param unite_threads,unite_max_target_seqs,unite_evalue,unite_task,unite_word_size,unite_strand,unite_dust,unite_perc_identity,unite_extra_args,unite_species_identity,unite_genus_identity,unite_family_identity,unite_min_query_coverage,unite_min_reference_coverage,unite_novel_identity,unite_blastn,unite_makeblastdb,unite_conda_env,conda
#'   Parameters passed to [annotate_consensus_blast()] for UNITE annotation.
#'   Defaults use `blastn`, `makeblastdb`, `threads = 10`,
#'   `max_target_seqs = 20`, `evalue = "1e-20"`, species/genus/family identity
#'   thresholds of 98.5/95/90 percent, minimum query/reference coverage of
#'   80/50 percent, and `novel_identity = 97`. `unite_conda_env = NULL` means
#'   BLAST tools are called from the current environment.
#' @param readme_name,chinese_readme_name README filenames written under
#'   `path_delivery`. Defaults are `"README.txt"` and `"README.zh-CN.txt"`. Set
#'   `chinese_readme_name = NULL` to skip the Chinese README.
#' @param overwrite Logical. If `TRUE`, replace an existing `path_delivery`.
#'   Default is `FALSE` to avoid accidentally deleting previous deliveries.
#' @param dry_run Logical. If `TRUE`, return a plan without copying files or
#'   running external tools. Default is `FALSE`.
#' @param echo Logical. If `TRUE`, print commands from wrapped runners.
#'   Default is `TRUE`.
#' @param wait,stdout,stderr Passed to wrapped command runners. Defaults are
#'   `wait = TRUE`, `stdout = ""`, and `stderr = ""`.
#' @param ... Additional arguments passed to [make_consensus_delivery()] when a
#'   new consensus delivery must be generated.
#'
#' @return Invisibly returns a list describing generated/copied outputs.
#'
#' @details
#' The function performs four main tasks:
#'
#' 1. Prepare consensus-delivery content. If `consensus_delivery_path` exists,
#'    its barcode folders are copied into the final delivery. Otherwise,
#'    `make_consensus_delivery()` is run first.
#' 2. Prepare wf-16s ITS content. `move_ITS()` is called in a temporary
#'    directory to generate abundance plots and normalize wf-16s outputs. If
#'    `path_ITS_result` was not supplied, `run_ITS()` is called first.
#' 3. Reorganize files. When `path_sampleInfo_file_list` is supplied, one
#'    delivery directory is created for each named sample-info file:
#'    `path_delivery/<group>/ITS_results/`. Each group keeps only the barcode
#'    samples listed in its sample-info table. Without `path_sampleInfo_file_list`,
#'    `path_delivery` itself is treated as the final ITS delivery directory.
#'    Per-barcode consensus outputs are copied to
#'    `samples/barcode*/consensus_results/`, per-barcode wf-16s alignment
#'    tables are copied to `samples/barcode*/identification_tables/`, and
#'    project-level abundance files and figures are copied to `path_delivery`.
#' 4. Optionally run consensus-vs-UNITE review. Consensus FASTA files are
#'    collected into `consensus_for_unite.fasta`, BLAST output is written to
#'    `unite_consensus_annotation/consensus.blast.tsv`, and the top-hit summary
#'    is written both to
#'    `unite_consensus_annotation/consensus.top_hits.tsv` and to the root-level
#'    `unite_consensus_top_hits.tsv`.
#'
#' @examples
#' \dontrun{
#' make_ITS_delivery(
#'   path_delivery = "delivery/ITS",
#'   path_ITS_result = "results/wf_its",
#'   consensus_delivery_path = "delivery/consensus",
#'   unite_db = "/data/reference/UNITE/unite_eukaryotes",
#'   unite_threads = 10
#' )
#'
#' make_ITS_delivery(
#'   path_delivery = "delivery/ITS",
#'   path_ITS_result = NULL,
#'   ITS_fastq = "fastq_pass_trim",
#'   consensus_delivery_path = NULL,
#'   path_proj = "/data/minknow/project/run",
#'   path_sampleInfo_file_list = "SampleInfo.csv",
#'   unite_db_fasta = "/data/reference/UNITE/UNITE_eukaryotes_all.fasta"
#' )
#' }
#'
#' @export
make_ITS_delivery <- function(path_ITS_result = NULL,
                              path_delivery = "/data/project_delivery",
                              consensus_delivery_path = NULL,
                              path_proj = NULL,
                              path_sampleInfo_file_list = NULL,
                              consensus_delivery_output = NULL,
                              ITS_fastq = "./fastq_pass_trim",
                              ITS_out_dir = NULL,
                              ITS_work_dir = NULL,
                              ITS_profile = "standard",
                              ITS_resume = TRUE,
                              ITS_database_set = "ncbi_16s_18s_28s_ITS",
                              ITS_min_len = 300,
                              ITS_max_len = 2000,
                              ITS_workflow = "epi2me-labs/wf-16s",
                              ITS_nextflow = "nextflow",
                              ITS_quiet = FALSE,
                              ITS_extra_args = "--minimap2_by_reference",
                              ITS_syntax_parser = "v1",
                              ITS_ansi_log = FALSE,
                              ITS_nextflow_env = NULL,
                              group_delivery_dir = "ITS_results",
                              barcode_col = "Barcode_ID",
                              samples_dir = "samples",
                              abundance_table = "abundance_table_genus.tsv",
                              alignment_tables_dir = "alignment_tables",
                              figure_dir = "figures",
                              identification_dir = "identification_tables",
                              tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                              cutoff = 0.01,
                              width = 12,
                              height = 6,
                              consensus_file = "all-consensus-seqs.fasta",
                              trimmed_consensus_file = "all-consensus-seqs_trimmed.fasta",
                              barcode_pattern = "^barcode[0-9]+$",
                              run_unite_annotation = TRUE,
                              unite_db = NULL,
                              unite_db_fasta = NULL,
                              unite_dir = "unite_consensus_annotation",
                              unite_top_hits_name = "unite_consensus_top_hits.tsv",
                              unite_threads = 10,
                              unite_max_target_seqs = 20,
                              unite_evalue = "1e-20",
                              unite_task = NULL,
                              unite_word_size = NULL,
                              unite_strand = NULL,
                              unite_dust = NULL,
                              unite_perc_identity = NULL,
                              unite_extra_args = NULL,
                              unite_species_identity = 98.5,
                              unite_genus_identity = 95,
                              unite_family_identity = 90,
                              unite_min_query_coverage = 80,
                              unite_min_reference_coverage = 50,
                              unite_novel_identity = 97,
                              unite_blastn = "blastn",
                              unite_makeblastdb = "makeblastdb",
                              unite_conda_env = NULL,
                              conda = "conda",
                              readme_name = "README.txt",
                              chinese_readme_name = "README.zh-CN.txt",
                              overwrite = FALSE,
                              dry_run = FALSE,
                              echo = TRUE,
                              wait = TRUE,
                              stdout = "",
                              stderr = "",
                              ...) {
  dots <- list(...)
  if ("output_dir" %in% names(dots)) {
    stop(
      "`output_dir` has been removed from `make_ITS_delivery()`; ",
      "use `path_delivery` instead.",
      call. = FALSE
    )
  }
  if (!is.null(path_ITS_result)) {
    check_scalar_character(path_ITS_result, "path_ITS_result")
  }
  check_scalar_character(path_delivery, "path_delivery")
  output_dir <- path_delivery
  if (!is.null(consensus_delivery_path)) {
    check_scalar_character(consensus_delivery_path, "consensus_delivery_path")
  }
  if (!is.null(path_proj)) check_scalar_character(path_proj, "path_proj")
  check_scalar_character(ITS_fastq, "ITS_fastq")
  if (!is.null(ITS_out_dir)) check_scalar_character(ITS_out_dir, "ITS_out_dir")
  if (!is.null(ITS_work_dir)) check_scalar_character(ITS_work_dir, "ITS_work_dir")
  check_scalar_character(ITS_profile, "ITS_profile")
  check_logical_scalar(ITS_resume, "ITS_resume")
  check_scalar_character(ITS_database_set, "ITS_database_set")
  ITS_min_len <- validate_positive_integer(ITS_min_len, "ITS_min_len")
  ITS_max_len <- validate_positive_integer(ITS_max_len, "ITS_max_len")
  check_scalar_character(ITS_workflow, "ITS_workflow")
  check_scalar_character(ITS_nextflow, "ITS_nextflow")
  check_logical_scalar(ITS_quiet, "ITS_quiet")
  if (!is.null(ITS_extra_args)) check_scalar_character(ITS_extra_args, "ITS_extra_args")
  if (!is.null(ITS_syntax_parser)) check_scalar_character(ITS_syntax_parser, "ITS_syntax_parser")
  check_logical_scalar(ITS_ansi_log, "ITS_ansi_log")
  check_scalar_character(group_delivery_dir, "group_delivery_dir")
  check_scalar_character(barcode_col, "barcode_col")
  if (!is.null(samples_dir)) check_scalar_character(samples_dir, "samples_dir")
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(alignment_tables_dir, "alignment_tables_dir")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  check_scalar_character(consensus_file, "consensus_file")
  check_scalar_character(trimmed_consensus_file, "trimmed_consensus_file")
  check_scalar_character(barcode_pattern, "barcode_pattern")
  check_logical_scalar(run_unite_annotation, "run_unite_annotation")
  check_scalar_character(unite_dir, "unite_dir")
  check_scalar_character(unite_top_hits_name, "unite_top_hits_name")
  check_scalar_character(unite_evalue, "unite_evalue")
  check_scalar_character(unite_blastn, "unite_blastn")
  check_scalar_character(unite_makeblastdb, "unite_makeblastdb")
  check_scalar_character(conda, "conda")
  check_scalar_character(readme_name, "readme_name")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }
  if (!is.null(unite_db)) check_scalar_character(unite_db, "unite_db")
  if (!is.null(unite_db_fasta)) check_file_arg(unite_db_fasta, "unite_db_fasta")
  if (!is.null(unite_conda_env)) check_scalar_character(unite_conda_env, "unite_conda_env")

  if (!is.character(tax_levels) || length(tax_levels) == 0L || anyNA(tax_levels)) {
    stop("`tax_levels` must be a non-empty character vector.", call. = FALSE)
  }
  grouped_delivery <- !is.null(path_sampleInfo_file_list)
  sample_info_groups <- NULL
  if (isTRUE(grouped_delivery)) {
    path_sampleInfo_file_list <- validate_sample_info_file_list(
      path_sampleInfo_file_list
    )
    sample_info_groups <- make_ITS_sample_info_groups(
      path_sampleInfo_file_list,
      barcode_col = barcode_col
    )
  }

  need_consensus <- is.null(consensus_delivery_path) ||
    !dir.exists(consensus_delivery_path)
  need_ITS <- is.null(path_ITS_result) || !dir.exists(path_ITS_result)
  if (!isTRUE(dry_run)) {
    if (!isTRUE(need_ITS)) {
      check_dir_arg(path_ITS_result, "path_ITS_result")
    } else if (!isTRUE(wait)) {
      stop(
        "`path_ITS_result` is NULL or missing, so `run_ITS()` must be run ",
        "before delivery assembly. Use `wait = TRUE` so the ITS workflow can ",
        "finish before the function continues.",
        call. = FALSE
      )
    }
    if (isTRUE(need_consensus)) {
      if (is.null(path_proj) || is.null(path_sampleInfo_file_list)) {
        stop(
          "When `consensus_delivery_path` is NULL or missing, supply ",
          "`path_proj` and `path_sampleInfo_file_list` so ",
          "`make_consensus_delivery()` can be run.",
          call. = FALSE
        )
      }
      check_dir_arg(path_proj, "path_proj")
    }
  }

  if (is.null(consensus_delivery_output)) {
    consensus_delivery_output <- file.path(
      dirname(output_dir),
      paste0(basename(output_dir), "_consensus_work")
    )
  } else {
    check_scalar_character(consensus_delivery_output, "consensus_delivery_output")
  }
  if (is.null(ITS_out_dir)) {
    ITS_out_dir <- file.path(
      dirname(output_dir),
      paste0(basename(output_dir), "_wf_ITS")
    )
  }
  if (is.null(ITS_work_dir)) {
    ITS_work_dir <- file.path(
      dirname(output_dir),
      paste0(basename(output_dir), "_wf_ITS_work")
    )
  }

  if (isTRUE(dry_run)) {
    ITS_plan <- NULL
    if (isTRUE(need_ITS)) {
      ITS_plan <- run_ITS(
        fastq = ITS_fastq,
        out_dir = ITS_out_dir,
        work_dir = ITS_work_dir,
        profile = ITS_profile,
        resume = ITS_resume,
        database_set = ITS_database_set,
        min_len = ITS_min_len,
        max_len = ITS_max_len,
        workflow = ITS_workflow,
        nextflow = ITS_nextflow,
        quiet = ITS_quiet,
        extra_args = ITS_extra_args,
        syntax_parser = ITS_syntax_parser,
        ansi_log = ITS_ansi_log,
        nextflow_env = ITS_nextflow_env,
        dry_run = TRUE,
        echo = FALSE,
        wait = TRUE,
        stdout = stdout,
        stderr = stderr
      )
    }
    return(invisible(make_ITS_delivery_dry_plan(
      path_ITS_result = path_ITS_result,
      ITS_out_dir = ITS_out_dir,
      need_ITS = need_ITS,
      ITS_plan = ITS_plan,
      grouped_delivery = grouped_delivery,
      sample_info_groups = sample_info_groups,
      group_delivery_dir = group_delivery_dir,
      output_dir = output_dir,
      consensus_delivery_path = consensus_delivery_path,
      consensus_delivery_output = consensus_delivery_output,
      need_consensus = need_consensus,
      samples_dir = samples_dir,
      abundance_table = abundance_table,
      figure_dir = figure_dir,
      identification_dir = identification_dir,
      unite_dir = unite_dir,
      unite_top_hits_name = unite_top_hits_name,
      run_unite_annotation = run_unite_annotation,
      readme_name = readme_name,
      chinese_readme_name = chinese_readme_name
    )))
  }

  if (dir.exists(output_dir) && !isTRUE(grouped_delivery)) {
    if (!isTRUE(overwrite)) {
      stop("`path_delivery` already exists. Use `overwrite = TRUE` to replace it: ",
           output_dir, call. = FALSE)
    }
    unlink(output_dir, recursive = TRUE)
  }
  if (isTRUE(grouped_delivery)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
  }
  ITS_result <- NULL
  if (isTRUE(need_ITS)) {
    ITS_result <- run_ITS(
      fastq = ITS_fastq,
      out_dir = ITS_out_dir,
      work_dir = ITS_work_dir,
      profile = ITS_profile,
      resume = ITS_resume,
      database_set = ITS_database_set,
      min_len = ITS_min_len,
      max_len = ITS_max_len,
      workflow = ITS_workflow,
      nextflow = ITS_nextflow,
      quiet = ITS_quiet,
      extra_args = ITS_extra_args,
      syntax_parser = ITS_syntax_parser,
      ansi_log = ITS_ansi_log,
      nextflow_env = ITS_nextflow_env,
      dry_run = FALSE,
      echo = echo,
      wait = wait,
      stdout = stdout,
      stderr = stderr
    )
    path_ITS_result <- ITS_out_dir
  }
  path_ITS_result <- normalizePath(path_ITS_result, mustWork = TRUE)

  consensus_result <- NULL
  if (isTRUE(need_consensus)) {
    consensus_result <- make_consensus_delivery(
      path_proj = path_proj,
      path_sampleInfo_file_list = path_sampleInfo_file_list,
      path_delivery = consensus_delivery_output,
      dry_run = FALSE,
      echo = echo,
      wait = wait,
      stdout = stdout,
      stderr = stderr,
      ...
    )
    consensus_delivery_path <- consensus_delivery_output
  }
  consensus_delivery_path <- normalizePath(consensus_delivery_path, mustWork = TRUE)

  targets <- make_ITS_delivery_targets(
    output_dir = output_dir,
    grouped_delivery = grouped_delivery,
    sample_info_groups = sample_info_groups,
    group_delivery_dir = group_delivery_dir
  )
  deliveries <- list()
  for (target_name in names(targets)) {
    target <- targets[[target_name]]
    if (isTRUE(grouped_delivery)) message(target_name)
    deliveries[[target_name]] <- assemble_single_ITS_delivery(
      path_ITS_result = path_ITS_result,
      output_dir = target$output_dir,
      consensus_delivery_path = resolve_ITS_group_consensus_root(
        consensus_delivery_path,
        target$group
      ),
      barcodes = target$barcodes,
      samples_dir = samples_dir,
      abundance_table = abundance_table,
      alignment_tables_dir = alignment_tables_dir,
      figure_dir = figure_dir,
      identification_dir = identification_dir,
      tax_levels = tax_levels,
      cutoff = cutoff,
      width = width,
      height = height,
      consensus_file = consensus_file,
      trimmed_consensus_file = trimmed_consensus_file,
      barcode_pattern = barcode_pattern,
      run_unite_annotation = run_unite_annotation,
      unite_db = unite_db,
      unite_db_fasta = unite_db_fasta,
      unite_dir = unite_dir,
      unite_top_hits_name = unite_top_hits_name,
      unite_threads = unite_threads,
      unite_max_target_seqs = unite_max_target_seqs,
      unite_evalue = unite_evalue,
      unite_task = unite_task,
      unite_word_size = unite_word_size,
      unite_strand = unite_strand,
      unite_dust = unite_dust,
      unite_perc_identity = unite_perc_identity,
      unite_extra_args = unite_extra_args,
      unite_species_identity = unite_species_identity,
      unite_genus_identity = unite_genus_identity,
      unite_family_identity = unite_family_identity,
      unite_min_query_coverage = unite_min_query_coverage,
      unite_min_reference_coverage = unite_min_reference_coverage,
      unite_novel_identity = unite_novel_identity,
      unite_blastn = unite_blastn,
      unite_makeblastdb = unite_makeblastdb,
      unite_conda_env = unite_conda_env,
      conda = conda,
      readme_name = readme_name,
      chinese_readme_name = chinese_readme_name,
      overwrite = overwrite,
      echo = echo,
      stdout = stdout,
      stderr = stderr
    )
  }

  if (!isTRUE(grouped_delivery)) {
    single <- deliveries[[1L]]
    return(invisible(c(
      single,
      list(
        consensus_delivery_path = consensus_delivery_path,
        consensus_generated = need_consensus,
        consensus_result = consensus_result,
        ITS_generated = need_ITS,
        ITS_result = ITS_result,
        grouped_delivery = FALSE
      )
    )))
  }

  invisible(list(
    path_delivery = output_dir,
    grouped_delivery = TRUE,
    group_delivery_dir = group_delivery_dir,
    sample_info_groups = sample_info_groups,
    consensus_delivery_path = consensus_delivery_path,
    consensus_generated = need_consensus,
    consensus_result = consensus_result,
    ITS_generated = need_ITS,
    ITS_result = ITS_result,
    delivery = deliveries
  ))
}

#' Write README files for an integrated ITS delivery folder
#'
#' @param output_dir Final ITS delivery directory.
#' @param samples_dir,abundance_table,figure_dir,identification_dir,unite_dir,unite_top_hits_name
#'   Directory and file names used in the delivery layout.
#' @param has_unite_annotation Logical. Whether UNITE annotation files were
#'   generated.
#' @param readme_name,chinese_readme_name README filenames.
#'
#' @return Invisibly returns generated README file paths.
#'
#' @export
write_ITS_delivery_readme <- function(output_dir,
                                      samples_dir = "samples",
                                      abundance_table = "abundance_table_genus.tsv",
                                      figure_dir = "figures",
                                      identification_dir = "identification_tables",
                                      unite_dir = "unite_consensus_annotation",
                                      unite_top_hits_name = "unite_consensus_top_hits.tsv",
                                      has_unite_annotation = TRUE,
                                      readme_name = "README.txt",
                                      chinese_readme_name = "README.zh-CN.txt") {
  check_dir_arg(output_dir, "output_dir")
  if (!is.null(samples_dir)) check_scalar_character(samples_dir, "samples_dir")
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  check_scalar_character(unite_dir, "unite_dir")
  check_scalar_character(unite_top_hits_name, "unite_top_hits_name")
  check_logical_scalar(has_unite_annotation, "has_unite_annotation")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }

  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  readme <- file.path(output_dir, readme_name)
  writeLines(
    ITS_delivery_readme(
      samples_dir = samples_dir,
      abundance_table = abundance_table,
      figure_dir = figure_dir,
      identification_dir = identification_dir,
      unite_dir = unite_dir,
      unite_top_hits_name = unite_top_hits_name,
      has_unite_annotation = has_unite_annotation
    ),
    readme,
    useBytes = TRUE
  )

  files <- c(README = readme)
  if (!is.null(chinese_readme_name)) {
    readme_zh <- file.path(output_dir, chinese_readme_name)
    writeLines(
      ITS_delivery_readme_zh(
        samples_dir = samples_dir,
        abundance_table = abundance_table,
        figure_dir = figure_dir,
        identification_dir = identification_dir,
        unite_dir = unite_dir,
        unite_top_hits_name = unite_top_hits_name,
        has_unite_annotation = has_unite_annotation
      ),
      readme_zh,
      useBytes = TRUE
    )
    files <- c(files, README_zh_CN = readme_zh)
  }

  invisible(files)
}

make_ITS_sample_info_groups <- function(path_sampleInfo_file_list,
                                        barcode_col) {
  groups <- list()
  for (group in names(path_sampleInfo_file_list)) {
    sample_info <- read_sample_info_table(path_sampleInfo_file_list[[group]])
    sample_barcode_col <- resolve_barcode_col(sample_info, barcode_col)
    validate_amplicon_delivery_sample_info(
      sample_info = sample_info,
      required_cols = sample_barcode_col
    )
    groups[[group]] <- list(
      sample_info = sample_info,
      barcode_col = sample_barcode_col,
      barcodes = sort(unique(delivery_barcode_names(sample_info, sample_barcode_col))),
      sample_info_file = path_sampleInfo_file_list[[group]]
    )
  }
  groups
}

make_ITS_delivery_targets <- function(output_dir,
                                      grouped_delivery,
                                      sample_info_groups,
                                      group_delivery_dir) {
  if (!isTRUE(grouped_delivery)) {
    return(list(.single = list(
      group = NA_character_,
      output_dir = output_dir,
      barcodes = NULL
    )))
  }

  targets <- list()
  for (group in names(sample_info_groups)) {
    targets[[group]] <- list(
      group = group,
      output_dir = file.path(output_dir, group, group_delivery_dir),
      barcodes = sample_info_groups[[group]]$barcodes
    )
  }
  targets
}

resolve_ITS_group_consensus_root <- function(consensus_delivery_path, group) {
  if (is.na(group) || !nzchar(group)) return(consensus_delivery_path)

  grouped_consensus <- file.path(consensus_delivery_path, group, "consensus_results")
  if (dir.exists(grouped_consensus)) return(grouped_consensus)

  grouped_root <- file.path(consensus_delivery_path, group)
  if (dir.exists(grouped_root)) return(grouped_root)

  consensus_delivery_path
}

assemble_single_ITS_delivery <- function(path_ITS_result,
                                         output_dir,
                                         consensus_delivery_path,
                                         barcodes,
                                         samples_dir,
                                         abundance_table,
                                         alignment_tables_dir,
                                         figure_dir,
                                         identification_dir,
                                         tax_levels,
                                         cutoff,
                                         width,
                                         height,
                                         consensus_file,
                                         trimmed_consensus_file,
                                         barcode_pattern,
                                         run_unite_annotation,
                                         unite_db,
                                         unite_db_fasta,
                                         unite_dir,
                                         unite_top_hits_name,
                                         unite_threads,
                                         unite_max_target_seqs,
                                         unite_evalue,
                                         unite_task,
                                         unite_word_size,
                                         unite_strand,
                                         unite_dust,
                                         unite_perc_identity,
                                         unite_extra_args,
                                         unite_species_identity,
                                         unite_genus_identity,
                                         unite_family_identity,
                                         unite_min_query_coverage,
                                         unite_min_reference_coverage,
                                         unite_novel_identity,
                                         unite_blastn,
                                         unite_makeblastdb,
                                         unite_conda_env,
                                         conda,
                                         readme_name,
                                         chinese_readme_name,
                                         overwrite,
                                         echo,
                                         stdout,
                                         stderr) {
  if (dir.exists(output_dir)) {
    if (!isTRUE(overwrite)) {
      stop("ITS delivery directory already exists. Use `overwrite = TRUE` to replace it: ",
           output_dir, call. = FALSE)
    }
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  sample_root <- if (is.null(samples_dir)) output_dir else file.path(output_dir, samples_dir)
  dir.create(sample_root, recursive = TRUE, showWarnings = FALSE)

  consensus_copy <- copy_consensus_delivery_to_ITS_samples(
    consensus_delivery_path = consensus_delivery_path,
    sample_root = sample_root,
    barcode_pattern = barcode_pattern,
    barcodes = barcodes,
    overwrite = TRUE
  )

  prepared_ITS <- prepare_ITS_result_for_barcodes(
    path_ITS_result = path_ITS_result,
    barcodes = barcodes,
    abundance_table = abundance_table,
    alignment_tables_dir = alignment_tables_dir
  )
  on.exit(if (!is.null(prepared_ITS$tmpdir)) unlink(prepared_ITS$tmpdir, recursive = TRUE),
          add = TRUE)

  tmp_delivery <- tempfile("its-move-")
  on.exit(unlink(tmp_delivery, recursive = TRUE), add = TRUE)
  its_moved <- move_ITS(
    path_result = prepared_ITS$path_result,
    path_delivery = tmp_delivery,
    overwrite = TRUE,
    tax_levels = tax_levels,
    abundance_table = abundance_table,
    alignment_tables_dir = alignment_tables_dir,
    figure_dir = figure_dir,
    identification_dir = identification_dir,
    cutoff = cutoff,
    width = width,
    height = height,
    readme_name = "README.wf-ITS.txt",
    chinese_readme_name = NULL
  )

  moved_its_dir <- its_moved$path_ITS
  copied_abundance <- copy_file_required(
    file.path(moved_its_dir, abundance_table),
    file.path(output_dir, abundance_table),
    overwrite = TRUE
  )
  copied_figures <- copy_directory_contents(
    file.path(moved_its_dir, figure_dir),
    file.path(output_dir, figure_dir),
    overwrite = TRUE
  )
  identification_copy <- distribute_ITS_identification_tables(
    identification_dir = file.path(moved_its_dir, identification_dir),
    sample_root = sample_root,
    destination_dir_name = identification_dir,
    barcode_pattern = barcode_pattern,
    barcodes = barcodes,
    overwrite = TRUE
  )

  consensus_fasta <- collect_ITS_consensus_fasta(
    consensus_delivery_path = consensus_delivery_path,
    output_dir = output_dir,
    trimmed_consensus_file = trimmed_consensus_file,
    consensus_file = consensus_file
  )

  unite_annotation <- NULL
  unite_top_hits <- NA_character_
  if (isTRUE(run_unite_annotation)) {
    if (is.null(unite_db) && is.null(unite_db_fasta)) {
      warning(
        "Skipping UNITE annotation because neither `unite_db` nor ",
        "`unite_db_fasta` was supplied.",
        call. = FALSE
      )
    } else if (is.na(consensus_fasta) || !file.exists(consensus_fasta)) {
      warning(
        "Skipping UNITE annotation because no consensus FASTA was found.",
        call. = FALSE
      )
    } else {
      unite_output_dir <- file.path(output_dir, unite_dir)
      unite_annotation <- annotate_consensus_blast(
        consensus_fasta = consensus_fasta,
        db = unite_db,
        db_fasta = unite_db_fasta,
        out_dir = unite_output_dir,
        output_tsv = file.path(unite_output_dir, "consensus.blast.tsv"),
        prefix = "consensus",
        threads = unite_threads,
        max_target_seqs = unite_max_target_seqs,
        evalue = unite_evalue,
        task = unite_task,
        word_size = unite_word_size,
        strand = unite_strand,
        dust = unite_dust,
        perc_identity = unite_perc_identity,
        extra_args = unite_extra_args,
        species_identity = unite_species_identity,
        genus_identity = unite_genus_identity,
        family_identity = unite_family_identity,
        min_query_coverage = unite_min_query_coverage,
        min_reference_coverage = unite_min_reference_coverage,
        novel_identity = unite_novel_identity,
        blastn = unite_blastn,
        makeblastdb = unite_makeblastdb,
        conda_env = unite_conda_env,
        conda = conda,
        dry_run = FALSE,
        echo = echo,
        stdout = stdout,
        stderr = stderr
      )
      detailed_top_hits <- file.path(unite_output_dir, "consensus.top_hits.tsv")
      utils::write.table(
        unite_annotation$top_hits,
        detailed_top_hits,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      unite_top_hits <- file.path(output_dir, unite_top_hits_name)
      copy_file_required(detailed_top_hits, unite_top_hits, overwrite = TRUE)
    }
  }

  readme_files <- write_ITS_delivery_readme(
    output_dir = output_dir,
    samples_dir = samples_dir,
    abundance_table = abundance_table,
    figure_dir = figure_dir,
    identification_dir = identification_dir,
    unite_dir = unite_dir,
    unite_top_hits_name = unite_top_hits_name,
    has_unite_annotation = !is.null(unite_annotation),
    readme_name = readme_name,
    chinese_readme_name = chinese_readme_name
  )

  list(
    path_delivery = output_dir,
    samples_dir = sample_root,
    barcodes = barcodes,
    consensus_copy = consensus_copy,
    its_result = its_moved,
    prepared_ITS_result = prepared_ITS,
    abundance_table = copied_abundance,
    figures = copied_figures,
    identification_tables = identification_copy,
    consensus_fasta = consensus_fasta,
    unite_annotation = unite_annotation,
    unite_top_hits = unite_top_hits,
    readme_files = readme_files
  )
}

make_ITS_delivery_dry_plan <- function(path_ITS_result,
                                       ITS_out_dir,
                                       need_ITS,
                                       ITS_plan,
                                       grouped_delivery,
                                       sample_info_groups,
                                       group_delivery_dir,
                                       output_dir,
                                       consensus_delivery_path,
                                       consensus_delivery_output,
                                       need_consensus,
                                       samples_dir,
                                       abundance_table,
                                       figure_dir,
                                       identification_dir,
                                       unite_dir,
                                       unite_top_hits_name,
                                       run_unite_annotation,
                                       readme_name,
                                       chinese_readme_name) {
  readme_paths <- list(readme = file.path(output_dir, readme_name))
  if (!is.null(chinese_readme_name)) {
    readme_paths$readme_zh <- file.path(output_dir, chinese_readme_name)
  }
  list(
    status = "dry_run",
    path_delivery = output_dir,
    need_ITS = need_ITS,
    ITS_out_dir = ITS_out_dir,
    grouped_delivery = grouped_delivery,
    group_delivery_dir = group_delivery_dir,
    groups = if (is.null(sample_info_groups)) NULL else stats::setNames(
      lapply(names(sample_info_groups), function(group) {
        list(
          path_delivery = file.path(output_dir, group, group_delivery_dir),
          barcodes = sample_info_groups[[group]]$barcodes
        )
      }),
      names(sample_info_groups)
    ),
    need_consensus = need_consensus,
    consensus_delivery_path = consensus_delivery_path,
    consensus_delivery_output = consensus_delivery_output,
    paths = list(
      path_ITS_result = path_ITS_result,
      effective_ITS_result = if (isTRUE(need_ITS)) ITS_out_dir else path_ITS_result,
      samples = if (is.null(samples_dir)) output_dir else file.path(output_dir, samples_dir),
      abundance_table = file.path(output_dir, abundance_table),
      figures = file.path(output_dir, figure_dir),
      identification_tables = file.path(
        if (is.null(samples_dir)) output_dir else file.path(output_dir, samples_dir),
        "<barcode>",
        identification_dir
      ),
      unite_dir = file.path(output_dir, unite_dir),
      unite_top_hits = file.path(output_dir, unite_top_hits_name),
      readme = readme_paths
    ),
    ITS_plan = ITS_plan,
    run_unite_annotation = run_unite_annotation
  )
}

copy_consensus_delivery_to_ITS_samples <- function(consensus_delivery_path,
                                                  sample_root,
                                                  barcode_pattern,
                                                  barcodes = NULL,
                                                  overwrite) {
  barcode_dirs <- list.dirs(
    consensus_delivery_path,
    recursive = TRUE,
    full.names = TRUE
  )
  barcode_dirs <- barcode_dirs[grepl(barcode_pattern, basename(barcode_dirs))]
  if (!is.null(barcodes)) {
    barcode_dirs <- barcode_dirs[basename(barcode_dirs) %in% barcodes]
  }
  barcode_dirs <- sort(unique(barcode_dirs))

  items <- data.frame(
    barcode = character(),
    source = character(),
    destination = character(),
    copied = logical(),
    stringsAsFactors = FALSE
  )
  for (barcode_dir in barcode_dirs) {
    barcode <- basename(barcode_dir)
    dest <- file.path(sample_root, barcode, "consensus_results")
    copied <- copy_directory_contents(barcode_dir, dest, overwrite = overwrite)
    items <- rbind(
      items,
      data.frame(
        barcode = barcode,
        source = barcode_dir,
        destination = dest,
        copied = length(copied) > 0L || dir.exists(dest),
        stringsAsFactors = FALSE
      )
    )
  }
  items
}

distribute_ITS_identification_tables <- function(identification_dir,
                                                sample_root,
                                                destination_dir_name,
                                                barcode_pattern,
                                                barcodes = NULL,
                                                overwrite) {
  files <- list.files(
    identification_dir,
    pattern = "-alignment-stats[.]tsv$",
    full.names = TRUE
  )
  files <- sort(files)
  items <- data.frame(
    barcode = character(),
    source = character(),
    destination = character(),
    copied = logical(),
    stringsAsFactors = FALSE
  )
  for (file in files) {
    barcode <- sub("-alignment-stats[.]tsv$", "", basename(file))
    if (!grepl(barcode_pattern, barcode)) {
      warning("Could not infer barcode from alignment table: ", file, call. = FALSE)
      next
    }
    if (!is.null(barcodes) && !(barcode %in% barcodes)) {
      next
    }
    dest_dir <- file.path(sample_root, barcode, destination_dir_name)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(dest_dir, basename(file))
    copied <- file.copy(file, dest, overwrite = overwrite)
    items <- rbind(
      items,
      data.frame(
        barcode = barcode,
        source = file,
        destination = dest,
        copied = isTRUE(copied),
        stringsAsFactors = FALSE
      )
    )
  }
  items
}

prepare_ITS_result_for_barcodes <- function(path_ITS_result,
                                            barcodes,
                                            abundance_table,
                                            alignment_tables_dir) {
  if (is.null(barcodes)) {
    return(list(path_result = path_ITS_result, tmpdir = NULL, barcodes = NULL))
  }

  tmpdir <- tempfile("its-result-subset-")
  dir.create(file.path(tmpdir, alignment_tables_dir), recursive = TRUE)

  subset_ITS_abundance_table(
    from = file.path(path_ITS_result, abundance_table),
    to = file.path(tmpdir, abundance_table),
    barcodes = barcodes
  )
  subset_ITS_alignment_tables(
    from_dir = file.path(path_ITS_result, alignment_tables_dir),
    to_dir = file.path(tmpdir, alignment_tables_dir),
    barcodes = barcodes
  )

  list(path_result = tmpdir, tmpdir = tmpdir, barcodes = barcodes)
}

subset_ITS_abundance_table <- function(from, to, barcodes) {
  check_file_arg(from, "abundance_table")
  abundance <- utils::read.delim(
    from,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!("tax" %in% names(abundance))) {
    stop("ITS abundance table must contain a `tax` column: ", from, call. = FALSE)
  }
  keep_barcodes <- intersect(barcodes, names(abundance))
  missing_barcodes <- setdiff(barcodes, names(abundance))
  if (length(missing_barcodes) > 0L) {
    warning(
      "Barcode column(s) were not found in ITS abundance table: ",
      paste(missing_barcodes, collapse = ", "),
      call. = FALSE
    )
  }
  out <- abundance[, c("tax", keep_barcodes), drop = FALSE]
  if (length(keep_barcodes) > 0L) {
    counts <- as.data.frame(lapply(out[keep_barcodes], function(x) {
      suppressWarnings(as.numeric(x))
    }), check.names = FALSE)
    out$total <- rowSums(counts, na.rm = TRUE)
  } else {
    out$total <- 0
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, to, sep = "\t", quote = FALSE, row.names = FALSE)
  normalizePath(to, mustWork = TRUE)
}

subset_ITS_alignment_tables <- function(from_dir, to_dir, barcodes) {
  check_dir_arg(from_dir, "alignment_tables_dir")
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- logical(length(barcodes))
  for (i in seq_along(barcodes)) {
    src <- file.path(from_dir, paste0(barcodes[[i]], "-alignment-stats.tsv"))
    if (!file.exists(src)) {
      warning("ITS alignment table not found for ", barcodes[[i]], ": ", src,
              call. = FALSE)
      copied[[i]] <- FALSE
      next
    }
    copied[[i]] <- file.copy(src, file.path(to_dir, basename(src)), overwrite = TRUE)
  }
  invisible(copied)
}

collect_ITS_consensus_fasta <- function(consensus_delivery_path,
                                        output_dir,
                                        trimmed_consensus_file,
                                        consensus_file) {
  candidates <- list.files(
    consensus_delivery_path,
    pattern = paste0("^", gsub("([.])", "\\\\\\1", trimmed_consensus_file), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(candidates) == 0L) {
    candidates <- list.files(
      consensus_delivery_path,
      pattern = paste0("^", gsub("([.])", "\\\\\\1", consensus_file), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
  }
  candidates <- sort(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) return(NA_character_)

  out <- file.path(output_dir, "consensus_for_unite.fasta")
  lines <- unlist(lapply(candidates, readLines, warn = FALSE), use.names = FALSE)
  writeLines(lines, out, useBytes = TRUE)
  out
}

copy_file_required <- function(from, to, overwrite) {
  if (!file.exists(from)) {
    stop("Required file does not exist: ", from, call. = FALSE)
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(from, to, overwrite = overwrite)
  if (!isTRUE(copied)) {
    stop("Failed to copy: ", from, " -> ", to, call. = FALSE)
  }
  normalizePath(to, mustWork = TRUE)
}

ITS_delivery_readme <- function(samples_dir,
                                abundance_table,
                                figure_dir,
                                identification_dir,
                                unite_dir,
                                unite_top_hits_name,
                                has_unite_annotation) {
  sample_prefix <- if (is.null(samples_dir)) "barcode*/" else paste0(samples_dir, "/barcode*/")
  c(
    "Integrated ITS Delivery",
    "",
    "Overview",
    "This delivery combines three complementary ITS result types: amplicon consensus results, wf-16s ITS taxonomic profiling results, and consensus-vs-UNITE BLAST review results when a UNITE database was supplied.",
    "",
    "Directory Contents",
    paste0("- ", sample_prefix, "consensus_results/: per-barcode consensus-analysis outputs."),
    paste0("- ", sample_prefix, identification_dir, "/: per-barcode wf-16s alignment-stat tables."),
    paste0("- ", abundance_table, ": genus-level abundance table from wf-16s ITS profiling."),
    paste0("- ", figure_dir, "/: abundance bar plots generated from the abundance table."),
    paste0("- ", unite_top_hits_name, ": root-level copy of the recommended UNITE top-hit summary."),
    paste0("- ", unite_dir, "/: detailed consensus-vs-UNITE BLAST results, including `consensus.blast.tsv` and `consensus.top_hits.tsv`."),
    "",
    "How To Interpret",
    "Use the abundance table and figures for sample-level composition summaries. Use the per-barcode alignment-stat tables to review reference-level evidence from wf-16s. Use the UNITE BLAST top-hit table as an independent consensus-sequence review.",
    "The wf-16s ITS results are suitable for routine composition analysis and candidate taxon screening. Strict species confirmation or novel-species assessment should not rely on a single database top hit.",
    "For potential novel species, review consensus quality, percent identity, query coverage, reference coverage, top-hit versus second-hit separation, Species Hypothesis information when available, and whether the closest matches are well-curated or type-material records.",
    "",
    if (isTRUE(has_unite_annotation)) {
      "UNITE annotation was generated for this delivery."
    } else {
      "UNITE annotation was not generated for this delivery because no UNITE database was supplied or no consensus FASTA was found."
    }
  )
}

ITS_delivery_readme_zh <- function(samples_dir,
                                   abundance_table,
                                   figure_dir,
                                   identification_dir,
                                   unite_dir,
                                   unite_top_hits_name,
                                   has_unite_annotation) {
  sample_prefix <- if (is.null(samples_dir)) "barcode*/" else paste0(samples_dir, "/barcode*/")
  c(
    "ITS 综合交付结果说明",
    "",
    "概述",
    "本交付目录整合三类 ITS 结果：扩增子共识序列结果、wf-16s ITS 分类丰度结果，以及在提供 UNITE 数据库时生成的 consensus vs UNITE BLAST 复核结果。",
    "",
    "目录内容",
    paste0("- ", sample_prefix, "consensus_results/：每个 barcode 的共识序列分析结果。"),
    paste0("- ", sample_prefix, identification_dir, "/：每个 barcode 的 wf-16s 比对统计表。"),
    paste0("- ", abundance_table, "：wf-16s ITS 分析得到的属水平丰度表。"),
    paste0("- ", figure_dir, "/：基于丰度表生成的丰度柱状图。"),
    paste0("- ", unite_top_hits_name, "：放在根目录下的 UNITE top-hit 推荐查看汇总表。"),
    paste0("- ", unite_dir, "/：完整的 consensus vs UNITE BLAST 复核结果，包括 `consensus.blast.tsv` 和 `consensus.top_hits.tsv`。"),
    "",
    "结果解读",
    "丰度表和图片适合用于样本整体组成展示；每个 barcode 的比对统计表适合复核 wf-16s 的参考序列命中证据；UNITE BLAST top-hit 表适合作为共识序列层面的独立复核结果。",
    "wf-16s ITS 结果适合常规组成分析和候选分类筛查，但严格种水平确认或新物种判断不应只依赖单个数据库 best hit。",
    "如果需要评估潜在新物种，建议重点查看共识序列质量、percent identity、query coverage、reference coverage、最佳命中与次佳命中的差距、Species Hypothesis 信息，以及最近命中是否来自高质量 curated record 或 type material。",
    "",
    if (isTRUE(has_unite_annotation)) {
      "本次交付已生成 UNITE 共识序列复核结果。"
    } else {
      "本次交付未生成 UNITE 复核结果，原因可能是未提供 UNITE 数据库，或没有找到可用于复核的共识序列 FASTA。"
    }
  )
}
