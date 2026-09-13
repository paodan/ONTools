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
#'   deliveries are created under `path_delivery/<group>/ITS/`.
#' @param consensus_delivery_path Existing consensus delivery directory. If
#'   `NULL` or if the directory does not exist, [make_consensus_delivery()] is
#'   run first using `path_proj`, `path_sampleInfo_file_list`, and the
#'   consensus-related parameters in this function.
#' @param path_proj,path_sampleInfo_file_list Arguments passed to
#'   [make_consensus_delivery()] when `consensus_delivery_path` is missing.
#'   When `path_sampleInfo_file_list` is supplied, its names are also used as
#'   delivery group names, matching [make_consensus_delivery()]. Each sample
#'   information file is copied to the corresponding project root.
#' @param consensus_delivery_output Directory used to store newly generated
#'   consensus delivery results. The default `NULL` creates a sibling work
#'   directory named `basename(path_delivery)_consensus_work`.
#' @param kit_name,model,demux_out,barcode_both_ends Dorado basecalling and
#'   demultiplexing parameters passed to [make_consensus_delivery()] when
#'   consensus results need to be generated. Defaults match
#'   [make_consensus_delivery()].
#' @param fastq_out FASTQ output directory name created by Dorado conversion.
#'   This matches the `fastq_out` argument in [make_consensus_delivery()].
#'   When `make_consensus_delivery()` is run in this function, the actual
#'   per-project FASTQ path passed to [run_ITS()] is taken from the consensus
#'   workflow result, typically `.../fastq_pass_trim/<group>`. If consensus
#'   results are not generated in this call, `fastq_out` may also be an existing
#'   FASTQ directory path. Default is `"fastq_pass_trim"`.
#' @param run_basecalling_demux_step,run_dorado_basecall_step,run_dorado_demux_step,run_dorado_fastq_step,run_QC_step,move_fastq_step,move_fastq_mode,run_amplicon_step,trim_consensus_step,run_filtered_QC_step,run_igv_step,collect_results_step,make_ab1
#'   Consensus-delivery step controls passed to [make_consensus_delivery()] when
#'   `consensus_delivery_path` is `NULL` or missing. Defaults match
#'   [make_consensus_delivery()]. They have no effect when an existing
#'   `consensus_delivery_path` is supplied.
#' @param ab1_name_template,ab1_samtools,consensus_index_file,sample_length_plot_pattern,project_col,amplicon_size_col,min_read_length_col,max_read_length_col,f_primer_col,r_primer_col,min_read_qual,min_n_reads,force_spoa_length_threshold,override_basecaller_cfg,amplicon_extra_args,barcode_digits,dorado_threads,dorado_fastq_write_md5,dorado,samtools,gzip,dorado_conda_env,igv,overwrite_fastq,overwrite_delivery,include_execution
#'   Additional [make_consensus_delivery()] parameters passed through when
#'   consensus results need to be generated. Defaults match
#'   [make_consensus_delivery()].
#' @param path_work Shared work root used to derive default ITS workflow output
#'   and work directories. If `NULL`, sibling directories of `path_delivery` are
#'   used.
#' @param out_dir,work_dir,profile,resume,database_set,min_len,max_len,workflow,nextflow,quiet,extra_args,syntax_parser,ansi_log,nextflow_env
#'   Parameters passed to [run_ITS()] when `path_ITS_result` is `NULL` or
#'   missing. Defaults mirror [run_ITS()] for ITS use: `database_set =
#'   "ncbi_16s_18s_28s_ITS"`, `min_len = 300`, `max_len = 2000`, and
#'   `extra_args = "--minimap2_by_reference"`.
#' @param run_ITS_step Logical. If `TRUE`, run [run_ITS()] when
#'   `path_ITS_result` is `NULL` or missing. Default is `TRUE`.
#' @param move_ITS_step Logical. If `TRUE`, organize wf-16s ITS results into
#'   the delivery structure using [move_ITS()] outputs. Default is `TRUE`.
#' @param ITS_dir Directory name created under each barcode sample directory
#'   for ITS-specific per-sample evidence. Default is `"ITS_results"`,
#'   producing `path_delivery/samples/barcode*/ITS_results/` in single-project
#'   mode and `path_delivery/<group>/samples/barcode*/ITS_results/` in grouped
#'   mode.
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
#'   `abundance_table` and `figure_dir` are written at the project root, while
#'   each barcode alignment table is copied to
#'   `samples/barcode*/ITS_results/identification_tables/`.
#' @param tax_levels,cutoff,width,height Plotting parameters passed to
#'   [move_ITS()]. By default, plots are generated for Kingdom, Phylum, Class,
#'   Order, Family, and Genus; taxa below `cutoff = 0.01` relative abundance
#'   are grouped by the plotting helper; figures are saved at 12 x 6 inches.
#' @param consensus_file,trimmed_consensus_file Consensus FASTA names searched
#'   under the consensus delivery. The default prefers
#'   `"all-consensus-seqs_trimmed.fasta"` and falls back to
#'   `"all-consensus-seqs.fasta"` when the trimmed file is absent. The selected
#'   consensus sequences are copied or collected at the project root as
#'   `consensus_file` and used for optional UNITE review. When trimmed consensus
#'   is available, it is also preserved under `trimmed_consensus_file`; matching
#'   FASTA indexes are copied when available.
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
#'
#' @return Invisibly returns a list describing generated/copied outputs.
#'
#' @details
#' The function performs four main tasks:
#'
#' 1. Prepare consensus-delivery content. If `consensus_delivery_path` exists,
#'    its barcode folders are copied into the final delivery. Otherwise,
#'    `make_consensus_delivery()` is run first. This step also creates or
#'    locates the grouped FASTQ directory used as input for ITS profiling.
#' 2. Prepare wf-16s ITS content. If `path_ITS_result` was not supplied,
#'    `run_ITS()` is called after consensus preparation, using the grouped FASTQ
#'    path returned by [make_consensus_delivery()]. `move_ITS()` is then called
#'    in a temporary directory to generate abundance plots and normalize wf-16s
#'    outputs.
#' 3. Reorganize files. When `path_sampleInfo_file_list` is supplied, one
#'    delivery directory is created for each named sample-info file:
#'    `path_delivery/<group>/ITS/`. Each group keeps only the barcode samples
#'    listed in its sample-info table. Without `path_sampleInfo_file_list`,
#'    `path_delivery` itself is treated as the final ITS delivery directory.
#'    Per-barcode consensus outputs are copied to
#'    `samples/barcode*/consensus_results/`, per-barcode wf-16s alignment
#'    tables are copied to
#'    `samples/barcode*/ITS_results/identification_tables/`, and project-level
#'    abundance files, figures, sample information, and consensus FASTA are
#'    copied to the project root.
#' 4. Optionally run consensus-vs-UNITE review. Consensus FASTA files are
#'    collected into the project-root `all-consensus-seqs.fasta`, BLAST output
#'    is written to `unite_consensus_annotation/consensus.blast.tsv`, and the
#'    top-hit summary is written to the root-level
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
#'   fastq_out = "fastq_pass_trim",
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
                              kit_name = "YS-NB576",
                              model = "sup",
                              demux_out = NULL,
                              fastq_out = "fastq_pass_trim",
                              barcode_both_ends = FALSE,
                              run_basecalling_demux_step = TRUE,
                              run_dorado_basecall_step = run_basecalling_demux_step,
                              run_dorado_demux_step = run_basecalling_demux_step,
                              run_dorado_fastq_step = run_basecalling_demux_step,
                              run_QC_step = TRUE,
                              move_fastq_step = TRUE,
                              move_fastq_mode = c("move", "reuse", "auto"),
                              run_amplicon_step = TRUE,
                              trim_consensus_step = TRUE,
                              run_filtered_QC_step = TRUE,
                              run_igv_step = TRUE,
                              collect_results_step = TRUE,
                              make_ab1 = TRUE,
                              ab1_name_template = "{barcode}.synthetic.ab1",
                              ab1_samtools = "samtools",
                              consensus_file = "all-consensus-seqs.fasta",
                              consensus_index_file = paste0(consensus_file, ".fai"),
                              trimmed_consensus_file = "all-consensus-seqs_trimmed.fasta",
                              barcode_pattern = "^barcode[0-9]+$",
                              sample_length_plot_pattern = "^Distribution_seqLength__.*\\.png$",
                              readme_name = "README.txt",
                              chinese_readme_name = "README.zh-CN.txt",
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
                              path_work = NULL,
                              out_dir = NULL,
                              work_dir = NULL,
                              profile = "standard",
                              resume = TRUE,
                              amplicon_extra_args = NULL,
                              barcode_digits = 3,
                              dorado_threads = NULL,
                              dorado_fastq_write_md5 = TRUE,
                              dorado = "dorado",
                              samtools = "samtools",
                              gzip = "gzip",
                              dorado_conda_env = NULL,
                              conda = "conda",
                              igv = "/usr/local/bin/IGV_Linux_2.19.8/igv.sh",
                              overwrite_fastq = FALSE,
                              overwrite_delivery = TRUE,
                              include_execution = FALSE,
                              run_ITS_step = TRUE,
                              move_ITS_step = TRUE,
                              database_set = "ncbi_16s_18s_28s_ITS",
                              min_len = 300,
                              max_len = 900,
                              workflow = "epi2me-labs/wf-16s",
                              nextflow = "nextflow",
                              quiet = FALSE,
                              extra_args = "--minimap2_by_reference",
                              syntax_parser = "v1",
                              ansi_log = FALSE,
                              nextflow_env = NULL,
                              ITS_dir = "ITS_results",
                              samples_dir = "samples",
                              abundance_table = "abundance_table_genus.tsv",
                              alignment_tables_dir = "alignment_tables",
                              figure_dir = "figures",
                              identification_dir = "identification_tables",
                              tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                              cutoff = 0.01,
                              width = 12,
                              height = 6,
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
                              overwrite = FALSE,
                              dry_run = FALSE,
                              echo = TRUE,
                              wait = TRUE,
                              stdout = "",
                              stderr = "") {
  if (!is.null(path_ITS_result)) {
    check_scalar_character(path_ITS_result, "path_ITS_result")
  }
  check_scalar_character(path_delivery, "path_delivery")
  output_dir <- path_delivery
  if (!is.null(consensus_delivery_path)) {
    check_scalar_character(consensus_delivery_path, "consensus_delivery_path")
  }
  if (!is.null(path_proj)) check_scalar_character(path_proj, "path_proj")
  check_scalar_character(kit_name, "kit_name")
  check_scalar_character(model, "model")
  if (!is.null(demux_out)) check_scalar_character(demux_out, "demux_out")
  check_scalar_character(fastq_out, "fastq_out")
  check_logical_scalar(barcode_both_ends, "barcode_both_ends")
  check_logical_scalar(run_basecalling_demux_step, "run_basecalling_demux_step")
  check_logical_scalar(run_dorado_basecall_step, "run_dorado_basecall_step")
  check_logical_scalar(run_dorado_demux_step, "run_dorado_demux_step")
  check_logical_scalar(run_dorado_fastq_step, "run_dorado_fastq_step")
  check_logical_scalar(run_QC_step, "run_QC_step")
  check_logical_scalar(move_fastq_step, "move_fastq_step")
  move_fastq_mode <- match.arg(move_fastq_mode)
  check_logical_scalar(run_amplicon_step, "run_amplicon_step")
  check_logical_scalar(trim_consensus_step, "trim_consensus_step")
  check_logical_scalar(run_filtered_QC_step, "run_filtered_QC_step")
  check_logical_scalar(run_igv_step, "run_igv_step")
  check_logical_scalar(collect_results_step, "collect_results_step")
  check_logical_scalar(make_ab1, "make_ab1")
  check_scalar_character(ab1_name_template, "ab1_name_template")
  check_scalar_character(ab1_samtools, "ab1_samtools")
  check_scalar_character(consensus_file, "consensus_file")
  check_scalar_character(consensus_index_file, "consensus_index_file")
  check_scalar_character(trimmed_consensus_file, "trimmed_consensus_file")
  check_scalar_character(barcode_pattern, "barcode_pattern")
  check_scalar_character(sample_length_plot_pattern, "sample_length_plot_pattern")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }
  check_scalar_character(project_col, "project_col")
  check_scalar_character(amplicon_size_col, "amplicon_size_col")
  check_scalar_character(barcode_col, "barcode_col")
  check_scalar_character(min_read_length_col, "min_read_length_col")
  check_scalar_character(max_read_length_col, "max_read_length_col")
  check_scalar_character(f_primer_col, "f_primer_col")
  check_scalar_character(r_primer_col, "r_primer_col")
  min_read_qual <- validate_nonnegative_number(min_read_qual, "min_read_qual")
  min_n_reads <- validate_positive_integer(min_n_reads, "min_n_reads")
  force_spoa_length_threshold <- validate_positive_integer(
    force_spoa_length_threshold,
    "force_spoa_length_threshold"
  )
  check_scalar_character(override_basecaller_cfg, "override_basecaller_cfg")
  if (!is.null(path_work)) check_scalar_character(path_work, "path_work")
  if (!is.null(out_dir)) check_scalar_character(out_dir, "out_dir")
  if (!is.null(work_dir)) check_scalar_character(work_dir, "work_dir")
  check_scalar_character(profile, "profile")
  check_logical_scalar(resume, "resume")
  check_scalar_character(database_set, "database_set")
  min_len <- validate_positive_integer(min_len, "min_len")
  max_len <- validate_positive_integer(max_len, "max_len")
  check_scalar_character(workflow, "workflow")
  check_scalar_character(nextflow, "nextflow")
  check_logical_scalar(quiet, "quiet")
  if (!is.null(extra_args)) check_scalar_character(extra_args, "extra_args")
  if (!is.null(syntax_parser)) check_scalar_character(syntax_parser, "syntax_parser")
  check_logical_scalar(ansi_log, "ansi_log")
  if (!is.null(amplicon_extra_args)) check_scalar_character(amplicon_extra_args, "amplicon_extra_args")
  barcode_digits <- validate_positive_integer(barcode_digits, "barcode_digits")
  dorado_threads <- validate_optional_positive_integer(dorado_threads, "dorado_threads")
  check_logical_scalar(dorado_fastq_write_md5, "dorado_fastq_write_md5")
  check_scalar_character(dorado, "dorado")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(gzip, "gzip")
  if (!is.null(dorado_conda_env)) check_scalar_character(dorado_conda_env, "dorado_conda_env")
  check_scalar_character(conda, "conda")
  check_scalar_character(igv, "igv")
  check_logical_scalar(overwrite_fastq, "overwrite_fastq")
  check_logical_scalar(overwrite_delivery, "overwrite_delivery")
  check_logical_scalar(include_execution, "include_execution")
  check_logical_scalar(run_ITS_step, "run_ITS_step")
  check_logical_scalar(move_ITS_step, "move_ITS_step")
  check_scalar_character(ITS_dir, "ITS_dir")
  if (!is.null(samples_dir)) check_scalar_character(samples_dir, "samples_dir")
  check_scalar_character(abundance_table, "abundance_table")
  check_scalar_character(alignment_tables_dir, "alignment_tables_dir")
  check_scalar_character(figure_dir, "figure_dir")
  check_scalar_character(identification_dir, "identification_dir")
  check_logical_scalar(run_unite_annotation, "run_unite_annotation")
  check_scalar_character(unite_dir, "unite_dir")
  check_scalar_character(unite_top_hits_name, "unite_top_hits_name")
  check_scalar_character(unite_evalue, "unite_evalue")
  check_scalar_character(unite_blastn, "unite_blastn")
  check_scalar_character(unite_makeblastdb, "unite_makeblastdb")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  check_logical_scalar(wait, "wait")
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
    } else if (!isTRUE(run_ITS_step)) {
      stop(
        "`path_ITS_result` is NULL or missing, but `run_ITS_step = FALSE`. ",
        "Supply an existing `path_ITS_result` or set `run_ITS_step = TRUE`.",
        call. = FALSE
      )
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
  default_work_root <- if (is.null(path_work)) dirname(output_dir) else path_work
  if (is.null(out_dir)) {
    out_dir <- file.path(
      default_work_root,
      paste0(basename(output_dir), "_wf_ITS")
    )
  }
  if (is.null(work_dir)) {
    work_dir <- file.path(
      default_work_root,
      paste0(basename(output_dir), "_wf_ITS_work")
    )
  }
  targets <- make_ITS_delivery_targets(
    output_dir = output_dir,
    grouped_delivery = grouped_delivery,
    sample_info_groups = sample_info_groups
  )

  if (!isTRUE(dry_run)) {
    check_ITS_delivery_targets_available(targets, overwrite = overwrite)
  }

  if (isTRUE(dry_run)) {
    ITS_plan <- NULL
    if (isTRUE(need_ITS) && isTRUE(run_ITS_step)) {
      ITS_plan <- run_ITS(
        fastq = fastq_out,
        out_dir = out_dir,
        work_dir = work_dir,
        profile = profile,
        resume = resume,
        database_set = database_set,
        min_len = min_len,
        max_len = max_len,
        workflow = workflow,
        nextflow = nextflow,
        quiet = quiet,
        extra_args = extra_args,
        syntax_parser = syntax_parser,
        ansi_log = ansi_log,
        nextflow_env = nextflow_env,
        dry_run = TRUE,
        echo = FALSE,
        wait = TRUE,
        stdout = stdout,
        stderr = stderr
      )
    }
    return(invisible(make_ITS_delivery_dry_plan(
      path_ITS_result = path_ITS_result,
      out_dir = out_dir,
      need_ITS = need_ITS,
      ITS_plan = ITS_plan,
      grouped_delivery = grouped_delivery,
      sample_info_groups = sample_info_groups,
      ITS_dir = ITS_dir,
      output_dir = output_dir,
      consensus_delivery_path = consensus_delivery_path,
      consensus_delivery_output = consensus_delivery_output,
      need_consensus = need_consensus,
      run_basecalling_demux_step = run_basecalling_demux_step,
      run_dorado_basecall_step = run_dorado_basecall_step,
      run_dorado_demux_step = run_dorado_demux_step,
      run_dorado_fastq_step = run_dorado_fastq_step,
      run_QC_step = run_QC_step,
      move_fastq_step = move_fastq_step,
      move_fastq_mode = move_fastq_mode,
      run_amplicon_step = run_amplicon_step,
      trim_consensus_step = trim_consensus_step,
      run_filtered_QC_step = run_filtered_QC_step,
      run_igv_step = run_igv_step,
      collect_results_step = collect_results_step,
      make_ab1 = make_ab1,
      run_ITS_step = run_ITS_step,
      move_ITS_step = move_ITS_step,
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
  consensus_result <- NULL
  if (isTRUE(need_consensus)) {
    consensus_result <- make_consensus_delivery(
      path_proj = path_proj,
      path_sampleInfo_file_list = path_sampleInfo_file_list,
      path_delivery = consensus_delivery_output,
      kit_name = kit_name,
      model = model,
      demux_out = demux_out,
      fastq_out = fastq_out,
      barcode_both_ends = barcode_both_ends,
      run_basecalling_demux_step = run_basecalling_demux_step,
      run_dorado_basecall_step = run_dorado_basecall_step,
      run_dorado_demux_step = run_dorado_demux_step,
      run_dorado_fastq_step = run_dorado_fastq_step,
      run_QC_step = run_QC_step,
      move_fastq_step = move_fastq_step,
      move_fastq_mode = move_fastq_mode,
      run_amplicon_step = run_amplicon_step,
      trim_consensus_step = trim_consensus_step,
      run_filtered_QC_step = run_filtered_QC_step,
      run_igv_step = run_igv_step,
      collect_results_step = collect_results_step,
      make_ab1 = make_ab1,
      ab1_name_template = ab1_name_template,
      ab1_samtools = ab1_samtools,
      consensus_file = consensus_file,
      consensus_index_file = consensus_index_file,
      trimmed_consensus_file = trimmed_consensus_file,
      barcode_pattern = barcode_pattern,
      sample_length_plot_pattern = sample_length_plot_pattern,
      readme_name = readme_name,
      chinese_readme_name = chinese_readme_name,
      project_col = project_col,
      amplicon_size_col = amplicon_size_col,
      barcode_col = barcode_col,
      min_read_length_col = min_read_length_col,
      max_read_length_col = max_read_length_col,
      f_primer_col = f_primer_col,
      r_primer_col = r_primer_col,
      min_read_qual = min_read_qual,
      min_n_reads = min_n_reads,
      force_spoa_length_threshold = force_spoa_length_threshold,
      override_basecaller_cfg = override_basecaller_cfg,
      profile = profile,
      resume = resume,
      amplicon_extra_args = amplicon_extra_args,
      barcode_digits = barcode_digits,
      dorado_threads = dorado_threads,
      dorado_fastq_write_md5 = dorado_fastq_write_md5,
      dorado = dorado,
      samtools = samtools,
      gzip = gzip,
      dorado_conda_env = dorado_conda_env,
      conda = conda,
      igv = igv,
      overwrite_fastq = overwrite_fastq,
      overwrite_delivery = overwrite_delivery,
      include_execution = include_execution,
      dry_run = FALSE,
      echo = echo,
      wait = wait,
      stdout = stdout,
      stderr = stderr
    )
    consensus_delivery_path <- consensus_delivery_output
  }
  consensus_delivery_path <- normalizePath(consensus_delivery_path, mustWork = TRUE)

  ITS_results <- list()
  path_ITS_results <- list()
  for (target_name in names(targets)) {
    target <- targets[[target_name]]
    if (isTRUE(need_ITS) && isTRUE(run_ITS_step)) {
      group_fastq <- resolve_ITS_fastq_input(
        target = target,
        consensus_result = consensus_result,
        fastq_out = fastq_out
      )
      target_out_dir <- if (isTRUE(grouped_delivery)) {
        file.path(out_dir, target_name)
      } else {
        out_dir
      }
      target_work_dir <- if (isTRUE(grouped_delivery)) {
        file.path(work_dir, target_name)
      } else {
        work_dir
      }
      ITS_results[[target_name]] <- run_ITS(
        fastq = group_fastq,
        out_dir = target_out_dir,
        work_dir = target_work_dir,
        profile = profile,
        resume = resume,
        database_set = database_set,
        min_len = min_len,
        max_len = max_len,
        workflow = workflow,
        nextflow = nextflow,
        quiet = quiet,
        extra_args = extra_args,
        syntax_parser = syntax_parser,
        ansi_log = ansi_log,
        nextflow_env = nextflow_env,
        dry_run = FALSE,
        echo = echo,
        wait = wait,
        stdout = stdout,
        stderr = stderr
      )
      path_ITS_results[[target_name]] <- normalizePath(target_out_dir, mustWork = TRUE)
    } else {
      path_ITS_results[[target_name]] <- normalizePath(path_ITS_result, mustWork = TRUE)
    }
  }

  deliveries <- list()
  for (target_name in names(targets)) {
    target <- targets[[target_name]]
    if (isTRUE(grouped_delivery)) message(target_name)
    deliveries[[target_name]] <- assemble_single_ITS_delivery(
      path_ITS_result = path_ITS_results[[target_name]],
      output_dir = target$output_dir,
      sample_info_file = target$sample_info_file,
      consensus_delivery_path = resolve_ITS_group_consensus_root(
        consensus_delivery_path,
        target$group
      ),
      barcodes = target$barcodes,
      samples_dir = samples_dir,
      ITS_dir = ITS_dir,
      move_ITS_step = move_ITS_step,
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
        ITS_result = if (length(ITS_results) == 0L) NULL else ITS_results[[1L]],
        grouped_delivery = FALSE
      )
    )))
  }

  invisible(list(
    path_delivery = output_dir,
    grouped_delivery = TRUE,
    ITS_dir = ITS_dir,
    sample_info_groups = sample_info_groups,
    consensus_delivery_path = consensus_delivery_path,
    consensus_generated = need_consensus,
    consensus_result = consensus_result,
    ITS_generated = need_ITS,
    ITS_result = ITS_results,
    path_ITS_results = path_ITS_results,
    delivery = deliveries
  ))
}

#' Write README files for an integrated ITS delivery folder
#'
#' @param output_dir Final ITS delivery directory.
#' @param samples_dir,ITS_dir,abundance_table,figure_dir,identification_dir,unite_dir,unite_top_hits_name
#'   Directory and file names used in the delivery layout. `ITS_dir` is the
#'   ITS-specific subdirectory below each barcode sample directory.
#' @param has_unite_annotation Logical. Whether UNITE annotation files were
#'   generated.
#' @param readme_name,chinese_readme_name README filenames.
#'
#' @return Invisibly returns generated README file paths.
#'
#' @export
write_ITS_delivery_readme <- function(output_dir,
                                      samples_dir = "samples",
                                      ITS_dir = "ITS_results",
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
  check_scalar_character(ITS_dir, "ITS_dir")
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
      ITS_dir = ITS_dir,
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
        ITS_dir = ITS_dir,
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
                                      sample_info_groups) {
  if (!isTRUE(grouped_delivery)) {
    return(list(.single = list(
      group = NA_character_,
      output_dir = output_dir,
      barcodes = NULL,
      sample_info_file = NULL
    )))
  }

  targets <- list()
  for (group in names(sample_info_groups)) {
    targets[[group]] <- list(
      group = group,
      output_dir = file.path(output_dir, group, "ITS"),
      barcodes = sample_info_groups[[group]]$barcodes,
      sample_info_file = sample_info_groups[[group]]$sample_info_file
    )
  }
  targets
}

check_ITS_delivery_targets_available <- function(targets, overwrite) {
  existing <- vapply(targets, function(target) {
    dir.exists(target$output_dir)
  }, logical(1L))
  if (any(existing) && !isTRUE(overwrite)) {
    stop(
      "ITS delivery directory already exists. Use `overwrite = TRUE` to replace it: ",
      paste(vapply(targets[existing], `[[`, character(1L), "output_dir"), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

resolve_ITS_group_consensus_root <- function(consensus_delivery_path, group) {
  if (is.na(group) || !nzchar(group)) return(consensus_delivery_path)

  grouped_consensus <- file.path(consensus_delivery_path, group, "consensus_results")
  if (dir.exists(grouped_consensus)) return(grouped_consensus)

  grouped_root <- file.path(consensus_delivery_path, group)
  if (dir.exists(grouped_root)) return(grouped_root)

  consensus_delivery_path
}

resolve_ITS_fastq_input <- function(target, consensus_result, fastq_out) {
  if (!is.null(consensus_result) && !is.null(consensus_result$workflow)) {
    workflow <- consensus_result$workflow
    group <- target$group
    if (!is.na(group) && nzchar(group) && group %in% names(workflow)) {
      group_fastq <- workflow[[group]]$paths$fastq
      if (!is.null(group_fastq) && length(group_fastq) == 1L && !is.na(group_fastq)) {
        return(group_fastq)
      }
    }
    if (length(workflow) == 1L) {
      group_fastq <- workflow[[1L]]$paths$fastq
      if (!is.null(group_fastq) && length(group_fastq) == 1L && !is.na(group_fastq)) {
        return(group_fastq)
      }
    }
  }

  if (!is.na(target$group) && nzchar(target$group)) {
    grouped_fastq <- file.path(fastq_out, target$group)
    if (dir.exists(grouped_fastq)) return(normalizePath(grouped_fastq, mustWork = TRUE))
  }
  fastq_out
}

assemble_single_ITS_delivery <- function(path_ITS_result,
                                         output_dir,
                                         sample_info_file,
                                         consensus_delivery_path,
                                         barcodes,
                                         samples_dir,
                                         ITS_dir,
                                         move_ITS_step,
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

  copied_sample_info <- copy_optional_sample_info_file(
    sample_info_file = sample_info_file,
    output_dir = output_dir,
    overwrite = TRUE
  )

  consensus_copy <- copy_consensus_delivery_to_ITS_samples(
    consensus_delivery_path = consensus_delivery_path,
    sample_root = sample_root,
    barcode_pattern = barcode_pattern,
    barcodes = barcodes,
    overwrite = TRUE
  )

  prepared_ITS <- NULL
  its_moved <- NULL
  copied_abundance <- NA_character_
  copied_figures <- character()
  identification_copy <- data.frame(
    barcode = character(),
    source = character(),
    destination = character(),
    copied = logical(),
    stringsAsFactors = FALSE
  )
  if (isTRUE(move_ITS_step)) {
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
      ITS_dir = ITS_dir,
      destination_dir_name = identification_dir,
      barcode_pattern = barcode_pattern,
      barcodes = barcodes,
      overwrite = TRUE
    )
  }

  consensus_fasta <- collect_ITS_consensus_fasta(
    consensus_delivery_path = consensus_delivery_path,
    output_dir = output_dir,
    trimmed_consensus_file = trimmed_consensus_file,
    consensus_file = consensus_file,
    output_name = consensus_file
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
      utils::write.table(
        format_unite_top_hits(unite_annotation$top_hits),
        file.path(output_dir, unite_top_hits_name),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      unite_top_hits <- file.path(output_dir, unite_top_hits_name)
    }
  }

  readme_files <- write_ITS_delivery_readme(
    output_dir = output_dir,
    samples_dir = samples_dir,
    ITS_dir = ITS_dir,
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
    ITS_dir = ITS_dir,
    sample_info_file = copied_sample_info,
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

format_unite_top_hits <- function(top_hits) {
  key_cols <- c(
    "qseqid", "qlen", "qstart", "qend",
    "sseqid", "slen", "sstart", "send",
    "length", "pident", "qcovs", "query_coverage", "reference_coverage",
    "mismatch", "gapopen", "evalue", "bitscore",
    "salltitles", "staxids", "taxonomy_path",
    "kingdom", "phylum", "class", "order", "family", "genus", "species",
    "rank", "annotation_level", "novel_candidate"
  )
  ordered <- intersect(key_cols, names(top_hits))
  remaining <- setdiff(names(top_hits), ordered)
  top_hits[, c(ordered, remaining), drop = FALSE]
}

make_ITS_delivery_dry_plan <- function(path_ITS_result,
                                       out_dir,
                                       need_ITS,
                                       ITS_plan,
                                       grouped_delivery,
                                       sample_info_groups,
                                       ITS_dir,
                                       output_dir,
                                       consensus_delivery_path,
                                       consensus_delivery_output,
                                       need_consensus,
                                       run_basecalling_demux_step,
                                       run_dorado_basecall_step,
                                       run_dorado_demux_step,
                                       run_dorado_fastq_step,
                                       run_QC_step,
                                       move_fastq_step,
                                       move_fastq_mode,
                                       run_amplicon_step,
                                       trim_consensus_step,
                                       run_filtered_QC_step,
                                       run_igv_step,
                                       collect_results_step,
                                       make_ab1,
                                       run_ITS_step,
                                       move_ITS_step,
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
    out_dir = out_dir,
    grouped_delivery = grouped_delivery,
    ITS_dir = ITS_dir,
    groups = if (is.null(sample_info_groups)) NULL else stats::setNames(
      lapply(names(sample_info_groups), function(group) {
        list(
          path_delivery = file.path(output_dir, group),
          ITS_delivery = file.path(output_dir, group, "ITS"),
          barcodes = sample_info_groups[[group]]$barcodes
        )
      }),
      names(sample_info_groups)
    ),
    need_consensus = need_consensus,
    consensus_delivery_path = consensus_delivery_path,
    consensus_delivery_output = consensus_delivery_output,
    consensus_steps = list(
      run_basecalling_demux_step = run_basecalling_demux_step,
      run_dorado_basecall_step = run_dorado_basecall_step,
      run_dorado_demux_step = run_dorado_demux_step,
      run_dorado_fastq_step = run_dorado_fastq_step,
      run_QC_step = run_QC_step,
      move_fastq_step = move_fastq_step,
      move_fastq_mode = move_fastq_mode,
      run_amplicon_step = run_amplicon_step,
      trim_consensus_step = trim_consensus_step,
      run_filtered_QC_step = run_filtered_QC_step,
      run_igv_step = run_igv_step,
      collect_results_step = collect_results_step,
      make_ab1 = make_ab1
    ),
    ITS_steps = list(
      run_ITS_step = run_ITS_step,
      move_ITS_step = move_ITS_step
    ),
    paths = list(
      path_ITS_result = path_ITS_result,
      effective_ITS_result = if (isTRUE(need_ITS)) out_dir else path_ITS_result,
      samples = if (is.null(samples_dir)) output_dir else file.path(output_dir, samples_dir),
      abundance_table = file.path(output_dir, abundance_table),
      figures = file.path(output_dir, figure_dir),
      identification_tables = file.path(
        if (is.null(samples_dir)) output_dir else file.path(output_dir, samples_dir),
        "<barcode>",
        ITS_dir,
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
                                                ITS_dir,
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
    dest_dir <- file.path(sample_root, barcode, ITS_dir, destination_dir_name)
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
                                        consensus_file,
                                        output_name = consensus_file) {
  trimmed_candidates <- list.files(
    consensus_delivery_path,
    pattern = paste0("^", gsub("([.])", "\\\\\\1", trimmed_consensus_file), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  trimmed_candidates <- sort(trimmed_candidates[file.exists(trimmed_candidates)])
  using_trimmed <- length(trimmed_candidates) > 0L
  candidates <- trimmed_candidates
  if (!isTRUE(using_trimmed)) {
    candidates <- list.files(
      consensus_delivery_path,
      pattern = paste0("^", gsub("([.])", "\\\\\\1", consensus_file), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    candidates <- sort(candidates[file.exists(candidates)])
  }
  if (length(candidates) == 0L) return(NA_character_)

  out <- file.path(output_dir, output_name)
  write_or_copy_consensus_candidates(candidates, out)

  if (isTRUE(using_trimmed)) {
    trimmed_out <- file.path(output_dir, trimmed_consensus_file)
    if (!identical(normalizePath(out, mustWork = FALSE),
                   normalizePath(trimmed_out, mustWork = FALSE))) {
      write_or_copy_consensus_candidates(candidates, trimmed_out)
    }
  }
  normalizePath(out, mustWork = TRUE)
}

write_or_copy_consensus_candidates <- function(candidates, out) {
  if (length(candidates) == 1L) {
    copy_file_required(candidates[[1L]], out, overwrite = TRUE)
    index <- paste0(candidates[[1L]], ".fai")
    if (file.exists(index)) {
      copy_file_required(index, paste0(out, ".fai"), overwrite = TRUE)
    }
    return(invisible(out))
  }

  lines <- unlist(lapply(candidates, readLines, warn = FALSE), use.names = FALSE)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, out, useBytes = TRUE)
  invisible(out)
}

copy_optional_sample_info_file <- function(sample_info_file,
                                           output_dir,
                                           overwrite) {
  if (is.null(sample_info_file) || is.na(sample_info_file) || !nzchar(sample_info_file)) {
    return(NA_character_)
  }
  check_file_arg(sample_info_file, "sample_info_file")
  copy_file_required(
    sample_info_file,
    file.path(output_dir, basename(sample_info_file)),
    overwrite = overwrite
  )
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
                                ITS_dir,
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
    "This delivery combines three complementary ITS result types:",
    "1. Amplicon consensus results generated for each barcode/sample.",
    "2. wf-16s-style ITS taxonomic profiling results generated from the demultiplexed FASTQ reads.",
    "3. Optional consensus-vs-UNITE BLAST review results generated from the project-level consensus FASTA when a UNITE database was supplied. UNITE reference data can be obtained from https://unite.ut.ee/repository.php.",
    "",
    "These result types answer related but different questions. The abundance table and figures summarize read-level taxonomic composition. The per-barcode alignment-stat tables show which reference records were supported by reads. The UNITE BLAST tables review the final consensus sequences against a curated ITS database.",
    "",
    "Directory Contents",
    paste0("- ", sample_prefix, "consensus_results/: per-barcode consensus-analysis outputs copied from the amplicon consensus delivery."),
    paste0("- ", sample_prefix, "consensus_results/consensus/consensus.fastq: consensus sequence for this barcode/sample, when available."),
    paste0("- ", sample_prefix, "consensus_results/alignments/: read-to-consensus alignment files, BAM indexes, and IGV snapshot PNGs, when generated."),
    paste0("- ", sample_prefix, "consensus_results/*.synthetic.ab1: synthetic Sanger-like trace generated from the read pileup against the consensus, when AB1 generation was enabled."),
    paste0("- ", sample_prefix, "consensus_results/Distribution_seqLength__*.png: read-length distribution plots before and/or after filtering."),
    paste0("- ", sample_prefix, ITS_dir, "/: per-barcode ITS profiling evidence."),
    paste0("- ", sample_prefix, ITS_dir, "/", identification_dir, "/: per-barcode wf-16s alignment-stat tables. Each `barcode*-alignment-stats.tsv` file summarizes database reference hits for one barcode."),
    paste0("- ", abundance_table, ": genus-level abundance table from wf-16s ITS profiling. This is the main table for sample-level composition summaries."),
    paste0("- ", figure_dir, "/: abundance bar plots generated from the abundance table. Files ending in `_percentage.png` show relative abundance; files ending in `_count.png` show read counts."),
    "- all-consensus-seqs.fasta: project-level consensus FASTA used for optional UNITE BLAST review. If a trimmed consensus FASTA was available, the trimmed sequences are preferred and also written under this filename.",
    "- all-consensus-seqs.fasta.fai: FASTA index for all-consensus-seqs.fasta, when available.",
    "- all-consensus-seqs_trimmed.fasta: trimmed project-level consensus FASTA, included when it is available in the consensus results.",
    "- all-consensus-seqs_trimmed.fasta.fai: FASTA index for the trimmed consensus FASTA, when available.",
    "- <project>.csv: sample information table used to define barcode/sample grouping for this delivery.",
    paste0("- ", unite_top_hits_name, ": root-level copy of the top UNITE BLAST hit per consensus sequence. This is the most convenient UNITE review table."),
    paste0("- ", unite_dir, "/consensus.blast.tsv: full consensus-vs-UNITE BLAST table with column headers and ONTools-derived coverage/taxonomy fields."),
    "- README.txt and README.zh-CN.txt: English and Chinese descriptions of this delivery.",
    "",
    "How To Interpret The Main Results",
    "Use the abundance table and figures for routine sample-level composition summaries. Use the per-barcode alignment-stat tables to review reference-level evidence from wf-16s. Use the UNITE BLAST top-hit table as an independent consensus-sequence review.",
    "The wf-16s ITS results are suitable for routine composition analysis and candidate taxon screening. Strict species confirmation or novel-species assessment should not rely on a single database top hit.",
    "For potential novel species, review consensus quality, percent identity, query coverage, reference coverage, top-hit versus second-hit separation, Species Hypothesis information when available, and whether the closest matches are well-curated or type-material records.",
    "",
    "Why wf-16s and UNITE Results May Differ",
    paste0("- ", abundance_table, " and ", figure_dir, "/ are generated from the wf-16s-style ITS workflow. In this delivery, that workflow commonly uses the database set configured by `database_set`, for example `ncbi_16s_18s_28s_ITS`. These results are read-level profiling summaries after the workflow's filtering, classification, and abundance aggregation steps."),
    paste0("- ", unite_dir, "/consensus.blast.tsv and ", unite_top_hits_name, " are generated by BLASTN of the final consensus sequences against a UNITE database from https://unite.ut.ee/repository.php, for example `unite_eukaryotes` built from a UNITE Species Hypothesis FASTA."),
    "- Because these analyses use different input units, databases, and decision rules, their taxonomic labels can differ. wf-16s summarizes many reads against its configured database; UNITE BLAST reviews one consensus sequence per barcode/sample against the selected UNITE database.",
    "- If the abundance table/figures and the UNITE consensus review disagree substantially, first check which database was used for each analysis, then compare coverage, identity, mapping quality, read counts, and whether the closest database records are well curated. The more reliable result is usually the one supported by high-quality consensus sequence, high coverage, high identity, clear separation from the second-best hit, and an appropriate curated database for the target organism group.",
    "",
    paste0("Table: ", abundance_table),
    "- tax: semicolon-separated taxonomic path. For this ITS delivery the expected order is usually superkingdom; kingdom; phylum; class; order; family; genus.",
    "- barcode/sample columns: read counts assigned to each taxonomic path for each barcode or sample.",
    "- total: total read count for that taxonomic path across all barcode/sample columns, when present.",
    "- Relative abundance within a sample can be calculated as: reads for one taxon in that sample / total classified and unclassified reads in that sample.",
    "- `Unclassified;Unknown;Unknown;...` means reads were not assigned to the reported taxonomic level in the abundance aggregation step. Common reasons include insufficient identity, insufficient reference coverage, off-target sequences, low read quality, chimeric reads, or incomplete database representation.",
    "",
    "Table: barcode*-alignment-stats.tsv",
    "These tables are copied from wf-16s `alignment_tables/`. They are useful for manual review of candidate database references, but they are not the same as the abundance table and should not be interpreted as final per-read abundance calls.",
    "- reference: database reference sequence identifier.",
    "- startpos: first covered/reference start position reported for this hit.",
    "- ref length: length of the database reference sequence.",
    "- number of reads: number of reads aligned to this reference record.",
    "- covbases: number of reference bases covered by aligned reads.",
    "- % coverage: percentage of the reference covered by aligned reads, calculated from covbases / ref length.",
    "- meandepth: mean read depth across the reference sequence.",
    "- meanbaseq: average base quality of reads contributing to this reference hit.",
    "- meanmapq: average mapping quality. Low values can indicate ambiguous mapping among similar references.",
    "- taxid: NCBI taxonomy identifier assigned to the reference record, when available.",
    "- superkingdom, kingdom, phylum, class, order, family, genus, species: taxonomy attached to the reference record.",
    "- genuspecies: older/combined genus-species field in some wf-16s outputs. Prefer separate `genus` and `species` columns when present.",
    "- mean, sd, Coefficient of Variance: summary statistics of depth distribution across the reference. A high coefficient of variation may indicate uneven coverage.",
    "- pcreads: percentage of reads in the barcode/sample represented by this reference hit.",
    "- Important limitation: this table does not report per-read percent identity. A named reference hit with low coverage, low mapping quality, or few reads should be treated as weak evidence.",
    "",
    paste0("Table: ", unite_top_hits_name),
    "This table contains the rank-1 UNITE BLAST hit for each consensus sequence. It is intended as an independent review of the final consensus sequence, not as a replacement for read-level abundance profiling. Original BLAST key metrics are retained in the table so that the suggested annotation can be audited.",
    "- qseqid: query consensus sequence ID.",
    "- qlen: query consensus length.",
    "- qstart, qend: aligned start and end positions on the query.",
    "- sseqid: matched UNITE/database sequence ID.",
    "- slen: matched reference sequence length.",
    "- sstart, send: aligned start and end positions on the matched reference.",
    "- length: BLAST alignment length.",
    "- pident: percent identity across the aligned region.",
    "- qcovs: BLAST-reported query coverage percentage.",
    "- mismatch: number of mismatches in the alignment.",
    "- gapopen: number of gap openings.",
    "- evalue: BLAST E-value; smaller values indicate stronger sequence similarity.",
    "- bitscore: BLAST bit score; larger values indicate stronger alignments.",
    "- salltitles: full matched reference title. For UNITE, this often contains Species Hypothesis and taxonomy information.",
    "- staxids: taxonomy ID from BLAST, when present.",
    "- query_coverage: coverage of the consensus sequence calculated from qstart/qend/qlen.",
    "- reference_coverage: coverage of the matched reference calculated from sstart/send/slen.",
    "- taxonomy_path, kingdom, phylum, class, order, family, genus, species: taxonomy parsed from UNITE-style headers when available.",
    "- rank: hit rank within each query after sorting by bitscore, E-value, identity, and coverage. In this root table rank is normally 1.",
    "- annotation_level: conservative label assigned by ONTools from identity and coverage thresholds: species, genus, family, or low_confidence.",
    "- novel_candidate: TRUE only for the rank-1 hit when pident is below the configured novel-candidate identity threshold and coverage is sufficient. With the default parameters, this means pident < 97, query_coverage >= 80, and reference_coverage >= 50. This is a screening flag, not a formal new-species conclusion.",
    "",
    "Suggested Review Order",
    paste0("1. Start with ", abundance_table, " and ", figure_dir, "/ to understand the overall sample composition."),
    "2. For each important sample, inspect the corresponding barcode alignment-stat table to see which reference records support the call.",
    paste0("3. If consensus sequences were generated, review ", unite_top_hits_name, " for an independent UNITE-based consensus check."),
    "4. For possible novel species or unexpected taxa, inspect the full BLAST hits, the consensus sequence, alignment coverage, read depth, and database metadata rather than relying only on the best hit.",
    "",
    if (isTRUE(has_unite_annotation)) {
      "UNITE annotation was generated for this delivery."
    } else {
      "UNITE annotation was not generated for this delivery because no UNITE database was supplied or no consensus FASTA was found."
    }
  )
}

ITS_delivery_readme_zh <- function(samples_dir,
                                   ITS_dir,
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
    "本交付目录整合三类 ITS 结果：",
    "1. 每个 barcode/样本的扩增子共识序列结果。",
    "2. 基于拆分后 FASTQ reads 生成的 wf-16s 风格 ITS 分类丰度结果。",
    "3. 如果提供了 UNITE 数据库，则额外提供基于项目共识序列的 consensus vs UNITE BLAST 复核结果。UNITE 参考数据可从 https://unite.ut.ee/repository.php 获取。",
    "",
    "这三类结果回答的问题不同。丰度表和图片用于概括 read 层面的样本组成；每个 barcode 的 alignment-stats 表用于查看 reads 支持了哪些数据库参考序列；UNITE BLAST 表用于对最终共识序列进行独立复核。",
    "",
    "目录内容",
    paste0("- ", sample_prefix, "consensus_results/：从扩增子共识交付结果复制来的每个 barcode 的共识序列分析结果。"),
    paste0("- ", sample_prefix, "consensus_results/consensus/consensus.fastq：该 barcode/样本的共识序列文件，如果流程成功生成则会存在。"),
    paste0("- ", sample_prefix, "consensus_results/alignments/：reads 回比对到 consensus 的 BAM、BAM index 以及 IGV 截图，如果对应步骤已运行则会存在。"),
    paste0("- ", sample_prefix, "consensus_results/*.synthetic.ab1：根据 reads pileup 合成的类 Sanger AB1 峰图文件，如果启用了 AB1 生成则会存在。"),
    paste0("- ", sample_prefix, "consensus_results/Distribution_seqLength__*.png：过滤前和/或过滤后的 reads 长度分布图。"),
    paste0("- ", sample_prefix, ITS_dir, "/：每个 barcode 的 ITS 分析证据。"),
    paste0("- ", sample_prefix, ITS_dir, "/", identification_dir, "/：每个 barcode 的 wf-16s 比对统计表。每个 `barcode*-alignment-stats.tsv` 文件汇总一个 barcode 的数据库参考序列命中情况。"),
    paste0("- ", abundance_table, "：wf-16s ITS 分析得到的属水平丰度表，是查看样本组成的主表。"),
    paste0("- ", figure_dir, "/：基于丰度表生成的丰度柱状图。`_percentage.png` 表示相对丰度，`_count.png` 表示 reads 数。"),
    "- all-consensus-seqs.fasta：用于 UNITE BLAST 复核的项目级共识序列 FASTA。如果存在 trimmed consensus，默认优先使用 trimmed 序列并同时保存为该文件名。",
    "- all-consensus-seqs.fasta.fai：all-consensus-seqs.fasta 的 FASTA index，如果可用则会存在。",
    "- all-consensus-seqs_trimmed.fasta：如果共识结果中存在 trimmed consensus FASTA，则按原文件名额外保留一份。",
    "- all-consensus-seqs_trimmed.fasta.fai：trimmed consensus FASTA 的 FASTA index，如果可用则会存在。",
    "- <project>.csv：本项目对应的样本信息表，用于定义 barcode/样本分组。",
    paste0("- ", unite_top_hits_name, "：放在根目录下的 UNITE top-hit 汇总表。每条 consensus 序列保留 rank 1 的最佳命中，是最方便查看 UNITE 复核结果的表。"),
    paste0("- ", unite_dir, "/consensus.blast.tsv：带有表头的完整 consensus vs UNITE BLAST 表格结果，并包含 ONTools 计算得到的 coverage 和分类解析字段。"),
    "- README.txt 和 README.zh-CN.txt：英文和中文交付说明。",
    "",
    "主要结果如何解读",
    "丰度表和图片适合用于样本整体组成展示；每个 barcode 的比对统计表适合复核 wf-16s 的参考序列命中证据；UNITE BLAST top-hit 表适合作为共识序列层面的独立复核结果。",
    "wf-16s ITS 结果适合常规组成分析和候选分类筛查，但严格种水平确认或新物种判断不应只依赖单个数据库 best hit。",
    "如果需要评估潜在新物种，建议重点查看共识序列质量、percent identity、query coverage、reference coverage、最佳命中与次佳命中的差距、Species Hypothesis 信息，以及最近命中是否来自高质量 curated record 或 type material。",
    "",
    "为什么 wf-16s 结果和 UNITE 复核结果可能不一致",
    paste0("- ", abundance_table, " 和 ", figure_dir, "/ 来自 wf-16s 风格 ITS 流程。该流程使用 `database_set` 指定的数据库，例如 `ncbi_16s_18s_28s_ITS`。这些结果是基于 reads 的分类、过滤和丰度汇总结果。"),
    paste0("- ", unite_dir, "/consensus.blast.tsv 和 ", unite_top_hits_name, " 来自最终 consensus 序列与 UNITE 数据库之间的 BLASTN 比对。UNITE 数据库来源为 https://unite.ut.ee/repository.php，例如可用 UNITE Species Hypothesis FASTA 构建 `unite_eukaryotes` 数据库。"),
    "- 两者使用的输入对象、数据库和判定规则都不同，因此分类结果可能存在差异。wf-16s 是对大量 reads 进行分类和丰度汇总；UNITE BLAST 是对每个 barcode/样本的一条或多条 consensus 序列与所选 UNITE 数据库进行复核。",
    "- 如果丰度表/图片与 UNITE consensus 复核结果差别很大，建议先确认两个分析分别使用了哪个数据库，再比较 coverage、identity、mapping quality、reads 数、最佳命中和次佳命中的差距，以及最近命中的数据库记录是否为高质量 curated record。通常更可靠的结果应同时具备高质量 consensus、高覆盖度、高 identity、与次佳命中有清晰差距，并且使用了适合目标类群的 curated 数据库。",
    "",
    paste0("表格说明：", abundance_table),
    "- tax：分号分隔的分类路径。ITS 交付结果中通常为 superkingdom; kingdom; phylum; class; order; family; genus。",
    "- barcode/样本列：每个 barcode 或样本中，被分配到该分类路径的 reads 数。",
    "- total：如果存在该列，表示该分类路径在所有 barcode/样本中的 reads 合计数。",
    "- 样本内相对丰度可按以下方式计算：该样本中某分类 reads 数 / 该样本总 reads 数。",
    "- `Unclassified;Unknown;Unknown;...` 表示这些 reads 在丰度汇总阶段没有被注释到报告层级。常见原因包括 identity 不足、reference coverage 不足、非目标扩增、reads 质量较低、嵌合序列，或数据库中缺少足够接近的参考序列。",
    "",
    "表格说明：barcode*-alignment-stats.tsv",
    "这些表来自 wf-16s 的 `alignment_tables/`，适合人工复核候选数据库参考序列，但它不是最终丰度表，也不能直接等同于逐条 read 的分类结果。",
    "- reference：数据库参考序列 ID。",
    "- startpos：该命中在参考序列上的起始位置。",
    "- ref length：数据库参考序列长度。",
    "- number of reads：比对到该参考序列的 reads 数。",
    "- covbases：参考序列上被 reads 覆盖到的碱基数。",
    "- % coverage：参考序列被覆盖的比例，通常由 covbases / ref length 计算得到。",
    "- meandepth：参考序列范围内的平均测序深度。",
    "- meanbaseq：支持该参考命中的 reads 的平均碱基质量。",
    "- meanmapq：平均比对质量。较低的 meanmapq 可能提示 reads 在多个相似参考之间存在模糊比对。",
    "- taxid：该参考序列对应的 NCBI taxonomy ID，如果数据库中提供则会显示。",
    "- superkingdom、kingdom、phylum、class、order、family、genus、species：数据库参考序列携带的分类注释。",
    "- genuspecies：部分 wf-16s 输出中的旧式或合并属种字段。如果同时存在 genus 和 species，优先看分开的两列。",
    "- mean、sd、Coefficient of Variance：参考序列深度分布的统计量。Coefficient of Variance 越高，通常表示覆盖越不均一。",
    "- pcreads：该参考命中的 reads 占该 barcode/样本 reads 的比例。",
    "- 重要限制：该表不包含逐条 read 的 percent identity。对于 reads 数少、coverage 低、meanmapq 低的命中，即使有明确物种名，也应作为弱证据谨慎解释。",
    "",
    paste0("表格说明：", unite_top_hits_name),
    "该表为每条 consensus 序列保留 rank 1 的 UNITE BLAST 最佳命中，用于对最终共识序列进行独立复核，不能替代 read 层面的丰度分析。表中保留了原始 BLAST 的关键比对指标，便于复核该注释建议是否可靠。",
    "- qseqid：查询序列 ID，即 consensus 序列 ID。",
    "- qlen：查询 consensus 序列长度。",
    "- qstart、qend：比对区域在查询序列上的起止位置。",
    "- sseqid：命中的 UNITE/数据库参考序列 ID。",
    "- slen：命中参考序列长度。",
    "- sstart、send：比对区域在命中参考序列上的起止位置。",
    "- length：BLAST 比对长度。",
    "- pident：比对区域内的百分比一致性。",
    "- qcovs：BLAST 输出的查询序列覆盖比例。",
    "- mismatch：错配数量。",
    "- gapopen：gap open 数量。",
    "- evalue：BLAST E-value，越小表示序列相似性越显著。",
    "- bitscore：BLAST bit score，越大表示比对越强。",
    "- salltitles：命中参考序列的完整标题。对于 UNITE 数据库，通常包含 Species Hypothesis 和分类信息。",
    "- staxids：BLAST 给出的 taxonomy ID，如果可用则会显示。",
    "- query_coverage：ONTools 根据 qstart、qend 和 qlen 计算的 consensus 覆盖比例。",
    "- reference_coverage：ONTools 根据 sstart、send 和 slen 计算的参考序列覆盖比例。",
    "- taxonomy_path、kingdom、phylum、class、order、family、genus、species：从 UNITE 风格标题中解析出的分类信息，如果标题中存在则会显示。",
    "- rank：同一条查询序列内的命中排名，排序依据包括 bitscore、E-value、identity 和 coverage。根目录的 top-hit 表中通常为 1。",
    "- annotation_level：ONTools 根据 identity 和 coverage 阈值给出的保守注释层级，包括 species、genus、family 或 low_confidence。",
    "- novel_candidate：只有 rank 1 的最佳命中在覆盖度足够、但 pident 低于设定的新物种候选 identity 阈值时才为 TRUE。在默认参数下，判断标准为 pident < 97、query_coverage >= 80 且 reference_coverage >= 50。该列只能作为筛查提示，不能单独作为新物种结论。",
    "",
    "推荐查看顺序",
    paste0("1. 先看 ", abundance_table, " 和 ", figure_dir, "/，了解样本整体组成。"),
    "2. 对重点样本，查看对应 barcode 的 alignment-stats 表，确认哪些数据库参考序列支持该分类结果。",
    paste0("3. 如果生成了 consensus 序列，再查看 ", unite_top_hits_name, "，用 UNITE 数据库对共识序列进行独立复核。"),
    "4. 如果涉及潜在新物种或意外分类，需要进一步查看完整 BLAST hits、consensus 序列、coverage、depth、identity、top hit 与 second hit 差距以及数据库记录来源。",
    "",
    if (isTRUE(has_unite_annotation)) {
      "本次交付已生成 UNITE 共识序列复核结果。"
    } else {
      "本次交付未生成 UNITE 复核结果，原因可能是未提供 UNITE 数据库，或没有找到可用于复核的共识序列 FASTA。"
    }
  )
}
