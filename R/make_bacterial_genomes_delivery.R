#' Build a bacterial genomes delivery package
#'
#' `make_bacterial_genomes_delivery()` collects bacterial genome analysis
#' deliverables into a stable customer-facing directory layout. Existing sample
#' filenames are preserved; files are only sorted into `Bam/`, `QC/`,
#' `Sequence/`, and `Var/` directories. When variant tables are not already
#' present, the function can generate simple `*.var.xls` and `*.filt.var.xls`
#' tables from VCF files with [read_vcf()].
#'
#' @param input_dir Directory containing files to deliver. Files may already be
#'   organized in subdirectories; scanning is recursive. If `fastq` is supplied,
#'   this can be `NULL` and an intermediate delivery-input directory is prepared
#'   from workflow outputs.
#' @param fastq Optional demultiplexed FASTQ directory. When supplied, the
#'   function can first run [run_wf_bacterial_genomes()], collect workflow
#'   outputs, optionally map reads and generate synthetic AB1 files, then build
#'   the delivery package.
#' @param output Output root directory.
#' @param project Project identifier used to create the delivery directory and,
#'   when requested, the archive.
#' @param delivery_suffix Suffix appended to `project` for the final delivery
#'   directory and archive names. Defaults to `""`, so the delivery directory is
#'   `<output>/<project>/`. Use `"_delivery"` to keep the older naming style.
#' @param sample_sheet Optional sample sheet or sample metadata file to copy
#'   into `Metadata/` in the delivery directory.
#' @param overwrite Logical. If `TRUE`, replace an existing delivery directory
#'   and archive.
#' @param make_archive Logical. If `TRUE`, create a `.tar.gz` archive next to
#'   the delivery directory.
#' @param run_wf Logical. If `TRUE` and `fastq` is supplied, run
#'   [run_wf_bacterial_genomes()] before collecting outputs.
#' @param wf_out_dir,work_dir Output and work directories used by
#'   [run_wf_bacterial_genomes()] when `fastq` is supplied.
#' @param staging_dir Intermediate directory used to collect workflow outputs
#'   before final packaging when `fastq` is supplied.
#' @param wf_min_len,wf_max_len,threads,profile,run_plasmid_id,resume,workflow,nextflow,seqkit,gzip,extra_args
#'   Parameters passed to [run_wf_bacterial_genomes()] when `run_wf = TRUE`.
#' @param run_mapping Logical. If `TRUE`, run [map_reads_to_assembly()] for
#'   consensus FASTA files that can be matched to FASTQ files by basename.
#' @param run_variant_calling Logical. If `TRUE`, run `bcftools mpileup` and
#'   `bcftools call` for consensus/BAM pairs to generate VCF files before
#'   creating variant tables.
#' @param make_ab1 Logical. If `TRUE`, run [synthetic_ab1_from_bam()] for
#'   consensus/BAM pairs that can be matched by basename.
#' @param sequencing_summary Optional sequencing summary file passed to
#'   [plot_seqQC()] to generate read-length QC PNGs.
#' @param mapping_threads Thread count passed to [map_reads_to_assembly()].
#' @param minimap2,samtools,bcftools Command names or paths used by mapping,
#'   variant calling, and AB1 generation.
#' @param variant_conda_env Optional conda environment name used for the
#'   `samtools faidx` and `bcftools` commands in variant calling. For example,
#'   use `"variant_qc"` to run `conda run -n variant_qc bcftools ...`.
#' @param conda Conda executable name or path used when `variant_conda_env` is
#'   supplied.
#' @param make_variant_tables Logical. If `TRUE`, generate missing variant
#'   tables from VCF files found under `input_dir`.
#' @param min_variant_percent Minimum value used for the generated filtered
#'   variant table when `variant_filter_column` is present. Use `NULL` to write
#'   the unfiltered table as the filtered table. Defaults to `0.05`, matching
#'   the fractional `Freq` column generated for bacterial variant tables.
#' @param variant_filter_column Column name used for filtering generated variant
#'   tables. Defaults to `"Freq"`.
#' @param variant_vcf_pattern Regular expression used to find VCF files for
#'   generated variant tables.
#' @param readme_name,chinese_readme_name README filenames. Set
#'   `chinese_readme_name = NULL` to skip the Chinese README.
#' @param include_manifest Logical. If `TRUE`, write `manifest.tsv` with the
#'   files included in the delivery and their source paths. Defaults to `FALSE`
#'   for cleaner customer-facing delivery packages.
#' @param dry_run Logical. If `TRUE`, return the planned copies and generated
#'   files without writing anything.
#' @param echo Logical. If `TRUE`, print a short summary.
#'
#' @return Invisibly returns a list with `paths`, `plan`, `manifest`, `archive`,
#'   and `status`.
#'
#' @examples
#' input <- tempfile("bg-input-")
#' output <- tempfile("bg-output-")
#' dir.create(input)
#' dir.create(output)
#' writeLines(c(">consensus", "ACGT"), file.path(input, "sample.consensus.fasta"))
#' res <- make_bacterial_genomes_delivery(
#'   input_dir = input,
#'   output = output,
#'   project = "PROJECT001",
#'   dry_run = TRUE
#' )
#' res$paths$delivery_dir
#'
#' @export
make_bacterial_genomes_delivery <- function(input_dir = NULL,
                                            fastq = NULL,
                                            output,
                                            project,
                                            delivery_suffix = "",
                                            sample_sheet = NULL,
                                            overwrite = FALSE,
                                            make_archive = TRUE,
                                            run_wf = !is.null(fastq),
                                            wf_out_dir = file.path(output, paste0(project, "_wf_bacterial_genomes")),
                                            work_dir = file.path(output, paste0(project, "_work_wf_bacterial_genomes")),
                                            staging_dir = file.path(output, paste0(project, "_delivery_input")),
                                            wf_min_len = 3000,
                                            wf_max_len = 12000,
                                            threads = 22,
                                            profile = "standard",
                                            run_plasmid_id = TRUE,
                                            resume = TRUE,
                                            workflow = "epi2me-labs/wf-bacterial-genomes",
                                            nextflow = "nextflow",
                                            seqkit = "seqkit",
                                            gzip = "gzip",
                                            extra_args = "-offline",
                                            run_mapping = FALSE,
                                            run_variant_calling = FALSE,
                                            make_ab1 = FALSE,
                                            sequencing_summary = NULL,
                                            mapping_threads = threads,
                                            minimap2 = "minimap2",
                                            samtools = "samtools",
                                            bcftools = "bcftools",
                                            variant_conda_env = NULL,
                                            conda = "conda",
                                            make_variant_tables = TRUE,
                                            min_variant_percent = 0.05,
                                            variant_filter_column = "Freq",
                                            variant_vcf_pattern = "[.]vcf([.]gz)?$",
                                            readme_name = "README.txt",
                                            chinese_readme_name = "README.zh-CN.txt",
                                            include_manifest = FALSE,
                                            dry_run = FALSE,
                                            echo = TRUE) {
  if (is.null(input_dir) && is.null(fastq)) {
    stop("Supply either `input_dir` or `fastq`.", call. = FALSE)
  }
  if (!is.null(input_dir)) check_dir_arg(input_dir, "input_dir")
  if (!is.null(fastq)) check_dir_arg(fastq, "fastq")
  if (!is.null(sample_sheet)) check_file_arg(sample_sheet, "sample_sheet")
  check_scalar_character(output, "output")
  check_scalar_character(project, "project")
  if (!is.character(delivery_suffix) || length(delivery_suffix) != 1L || is.na(delivery_suffix)) {
    stop("`delivery_suffix` must be a single character string.", call. = FALSE)
  }
  check_scalar_character(wf_out_dir, "wf_out_dir")
  check_scalar_character(work_dir, "work_dir")
  check_scalar_character(staging_dir, "staging_dir")
  check_scalar_character(profile, "profile")
  check_scalar_character(workflow, "workflow")
  check_scalar_character(nextflow, "nextflow")
  check_scalar_character(seqkit, "seqkit")
  check_scalar_character(gzip, "gzip")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(bcftools, "bcftools")
  check_scalar_character(conda, "conda")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(make_archive, "make_archive")
  check_logical_scalar(run_wf, "run_wf")
  check_logical_scalar(run_plasmid_id, "run_plasmid_id")
  check_logical_scalar(resume, "resume")
  check_logical_scalar(run_mapping, "run_mapping")
  check_logical_scalar(run_variant_calling, "run_variant_calling")
  check_logical_scalar(make_ab1, "make_ab1")
  check_logical_scalar(make_variant_tables, "make_variant_tables")
  check_scalar_character(variant_filter_column, "variant_filter_column")
  check_scalar_character(variant_vcf_pattern, "variant_vcf_pattern")
  check_scalar_character(readme_name, "readme_name")
  if (!is.null(chinese_readme_name)) {
    check_scalar_character(chinese_readme_name, "chinese_readme_name")
  }
  if (!is.null(variant_conda_env)) {
    check_scalar_character(variant_conda_env, "variant_conda_env")
  }
  check_logical_scalar(include_manifest, "include_manifest")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  if (!is.null(sequencing_summary)) check_file_arg(sequencing_summary, "sequencing_summary")
  if (!is.null(extra_args)) check_scalar_character(extra_args, "extra_args")
  wf_min_len <- validate_positive_integer(wf_min_len, "wf_min_len")
  wf_max_len <- validate_optional_positive_integer(wf_max_len, "wf_max_len")
  threads <- validate_positive_integer(threads, "threads")
  mapping_threads <- validate_positive_integer(mapping_threads, "mapping_threads")
  if (!is.null(min_variant_percent)) {
    min_variant_percent <- validate_nonnegative_number(
      min_variant_percent,
      "min_variant_percent"
    )
  }

  if (!is.null(input_dir)) input_dir <- normalizePath(input_dir, mustWork = TRUE)
  if (!is.null(fastq)) fastq <- normalizePath(fastq, mustWork = TRUE)
  if (!is.null(sequencing_summary)) {
    sequencing_summary <- normalizePath(sequencing_summary, mustWork = TRUE)
  }
  if (!is.null(sample_sheet)) {
    sample_sheet <- normalizePath(sample_sheet, mustWork = TRUE)
  }
  output <- normalizePath(output, mustWork = FALSE)
  wf_out_dir <- normalizePath(wf_out_dir, mustWork = FALSE)
  work_dir <- normalizePath(work_dir, mustWork = FALSE)
  staging_dir <- normalizePath(staging_dir, mustWork = FALSE)
  delivery_name <- paste0(project, delivery_suffix)
  delivery_dir <- file.path(output, delivery_name)
  archive <- file.path(output, paste0(delivery_name, ".tar.gz"))
  paths <- list(
    input_dir = input_dir,
    fastq = fastq,
    sample_sheet = sample_sheet,
    wf_out_dir = wf_out_dir,
    work_dir = work_dir,
    staging_dir = if (!is.null(fastq)) staging_dir else NULL,
    output_root = output,
    delivery_dir = delivery_dir,
    archive = if (isTRUE(make_archive)) archive else NULL
  )

  wf_result <- NULL
  collection <- NULL
  input_for_packaging <- input_dir

  if (!is.null(fastq)) {
    if (isTRUE(run_wf)) {
      wf_result <- run_wf_bacterial_genomes(
        fastq = fastq,
        out_dir = wf_out_dir,
        work_dir = work_dir,
        threads = threads,
        profile = profile,
        run_plasmid_id = run_plasmid_id,
        resume = resume,
        min_len = wf_min_len,
        max_len = wf_max_len,
        workflow = workflow,
        nextflow = nextflow,
        seqkit = seqkit,
        gzip = gzip,
        extra_args = extra_args,
        dry_run = dry_run,
        echo = echo
      )
    }

    if (isTRUE(dry_run) && !dir.exists(wf_out_dir)) {
      input_for_packaging <- if (!is.null(input_dir)) input_dir else staging_dir
    } else {
      collection <- collect_bacterial_genomes_workflow_outputs(
        wf_out_dir = wf_out_dir,
        fastq = fastq,
        staging_dir = staging_dir,
        sequencing_summary = sequencing_summary,
        run_mapping = run_mapping,
        run_variant_calling = run_variant_calling,
        make_ab1 = make_ab1,
        mapping_threads = mapping_threads,
        minimap2 = minimap2,
        samtools = samtools,
        bcftools = bcftools,
        variant_conda_env = variant_conda_env,
        conda = conda,
        overwrite = overwrite,
        dry_run = dry_run,
        echo = echo
      )
      input_for_packaging <- staging_dir
    }
  }

  files <- if (dir.exists(input_for_packaging)) {
    list.files(input_for_packaging, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  } else {
    character()
  }
  files <- files[file.exists(files) & !dir.exists(files)]
  copy_plan <- bacterial_delivery_copy_plan(files, delivery_dir)
  generated_plan <- bacterial_delivery_variant_generation_plan(
    files = files,
    delivery_dir = delivery_dir,
    copy_plan = copy_plan,
    enabled = make_variant_tables,
    variant_vcf_pattern = variant_vcf_pattern
  )
  plan <- rbind(copy_plan, generated_plan)
  check_bacterial_delivery_duplicate_destinations(plan)

  if (isTRUE(dry_run)) {
    if (isTRUE(echo)) {
      message("Planned bacterial genomes delivery: ", delivery_dir)
      message("Files to copy: ", sum(plan$action == "copy"))
      message("Variant tables to generate: ", sum(plan$action == "generate_variant"))
    }
    return(invisible(list(
      paths = paths,
      wf = wf_result,
      collection = collection,
      plan = plan,
      manifest = NULL,
      archive = if (isTRUE(make_archive)) archive else NULL,
      status = NA_integer_
    )))
  }

  if (dir.exists(delivery_dir) || file.exists(archive)) {
    if (!isTRUE(overwrite)) {
      stop(
        "Delivery output already exists. Use `overwrite = TRUE` to replace: ",
        delivery_dir,
        call. = FALSE
      )
    }
    unlink(delivery_dir, recursive = TRUE, force = TRUE)
    unlink(archive, force = TRUE)
  }

  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(delivery_dir, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c("Bam", "QC", "Sequence", "Var", "md5")) {
    dir.create(file.path(delivery_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(sample_sheet)) {
    dir.create(file.path(delivery_dir, "Metadata"), recursive = TRUE, showWarnings = FALSE)
  }

  manifest_rows <- bacterial_delivery_execute_copy_plan(plan)
  variant_rows <- bacterial_delivery_execute_variant_plan(
    generated_plan,
    min_variant_percent = min_variant_percent,
    variant_filter_column = variant_filter_column
  )
  manifest_rows <- rbind(manifest_rows, variant_rows)

  merged_rows <- bacterial_delivery_write_merged_variant_tables(delivery_dir)
  manifest_rows <- rbind(manifest_rows, merged_rows)

  readme_rows <- bacterial_delivery_write_readmes(
    delivery_dir = delivery_dir,
    project = project,
    readme_name = readme_name,
    chinese_readme_name = chinese_readme_name
  )
  manifest_rows <- rbind(manifest_rows, readme_rows)

  metadata_rows <- bacterial_delivery_copy_sample_sheet(
    sample_sheet = sample_sheet,
    delivery_dir = delivery_dir
  )
  manifest_rows <- rbind(manifest_rows, metadata_rows)

  manifest_path <- NULL
  if (isTRUE(include_manifest)) {
    manifest_path <- file.path(delivery_dir, "manifest.tsv")
    utils::write.table(
      manifest_rows,
      manifest_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  bacterial_delivery_write_md5(delivery_dir)

  archive_status <- NA_integer_
  if (isTRUE(make_archive)) {
    archive_status <- bacterial_delivery_make_archive(delivery_dir, archive)
  }

  if (isTRUE(echo)) {
    message("Bacterial genomes delivery written: ", normalizePath(delivery_dir, mustWork = TRUE))
    if (isTRUE(make_archive)) {
      message("Archive written: ", normalizePath(archive, mustWork = TRUE))
    }
  }

  invisible(list(
    paths = list(
      input_dir = input_for_packaging,
      fastq = fastq,
      sample_sheet = sample_sheet,
      wf_out_dir = wf_out_dir,
      work_dir = work_dir,
      staging_dir = if (!is.null(fastq)) staging_dir else NULL,
      output_root = normalizePath(output, mustWork = TRUE),
      delivery_dir = normalizePath(delivery_dir, mustWork = TRUE),
      archive = if (isTRUE(make_archive)) normalizePath(archive, mustWork = TRUE) else NULL
    ),
    wf = wf_result,
    collection = collection,
    plan = plan,
    manifest = manifest_path,
    archive = if (isTRUE(make_archive)) archive else NULL,
    status = if (isTRUE(make_archive)) archive_status else 0L
  ))
}

bacterial_delivery_copy_plan <- function(files, delivery_dir) {
  rows <- lapply(files, function(path) {
    category <- bacterial_delivery_category(path)
    if (is.null(category)) {
      return(NULL)
    }
    destination <- if (identical(category, "Root")) {
      file.path(delivery_dir, basename(path))
    } else {
      file.path(delivery_dir, category, basename(path))
    }
    data.frame(
      action = "copy",
      category = category,
      source = path,
      destination = destination,
      label = bacterial_delivery_label(category, path),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(bacterial_delivery_empty_plan())
  }
  do.call(rbind, rows)
}

collect_bacterial_genomes_workflow_outputs <- function(wf_out_dir,
                                                       fastq,
                                                       staging_dir,
                                                       sequencing_summary,
                                                       run_mapping,
                                                       run_variant_calling,
                                                       make_ab1,
                                                       mapping_threads,
                                                       minimap2,
                                                       samtools,
                                                       bcftools,
                                                       variant_conda_env,
                                                       conda,
                                                       overwrite,
                                                       dry_run,
                                                       echo) {
  if (!dir.exists(wf_out_dir)) {
    stop("`wf_out_dir` does not exist: ", wf_out_dir, call. = FALSE)
  }

  wf_out_dir <- normalizePath(wf_out_dir, mustWork = TRUE)
  fastq <- normalizePath(fastq, mustWork = TRUE)
  plan <- bacterial_genomes_workflow_collect_plan(wf_out_dir, staging_dir)
  if (isTRUE(dry_run)) {
    return(list(
      staging_dir = staging_dir,
      plan = plan,
      mapping = list(),
      variants = list(),
      ab1 = list(),
      qc = NULL,
      status = NA_integer_
    ))
  }

  if (dir.exists(staging_dir) && isTRUE(overwrite)) {
    unlink(staging_dir, recursive = TRUE, force = TRUE)
  }
  for (subdir in c("Bam", "QC", "Sequence", "Var")) {
    dir.create(file.path(staging_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }

  manifest <- bacterial_genomes_execute_collect_plan(plan)

  qc <- NULL
  if (!is.null(sequencing_summary)) {
    qc <- plot_seqQC(
      sequencing_summary,
      runName = basename(staging_dir),
      device = "png",
      out_dir = file.path(staging_dir, "QC"),
      min_read_length = NULL,
      max_read_length = NULL,
      barcode_digits = 3
    )
  }

  mapping <- list()
  if (isTRUE(run_mapping)) {
    mapping <- bacterial_genomes_generate_missing_bams(
      staging_dir = staging_dir,
      fastq = fastq,
      threads = mapping_threads,
      minimap2 = minimap2,
      samtools = samtools,
      echo = echo
    )
    coverage <- bacterial_genomes_generate_missing_coverage_plots(
      staging_dir = staging_dir,
      samtools = samtools
    )
    mapping$coverage <- coverage
  }

  variants <- list()
  if (isTRUE(run_variant_calling)) {
    variants <- bacterial_genomes_call_variants(
      staging_dir = staging_dir,
      samtools = samtools,
      bcftools = bcftools,
      conda_env = variant_conda_env,
      conda = conda,
      echo = echo
    )
  }

  ab1 <- list()
  if (isTRUE(make_ab1)) {
    ab1 <- bacterial_genomes_generate_missing_ab1(
      staging_dir = staging_dir,
      samtools = samtools,
      echo = echo
    )
  }

  list(
    staging_dir = normalizePath(staging_dir, mustWork = TRUE),
    plan = plan,
    manifest = manifest,
    mapping = mapping,
    variants = variants,
    ab1 = ab1,
    qc = qc,
    status = 0L
  )
}

bacterial_genomes_workflow_collect_plan <- function(wf_out_dir, staging_dir) {
  files <- list.files(wf_out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.exists(files) & !dir.exists(files)]
  rows <- list()

  consensus <- bacterial_genomes_consensus_files(files)
  for (path in consensus) {
    stem <- bacterial_genomes_consensus_stem(path)
    rows[[length(rows) + 1L]] <- data.frame(
      action = if (grepl("[.]gz$", path, ignore.case = TRUE)) "gunzip_copy" else "copy",
      category = "Sequence",
      source = path,
      destination = file.path(staging_dir, "Sequence", paste0(stem, ".consensus.fasta")),
      label = "consensus_sequence",
      stringsAsFactors = FALSE
    )
  }

  bam_files <- files[grepl("[.]bam([.]bai)?$", basename(files), ignore.case = TRUE) |
                       grepl("[.]bai$", basename(files), ignore.case = TRUE)]
  for (path in bam_files) {
    rows[[length(rows) + 1L]] <- data.frame(
      action = "copy",
      category = "Bam",
      source = path,
      destination = file.path(staging_dir, "Bam", basename(path)),
      label = "bam",
      stringsAsFactors = FALSE
    )
  }

  png_files <- files[grepl("[.]png$", basename(files), ignore.case = TRUE)]
  for (path in png_files) {
    rows[[length(rows) + 1L]] <- data.frame(
      action = "copy",
      category = "QC",
      source = path,
      destination = file.path(staging_dir, "QC", basename(path)),
      label = "qc_plot",
      stringsAsFactors = FALSE
    )
  }

  vcf_files <- files[grepl("[.]vcf([.]gz)?$", basename(files), ignore.case = TRUE)]
  for (path in vcf_files) {
    rows[[length(rows) + 1L]] <- data.frame(
      action = "copy",
      category = "Var",
      source = path,
      destination = file.path(staging_dir, "Var", basename(path)),
      label = "variant_vcf",
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) return(bacterial_delivery_empty_plan())
  plan <- do.call(rbind, rows)
  check_bacterial_delivery_duplicate_destinations(plan)
  plan
}

bacterial_genomes_execute_collect_plan <- function(plan) {
  rows <- bacterial_delivery_empty_manifest()
  for (i in seq_len(nrow(plan))) {
    dir.create(dirname(plan$destination[[i]]), recursive = TRUE, showWarnings = FALSE)
    if (identical(plan$action[[i]], "gunzip_copy")) {
      input <- gzfile(plan$source[[i]], open = "rt")
      output <- file(plan$destination[[i]], open = "wt")
      repeat {
        chunk <- readLines(input, n = 10000)
        if (length(chunk) == 0L) break
        writeLines(chunk, output)
      }
      close(input)
      close(output)
    } else {
      ok <- file.copy(plan$source[[i]], plan$destination[[i]], overwrite = TRUE)
      if (!isTRUE(ok)) {
        stop("Failed to copy workflow output: ", plan$source[[i]], call. = FALSE)
      }
    }
    rows <- rbind(rows, bacterial_delivery_manifest_row(
      label = plan$label[[i]],
      source = plan$source[[i]],
      destination = plan$destination[[i]],
      status = if (identical(plan$action[[i]], "gunzip_copy")) "decompressed" else "copied"
    ))
  }
  rows
}

bacterial_genomes_consensus_files <- function(files) {
  preferred <- files[grepl("[.]medaka[.]fasta([.]gz)?$", basename(files), ignore.case = TRUE)]
  if (length(preferred) > 0L) return(preferred)
  files[grepl("[.](fa|fasta|fna)([.]gz)?$", basename(files), ignore.case = TRUE)]
}

bacterial_genomes_consensus_stem <- function(path) {
  stem <- basename(path)
  stem <- sub("[.]gz$", "", stem, ignore.case = TRUE)
  stem <- sub("[.]medaka[.]fasta$", "", stem, ignore.case = TRUE)
  stem <- sub("[.](fa|fasta|fna)$", "", stem, ignore.case = TRUE)
  stem
}

bacterial_genomes_fastq_files <- function(fastq) {
  files <- list.files(fastq, recursive = TRUE, full.names = TRUE)
  files[file.exists(files) & grepl("[.]fastq([.]gz)?$", basename(files), ignore.case = TRUE)]
}

bacterial_genomes_match_by_stem <- function(query_stem, candidates) {
  if (length(candidates) == 0L) return(NA_character_)
  stems <- basename(candidates)
  stems <- sub("[.]gz$", "", stems, ignore.case = TRUE)
  stems <- sub("[.]fastq$", "", stems, ignore.case = TRUE)
  exact <- match(query_stem, stems)
  if (!is.na(exact)) return(candidates[[exact]])
  starts <- which(startsWith(query_stem, stems) | startsWith(stems, query_stem))
  if (length(starts) == 1L) return(candidates[[starts]])
  if (length(candidates) == 1L) return(candidates[[1L]])
  NA_character_
}

bacterial_genomes_generate_missing_bams <- function(staging_dir,
                                                    fastq,
                                                    threads,
                                                    minimap2,
                                                    samtools,
                                                    echo) {
  consensus_files <- list.files(file.path(staging_dir, "Sequence"),
                                pattern = "[.]consensus[.]fasta$",
                                full.names = TRUE)
  fastq_files <- bacterial_genomes_fastq_files(fastq)
  results <- list()
  for (consensus in consensus_files) {
    stem <- sub("[.]consensus[.]fasta$", "", basename(consensus), ignore.case = TRUE)
    bam <- file.path(staging_dir, "Bam", paste0(stem, ".bam"))
    if (file.exists(bam)) next
    reads <- bacterial_genomes_match_by_stem(stem, fastq_files)
    if (is.na(reads)) {
      warning("No matching FASTQ found for consensus: ", basename(consensus), call. = FALSE)
      next
    }
    results[[stem]] <- map_reads_to_assembly(
      assembly = consensus,
      reads = reads,
      align_bam = bam,
      depth_file = file.path(staging_dir, "Bam", paste0(stem, ".depth.txt")),
      depth_plot = file.path(staging_dir, "QC", paste0(stem, ".coverage.png")),
      threads = threads,
      minimap2 = minimap2,
      samtools = samtools,
      echo = echo
    )
  }
  results
}

bacterial_genomes_generate_missing_ab1 <- function(staging_dir, samtools, echo) {
  consensus_files <- list.files(file.path(staging_dir, "Sequence"),
                                pattern = "[.]consensus[.]fasta$",
                                full.names = TRUE)
  results <- list()
  for (consensus in consensus_files) {
    stem <- sub("[.]consensus[.]fasta$", "", basename(consensus), ignore.case = TRUE)
    ab1 <- file.path(staging_dir, "Sequence", paste0(stem, ".ab1"))
    if (file.exists(ab1)) next
    bam <- file.path(staging_dir, "Bam", paste0(stem, ".bam"))
    if (!file.exists(bam)) {
      warning("No matching BAM found for AB1 generation: ", basename(consensus), call. = FALSE)
      next
    }
    results[[stem]] <- synthetic_ab1_from_bam(
      consensus = consensus,
      bam = bam,
      output_ab1 = ab1,
      sample = stem,
      samtools = samtools,
      overwrite = TRUE,
      echo = echo
    )
  }
  results
}

bacterial_genomes_generate_missing_coverage_plots <- function(staging_dir, samtools) {
  consensus_files <- list.files(file.path(staging_dir, "Sequence"),
                                pattern = "[.]consensus[.]fasta$",
                                full.names = TRUE)
  results <- list()
  for (consensus in consensus_files) {
    stem <- sub("[.]consensus[.]fasta$", "", basename(consensus), ignore.case = TRUE)
    bam <- file.path(staging_dir, "Bam", paste0(stem, ".bam"))
    if (!file.exists(bam)) next

    depth_file <- file.path(staging_dir, "Bam", paste0(stem, ".depth.txt"))
    coverage_png <- file.path(staging_dir, "QC", paste0(stem, ".coverage.png"))
    if (file.exists(coverage_png)) next

    dir.create(dirname(depth_file), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(coverage_png), recursive = TRUE, showWarnings = FALSE)

    status <- system2(
      samtools,
      args = c("depth", "-aa", bam),
      stdout = depth_file,
      stderr = ""
    )
    if (!identical(status, 0L)) {
      warning("samtools depth failed for BAM: ", bam, call. = FALSE)
      next
    }
    if (!file.exists(depth_file) || file.info(depth_file)$size == 0L) {
      warning("Depth file is empty for BAM: ", bam, call. = FALSE)
      next
    }

    results[[stem]] <- list(
      status = status,
      paths = list(
        consensus = consensus,
        bam = bam,
        depth = depth_file,
        coverage_plot = coverage_png
      ),
      plot = plot_read_depth(
        depth_file = depth_file,
        depth_plot = coverage_png,
        width = 8,
        height = 5,
        facet_nrow = 3
      )
    )
  }
  results
}

bacterial_genomes_call_variants <- function(staging_dir,
                                            samtools,
                                            bcftools,
                                            conda_env,
                                            conda,
                                            echo) {
  consensus_files <- list.files(file.path(staging_dir, "Sequence"),
                                pattern = "[.]consensus[.]fasta$",
                                full.names = TRUE)
  results <- list()
  for (consensus in consensus_files) {
    stem <- sub("[.]consensus[.]fasta$", "", basename(consensus), ignore.case = TRUE)
    bam <- file.path(staging_dir, "Bam", paste0(stem, ".bam"))
    if (!file.exists(bam)) {
      warning("No matching BAM found for variant calling: ", basename(consensus), call. = FALSE)
      next
    }

    vcf <- file.path(staging_dir, "Var", paste0(stem, ".vcf.gz"))
    vcf_index <- paste0(vcf, ".csi")
    if (file.exists(vcf)) {
      results[[stem]] <- list(
        status = 0L,
        skipped = TRUE,
        paths = list(consensus = consensus, bam = bam, vcf = vcf, index = vcf_index)
      )
      next
    }

    dir.create(dirname(vcf), recursive = TRUE, showWarnings = FALSE)
    fai <- paste0(consensus, ".fai")
    if (!file.exists(fai)) {
      faidx_call <- dehost_fastq_external_call(
        command = samtools,
        args = c("faidx", consensus),
        conda_env = conda_env,
        conda = conda
      )
      if (isTRUE(echo)) {
        message(paste(c(shQuote(faidx_call$command), shQuote(faidx_call$args)), collapse = " "))
      }
      faidx_status <- system2(faidx_call$command, args = faidx_call$args, stderr = "")
      if (!identical(faidx_status, 0L)) {
        warning("samtools faidx failed for consensus: ", consensus, call. = FALSE)
        next
      }
    }

    tmp_bcf <- tempfile("mpileup_", fileext = ".bcf")
    on.exit(unlink(tmp_bcf), add = TRUE)
    mpileup_args <- c("mpileup", "-Ou", "-f", consensus, "-o", tmp_bcf, bam)
    call_args <- c("call", "-mv", "-Oz", "-o", vcf, tmp_bcf)
    index_args <- c("index", vcf)
    mpileup_call <- dehost_fastq_external_call(
      command = bcftools,
      args = mpileup_args,
      conda_env = conda_env,
      conda = conda
    )
    call_call <- dehost_fastq_external_call(
      command = bcftools,
      args = call_args,
      conda_env = conda_env,
      conda = conda
    )
    index_call <- dehost_fastq_external_call(
      command = bcftools,
      args = index_args,
      conda_env = conda_env,
      conda = conda
    )
    if (isTRUE(echo)) {
      message(paste(c(shQuote(mpileup_call$command), shQuote(mpileup_call$args)), collapse = " "))
      message(paste(c(shQuote(call_call$command), shQuote(call_call$args)), collapse = " "))
      message(paste(c(shQuote(index_call$command), shQuote(index_call$args)), collapse = " "))
    }

    mpileup_status <- system2(mpileup_call$command, args = mpileup_call$args, stderr = "")
    if (!identical(mpileup_status, 0L)) {
      warning("bcftools mpileup failed for BAM: ", bam, call. = FALSE)
      next
    }
    call_status <- system2(call_call$command, args = call_call$args, stderr = "")
    if (!identical(call_status, 0L)) {
      warning("bcftools call failed for BAM: ", bam, call. = FALSE)
      next
    }
    index_status <- system2(index_call$command, args = index_call$args, stderr = "")
    if (!identical(index_status, 0L)) {
      warning("bcftools index failed for VCF: ", vcf, call. = FALSE)
    }

    results[[stem]] <- list(
      status = if (identical(index_status, 0L)) 0L else index_status,
      skipped = FALSE,
      commands = list(
        mpileup = paste(c(shQuote(mpileup_call$command), shQuote(mpileup_call$args)), collapse = " "),
        call = paste(c(shQuote(call_call$command), shQuote(call_call$args)), collapse = " "),
        index = paste(c(shQuote(index_call$command), shQuote(index_call$args)), collapse = " ")
      ),
      paths = list(consensus = consensus, bam = bam, vcf = vcf, index = vcf_index)
    )
  }

  results
}

bacterial_delivery_variant_generation_plan <- function(files,
                                                       delivery_dir,
                                                       copy_plan,
                                                       enabled,
                                                       variant_vcf_pattern) {
  if (!isTRUE(enabled) || length(files) == 0L) {
    return(bacterial_delivery_empty_plan())
  }

  vcfs <- files[grepl(variant_vcf_pattern, basename(files), ignore.case = TRUE)]
  if (length(vcfs) == 0L) {
    return(bacterial_delivery_empty_plan())
  }

  existing_dest <- copy_plan$destination
  rows <- list()
  for (vcf in vcfs) {
    stem <- sub("[.]vcf([.]gz)?$", "", basename(vcf), ignore.case = TRUE)
    var_dest <- file.path(delivery_dir, "Var", paste0(stem, ".var.xls"))
    filt_dest <- file.path(delivery_dir, "Var", paste0(stem, ".filt.var.xls"))
    if (!var_dest %in% existing_dest) {
      rows[[length(rows) + 1L]] <- data.frame(
        action = "generate_variant",
        category = "Var",
        source = vcf,
        destination = var_dest,
        label = "generated_variant_table",
        stringsAsFactors = FALSE
      )
    }
    if (!filt_dest %in% existing_dest) {
      rows[[length(rows) + 1L]] <- data.frame(
        action = "generate_variant",
        category = "Var",
        source = vcf,
        destination = filt_dest,
        label = "generated_filtered_variant_table",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(bacterial_delivery_empty_plan())
  }
  do.call(rbind, rows)
}

bacterial_delivery_category <- function(path) {
  name <- basename(path)
  lower <- tolower(name)

  if (lower %in% c("merged_data.var.xls", "merged_data.filt.var.xls")) {
    return("Root")
  }
  if (grepl("[.]bam([.]bai)?$", lower) || grepl("[.]bai$", lower)) {
    return("Bam")
  }
  if (grepl("[.]png$", lower)) {
    return("QC")
  }
  if (grepl("[.]ab1$", lower) ||
      grepl("[.](fa|fasta|fna)([.]gz)?$", lower)) {
    return("Sequence")
  }
  if (grepl("[.]var[.]xls$", lower)) {
    return("Var")
  }

  NULL
}

bacterial_delivery_label <- function(category, path) {
  lower <- tolower(basename(path))
  if (identical(category, "Root")) return("merged_variant_table")
  if (identical(category, "Bam")) return("bam")
  if (identical(category, "QC")) return("qc_plot")
  if (identical(category, "Sequence")) {
    if (grepl("[.]ab1$", lower)) return("ab1")
    return("consensus_sequence")
  }
  if (identical(category, "Var")) {
    if (grepl("[.]filt[.]var[.]xls$", lower)) return("filtered_variant_table")
    return("variant_table")
  }
  "file"
}

bacterial_delivery_empty_plan <- function() {
  data.frame(
    action = character(),
    category = character(),
    source = character(),
    destination = character(),
    label = character(),
    stringsAsFactors = FALSE
  )
}

check_bacterial_delivery_duplicate_destinations <- function(plan) {
  if (nrow(plan) == 0L) return(invisible(TRUE))
  duplicated_dest <- unique(plan$destination[duplicated(plan$destination)])
  if (length(duplicated_dest) > 0L) {
    stop(
      "Multiple source files map to the same delivery filename: ",
      paste(duplicated_dest, collapse = ", "),
      ". Rename one of the source files or separate the inputs.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

bacterial_delivery_execute_copy_plan <- function(plan) {
  copy_plan <- plan[plan$action == "copy", , drop = FALSE]
  rows <- bacterial_delivery_empty_manifest()
  for (i in seq_len(nrow(copy_plan))) {
    dir.create(dirname(copy_plan$destination[[i]]), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(copy_plan$source[[i]], copy_plan$destination[[i]], overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to copy file: ", copy_plan$source[[i]], call. = FALSE)
    }
    rows <- rbind(rows, bacterial_delivery_manifest_row(
      label = copy_plan$label[[i]],
      source = copy_plan$source[[i]],
      destination = copy_plan$destination[[i]],
      status = "copied"
    ))
  }
  rows
}

bacterial_delivery_execute_variant_plan <- function(plan,
                                                    min_variant_percent,
                                                    variant_filter_column) {
  variant_plan <- plan[plan$action == "generate_variant", , drop = FALSE]
  rows <- bacterial_delivery_empty_manifest()
  if (nrow(variant_plan) == 0L) return(rows)

  for (vcf in unique(variant_plan$source)) {
    table <- bacterial_delivery_variant_table_from_vcf(vcf)
    filtered <- bacterial_delivery_filter_variant_table(
      table,
      min_variant_percent = min_variant_percent,
      variant_filter_column = variant_filter_column
    )

    vcf_plan <- variant_plan[variant_plan$source == vcf, , drop = FALSE]
    for (i in seq_len(nrow(vcf_plan))) {
      dir.create(dirname(vcf_plan$destination[[i]]), recursive = TRUE, showWarnings = FALSE)
      out <- if (identical(vcf_plan$label[[i]], "generated_filtered_variant_table")) {
        filtered
      } else {
        table
      }
      utils::write.table(
        out,
        vcf_plan$destination[[i]],
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      rows <- rbind(rows, bacterial_delivery_manifest_row(
        label = vcf_plan$label[[i]],
        source = vcf,
        destination = vcf_plan$destination[[i]],
        status = "generated"
      ))
    }
  }

  rows
}

bacterial_delivery_empty_variant_table <- function() {
  data.frame(
    Chr = character(),
    Pos = integer(),
    Ref = character(),
    Alt = character(),
    DP = integer(),
    Ref_dp = integer(),
    Alt_dp = integer(),
    Freq = numeric(),
    DP4 = character(),
    Seq = character(),
    stringsAsFactors = FALSE
  )
}

bacterial_delivery_variant_table_from_vcf <- function(vcf) {
  raw <- read_vcf(
    vcf,
    parse_info = TRUE,
    parse_genotypes = TRUE,
    add_variant_type = FALSE,
    add_allele_depth = TRUE,
    keep_info = TRUE,
    keep_format = TRUE
  )
  if (is.null(raw)) return(bacterial_delivery_empty_variant_table())

  info <- if ("INFO" %in% names(raw)) {
    parse_vcf_info(raw$INFO, simplify = FALSE)
  } else {
    data.frame(row.names = seq_len(nrow(raw)))
  }

  dp4 <- bacterial_delivery_variant_dp4(raw, info)
  depth <- bacterial_delivery_variant_depths(raw, info, dp4)
  seq_context <- bacterial_delivery_variant_sequence_context(raw)

  data.frame(
    Chr = as.character(raw$CHROM),
    Pos = as.integer(raw$POS),
    Ref = bacterial_delivery_variant_ref(raw$REF),
    Alt = bacterial_delivery_variant_alt(raw$ALT),
    DP = as.integer(depth$dp),
    Ref_dp = as.integer(depth$ref_dp),
    Alt_dp = as.integer(depth$alt_dp),
    Freq = bacterial_delivery_round_freq(depth$freq),
    DP4 = dp4$label,
    Seq = seq_context,
    stringsAsFactors = FALSE
  )
}

bacterial_delivery_variant_dp4 <- function(raw, info) {
  value <- if ("DP4" %in% names(info)) as.character(info$DP4) else rep(NA_character_, nrow(raw))
  parsed <- lapply(value, function(x) {
    if (is.na(x) || !nzchar(x)) return(rep(NA_integer_, 4L))
    nums <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1L]]))
    nums <- nums[seq_len(min(length(nums), 4L))]
    c(nums, rep(NA_integer_, 4L - length(nums)))
  })
  matrix <- do.call(rbind, parsed)
  list(
    value = matrix,
    label = ifelse(is.na(value) | !nzchar(value), "", paste0("DP4=", value))
  )
}

bacterial_delivery_variant_depths <- function(raw, info, dp4) {
  ref_dp <- rep(NA_real_, nrow(raw))
  alt_dp <- rep(NA_real_, nrow(raw))

  has_dp4 <- rowSums(is.na(dp4$value)) == 0L
  if (any(has_dp4)) {
    ref_dp[has_dp4] <- dp4$value[has_dp4, 1L] + dp4$value[has_dp4, 2L]
    alt_dp[has_dp4] <- dp4$value[has_dp4, 3L] + dp4$value[has_dp4, 4L]
  }

  if ("ref_depth_downsampled" %in% names(raw)) {
    missing_ref <- is.na(ref_dp)
    ref_dp[missing_ref] <- suppressWarnings(as.numeric(raw$ref_depth_downsampled[missing_ref]))
  }
  if ("alt_depth_downsampled" %in% names(raw)) {
    missing_alt <- is.na(alt_dp)
    alt_dp[missing_alt] <- suppressWarnings(as.numeric(raw$alt_depth_downsampled[missing_alt]))
  }

  dp <- rep(NA_real_, nrow(raw))
  if ("DP" %in% names(info)) {
    dp <- suppressWarnings(as.numeric(info$DP))
  }
  if ("DP" %in% names(raw)) {
    missing_dp <- is.na(dp)
    dp[missing_dp] <- suppressWarnings(as.numeric(raw$DP[missing_dp]))
  }
  missing_dp <- is.na(dp) & (!is.na(ref_dp) | !is.na(alt_dp))
  dp[missing_dp] <- rowSums(cbind(ref_dp[missing_dp], alt_dp[missing_dp]), na.rm = TRUE)

  freq <- rep(NA_real_, nrow(raw))
  denom <- ref_dp + alt_dp
  freq[!is.na(denom) & denom > 0] <- alt_dp[!is.na(denom) & denom > 0] / denom[!is.na(denom) & denom > 0]
  if ("variant_percent" %in% names(raw)) {
    missing_freq <- is.na(freq)
    freq[missing_freq] <- suppressWarnings(as.numeric(raw$variant_percent[missing_freq])) / 100
  }

  list(dp = dp, ref_dp = ref_dp, alt_dp = alt_dp, freq = freq)
}

bacterial_delivery_variant_ref <- function(ref) {
  ref <- as.character(ref)
  ref[is.na(ref) | ref == "."] <- "-"
  ref
}

bacterial_delivery_variant_alt <- function(alt) {
  alt <- as.character(alt)
  alt[is.na(alt) | alt == "."] <- "-"
  alt
}

bacterial_delivery_round_freq <- function(freq) {
  ifelse(is.na(freq), NA_real_, round(freq, 5))
}

bacterial_delivery_variant_sequence_context <- function(raw) {
  upstream_cols <- grep("^ref_upstream_", names(raw), value = TRUE)
  downstream_cols <- grep("^ref_downstream_", names(raw), value = TRUE)
  upstream <- if (length(upstream_cols) > 0L) as.character(raw[[upstream_cols[[1L]]]]) else rep("", nrow(raw))
  downstream <- if (length(downstream_cols) > 0L) as.character(raw[[downstream_cols[[1L]]]]) else rep("", nrow(raw))
  upstream[is.na(upstream)] <- ""
  downstream[is.na(downstream)] <- ""
  paste0(
    upstream,
    "[",
    bacterial_delivery_variant_ref(raw$REF),
    "/",
    bacterial_delivery_variant_alt(raw$ALT),
    "]",
    downstream
  )
}

bacterial_delivery_filter_variant_table <- function(table,
                                                    min_variant_percent,
                                                    variant_filter_column) {
  if (is.null(min_variant_percent) || !variant_filter_column %in% names(table)) {
    return(table)
  }
  value <- suppressWarnings(as.numeric(table[[variant_filter_column]]))
  table[!is.na(value) & value >= min_variant_percent, , drop = FALSE]
}

bacterial_delivery_write_merged_variant_tables <- function(delivery_dir) {
  rows <- bacterial_delivery_empty_manifest()
  var_dir <- file.path(delivery_dir, "Var")
  if (!dir.exists(var_dir)) return(rows)

  merged_var <- file.path(delivery_dir, "merged_data.var.xls")
  merged_filt <- file.path(delivery_dir, "merged_data.filt.var.xls")

  if (!file.exists(merged_var)) {
    var_files <- list.files(var_dir, pattern = "[.]var[.]xls$", full.names = TRUE)
    var_files <- var_files[!grepl("[.]filt[.]var[.]xls$", basename(var_files), ignore.case = TRUE)]
    if (length(var_files) > 0L) {
      bacterial_delivery_write_merged_table(var_files, merged_var)
      rows <- rbind(rows, bacterial_delivery_manifest_row(
        label = "merged_variant_table",
        source = paste(var_files, collapse = ";"),
        destination = merged_var,
        status = "generated"
      ))
    }
  }

  if (!file.exists(merged_filt)) {
    filt_files <- list.files(var_dir, pattern = "[.]filt[.]var[.]xls$", full.names = TRUE)
    if (length(filt_files) > 0L) {
      bacterial_delivery_write_merged_table(filt_files, merged_filt)
      rows <- rbind(rows, bacterial_delivery_manifest_row(
        label = "merged_filtered_variant_table",
        source = paste(filt_files, collapse = ";"),
        destination = merged_filt,
        status = "generated"
      ))
    }
  }

  rows
}

bacterial_delivery_write_merged_table <- function(files, output) {
  tables <- lapply(files, function(path) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  })
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(table) {
    missing <- setdiff(columns, names(table))
    for (col in missing) table[[col]] <- NA
    table[, columns, drop = FALSE]
  })
  merged <- do.call(rbind, tables)
  utils::write.table(merged, output, sep = "\t", quote = FALSE, row.names = FALSE)
}

bacterial_delivery_write_readmes <- function(delivery_dir,
                                             project,
                                             readme_name,
                                             chinese_readme_name) {
  rows <- bacterial_delivery_empty_manifest()
  readme <- file.path(delivery_dir, readme_name)
  writeLines(bacterial_delivery_readme(project), readme)
  rows <- rbind(rows, bacterial_delivery_manifest_row(
    label = "README",
    source = NA_character_,
    destination = readme,
    status = "generated"
  ))

  if (!is.null(chinese_readme_name)) {
    readme_zh <- file.path(delivery_dir, chinese_readme_name)
    writeLines(bacterial_delivery_readme_zh(project), readme_zh)
    rows <- rbind(rows, bacterial_delivery_manifest_row(
      label = "README_zh_CN",
      source = NA_character_,
      destination = readme_zh,
      status = "generated"
    ))
  }

  rows
}

bacterial_delivery_copy_sample_sheet <- function(sample_sheet, delivery_dir) {
  rows <- bacterial_delivery_empty_manifest()
  if (is.null(sample_sheet)) return(rows)

  destination <- file.path(delivery_dir, "Metadata", basename(sample_sheet))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(sample_sheet, destination, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy sample sheet: ", sample_sheet, call. = FALSE)
  }

  rbind(rows, bacterial_delivery_manifest_row(
    label = "sample_sheet",
    source = sample_sheet,
    destination = destination,
    status = "copied"
  ))
}

bacterial_delivery_readme <- function(project) {
  c(
    paste0("Project: ", project),
    "Bacterial genome/plasmid results delivery package",
    "=================================================",
    "",
    "Recommended files to review:",
    "- Sequence/*.consensus.fasta: final consensus sequences for each sample.",
    "- Var/*.filt.var.xls: filtered per-sample variant tables for routine review.",
    "- merged_data.filt.var.xls: filtered variants merged across delivered samples.",
    "- QC/*.png: quality-control and coverage figures, when available.",
    "- Metadata/: sample sheet or sample metadata table, when provided.",
    "",
    "Directory contents:",
    "- Bam/: read alignments against the consensus sequence and BAM index files.",
    "- QC/: coverage, read-length, or other quality-control figures.",
    "- Sequence/: consensus FASTA files and synthetic AB1 traces, when present.",
    "- Var/: per-sample variant tables and source VCF files, when generated.",
    "- Metadata/: sample sheet or project metadata table, when provided.",
    "- merged_data.var.xls: merged unfiltered variant table, when variant tables are available.",
    "- merged_data.filt.var.xls: merged filtered variant table, when filtered variant tables are available.",
    "- manifest.tsv: optional internal file list with source paths, present only when requested.",
    "- md5/md5.txt: MD5 checksums for delivered files.",
    "",
    "Variant table columns:",
    "- Chr, Pos, Ref, Alt: reference sequence, position, reference allele, and alternate allele.",
    "- DP, Ref_dp, Alt_dp, Freq: total depth, reference-supporting depth, alternate-supporting depth, and alternate allele frequency.",
    "- DP4: strand-level support in the order ref-forward, ref-reverse, alt-forward, alt-reverse.",
    "- Seq: local sequence context around the variant, when available."
  )
}

bacterial_delivery_readme_zh <- function(project) {
  c(
    paste0("项目：", project),
    "细菌基因组/质粒结果交付包",
    "========================",
    "",
    "建议优先查看：",
    "- Sequence/*.consensus.fasta：每个样本的最终共识序列。",
    "- Var/*.filt.var.xls：单样本过滤后的变异表，适合常规查看。",
    "- merged_data.filt.var.xls：所有交付样本合并后的过滤变异表。",
    "- QC/*.png：质控图、覆盖度图或读长分布图（如有）。",
    "- Metadata/：样本信息表或项目 metadata 表（如提供）。",
    "",
    "目录内容：",
    "- Bam/：reads 回帖到共识序列后的 BAM 比对文件及其索引文件。",
    "- QC/：覆盖度、读长分布或其他质控图片。",
    "- Sequence/：共识序列 FASTA 文件和合成 AB1 文件（如存在）。",
    "- Var/：单样本变异表以及生成的 VCF 文件（如有）。",
    "- Metadata/：样本信息表或项目 metadata 表（如提供）。",
    "- merged_data.var.xls：合并后的未过滤变异表（如有变异表）。",
    "- merged_data.filt.var.xls：合并后的过滤变异表（如有过滤变异表）。",
    "- manifest.tsv：可选的内部文件清单及来源路径，仅在请求时提供。",
    "- md5/md5.txt：交付文件的 MD5 校验值。",
    "",
    "变异表字段说明：",
    "- Chr、Pos、Ref、Alt：参考序列、位置、参考等位基因和替代等位基因。",
    "- DP、Ref_dp、Alt_dp、Freq：总深度、支持参考的深度、支持突变的深度和突变频率。",
    "- DP4：链向支持数，顺序为 ref 正链、ref 负链、alt 正链、alt 负链。",
    "- Seq：变异位点附近的序列上下文（如可获得）。"
  )
}

bacterial_delivery_empty_manifest <- function() {
  data.frame(
    label = character(),
    source = character(),
    delivery_path = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

bacterial_delivery_manifest_row <- function(label, source, destination, status) {
  data.frame(
    label = label,
    source = source,
    delivery_path = bacterial_delivery_relative_path(destination),
    status = status,
    stringsAsFactors = FALSE
  )
}

bacterial_delivery_relative_path <- function(path) {
  parts <- strsplit(normalizePath(path, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
  root_index <- match(TRUE, grepl("_delivery$", parts))
  if (!is.na(root_index)) {
    return(paste(parts[seq.int(root_index + 1L, length(parts))], collapse = "/"))
  }

  top_dirs <- c("Bam", "QC", "Sequence", "Var", "Metadata", "md5")
  top_index <- match(TRUE, parts %in% top_dirs)
  if (!is.na(top_index)) {
    return(paste(parts[seq.int(top_index, length(parts))], collapse = "/"))
  }

  basename(path)
}

bacterial_delivery_write_md5 <- function(delivery_dir) {
  md5_dir <- file.path(delivery_dir, "md5")
  dir.create(md5_dir, recursive = TRUE, showWarnings = FALSE)
  md5_file <- file.path(md5_dir, "md5.txt")
  files <- list.files(delivery_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  files <- files[normalizePath(files, mustWork = TRUE) != normalizePath(md5_file, mustWork = FALSE)]
  sums <- tools::md5sum(files)
  rel <- vapply(files, bacterial_delivery_relative_path, character(1))
  writeLines(paste(unname(sums), rel), md5_file)
  invisible(md5_file)
}

bacterial_delivery_make_archive <- function(delivery_dir, archive) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(dirname(delivery_dir))
  utils::tar(
    tarfile = archive,
    files = basename(delivery_dir),
    compression = "gzip",
    tar = "internal"
  )
}
