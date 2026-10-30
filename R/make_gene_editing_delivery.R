#' Build a gene-editing analysis delivery
#'
#' Run the gene-editing delivery workflow for standardized sample metadata.
#' The workflow can subset source BAM files, call and summarize variants,
#' generate synthetic AB1 files and FASTQ files, draw read-length and IGV
#' quality-control figures, and collect the resulting files into a delivery
#' directory.
#'
#' `sample_info` is expected to contain one row per sample/barcode. A row may
#' specify multiple candidate references in `ref_name` by separating their
#' names with `"|"`, for example `"PLA3-A|PLA3-B|PLA3-D"`. The row's
#' `region_start` and `region_end` are applied to every reference named in that
#' row. The names in `ref_name` refer to the names used by `ref_all`; `re_names`
#' maps the sanitized names stored in the BAM and sanitized FASTA to those
#' reference names.
#'
#' When BAM subsetting is enabled, reads must be primary alignments with mapping
#' quality at least 20 and must completely cover the configured closed interval.
#' At most `max_reads` eligible alignments are retained across all references in
#' a row. The subset BAM is indexed and replaces an existing work BAM.
#'
#' Variant calling creates two workbooks. Both apply a minimum of three reads
#' per event and combination, `min_event_ratio`, at least 20 callable reads, and
#' mapping quality 20. The workbook at `paths$xlsx_file_full` uses a combination
#' ratio threshold of zero; the workbook at `paths$xlsx_file_filtered` uses
#' `min_comb_ratio`.
#'
#' One synthetic AB1 file is generated for every sanitized reference associated
#' with a sample. The reference name is inserted before the `.ab1` suffix of
#' `sample_info$ab1_file`. FASTQ and read-length QC outputs are generated once
#' per row. IGV snapshots are generated once per sanitized reference.
#'
#' @param sample_info A data frame in the format returned by
#'   [standard_sample_info_gene_editing()]. It must contain `sample_id`,
#'   `barcode`, `ref_name`, `region_start`, `region_end`, `bam_file_origin`,
#'   `bam_file`, `output_fastq`, `reference_file`, and `ab1_file`. Multiple
#'   references in `ref_name` must be separated by `"|"` and share the row's
#'   analysis interval.
#' @param paths A named path list, normally returned by
#'   [standard_paths_gene_editing()]. It must contain `ref_sanitized`,
#'   `xlsx_file_full`, `xlsx_file_filtered`, `sample_info_path`,
#'   `work_bam_path`, `work_qc_path`, `work_ab1_path`, and `delivery_path`.
#'   Variant calling also uses `event_cache_dir`.
#' @param ref_all Named reference sequences accepted by [run_editing_batch2()].
#'   Names must correspond to the reference names stored in `sample_info$ref_name`.
#' @param re_names A named character vector mapping sanitized BAM/FASTA reference
#'   names to the names in `ref_all`. For example,
#'   `c(PLA3_A_ = "PLA3-A")`. Every reference listed in `sample_info$ref_name`
#'   must have a matching value.
#' @param max_reads Positive number giving the maximum total number of eligible
#'   alignments retained in each sample's subset BAM across all references.
#' @param min_event_ratio Minimum event-support percentage passed to
#'   [run_editing_batch2()] for both output workbooks.
#' @param min_comb_ratio Minimum variant-combination percentage used for the
#'   filtered workbook. The full workbook uses zero for this threshold.
#' @param upstream = 20, upstream of the start editing site.
#' @param downstream = 20, downstream of the end editing site.
#' @param run_subset_bam_step Logical. If `TRUE`, generate indexed work BAMs
#'   from `bam_file_origin`; otherwise, existing paths in `bam_file` are used.
#' @param run_variant_calling Logical. If `TRUE`, run gene-editing variant
#'   calling and write the full and filtered Excel workbooks.
#' @param make_ab1 Logical. If `TRUE`, generate one synthetic AB1 file for every
#'   sanitized reference associated with each retained sample.
#' @param make_fastq Logical. If `TRUE`, convert each retained work BAM to the
#'   corresponding `output_fastq` file.
#' @param mc.cores Number of forked workers used by the per-row workflow stages
#'   and passed to [run_editing_batch2()]. If `NULL`, the smaller of the sample
#'   count and the number reported by [parallel::detectCores()] is used, with a
#'   minimum of one.
#' @param plot_qc_peaks Logical. If `TRUE`, draw a read-length peak plot for each
#'   existing output FASTQ and save it under `work_qc_path/Read_length`.
#' @param plot_qc_igv Logical. If `TRUE`, generate a primary-alignment IGV
#'   snapshot for each sanitized reference under `work_qc_path/IGV`.
#' @param run_delivery Logical. If `TRUE`, collect the sanitized reference and
#'   its index, both variant workbooks, the source sample-information file, and
#'   the work BAM, QC, and AB1 directories under `paths$delivery_path`.
#' @param igv Character scalar giving the path or command used to launch IGV.
#'
#' @return A list with two elements: `sample_info_valid`, containing the rows
#'   whose source and work BAM files were found at the relevant stages, and
#'   `paths`, containing the input path list. The function also writes the
#'   selected workflow outputs as side effects.
#'
#' @details
#' Source-BAM rows that do not exist are removed before processing. After the
#' optional subset step, rows whose work BAM does not exist are also removed.
#' An error is raised if no rows remain at either point.
#'
#' The per-row stages use fork-based [parallel::mclapply()]. On platforms where
#' forked processing is unavailable, the behavior and degree of parallelism are
#' governed by `parallel`'s platform-specific implementation.
#'
#' Existing subset BAM, AB1, and FASTQ outputs are overwritten. The subset BAM
#' sampling seed is fixed at 1234. AB1 generation uses the defaults of
#' [synthetic_ab1_from_bam()]. IGV output contains primary alignments only.
#'
#' @examples
#' \dontrun{
#' library(ONTools)
#'
#' paths <- standard_paths_gene_editing(
#'   input_path_root = "/analysis/wf_amplicon_variant/project01",
#'   delivery_path_root = "/delivery/project01",
#'   ref_sanitized_file = "reference_sanitized_seqIDs.fasta"
#' )
#' paths$sample_info_path <- "/metadata/sample_info_for_gene_editing.xlsx"
#' create_standard_directories(
#'   paths,
#'   create_work = TRUE,
#'   create_delivery = FALSE
#' )
#'
#' source_info <- as.data.frame(readxl::read_excel(paths$sample_info_path))
#' sample_info <- standard_sample_info_gene_editing(
#'   sample_info = source_info,
#'   input_path_root = paths$input_path_root,
#'   delivery_path_root = paths$delivery_path_root,
#'   delivery_folder = "gene_editing_delivery",
#'   work_path = paths$work_path,
#'   sample_id_column = "Sample_HiTOM",
#'   barcode_column = "barcode",
#'   ref_name_column = "Reference_name",
#'   ref_region_start_column = "Reference_Region_Start",
#'   ref_region_end_column = "Reference_Region_End",
#'   ref_sanitized_file = "reference_sanitized_seqIDs.fasta"
#' )
#'
#' re_names <- c(
#'   PLA3_A_ = "PLA3-A",
#'   PLA3_B_ = "PLA3-B",
#'   PLA3_D_ = "PLA3-D"
#' )
#' ref_all <- Biostrings::readDNAStringSet(
#'   "/references/HI-TOM-PLA3PLA10.txt"
#' )
#'
#' result <- make_gene_editing_delivery(
#'   sample_info = sample_info,
#'   paths = paths,
#'   ref_all = ref_all,
#'   re_names = re_names,
#'   max_reads = 2500,
#'   min_event_ratio = 0.5,
#'   min_comb_ratio = 0.5,
#'   mc.cores = 8,
#'   igv = "/usr/local/bin/IGV_Linux_2.19.8/igv.sh"
#' )
#' result$sample_info_valid
#' }
#'
#' @seealso [standard_paths_gene_editing()],
#'   [standard_sample_info_gene_editing()], [subset_bam_reads()],
#'   [run_editing_batch2()], [synthetic_ab1_from_bam()], [igv_snapshot()],
#'   [collect_items()]
#' @importFrom parallel mclapply
#' @export
make_gene_editing_delivery = function(sample_info, paths, ref_all, re_names,
                                      max_reads = 2500,
                                      min_event_ratio = 0.5, min_comb_ratio = 0.5,
                                      upstream = 20,
                                      downstream = 20,
                                      run_subset_bam_step = TRUE,
                                      run_variant_calling = TRUE,
                                      make_ab1 = TRUE,
                                      make_fastq = TRUE,
                                      mc.cores = NULL,
                                      plot_qc_peaks = TRUE,
                                      plot_qc_igv = TRUE,
                                      run_delivery = TRUE,
                                      igv = "/usr/local/bin/IGV_Linux_2.19.8/igv.sh"
                                      ){
  missing_col = c()
  for(mi in c("sample_id", "barcode", "ref_name", "region_start", "region_end",
              "bam_file_origin", "bam_file", "output_fastq", "reference_file", "ab1_file")){
    if(!mi %in% colnames(sample_info)){
      missing_col = c(missing_col, mi)
    }
  }

  if(length(missing_col) > 0){
    stop("`", paste(missing_col, collapse = "`, `"), "` column must exist in `sample_info`")
  }

  sample_info$region_start = max(1, sample_info$region_start - upstream)

  sample_info$region_end = unname(sapply(sample_info$ref_name, function(x) {
    nm = strsplit(x, split = "\\|")
    min(Biostrings::width(ref_all[nm]),
        sample_info$region_end + downstream)
  }))


  path_names = c("ref_sanitized", "xlsx_file_full", "xlsx_file_filtered", "sample_info_path",
                 "work_bam_path", "work_qc_path", "work_ab1_path", "delivery_path")
  not_exist_paths = path_names[!(path_names %in% names(paths))]
  if(length(not_exist_paths) > 0){
    stop("These path names do not exist:\n", paste(not_exist_paths, collapse = ", "))
  }

  sample_info = sample_info[file.exists(sample_info[["bam_file_origin"]]), , drop = FALSE]
  if(nrow(sample_info) < 1){
    stop("No bam file is found")
  }
  if(is.null(mc.cores)){
    mc.cores = max(1, min(nrow(sample_info), parallel::detectCores()))
  }

  sanitized_ref_name = function(x, re_names){
    names(re_names)[re_names %in% strsplit(x, "\\|")[[1]]]
  }

  if(run_subset_bam_step){
    message("run_subset_bam_step")
    tmp = mclapply(1:nrow(sample_info), function(mi){
      # if(file.exists(sample_info$bam_file_origin[mi])){}
      subset_bam_reads(bam = sample_info$bam_file_origin[mi],
                       output_bam = sample_info$bam_file[mi],
                       reference_name = sanitized_ref_name(sample_info$ref_name[mi], re_names),
                       min_mapq = 20, max_reads = max_reads,
                       region_start = sample_info$region_start[mi],
                       region_end = sample_info$region_end[mi],
                       coverage = "complete",
                       primary_only = TRUE,
                       seed = 1234, index_output = TRUE, overwrite = TRUE)

    }, mc.cores = mc.cores)
  }
  sample_info = sample_info[file.exists(sample_info[["bam_file"]]), , drop = FALSE]
  if(nrow(sample_info) < 1){
    stop("No bam_file is found")
  }

  if(run_variant_calling){
    message("run_variant_calling")
    # Generate the gene variant tables
    for(fi in c("xlsx_file_full", "xlsx_file_filtered")){
      # Remove the low quality events (100 seconds)
      res_all_rbind_0.50 = ONTools::run_editing_batch2(
        sample_info = sample_info[file.exists(sample_info$bam_file),,drop = FALSE],
        ref_all = ref_all,
        re_names = re_names,
        cache_dir = paths$event_cache_dir,
        min_anchor = 0,
        min_event_reads = 3,
        min_event_ratio = min_event_ratio,  # Important parameter
        min_comb_reads = 3,
        min_comb_ratio = c(xlsx_file_full = 0, xlsx_file_filtered = min_comb_ratio)[fi],     # Important parameter
        min_callable_reads = 20L,
        min_mapq = 20L,
        colored_xlsx = paths[[fi]], mc.cores = mc.cores)
    }
  }


  # Generate ab1 files
  if(make_ab1){
    message("make_ab1")

    tmp = mclapply(1:nrow(sample_info), function(mi){
      sample_info_mi = sample_info[mi, , drop = FALSE]
      if(file.exists(sample_info_mi$bam_file)){
        for(refi in sanitized_ref_name(sample_info_mi$ref_name, re_names)){
          synthetic_ab1_from_bam(sample_info_mi$reference_file,
                                 bam = sample_info_mi$bam_file,
                                 output_ab1 = gsub("\\.ab1$", paste0("_", refi, ".ab1"),
                                                   sample_info_mi$ab1_file),
                                 reference_name = refi,
                                 overwrite = TRUE)
        }
      }
    }, mc.cores = mc.cores)
  }

  # Convert bam to fastq files
  if(make_fastq){
    message("make_fastq")
    tmp = mclapply(1:nrow(sample_info), function(mi){
      if(file.exists(sample_info$bam_file[mi])){
        bam_to_fastq(sample_info$bam_file[mi], output_fastq = sample_info$output_fastq[mi],
                     threads = 1, copy_tags = TRUE, overwrite = TRUE)
      }
    }, mc.cores = mc.cores)
  }

  if(plot_qc_peaks){
    message("plot_qc_peaks")
    # Plot peaks
    tmp = mclapply(1:nrow(sample_info), function(mi){
      sample_info_mi = sample_info[mi,,drop = FALSE]
      fastq_file_mi = sample_info_mi$output_fastq
      if(file.exists(fastq_file_mi)){
        p = plot_fastq_peaks(fastq_file_mi,
                             min_distance = (sample_info_mi$region_end - sample_info_mi$region_start)/2,
                             left_offset = (sample_info_mi$region_end - sample_info_mi$region_start)/3,
                             right_offset = (sample_info_mi$region_end - sample_info_mi$region_start)/3,
                             print_plots = FALSE, subtitle = sample_info_mi$sample_id)

        file_mi = gsub(".fastq.gz", "_read_length.png", basename(fastq_file_mi))
        ggplot2::ggsave(filename = file.path(paths$work_qc_path,  "Read_length", file_mi),
                        plot = p[[1]], width = 7, height = 4.5, create.dir = TRUE)
        return(1)
      } else {
        return(0)
      }
    }, mc.cores = mc.cores)
  }

  # Save IGV figures
  if(plot_qc_igv){
    message("plot_qc_igv")
    tmp = mclapply(1:nrow(sample_info), function(i){
      if(file.exists(sample_info$bam_file[i])){
        ref_names = sanitized_ref_name(sample_info$ref_name[[i]], re_names)
        for(ci in ref_names){
          xxx = igv_snapshot(genome_fasta = sample_info$reference_file[i],
                             bam = sample_info$bam_file[i],
                             chr = ci,
                             out_dir = file.path(paths$work_qc_path, "IGV"),
                             igv = igv,
                             dry_run = F,
                             snapshot_name = paste0(sample_info$sample_id[i],"_", ci, "_primary_igv.png"),
                             igv_primary_only = T)
        }
      }
    }, mc.cores = mc.cores)
  }

  # Delivery
  if(run_delivery){
    message("run_delivery")
    dir.create(paths$delivery_path)
    collect_items(
      source = c(paths$ref_sanitized, paste0(paths$ref_sanitized, ".fai"),
                 unlist(paths[c("xlsx_file_full", "xlsx_file_filtered", "sample_info_path",
                                "work_bam_path", "work_qc_path", "work_ab1_path")])),
      destination_dir = paths$delivery_path)
  }
  return(list(sample_info_valid = sample_info, paths = paths))
}
