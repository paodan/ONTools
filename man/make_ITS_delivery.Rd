% Please edit documentation in R/make_ITS_delivery.R
\name{make_ITS_delivery}
\alias{make_ITS_delivery}
\title{Build an integrated ITS delivery folder}
\usage{
make_ITS_delivery(
  path_ITS_result = NULL,
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
  sample_length_plot_pattern = "^Distribution_seqLength__.*\\\\.png$",
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
  max_len = 2000,
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
  stderr = ""
)
}
\arguments{
\item{path_ITS_result}{Path to a completed wf-16s result directory for ITS data. Default is \code{NULL}. If \code{NULL} or if the directory does not exist, \code{run_ITS()} is run first. When an existing directory is supplied, the function expects it to contain \code{abundance_table_genus.tsv} and \code{alignment_tables/}.}
\item{path_delivery}{Final ITS delivery root. Default is \code{"/data/project_delivery"}. In grouped mode, per-group ITS deliveries are created under \code{path_delivery/<group>/ITS/}.}
\item{consensus_delivery_path}{Existing consensus delivery directory. If \code{NULL} or if the directory does not exist, \code{make_consensus_delivery()} is run first using \code{path_proj}, \code{path_sampleInfo_file_list}, and the consensus-related parameters in this function.}
\item{path_proj, path_sampleInfo_file_list}{Arguments passed to \code{make_consensus_delivery()} when \code{consensus_delivery_path} is missing. When \code{path_sampleInfo_file_list} is supplied, its names are also used as delivery group names, matching \code{make_consensus_delivery()}. Each sample information file is copied to the corresponding project root.}
\item{consensus_delivery_output}{Directory used to store newly generated consensus delivery results. The default \code{NULL} creates a sibling work directory named \code{basename(path_delivery)_consensus_work}.}
\item{kit_name, model, demux_out, barcode_both_ends}{Dorado basecalling and demultiplexing parameters passed to \code{make_consensus_delivery()} when consensus results need to be generated. Defaults match \code{make_consensus_delivery()}.}
\item{fastq_out}{FASTQ output directory name created by Dorado conversion. This matches the \code{fastq_out} argument in \code{make_consensus_delivery()}. When \code{make_consensus_delivery()} is run in this function, the actual per-project FASTQ path passed to \code{run_ITS()} is taken from the consensus workflow result, typically \code{.../fastq_pass_trim/<group>}. If consensus results are not generated in this call, \code{fastq_out} may also be an existing FASTQ directory path. Default is \code{"fastq_pass_trim"}.}
\item{run_basecalling_demux_step, run_dorado_basecall_step, run_dorado_demux_step, run_dorado_fastq_step, run_QC_step, move_fastq_step, move_fastq_mode, run_amplicon_step, trim_consensus_step, run_filtered_QC_step, run_igv_step, collect_results_step, make_ab1}{Consensus-delivery step controls passed to \code{make_consensus_delivery()} when \code{consensus_delivery_path} is \code{NULL} or missing. Defaults match \code{make_consensus_delivery()}. They have no effect when an existing \code{consensus_delivery_path} is supplied.}
\item{ab1_name_template, ab1_samtools, consensus_index_file, sample_length_plot_pattern, project_col, amplicon_size_col, min_read_length_col, max_read_length_col, f_primer_col, r_primer_col, min_read_qual, min_n_reads, force_spoa_length_threshold, override_basecaller_cfg, amplicon_extra_args, barcode_digits, dorado_threads, dorado_fastq_write_md5, dorado, samtools, gzip, dorado_conda_env, igv, overwrite_fastq, overwrite_delivery, include_execution}{Additional \code{make_consensus_delivery()} parameters passed through when consensus results need to be generated. Defaults match \code{make_consensus_delivery()}.}
\item{path_work}{Shared work root used to derive default ITS workflow output and work directories. If \code{NULL}, sibling directories of \code{path_delivery} are used.}
\item{out_dir}{Output directory passed to \code{run_ITS()}. Default \code{NULL} creates a sibling/work directory named \code{basename(path_delivery)_wf_ITS}; this directory then becomes the effective \code{path_ITS_result}.}
\item{work_dir}{Nextflow work directory passed to \code{run_ITS()}. Default \code{NULL} creates a sibling/work directory named \code{basename(path_delivery)_wf_ITS_work}.}
\item{profile}{Nextflow profile passed to \code{run_ITS()}. Default is \code{"standard"}.}
\item{resume}{Logical. If \code{TRUE} (default), pass \code{-resume} to Nextflow through \code{run_ITS()}.}
\item{database_set}{wf-16s database set used by \code{run_ITS()}. Default is \code{"ncbi_16s_18s_28s_ITS"}.}
\item{min_len, max_len}{Minimum and maximum read length passed to \code{run_ITS()}. Defaults are \code{300} and \code{2000}.}
\item{workflow}{Nextflow workflow name or path passed to \code{run_ITS()}. Default is \code{"epi2me-labs/wf-16s"}.}
\item{nextflow}{Nextflow executable name or path. Default is \code{"nextflow"}.}
\item{quiet}{Logical. If \code{TRUE}, pass \code{-q} to Nextflow. Default is \code{FALSE}.}
\item{extra_args}{Additional raw wf-16s command-line arguments passed through \code{run_ITS()}. Default is \code{"--minimap2_by_reference"}.}
\item{syntax_parser}{Nextflow syntax parser version passed as \code{NXF_SYNTAX_PARSER}. Default is \code{"v1"}. Use \code{NULL} to leave it unchanged.}
\item{ansi_log}{Logical. Passed as \code{NXF_ANSI_LOG}. Default is \code{FALSE}.}
\item{nextflow_env}{Optional additional environment variables passed to \code{run_ITS()}, formatted as \code{"NAME=value"}. Default is \code{NULL}.}
\item{run_ITS_step}{Logical. If \code{TRUE}, run \code{run_ITS()} when \code{path_ITS_result} is \code{NULL} or missing. Default is \code{TRUE}.}
\item{move_ITS_step}{Logical. If \code{TRUE}, organize wf-16s ITS results into the delivery structure using \code{move_ITS()} outputs. Default is \code{TRUE}.}
\item{ITS_dir}{Directory name created under each barcode sample directory for ITS-specific per-sample evidence. Default is \code{"ITS_results"}, producing \code{path_delivery/samples/barcode*/ITS_results/} in single-project mode and \code{path_delivery/<group>/samples/barcode*/ITS_results/} in grouped mode.}
\item{barcode_col}{Column name in each sample information file used to infer barcode directory names. Default is \code{"Barcode_ID"}, matching \code{make_consensus_delivery()}. Values like \code{PBC001-097} are converted to \code{barcode097}; otherwise values are used as-is after resolving the barcode column.}
\item{samples_dir}{Directory under \code{path_delivery} that stores per-barcode sample results. The default \code{"samples"} creates \code{path_delivery/samples/barcode*/}. Set to \code{NULL} to put \code{barcode*/} directly under \code{path_delivery}.}
\item{abundance_table}{Name of the wf-16s abundance table under \code{path_ITS_result}. Default is \code{"abundance_table_genus.tsv"}. In the final delivery this file is copied to the root of \code{path_delivery}.}
\item{alignment_tables_dir}{Name of the wf-16s per-barcode alignment table directory under \code{path_ITS_result}. Default is \code{"alignment_tables"}.}
\item{figure_dir}{Name of the figure directory created at the final delivery root. Default is \code{"figures"}.}
\item{identification_dir}{Name of the per-barcode directory that receives wf-16s alignment tables. Default is \code{"identification_tables"}, producing \code{samples/barcode*/ITS_results/identification_tables/}.}
\item{tax_levels}{Taxonomic levels plotted from the abundance table. Default is \code{c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")}.}
\item{cutoff}{Minimum relative abundance kept as a separate plotted taxon by \code{plot_abundance_bar()}. Default is \code{0.01}.}
\item{width, height}{Plot size in inches. Defaults are \code{width = 12} and \code{height = 6}.}
\item{consensus_file}{Fallback consensus FASTA filename searched under the consensus delivery. Default is \code{"all-consensus-seqs.fasta"}. The selected consensus sequences are copied or collected at the project root as this filename and used for optional UNITE review.}
\item{trimmed_consensus_file}{Preferred consensus FASTA filename used before falling back to \code{consensus_file}. Default is \code{"all-consensus-seqs_trimmed.fasta"}. When this file is available, it is also preserved under this original filename in the project root. Matching FASTA indexes are copied when available.}
\item{barcode_pattern}{Regular expression used to identify barcode directories and alignment-stat files. The default \code{"^barcode[0-9]+$"} matches names such as \code{barcode097} and \code{barcode303}.}
\item{run_unite_annotation}{Logical. If \code{TRUE}, run \code{annotate_consensus_blast()} on the collected consensus FASTA. Default is \code{TRUE}; if neither \code{unite_db} nor \code{unite_db_fasta} is supplied, this step is skipped with a warning.}
\item{unite_db}{Prefix of an already-built nucleotide BLAST database, for example \code{"/data/reference/UNITE/unite_eukaryotes"}. Default is \code{NULL}.}
\item{unite_db_fasta}{FASTA file used to build a BLAST database with \code{makeblastdb}. Default is \code{NULL}. Use this when the UNITE FASTA has been downloaded but a BLAST database has not yet been built.}
\item{unite_dir}{Directory under \code{path_delivery} for detailed UNITE BLAST results. Default is \code{"unite_consensus_annotation"}.}
\item{unite_top_hits_name}{Filename copied to \code{path_delivery} for the UNITE top-hit summary. Default is \code{"unite_consensus_top_hits.tsv"}, placed beside \code{abundance_table_genus.tsv} for quick review.}
\item{unite_threads}{Number of BLAST threads. Default is \code{10}.}
\item{unite_max_target_seqs}{Maximum number of target hits retained per consensus sequence. Default is \code{20}.}
\item{unite_evalue}{BLAST E-value cutoff. Default is \code{"1e-20"}.}
\item{unite_task}{Optional BLAST task, such as \code{"blastn"} or \code{"megablast"}. Default \code{NULL} lets BLAST choose its default.}
\item{unite_word_size}{Optional BLAST word size. Default is \code{NULL}.}
\item{unite_strand}{Optional BLAST strand argument, such as \code{"both"}, \code{"plus"}, or \code{"minus"}. Default is \code{NULL}.}
\item{unite_dust}{Optional BLAST low-complexity filtering setting, such as \code{"yes"} or \code{"no"}. Default is \code{NULL}.}
\item{unite_perc_identity}{Optional minimum BLAST percent identity. Default is \code{NULL}.}
\item{unite_extra_args}{Additional command-line arguments passed to \code{blastn}. Default is \code{NULL}.}
\item{unite_species_identity, unite_genus_identity, unite_family_identity}{Percent identity thresholds used by \code{annotate_consensus_blast()} to label the suggested annotation level. Defaults are \code{98.5}, \code{95}, and \code{90}.}
\item{unite_min_query_coverage}{Minimum query coverage percentage used for confident annotation-level assignment. Default is \code{80}.}
\item{unite_min_reference_coverage}{Minimum reference coverage percentage used for confident annotation-level assignment. Default is \code{50}.}
\item{unite_novel_identity}{Identity threshold below which a sufficiently covered top hit is flagged as a possible novel-candidate signal. Default is \code{97}.}
\item{unite_blastn}{Executable name or path for \code{blastn}. Default is \code{"blastn"}.}
\item{unite_makeblastdb}{Executable name or path for \code{makeblastdb}. Default is \code{"makeblastdb"}.}
\item{unite_conda_env}{Optional conda environment used for BLAST commands. Default \code{NULL} calls tools from the current environment.}
\item{conda}{Conda executable used when \code{unite_conda_env} is not \code{NULL}. Default is \code{"conda"}.}
\item{readme_name, chinese_readme_name}{README filenames written under \code{path_delivery}. Defaults are \code{"README.txt"} and \code{"README.zh-CN.txt"}. Set \code{chinese_readme_name = NULL} to skip the Chinese README.}
\item{overwrite}{Logical. If \code{TRUE}, replace an existing \code{path_delivery}. Default is \code{FALSE} to avoid accidentally deleting previous deliveries.}
\item{dry_run}{Logical. If \code{TRUE}, return a plan without copying files or running external tools. Default is \code{FALSE}.}
\item{echo}{Logical. If \code{TRUE}, print commands from wrapped runners. Default is \code{TRUE}.}
\item{wait, stdout, stderr}{Passed to wrapped command runners. Defaults are \code{wait = TRUE}, \code{stdout = ""}, and \code{stderr = ""}.}
}
\value{
Invisibly returns a list describing generated and copied outputs.
}
\description{
\code{make_ITS_delivery()} combines amplicon consensus delivery results,
wf-16s ITS taxonomic profiling results, and optional consensus-vs-UNITE BLAST
review results into one customer-facing delivery folder.

The default layout keeps per-sample files under \code{samples/barcode*/} and
keeps project-level summaries, figures, and UNITE top-hit results directly
under \code{path_delivery}. This makes the root directory useful for quick review
while keeping each barcode's detailed evidence in one place.
}
\details{
The function performs four main tasks:
\enumerate{
\item Prepare consensus-delivery content. If \code{consensus_delivery_path}
exists, its barcode folders are copied into the final delivery. Otherwise,
\code{make_consensus_delivery()} is run first. This step also creates or locates
the grouped FASTQ directory used as input for ITS profiling.
\item Prepare wf-16s ITS content. If \code{path_ITS_result} is \code{NULL} or
missing, \code{run_ITS()} is called after consensus preparation, using the
grouped FASTQ path returned by \code{make_consensus_delivery()}. \code{move_ITS()}
is then called in a temporary directory to generate abundance plots and
normalize wf-16s outputs. Because the delivery assembly needs the completed ITS
output files, \code{wait = TRUE} is required when ITS results must be generated.
\item Reorganize files. When \code{path_sampleInfo_file_list} is supplied, one
delivery directory is created for each named sample-info file:
\code{path_delivery/<group>/ITS/}. Each group keeps only the barcode samples
listed in its sample-info table. Without \code{path_sampleInfo_file_list},
\code{path_delivery} itself is treated as the final ITS delivery directory.
Per-barcode consensus outputs are copied to
\code{samples/barcode*/consensus_results/}, per-barcode wf-16s alignment tables
are copied to \code{samples/barcode*/ITS_results/identification_tables/}, and
project-level abundance files, figures, sample information, and consensus FASTA
are copied to the project root.
\item Optionally run consensus-vs-UNITE review. Consensus FASTA files are
collected into the project-root \code{all-consensus-seqs.fasta}, BLAST output is
written to \code{unite_consensus_annotation/consensus.blast.tsv}, and the
top-hit summary is written both to
\code{unite_consensus_annotation/consensus.top_hits.tsv} and to the root-level
\code{unite_consensus_top_hits.tsv}.
}
}
\examples{
\dontrun{
make_ITS_delivery(
  path_delivery = "delivery/ITS",
  path_ITS_result = "results/wf_its",
  consensus_delivery_path = "delivery/consensus",
  unite_db = "/data/reference/UNITE/unite_eukaryotes",
  unite_threads = 10
)

make_ITS_delivery(
  path_delivery = "delivery/ITS",
  path_ITS_result = NULL,
  fastq_out = "fastq_pass_trim",
  consensus_delivery_path = NULL,
  path_proj = "/data/minknow/project/run",
  path_sampleInfo_file_list = "SampleInfo.csv",
  unite_db_fasta = "/data/reference/UNITE/UNITE_eukaryotes_all.fasta"
)
}
}
