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
  ITS_dir = "ITS_results",
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
  ...
)
}
\arguments{
\item{path_ITS_result}{Path to a completed wf-16s result directory for ITS data. Default is \code{NULL}. If \code{NULL} or if the directory does not exist, \code{run_ITS()} is run first. When an existing directory is supplied, the function expects it to contain \code{abundance_table_genus.tsv} and \code{alignment_tables/}.}
\item{path_delivery}{Final ITS delivery root. Default is \code{"/data/project_delivery"}. In grouped mode, per-group ITS deliveries are created below this directory.}
\item{consensus_delivery_path}{Existing consensus delivery directory. If \code{NULL} or if the directory does not exist, \code{make_consensus_delivery()} is run first using \code{path_proj}, \code{path_sampleInfo_file_list}, and \code{...}.}
\item{path_proj, path_sampleInfo_file_list}{Arguments passed to \code{make_consensus_delivery()} when \code{consensus_delivery_path} is missing. When \code{path_sampleInfo_file_list} is supplied, its names are also used as delivery group names, matching \code{make_consensus_delivery()}.}
\item{consensus_delivery_output}{Directory used to store newly generated consensus delivery results. The default \code{NULL} creates a sibling work directory named \code{basename(path_delivery)_consensus_work}.}
\item{ITS_fastq}{FASTQ directory passed to \code{run_ITS()} when ITS results need to be generated. Default is \code{"./fastq_pass_trim"}.}
\item{ITS_out_dir}{Output directory passed to \code{run_ITS()}. Default \code{NULL} creates a sibling directory named \code{basename(path_delivery)_wf_ITS}; this directory then becomes the effective \code{path_ITS_result}.}
\item{ITS_work_dir}{Nextflow work directory passed to \code{run_ITS()}. Default \code{NULL} creates a sibling directory named \code{basename(path_delivery)_wf_ITS_work}.}
\item{ITS_profile}{Nextflow profile passed to \code{run_ITS()}. Default is \code{"standard"}.}
\item{ITS_resume}{Logical. If \code{TRUE} (default), pass \code{-resume} to Nextflow through \code{run_ITS()}.}
\item{ITS_database_set}{wf-16s database set used by \code{run_ITS()}. Default is \code{"ncbi_16s_18s_28s_ITS"}.}
\item{ITS_min_len, ITS_max_len}{Minimum and maximum read length passed to \code{run_ITS()}. Defaults are \code{300} and \code{2000}.}
\item{ITS_workflow}{Nextflow workflow name or path passed to \code{run_ITS()}. Default is \code{"epi2me-labs/wf-16s"}.}
\item{ITS_nextflow}{Nextflow executable name or path. Default is \code{"nextflow"}.}
\item{ITS_quiet}{Logical. If \code{TRUE}, pass \code{-q} to Nextflow. Default is \code{FALSE}.}
\item{ITS_extra_args}{Additional raw wf-16s command-line arguments passed through \code{run_ITS()}. Default is \code{"--minimap2_by_reference"}.}
\item{ITS_syntax_parser}{Nextflow syntax parser version passed as \code{NXF_SYNTAX_PARSER}. Default is \code{"v1"}. Use \code{NULL} to leave it unchanged.}
\item{ITS_ansi_log}{Logical. Passed as \code{NXF_ANSI_LOG}. Default is \code{FALSE}.}
\item{ITS_nextflow_env}{Optional additional environment variables passed to \code{run_ITS()}, formatted as \code{"NAME=value"}. Default is \code{NULL}.}
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
\item{consensus_file}{Fallback consensus FASTA filename searched under the consensus delivery. Default is \code{"all-consensus-seqs.fasta"}.}
\item{trimmed_consensus_file}{Preferred consensus FASTA filename for UNITE review. Default is \code{"all-consensus-seqs_trimmed.fasta"}.}
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
\item{...}{Additional arguments passed to \code{make_consensus_delivery()} when a new consensus delivery must be generated.}
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
\code{make_consensus_delivery()} is run first.
\item Prepare wf-16s ITS content. If \code{path_ITS_result} is \code{NULL} or
missing, \code{run_ITS()} is called first. \code{move_ITS()} is then called in a
temporary directory to generate abundance plots and normalize wf-16s outputs.
Because the delivery assembly needs the completed ITS output files,
\code{wait = TRUE} is required when ITS results must be generated.
\item Reorganize files. When \code{path_sampleInfo_file_list} is supplied, one
delivery directory is created for each named sample-info file:
\code{path_delivery/<group>/}. Each group keeps only the barcode samples listed
in its sample-info table. Without \code{path_sampleInfo_file_list},
\code{path_delivery} itself is treated as the final ITS delivery directory.
Per-barcode consensus outputs are copied to
\code{samples/barcode*/consensus_results/}, per-barcode wf-16s alignment tables
are copied to \code{samples/barcode*/ITS_results/identification_tables/}, and
project-level abundance files and figures are copied to the project root.
\item Optionally run consensus-vs-UNITE review. Consensus FASTA files are
collected into \code{consensus_for_unite.fasta}, BLAST output is written to
\code{unite_consensus_annotation/consensus.blast.tsv}, and the top-hit summary
is written both to
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
  ITS_fastq = "fastq_pass_trim",
  consensus_delivery_path = NULL,
  path_proj = "/data/minknow/project/run",
  path_sampleInfo_file_list = "SampleInfo.csv",
  unite_db_fasta = "/data/reference/UNITE/UNITE_eukaryotes_all.fasta"
)
}
}
