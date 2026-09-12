% Please edit documentation in R/make_ITS_delivery.R
\name{make_ITS_delivery}
\alias{make_ITS_delivery}
\title{Build an integrated ITS delivery folder}
\usage{
make_ITS_delivery(
  path_ITS_result,
  output_dir,
  consensus_delivery_path = NULL,
  path_proj = NULL,
  path_sampleInfo_file_list = NULL,
  consensus_delivery_output = NULL,
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
  overwrite = TRUE,
  dry_run = FALSE,
  echo = TRUE,
  wait = TRUE,
  stdout = "",
  stderr = "",
  ...
)
}
\arguments{
\item{path_ITS_result}{Path to a completed wf-16s result directory for ITS data.}
\item{output_dir}{Final ITS delivery directory.}
\item{consensus_delivery_path}{Existing consensus delivery directory. If \code{NULL} or missing, \code{make_consensus_delivery()} is run first.}
\item{path_proj, path_sampleInfo_file_list}{Arguments passed to \code{make_consensus_delivery()} when \code{consensus_delivery_path} is missing.}
\item{consensus_delivery_output}{Directory used to store newly generated consensus delivery results. If \code{NULL}, a sibling work directory is used.}
\item{samples_dir}{Directory under \code{output_dir} that stores per-barcode sample results.}
\item{abundance_table, alignment_tables_dir, figure_dir, identification_dir}{File or directory names used by \code{move_ITS()} and the final delivery layout.}
\item{tax_levels, cutoff, width, height}{Plotting parameters passed to \code{move_ITS()}.}
\item{consensus_file, trimmed_consensus_file}{Consensus FASTA names searched under the consensus delivery. Trimmed consensus files are preferred for UNITE annotation.}
\item{barcode_pattern}{Regular expression used to identify barcode directories and alignment-stat files.}
\item{run_unite_annotation}{Logical. If \code{TRUE}, run \code{annotate_consensus_blast()} on the collected consensus FASTA.}
\item{unite_db, unite_db_fasta}{BLAST database prefix or FASTA passed to \code{annotate_consensus_blast()}. When both are \code{NULL}, UNITE annotation is skipped with a warning.}
\item{unite_dir}{Directory under \code{output_dir} for detailed UNITE BLAST results.}
\item{unite_top_hits_name}{Filename copied to \code{output_dir} for the top-hit summary.}
\item{unite_threads, unite_max_target_seqs, unite_evalue, unite_task, unite_word_size, unite_strand, unite_dust, unite_perc_identity, unite_extra_args, unite_species_identity, unite_genus_identity, unite_family_identity, unite_min_query_coverage, unite_min_reference_coverage, unite_novel_identity, unite_blastn, unite_makeblastdb, unite_conda_env, conda}{Parameters passed to \code{annotate_consensus_blast()} for UNITE annotation.}
\item{readme_name, chinese_readme_name}{README filenames written under \code{output_dir}. Set \code{chinese_readme_name = NULL} to skip the Chinese README.}
\item{overwrite}{Logical. If \code{TRUE}, replace an existing \code{output_dir}.}
\item{dry_run}{Logical. If \code{TRUE}, return a plan without copying files or running external tools.}
\item{echo}{Logical. If \code{TRUE}, print commands from wrapped runners.}
\item{wait, stdout, stderr}{Passed to wrapped command runners.}
\item{...}{Additional arguments passed to \code{make_consensus_delivery()} when a new consensus delivery must be generated.}
}
\value{
Invisibly returns a list describing generated and copied outputs.
}
\description{
\code{make_ITS_delivery()} combines amplicon consensus delivery results,
wf-16s ITS taxonomic profiling results, and optional consensus-vs-UNITE BLAST
review results into one customer-facing delivery folder.
}
\details{
The final delivery puts per-barcode results under \code{samples/barcode*/}.
Consensus results are copied to \code{samples/barcode*/consensus_results/};
wf-16s alignment tables are copied to
\code{samples/barcode*/identification_tables/}. The abundance table and
figures are kept at the delivery root. When UNITE annotation is enabled, the
detailed BLAST table remains under \code{unite_consensus_annotation/}, and the
top-hit table is also copied to \code{unite_consensus_top_hits.tsv} at the root.
}
