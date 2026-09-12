% Please edit documentation in R/make_ITS_delivery.R
\name{write_ITS_delivery_readme}
\alias{write_ITS_delivery_readme}
\title{Write README files for an integrated ITS delivery folder}
\usage{
write_ITS_delivery_readme(
  output_dir,
  samples_dir = "samples",
  ITS_dir = "ITS_results",
  abundance_table = "abundance_table_genus.tsv",
  figure_dir = "figures",
  identification_dir = "identification_tables",
  unite_dir = "unite_consensus_annotation",
  unite_top_hits_name = "unite_consensus_top_hits.tsv",
  has_unite_annotation = TRUE,
  readme_name = "README.txt",
  chinese_readme_name = "README.zh-CN.txt"
)
}
\arguments{
\item{output_dir}{Final ITS delivery directory.}
\item{samples_dir, ITS_dir, abundance_table, figure_dir, identification_dir, unite_dir, unite_top_hits_name}{Directory and file names used in the delivery layout. \code{ITS_dir} is the ITS-specific subdirectory below each barcode sample directory.}
\item{has_unite_annotation}{Logical. Whether UNITE annotation files were generated.}
\item{readme_name, chinese_readme_name}{README filenames. Set \code{chinese_readme_name = NULL} to skip the Chinese README.}
}
\value{
Invisibly returns generated README file paths.
}
\description{
Writes English and optional Chinese README files describing an integrated ITS
delivery folder that contains consensus results, wf-16s ITS profiling results,
and optional UNITE consensus review results.
}
