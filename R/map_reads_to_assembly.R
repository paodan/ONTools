#' Map reads back to an assembly and plot read depth
#'
#' `map_reads_to_assembly()` maps Nanopore reads back to assembled contigs with
#' `minimap2`, sorts and indexes the BAM with `samtools`, writes per-base depth
#' with `samtools depth`, and optionally saves a depth plot.
#'
#' @param assembly Assembly FASTA file used as the mapping reference.
#' @param reads Input reads in FASTQ or FASTQ.GZ format.
#' @param align_bam Output coordinate-sorted BAM file.
#' @param depth_file Output depth table. If `NULL`, writes
#'   `<align_bam without .bam>.depth.txt`.
#' @param depth_plot Output depth plot image. If `NULL` and `plot_depth = TRUE`,
#'   writes `<depth_file>.png`.
#' @param preset Minimap2 preset passed to `minimap2 -x`. Default `"map-ont"`
#'   matches Nanopore reads.
#' @param threads Positive integer thread count passed to `minimap2 -t` and
#'   `samtools sort -@`.
#' @param secondary Logical. If `FALSE`, pass `--secondary=no` to minimap2.
#' @param depth_all_positions Logical. If `TRUE`, pass `samtools depth -a` so
#'   zero-depth positions are also reported.
#' @param depth_all_references Logical. If `TRUE`, pass `samtools depth -aa` so
#'   every position in every reference sequence is reported. This overrides
#'   `depth_all_positions`.
#' @param min_mapping_quality Optional minimum mapping quality passed to
#'   `samtools depth -Q`.
#' @param min_base_quality Optional minimum base quality passed to
#'   `samtools depth -q`.
#' @param plot_depth Logical. If `TRUE`, read the depth table and save a depth
#'   plot with ggplot2.
#' @param plot_width,plot_height Width and height in inches passed to
#'   [ggplot2::ggsave()] when `plot_depth = TRUE`.
#' @param facet_nrow Number of rows passed to [ggplot2::facet_wrap()] for contig
#'   panels in the depth plot.
#' @param minimap2,samtools Command names or executable paths.
#' @param conda_env Optional conda environment name. If supplied, external
#'   commands are run with `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used when `conda_env` is supplied.
#' @param dry_run Logical. If `TRUE`, return planned commands without running.
#' @param echo Logical. If `TRUE`, print planned commands before execution.
#' @param stderr Passed to [system2()]. Defaults stream errors to the R console.
#'
#' @return Invisibly returns a list with `status`, `commands`, `paths`, `plot`,
#'   `preset`, and `conda_env`.
#'
#' @examples
#' assembly <- tempfile(fileext = ".fasta")
#' reads <- tempfile(fileext = ".fastq")
#' bam <- tempfile(fileext = ".bam")
#' writeLines(c(">contig1", "ACGTACGT"), assembly)
#' writeLines(c("@read1", "ACGT", "+", "!!!!"), reads)
#' res <- map_reads_to_assembly(assembly, reads, bam, dry_run = TRUE)
#' res$commands$minimap2
#'
#' @export
map_reads_to_assembly <- function(assembly,
                                  reads,
                                  align_bam,
                                  depth_file = NULL,
                                  depth_plot = NULL,
                                  preset = "map-ont",
                                  threads = 16,
                                  secondary = TRUE,
                                  depth_all_positions = FALSE,
                                  depth_all_references = FALSE,
                                  min_mapping_quality = NULL,
                                  min_base_quality = NULL,
                                  plot_depth = TRUE,
                                  plot_width = 8,
                                  plot_height = 5,
                                  facet_nrow = 3,
                                  minimap2 = "minimap2",
                                  samtools = "samtools",
                                  conda_env = NULL,
                                  conda = "conda",
                                  dry_run = FALSE,
                                  echo = TRUE,
                                  stderr = "") {
  check_file_arg(assembly, "assembly")
  check_file_arg(reads, "reads")
  check_scalar_character(align_bam, "align_bam")
  check_scalar_character(preset, "preset")
  check_scalar_character(minimap2, "minimap2")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  check_logical_scalar(secondary, "secondary")
  check_logical_scalar(depth_all_positions, "depth_all_positions")
  check_logical_scalar(depth_all_references, "depth_all_references")
  check_logical_scalar(plot_depth, "plot_depth")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")

  if (!is.null(depth_file)) {
    check_scalar_character(depth_file, "depth_file")
  }
  if (!is.null(depth_plot)) {
    check_scalar_character(depth_plot, "depth_plot")
  }
  if (!is.null(conda_env)) {
    check_scalar_character(conda_env, "conda_env")
  }

  threads <- validate_positive_integer(threads, "threads")
  min_mapping_quality <- validate_optional_nonnegative_integer(
    min_mapping_quality,
    "min_mapping_quality"
  )
  min_base_quality <- validate_optional_nonnegative_integer(
    min_base_quality,
    "min_base_quality"
  )
  plot_width <- validate_positive_number(plot_width, "plot_width")
  plot_height <- validate_positive_number(plot_height, "plot_height")
  facet_nrow <- validate_positive_integer(facet_nrow, "facet_nrow")

  assembly <- normalizePath(assembly, mustWork = TRUE)
  reads <- normalizePath(reads, mustWork = TRUE)
  dir.create(dirname(align_bam), recursive = TRUE, showWarnings = FALSE)
  align_bam <- normalizePath(align_bam, mustWork = FALSE)
  align_bai <- paste0(align_bam, ".bai")
  sam <- tempfile(fileext = ".sam")

  if (is.null(depth_file)) {
    depth_file <- paste0(sub("[.]bam$", "", align_bam, ignore.case = TRUE), ".depth.txt")
  } else {
    dir.create(dirname(depth_file), recursive = TRUE, showWarnings = FALSE)
  }
  depth_file <- normalizePath(depth_file, mustWork = FALSE)

  if (isTRUE(plot_depth)) {
    if (is.null(depth_plot)) {
      depth_plot <- paste0(depth_file, ".png")
    } else {
      dir.create(dirname(depth_plot), recursive = TRUE, showWarnings = FALSE)
    }
    depth_plot <- normalizePath(depth_plot, mustWork = FALSE)
  } else {
    depth_plot <- NULL
  }

  minimap2_args <- c(
    "-a",
    "-x", preset,
    "-t", as.character(threads)
  )
  if (!isTRUE(secondary)) {
    minimap2_args <- c(minimap2_args, "--secondary=no")
  }
  minimap2_args <- c(minimap2_args, assembly, reads)

  sort_args <- c(
    "sort",
    "-@", as.character(threads),
    "-o", align_bam,
    sam
  )
  index_args <- c("index", align_bam)
  depth_args <- build_samtools_depth_args(
    bam = align_bam,
    depth_all_positions = depth_all_positions,
    depth_all_references = depth_all_references,
    min_mapping_quality = min_mapping_quality,
    min_base_quality = min_base_quality
  )

  minimap2_call <- dehost_fastq_external_call(
    command = minimap2,
    args = minimap2_args,
    conda_env = conda_env,
    conda = conda
  )
  sort_call <- dehost_fastq_external_call(
    command = samtools,
    args = sort_args,
    conda_env = conda_env,
    conda = conda
  )
  index_call <- dehost_fastq_external_call(
    command = samtools,
    args = index_args,
    conda_env = conda_env,
    conda = conda
  )
  depth_call <- dehost_fastq_external_call(
    command = samtools,
    args = depth_args,
    conda_env = conda_env,
    conda = conda
  )

  commands <- list(
    minimap2 = paste(c(shQuote(minimap2_call$command), shQuote(minimap2_call$args), ">", shQuote(sam)), collapse = " "),
    samtools_sort = paste(c(shQuote(sort_call$command), shQuote(sort_call$args)), collapse = " "),
    samtools_index = paste(c(shQuote(index_call$command), shQuote(index_call$args)), collapse = " "),
    samtools_depth = paste(c(shQuote(depth_call$command), shQuote(depth_call$args), ">", shQuote(depth_file)), collapse = " ")
  )
  if (isTRUE(plot_depth)) {
    commands$plot_depth <- paste("Plot depth in R:", shQuote(depth_file), ">", shQuote(depth_plot))
  }

  if (isTRUE(echo)) {
    message(commands$minimap2)
    message(commands$samtools_sort)
    message(commands$samtools_index)
    message(commands$samtools_depth)
    if (isTRUE(plot_depth)) message(commands$plot_depth)
  }

  paths <- list(
    assembly = assembly,
    reads = reads,
    sam = sam,
    bam = align_bam,
    bai = align_bai,
    depth = depth_file,
    depth_plot = depth_plot
  )

  if (isTRUE(dry_run)) {
    return(invisible(list(
      status = NA_integer_,
      commands = commands,
      paths = paths,
      plot = NULL,
      preset = preset,
      conda_env = conda_env
    )))
  }

  if (is.null(conda_env)) {
    require_external_command(minimap2)
    require_external_command(samtools)
  } else {
    require_external_command(conda)
  }

  on.exit(unlink(sam), add = TRUE)

  minimap2_status <- system2(
    minimap2_call$command,
    args = minimap2_call$args,
    stdout = sam,
    stderr = stderr
  )
  if (!identical(minimap2_status, 0L)) {
    stop("minimap2 failed with exit status: ", minimap2_status, call. = FALSE)
  }

  sort_status <- system2(
    sort_call$command,
    args = sort_call$args,
    stderr = stderr
  )
  if (!identical(sort_status, 0L)) {
    stop("samtools sort failed with exit status: ", sort_status, call. = FALSE)
  }

  index_status <- system2(
    index_call$command,
    args = index_call$args,
    stderr = stderr
  )
  if (!identical(index_status, 0L)) {
    stop("samtools index failed with exit status: ", index_status, call. = FALSE)
  }

  depth_status <- system2(
    depth_call$command,
    args = depth_call$args,
    stdout = depth_file,
    stderr = stderr
  )
  if (!identical(depth_status, 0L)) {
    stop("samtools depth failed with exit status: ", depth_status, call. = FALSE)
  }

  plot <- NULL
  if (isTRUE(plot_depth)) {
    plot <- plot_read_depth(
      depth_file = depth_file,
      depth_plot = depth_plot,
      width = plot_width,
      height = plot_height,
      facet_nrow = facet_nrow
    )
  }

  invisible(list(
    status = 0L,
    commands = commands,
    paths = paths,
    plot = plot,
    preset = preset,
    conda_env = conda_env
  ))
}

build_samtools_depth_args <- function(bam,
                                      depth_all_positions,
                                      depth_all_references,
                                      min_mapping_quality,
                                      min_base_quality) {
  args <- "depth"
  if (isTRUE(depth_all_references)) {
    args <- c(args, "-aa")
  } else if (isTRUE(depth_all_positions)) {
    args <- c(args, "-a")
  }
  if (!is.null(min_mapping_quality)) {
    args <- c(args, "-Q", as.character(min_mapping_quality))
  }
  if (!is.null(min_base_quality)) {
    args <- c(args, "-q", as.character(min_base_quality))
  }

  c(args, bam)
}

plot_read_depth <- function(depth_file, depth_plot, width, height, facet_nrow) {
  depth <- utils::read.table(
    depth_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("contig", "position", "depth")
  )

  plot <- ggplot2::ggplot(depth, ggplot2::aes(.data$position, .data$depth)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$contig), nrow = facet_nrow) +
    ggplot2::xlab("Position") +
    ggplot2::ylab("Depth") +
    ggplot2::theme_bw()

  ggplot2::ggsave(filename = depth_plot, plot = plot, width = width, height = height)
  plot
}

validate_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a single positive number.", call. = FALSE)
  }

  x
}
