#' Plot per-base sequencing coverage from a BAM file
#'
#' `plot_coverage()` runs `samtools depth` on a coordinate-aligned BAM file,
#' calculates the percentage of reference positions meeting each requested
#' depth threshold, and returns a publication-ready coverage plot.
#'
#' @param bam_file Path to an input BAM file. The BAM should be coordinate
#'   sorted. A neighboring `.bai` or `.csi` index is normally required by
#'   `samtools`; the function warns when neither standard index is found.
#' @param min_base_quality Non-negative integer passed to `samtools depth -q`.
#'   Only bases with base quality greater than or equal to this value are
#'   counted. Note that `samtools depth` uses lowercase `-q` for base quality.
#' @param min_mapping_quality Non-negative integer passed to `samtools depth
#'   -Q`. Only reads with mapping quality greater than or equal to this value
#'   are counted.
#' @param depth_thresholds Non-negative integer vector of coverage thresholds
#'   summarized in the plot annotation. A position is counted when its depth is
#'   greater than or equal to the threshold.
#' @param major_breaks,minor_breaks Positive numeric spacing, in reference
#'   bases, between major and minor x-axis grid lines.
#' @param facet Logical. If `TRUE`, BAMs containing multiple reference
#'   sequences are plotted in separate panels with independent x-axis ranges.
#' @param samtools Command name or executable path for `samtools`.
#' @param conda_env Optional conda environment name. When supplied, `samtools`
#'   is run through `conda run -n <conda_env>`.
#' @param conda Conda executable name or path used with `conda_env`.
#' @param echo Logical. If `TRUE`, print the command before execution.
#' @param stderr Passed to [system2()]. The default `""` streams standard error
#'   to the R console.
#' @param ... Compatibility arguments from the original function:
#'   `bamFile`, `minSeqQ`, `minMapQ`, `depth`, and `breaks`. New code should use
#'   the snake_case arguments above.
#'
#' @return A [ggplot2::ggplot()] object, or `NULL` with a warning if `samtools`
#'   returns no coverage rows. The plot's `data` member contains the per-base
#'   table with columns `reference`, `position`, and `depth`. The following
#'   attributes are also attached:
#'
#'   * `depth_summary`: a data frame containing each threshold, the number of
#'     positions meeting it, and the percentage of all reported positions.
#'   * `command`: the shell-readable `samtools depth` command.
#'
#' @details
#' The command uses `samtools depth -aa -H`, so zero-depth positions and
#' reference sequences with no aligned reads are included when they are present
#' in the BAM header. Consequently, threshold percentages use all reported
#' reference positions as their denominator, including positions with depth
#' zero.
#'
#' By default, `samtools depth` excludes unmapped, secondary, QC-fail, and
#' duplicate alignments. Deletions in read alignments are not counted as depth.
#' These are properties of `samtools depth`, not additional filtering performed
#' by ONTools.
#'
#' The returned values are coverage counts, not A/C/G/T allele counts. Use
#' `samtools mpileup` when per-base allele support is required.
#'
#' @examples
#' \dontrun{
#' p <- plot_coverage(
#'   bam_file = "alignments/sample.bam",
#'   min_base_quality = 10,
#'   min_mapping_quality = 20,
#'   depth_thresholds = c(10, 20, 50, 100),
#'   major_breaks = 500,
#'   minor_breaks = 50
#' )
#' p
#' attr(p, "depth_summary")
#' }
#'
#' @importFrom rlang .data
#' @export
plot_coverage <- function(bam_file = NULL,
                          min_base_quality = 0,
                          min_mapping_quality = 0,
                          depth_thresholds = c(10, 20, 50, 100, 200, 500),
                          major_breaks = 500,
                          minor_breaks = 50,
                          facet = TRUE,
                          samtools = "samtools",
                          conda_env = NULL,
                          conda = "conda",
                          echo = TRUE,
                          stderr = "",
                          ...) {
  bam_file_missing <- missing(bam_file)
  min_base_quality_missing <- missing(min_base_quality)
  min_mapping_quality_missing <- missing(min_mapping_quality)
  depth_thresholds_missing <- missing(depth_thresholds)
  major_breaks_missing <- missing(major_breaks)
  dots <- list(...)
  plot_coverage_check_legacy_args(dots)

  if (!is.null(dots$bamFile)) {
    if (!bam_file_missing) {
      stop("Supply only one of `bam_file` and `bamFile`.", call. = FALSE)
    }
    bam_file <- dots$bamFile
  }
  if (!is.null(dots$minSeqQ)) {
    if (!min_base_quality_missing) {
      stop("Supply only one of `min_base_quality` and `minSeqQ`.", call. = FALSE)
    }
    min_base_quality <- dots$minSeqQ
  }
  if (!is.null(dots$minMapQ)) {
    if (!min_mapping_quality_missing) {
      stop("Supply only one of `min_mapping_quality` and `minMapQ`.", call. = FALSE)
    }
    min_mapping_quality <- dots$minMapQ
  }
  if (!is.null(dots$depth)) {
    if (!depth_thresholds_missing) {
      stop("Supply only one of `depth_thresholds` and `depth`.", call. = FALSE)
    }
    depth_thresholds <- dots$depth
  }
  if (!is.null(dots$breaks)) {
    if (!major_breaks_missing) {
      stop("Supply only one of `major_breaks` and `breaks`.", call. = FALSE)
    }
    major_breaks <- dots$breaks
  }

  if (is.null(bam_file)) {
    stop("`bam_file` is required.", call. = FALSE)
  }
  check_file_arg(bam_file, "bam_file")
  check_scalar_character(samtools, "samtools")
  check_scalar_character(conda, "conda")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  check_logical_scalar(facet, "facet")
  check_logical_scalar(echo, "echo")

  min_base_quality <- validate_nonnegative_integer(
    min_base_quality,
    "min_base_quality"
  )
  min_mapping_quality <- validate_nonnegative_integer(
    min_mapping_quality,
    "min_mapping_quality"
  )
  depth_thresholds <- plot_coverage_validate_thresholds(depth_thresholds)
  major_breaks <- validate_positive_number(major_breaks, "major_breaks")
  minor_breaks <- validate_positive_number(minor_breaks, "minor_breaks")

  bam_file <- normalizePath(bam_file, mustWork = TRUE)
  if (!has_standard_bam_index(bam_file)) {
    warning(
      "No standard BAM index found next to `bam_file`; samtools depth may fail. ",
      "Create one with `samtools index` if needed.",
      call. = FALSE
    )
  }

  depth_call <- dehost_fastq_external_call(
    command = samtools,
    args = c(
      "depth", "-aa", "-H",
      "-q", as.character(min_base_quality),
      "-Q", as.character(min_mapping_quality),
      bam_file
    ),
    conda_env = conda_env,
    conda = conda
  )
  command <- paste(
    c(shQuote(depth_call$command), shQuote(depth_call$args)),
    collapse = " "
  )
  if (isTRUE(echo)) message(command)

  if (is.null(conda_env)) {
    require_external_command(samtools)
  } else {
    require_external_command(conda)
  }
  output <- suppressWarnings(system2(
    depth_call$command,
    args = depth_call$args,
    stdout = TRUE,
    stderr = stderr
  ))
  status <- attr(output, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("samtools depth failed with exit status: ", status, call. = FALSE)
  }

  coverage <- parse_samtools_depth_output(output)
  if (nrow(coverage) == 0L) {
    warning("No coverage data were generated.", call. = FALSE)
    return(NULL)
  }

  depth_summary <- data.frame(
    threshold = depth_thresholds,
    positions = vapply(
      depth_thresholds,
      function(value) sum(coverage$depth >= value),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  depth_summary$percent <- round(
    depth_summary$positions / nrow(coverage) * 100,
    digits = 1
  )
  annotation <- paste0(
    ">=", depth_summary$threshold, "X: ", depth_summary$percent, "%",
    collapse = "\n"
  )

  max_position <- max(coverage$position)
  max_depth <- max(coverage$depth)
  y_limit <- if (max_depth > 0L) max_depth else 1
  title <- sub("[.]bam$", "", basename(bam_file), ignore.case = TRUE)

  plot <- ggplot2::ggplot(
    coverage,
    ggplot2::aes(x = .data$position, y = .data$depth)
  ) +
    ggplot2::geom_col(fill = "lightblue", width = 1) +
    ggplot2::coord_cartesian(ylim = c(0, y_limit), clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max_position, by = major_breaks),
      minor_breaks = seq(0, max_position, by = minor_breaks)
    ) +
    ggplot2::labs(
      title = title,
      x = "Reference position (bp)",
      y = "Coverage depth"
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.05,
      vjust = 1.1
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "grey75",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_line(
        colour = "grey90",
        linewidth = 0.25
      ),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (isTRUE(facet) && length(unique(coverage$reference)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(
      ggplot2::vars(.data$reference),
      scales = "free_x",
      ncol = 1
    )
  }

  attr(plot, "depth_summary") <- depth_summary
  attr(plot, "command") <- command
  plot
}

plot_coverage_check_legacy_args <- function(dots) {
  if (length(dots) == 0L) return(invisible(NULL))
  if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
    stop("All arguments in `...` must be named.", call. = FALSE)
  }
  unknown <- setdiff(names(dots), c("bamFile", "minSeqQ", "minMapQ", "depth", "breaks"))
  if (length(unknown) > 0L) {
    stop("Unused argument `", unknown[[1L]], "`.", call. = FALSE)
  }
  invisible(NULL)
}

plot_coverage_validate_thresholds <- function(x) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0) || any(x != as.integer(x))) {
    stop("`depth_thresholds` must be a non-empty vector of non-negative integers.",
         call. = FALSE)
  }
  sort(unique(as.integer(x)))
}

parse_samtools_depth_output <- function(output) {
  rows <- output[nzchar(output) & !startsWith(output, "#")]
  if (length(rows) == 0L) {
    return(data.frame(
      reference = character(),
      position = integer(),
      depth = integer(),
      stringsAsFactors = FALSE
    ))
  }

  fields <- strsplit(rows, "\t", fixed = TRUE)
  if (any(lengths(fields) != 3L)) {
    stop("Unexpected output from samtools depth; expected three tab-separated columns.",
         call. = FALSE)
  }
  coverage <- data.frame(
    reference = vapply(fields, `[[`, character(1), 1L),
    position = suppressWarnings(as.integer(vapply(fields, `[[`, character(1), 2L))),
    depth = suppressWarnings(as.integer(vapply(fields, `[[`, character(1), 3L))),
    stringsAsFactors = FALSE
  )
  if (anyNA(coverage$position) || anyNA(coverage$depth) ||
      any(coverage$position < 1L) || any(coverage$depth < 0L)) {
    stop("samtools depth returned invalid position or depth values.", call. = FALSE)
  }
  coverage
}
