#' Split multiple FASTQ files by expected and detected read-length ranges
#'
#' `split_fastqs()` is the batch counterpart of [split_fastq()]. It validates a
#' vector of input files, plans their output directories, prevents predictable
#' output-name collisions, and calls [split_fastq()] once for every input. All
#' FASTQ parsing, peak detection, and `seqkit` filtering remain in the
#' single-file function.
#'
#' @param fastqs Non-empty character vector containing paths to `.fastq` or
#'   `.fastq.gz` input files.
#' @param out_fastq_path Output-directory specification. Use `NULL` to place
#'   each input's outputs in a `peaks` directory beside that input. Supply one
#'   path to place every output in a shared directory, or supply one path per
#'   input file to use separate output directories.
#' @inheritParams split_fastq
#'
#' @return Invisibly returns a named list with one element per input file. Each
#' non-`NULL` element is the complete result returned by [split_fastq()], with
#' components such as `output_fastq`, `peak_counts`, `filter_results`, `plot`,
#' and `read_count`. An empty input FASTQ produces a warning and a `NULL`
#' element without stopping later files. List names are unique versions of the
#' input base names.
#'
#' @details
#' Input validation and collision checks happen before the first file is split,
#' reducing the chance of a batch stopping after partially writing outputs.
#' Files with the same output stem and compression suffix cannot share an output
#' directory because they would generate identical `_peak_...` names. Choose
#' separate `out_fastq_path` values for such inputs.
#'
#' Errors other than an empty FASTQ stop the batch. This includes malformed
#' FASTQ records, invalid parameters, unavailable external commands, and failed
#' `seqkit` execution. Use `dry_run = TRUE` to validate inputs, calculate peaks,
#' inspect every planned command, and detect output collisions without running
#' `seqkit`.
#'
#' All filtering parameters are forwarded unchanged to [split_fastq()]. In
#' particular, `.fastq.gz` inputs produce genuinely gzip-compressed outputs
#' because [filter_fastq()] writes through `seqkit -o`.
#'
#' @examples
#' fastq_dir <- tempfile("split-fastqs-input-")
#' output_dir <- tempfile("split-fastqs-output-")
#' dir.create(fastq_dir)
#'
#' write_example_fastq <- function(path, lengths) {
#'   records <- unlist(lapply(seq_along(lengths), function(i) {
#'     sequence <- paste(rep("A", lengths[[i]]), collapse = "")
#'     quality <- paste(rep("!", lengths[[i]]), collapse = "")
#'     c(paste0("@read", i), sequence, "+", quality)
#'   }), use.names = FALSE)
#'   writeLines(records, path)
#' }
#'
#' fastqs <- file.path(fastq_dir, c("sample1.fastq", "sample2.fastq"))
#' write_example_fastq(fastqs[[1L]], c(4, 5, 6, 8))
#' write_example_fastq(fastqs[[2L]], c(4, 5, 5, 7))
#'
#' plans <- split_fastqs(
#'   fastqs = fastqs,
#'   out_fastq_path = output_dir,
#'   expected_size = 5,
#'   min_read_length = 4,
#'   max_read_length = 6,
#'   by_peak = FALSE,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' names(plans)
#' plans[[1L]]$peak_counts
#' plans[[1L]]$filter_results[[1L]]$command_string
#'
#' # One output directory per input is also supported.
#' separate_dirs <- file.path(output_dir, c("sample1", "sample2"))
#' separate_plans <- split_fastqs(
#'   fastqs,
#'   out_fastq_path = separate_dirs,
#'   expected_size = 5,
#'   min_read_length = 4,
#'   max_read_length = 6,
#'   by_peak = FALSE,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' vapply(separate_plans, `[[`, character(1), "output_dir")
#'
#' @seealso [split_fastq()], [read_fastq_lengths()], [filter_fastq()]
#' @export
split_fastqs <- function(fastqs,
                         out_fastq_path = NULL,
                         expected_size = 1600,
                         min_read_length = 1200,
                         max_read_length = 1800,
                         by_peak = TRUE,
                         min_relative_height = 0.05,
                         log_scale = FALSE,
                         min_distance = 400,
                         density_n = 5000,
                         x_breaks = c(1, seq(500, 3000, by = 500)),
                         x_scale = c("linear", "auto", "log10"),
                         seqkit = "seqkit",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stderr = "") {
  if (!is.character(fastqs) || length(fastqs) == 0L ||
      anyNA(fastqs) || any(!nzchar(fastqs))) {
    stop("`fastqs` must be a non-empty character vector of file paths.",
         call. = FALSE)
  }

  normalized_fastqs <- vapply(fastqs, function(fastq) {
    check_file_arg(fastq, "fastqs")
    split_fastq_suffix(fastq, name = "fastqs")
    normalizePath(fastq, mustWork = TRUE)
  }, character(1), USE.NAMES = FALSE)

  if (anyDuplicated(normalized_fastqs)) {
    duplicated_file <- normalized_fastqs[duplicated(normalized_fastqs)][[1L]]
    stop("`fastqs` contains a duplicated input file: ", duplicated_file,
         call. = FALSE)
  }

  output_dirs <- split_fastqs_output_dirs(
    fastqs = normalized_fastqs,
    out_fastq_path = out_fastq_path
  )
  split_fastqs_check_collisions(normalized_fastqs, output_dirs)

  results <- vector("list", length(normalized_fastqs))
  names(results) <- make.unique(basename(normalized_fastqs))
  for (i in seq_along(normalized_fastqs)) {
    result <- split_fastq(
      fastq_file = normalized_fastqs[[i]],
      out_fastq_path = output_dirs[[i]],
      expected_size = expected_size,
      min_read_length = min_read_length,
      max_read_length = max_read_length,
      by_peak = by_peak,
      min_relative_height = min_relative_height,
      log_scale = log_scale,
      min_distance = min_distance,
      density_n = density_n,
      x_breaks = x_breaks,
      x_scale = x_scale,
      seqkit = seqkit,
      conda_env = conda_env,
      conda = conda,
      dry_run = dry_run,
      echo = echo,
      stderr = stderr
    )
    results[i] <- list(result)
  }

  invisible(results)
}

split_fastqs_output_dirs <- function(fastqs, out_fastq_path) {
  if (is.null(out_fastq_path)) {
    return(file.path(dirname(fastqs), "peaks"))
  }
  if (!is.character(out_fastq_path) || anyNA(out_fastq_path) ||
      any(!nzchar(out_fastq_path)) ||
      !length(out_fastq_path) %in% c(1L, length(fastqs))) {
    stop(
      "`out_fastq_path` must be `NULL`, one non-empty path, or one path per ",
      "input file.",
      call. = FALSE
    )
  }
  if (length(out_fastq_path) == 1L) {
    out_fastq_path <- rep(out_fastq_path, length(fastqs))
  }
  vapply(out_fastq_path, normalizePath, character(1), mustWork = FALSE,
         USE.NAMES = FALSE)
}

split_fastqs_check_collisions <- function(fastqs, output_dirs) {
  output_keys <- vapply(seq_along(fastqs), function(i) {
    suffix <- split_fastq_suffix(fastqs[[i]])
    stem <- sub(
      "[.]fastq([.]gz)?$",
      "",
      basename(fastqs[[i]]),
      ignore.case = TRUE
    )
    file.path(output_dirs[[i]], paste0(stem, suffix))
  }, character(1))
  collision_key <- tolower(output_keys)
  if (anyDuplicated(collision_key)) {
    collision <- output_keys[duplicated(collision_key)][[1L]]
    stop(
      "Multiple inputs would create colliding outputs based on: ",
      collision,
      ". Supply separate `out_fastq_path` values.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
