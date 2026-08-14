#' Get primer sequences from a sample information table
#'
#' `get_primers()` extracts the forward and reverse primer sequences for a
#' project and/or expected amplicon size from a sample information table. It is
#' intended to pair with [trim_fasta_keep_primers()].
#'
#' @param sample_info Sample information table, either a data frame or a CSV
#'   file path.
#' @param size Optional expected amplicon size used to filter
#'   `amplicon_size_col`. Matching is performed after converting both values to
#'   character.
#' @param project Optional project label used to filter `project_col`. If
#'   supplied, `project_col` must also be supplied.
#' @param project_col Optional column in `sample_info` identifying the project.
#' @param amplicon_size_col Column in `sample_info` identifying expected
#'   amplicon size or length group.
#' @param f_primer_col Column containing the forward primer sequence.
#' @param r_primer_col Column containing the reverse primer sequence.
#' @param uppercase Logical. If `TRUE`, return primer sequences in upper case.
#' @param allow_multiple Logical. If `FALSE`, stop when the filtered table
#'   contains more than one unique primer pair.
#'
#' @return A list with `f_primer`, `r_primer`, `primer_table`, `n_rows`, and
#'   `filters`. For compatibility with the default column names, the list also
#'   contains `Primer_Sequence_5` and `Primer_Sequence_3`.
#'
#' @examples
#' sample_info <- data.frame(
#'   Project = c("P1", "P1", "P2"),
#'   Expected_Amplicon_Size_bp = c("1600bp", "1600bp", "3500bp"),
#'   Primer_Sequence_5 = c("AGAGTTTGATCMTGGCTCAG", "AGAGTTTGATCMTGGCTCAG", "AAA"),
#'   Primer_Sequence_3 = c("TACGGYTACCTTGTTACGACTT", "TACGGYTACCTTGTTACGACTT", "TTT")
#' )
#'
#' primers <- get_primers(
#'   sample_info,
#'   size = "1600bp",
#'   project = "P1",
#'   project_col = "Project"
#' )
#' primers$f_primer
#' primers$r_primer
#'
#' @export
get_primers <- function(sample_info,
                        size = NULL,
                        project = NULL,
                        project_col = NULL,
                        amplicon_size_col = "Expected_Amplicon_Size_bp",
                        f_primer_col = "Primer_Sequence_5",
                        r_primer_col = "Primer_Sequence_3",
                        uppercase = TRUE,
                        allow_multiple = FALSE) {
  check_scalar_character(amplicon_size_col, "amplicon_size_col")
  check_scalar_character(f_primer_col, "f_primer_col")
  check_scalar_character(r_primer_col, "r_primer_col")
  check_logical_scalar(uppercase, "uppercase")
  check_logical_scalar(allow_multiple, "allow_multiple")

  if (!is.null(project_col)) {
    check_scalar_character(project_col, "project_col")
  }
  if (!is.null(size)) {
    size <- check_scalar_filter_value(size, "size")
  }
  if (!is.null(project)) {
    project <- check_scalar_filter_value(project, "project")
    if (is.null(project_col)) {
      stop("`project_col` must be supplied when `project` is supplied.",
           call. = FALSE)
    }
  }

  sample_info <- read_sample_info_table(sample_info)
  required_cols <- c(f_primer_col, r_primer_col)
  if (!is.null(size)) {
    required_cols <- c(required_cols, amplicon_size_col)
  }
  if (!is.null(project)) {
    required_cols <- c(required_cols, project_col)
  }

  missing_cols <- setdiff(required_cols, names(sample_info))
  if (length(missing_cols) > 0L) {
    stop(
      "`sample_info` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  selected <- rep(TRUE, nrow(sample_info))
  if (!is.null(size)) {
    selected <- selected &
      as.character(sample_info[[amplicon_size_col]]) == size
  }
  if (!is.null(project)) {
    selected <- selected &
      as.character(sample_info[[project_col]]) == project
  }

  filtered <- sample_info[selected, , drop = FALSE]
  if (nrow(filtered) == 0L) {
    stop("No rows in `sample_info` match the requested filters.", call. = FALSE)
  }

  primer_table <- unique(data.frame(
    f_primer = trimws(as.character(filtered[[f_primer_col]])),
    r_primer = trimws(as.character(filtered[[r_primer_col]])),
    stringsAsFactors = FALSE
  ))

  complete_pair <- !is.na(primer_table$f_primer) &
    nzchar(primer_table$f_primer) &
    !is.na(primer_table$r_primer) &
    nzchar(primer_table$r_primer)
  if (any(!complete_pair)) {
    warning(
      "Dropping ",
      sum(!complete_pair),
      " row(s) with missing primer sequences.",
      call. = FALSE
    )
    primer_table <- primer_table[complete_pair, , drop = FALSE]
  }

  if (nrow(primer_table) == 0L) {
    stop("No complete primer pairs were found after filtering.", call. = FALSE)
  }

  if (isTRUE(uppercase)) {
    primer_table$f_primer <- toupper(primer_table$f_primer)
    primer_table$r_primer <- toupper(primer_table$r_primer)
    primer_table <- unique(primer_table)
  }

  if (nrow(primer_table) > 1L && !isTRUE(allow_multiple)) {
    stop(
      "More than one unique primer pair was found. ",
      "Use narrower filters or set `allow_multiple = TRUE`.",
      call. = FALSE
    )
  }

  filters <- list(
    size = size,
    project = project,
    project_col = project_col,
    amplicon_size_col = amplicon_size_col
  )

  list(
    f_primer = primer_table$f_primer,
    r_primer = primer_table$r_primer,
    Primer_Sequence_5 = primer_table$f_primer,
    Primer_Sequence_3 = primer_table$r_primer,
    primer_table = primer_table,
    n_rows = nrow(filtered),
    filters = filters
  )
}

check_scalar_filter_value <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) {
    stop("`", name, "` must be a single non-missing value.", call. = FALSE)
  }

  trimws(as.character(x))
}
