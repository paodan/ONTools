#' Read a PAF file into a data frame
#'
#' `read_paf()` reads a plain or gzipped PAF file and returns a data frame with
#' standard PAF column names. Optional SAM-style tags can be kept as raw fields
#' or parsed into separate columns.
#'
#' @param paf_file Path to a `.paf` or `.paf.gz` file.
#' @param parse_tags Logical. If `TRUE`, split optional `TAG:TYPE:VALUE` fields
#'   into separate columns named by `TAG`.
#' @param simplify Logical. If `TRUE`, convert standard numeric PAF columns and
#'   parsed integer/float tags to numeric columns when possible.
#'
#' @return A data frame, or `NULL` for empty PAF files. Standard PAF columns are
#'   named `query_name`, `query_length`, `query_start`, `query_end`, `strand`,
#'   `target_name`, `target_length`, `target_start`, `target_end`,
#'   `residue_matches`, `alignment_block_length`, and `mapping_quality`.
#'
#' @examples
#' paf <- tempfile(fileext = ".paf")
#' writeLines(
#'   "read1\t1000\t0\t900\t+\tchr1\t2000\t100\t1000\t890\t900\t60\tNM:i:10",
#'   paf
#' )
#'
#' read_paf(paf)
#'
#' @export
read_paf <- function(paf_file, parse_tags = TRUE, simplify = TRUE) {
  check_file_arg(paf_file, "paf_file")
  check_logical_scalar(parse_tags, "parse_tags")
  check_logical_scalar(simplify, "simplify")

  if (file.info(paf_file)$size == 0L) {
    return(NULL)
  }

  paf <- utils::read.table(
    paf_file,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )

  if (nrow(paf) == 0L) {
    return(NULL)
  }
  if (ncol(paf) < 12L) {
    stop("PAF file must contain at least 12 tab-delimited columns.",
         call. = FALSE)
  }

  standard_names <- c(
    "query_name", "query_length", "query_start", "query_end", "strand",
    "target_name", "target_length", "target_start", "target_end",
    "residue_matches", "alignment_block_length", "mapping_quality"
  )
  names(paf)[seq_along(standard_names)] <- standard_names

  tag_cols <- character()
  if (ncol(paf) > length(standard_names)) {
    tag_cols <- names(paf)[seq.int(length(standard_names) + 1L, ncol(paf))]
    names(paf)[match(tag_cols, names(paf))] <- paste0("tag_", seq_along(tag_cols))
    tag_cols <- paste0("tag_", seq_along(tag_cols))
  }

  if (isTRUE(simplify)) {
    numeric_cols <- c(
      "query_length", "query_start", "query_end", "target_length",
      "target_start", "target_end", "residue_matches",
      "alignment_block_length", "mapping_quality"
    )
    paf[numeric_cols] <- lapply(paf[numeric_cols], function(x) {
      suppressWarnings(as.numeric(x))
    })
  }

  if (isTRUE(parse_tags) && length(tag_cols) > 0L) {
    tags <- parse_paf_tags(paf[tag_cols], simplify = simplify)
    paf <- paf[setdiff(names(paf), tag_cols)]
    if (ncol(tags) > 0L) {
      paf <- cbind(paf, tags, stringsAsFactors = FALSE)
    }
  }

  paf
}

parse_paf_tags <- function(tag_data, simplify = TRUE) {
  rows <- lapply(seq_len(nrow(tag_data)), function(i) {
    fields <- unlist(tag_data[i, , drop = TRUE], use.names = FALSE)
    fields <- fields[nzchar(fields)]
    if (length(fields) == 0L) {
      return(list())
    }

    row <- list()
    for (field in fields) {
      parts <- strsplit(field, ":", fixed = TRUE)[[1L]]
      if (length(parts) < 3L || !nzchar(parts[[1L]])) {
        next
      }
      tag <- make.names(parts[[1L]])
      value <- paste(parts[-c(1L, 2L)], collapse = ":")
      row[[tag]] <- convert_paf_tag_value(value, parts[[2L]], simplify)
    }
    row
  })

  tag_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  if (length(tag_names) == 0L) {
    return(data.frame(row.names = seq_len(nrow(tag_data))))
  }

  out <- lapply(tag_names, function(tag) {
    values <- lapply(rows, function(row) {
      if (tag %in% names(row)) row[[tag]] else NA
    })
    unlist(values, use.names = FALSE)
  })
  names(out) <- tag_names
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

convert_paf_tag_value <- function(value, type, simplify = TRUE) {
  if (!isTRUE(simplify)) {
    return(value)
  }
  if (type %in% c("i", "f")) {
    return(suppressWarnings(as.numeric(value)))
  }
  value
}
