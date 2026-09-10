#' Read a VCF file into a data frame
#'
#' `read_vcf()` reads plain or gzipped VCF files and returns a data frame with
#' fixed VCF columns, optional REF/ALT variant types, optional medaka allele
#' depth summaries, parsed `INFO` fields, and optionally parsed sample genotype
#' fields from `FORMAT`.
#'
#' @param vcf_file Path to a `.vcf` or `.vcf.gz` file.
#' @param parse_info Logical. If `TRUE`, split the `INFO` column into separate
#'   columns.
#' @param parse_genotypes Logical. If `TRUE`, split sample columns according to
#'   the `FORMAT` column.
#' @param sample_format Output shape for parsed sample genotype fields.
#'   `"wide"` keeps one row per VCF record; `"long"` returns one row per
#'   VCF-record/sample pair.
#' @param samples Optional character vector of sample names to keep. If `NULL`,
#'   all VCF sample columns are used.
#' @param add_variant_type Logical. If `TRUE`, add a `Type` column inferred from
#'   the `REF` and `ALT` alleles. The rule mirrors the variant type shown in the
#'   wf-amplicon report with an extra normalized-indel guard: single-base
#'   substitutions are `SNP`, unequal-length alleles with a shared anchor base
#'   are `INDEL`, equal-length multi-base substitutions are `MNP`, unequal-length
#'   alleles without a shared anchor are `OTHER`, and reference/missing calls are
#'   `REF`/`NA`.
#' @param add_allele_depth Logical. If `TRUE`, add allele-depth summary columns
#'   from medaka/wf-amplicon `SR` and `AR` INFO fields when available. `SR` is
#'   interpreted as `ref_fwd,ref_rev,alt1_fwd,alt1_rev,...`; `AR` is interpreted
#'   as ambiguous spanning reads by strand, `ambiguous_fwd,ambiguous_rev`.
#' @param keep_info Logical. If `TRUE`, keep the original `INFO` column after
#'   parsing.
#' @param keep_format Logical. If `TRUE`, keep the original `FORMAT` and raw
#'   sample genotype columns after parsing.
#' @param simplify Logical. If `TRUE`, convert parsed scalar numeric fields to
#'   numeric columns when possible.
#' @param info_prefix Prefix added to parsed `INFO` column names.
#' @param genotype_prefix Prefix added to parsed genotype field names. When
#'   `NULL`, single-sample wide output uses bare names such as `GT` and `DP`;
#'   multi-sample wide output uses `<sample>_<field>`.
#'
#' @return A data frame. Empty VCF files return a zero-row data frame with the
#'   detected VCF columns. When `add_allele_depth = TRUE` and medaka/wf-amplicon
#'   INFO fields are available, additional columns are added:
#'   \itemize{
#'     \item `ref_fwd_depth`, `ref_rev_depth`: forward/reverse spanning reads
#'       best aligned to the reference allele.
#'     \item `alt_fwd_depth`, `alt_rev_depth`: forward/reverse spanning reads
#'       best aligned to ALT allele(s). For multi-allelic records, all ALT
#'       alleles are summed.
#'     \item `ref_depth`, `alt_depth`: strand-summed REF and ALT support.
#'     \item `variant_percent`: `alt_depth / (ref_depth + alt_depth) * 100`.
#'     \item `ambiguous_fwd_depth`, `ambiguous_rev_depth`, `ambiguous_depth`:
#'       ambiguous spanning reads from `AR`.
#'   }
#'
#' @examples
#' vcf <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
#'   "chr1\t10\t.\tA\tG\t60\tPASS\tDP=30;AF=0.4\tGT:DP:AD\t0/1:30:18,12"
#' ), vcf)
#'
#' read_vcf(vcf)
#'
#' @export
read_vcf <- function(vcf_file,
                     parse_info = TRUE,
                     parse_genotypes = TRUE,
                     sample_format = c("wide", "long"),
                     samples = NULL,
                     add_variant_type = TRUE,
                     add_allele_depth = TRUE,
                     keep_info = FALSE,
                     keep_format = FALSE,
                     simplify = TRUE,
                     info_prefix = "",
                     genotype_prefix = NULL) {
  check_file_arg(vcf_file, "vcf_file")
  check_logical_scalar(parse_info, "parse_info")
  check_logical_scalar(parse_genotypes, "parse_genotypes")
  check_logical_scalar(add_variant_type, "add_variant_type")
  check_logical_scalar(add_allele_depth, "add_allele_depth")
  check_logical_scalar(keep_info, "keep_info")
  check_logical_scalar(keep_format, "keep_format")
  check_logical_scalar(simplify, "simplify")
  check_optional_prefix(info_prefix, "info_prefix")
  if (!is.null(genotype_prefix)) check_scalar_character(genotype_prefix, "genotype_prefix")
  if (!is.null(samples) && (!is.character(samples) || anyNA(samples))) {
    stop("`samples` must be a character vector without missing values.",
         call. = FALSE)
  }
  sample_format <- match.arg(sample_format)

  vcf_file <- normalizePath(vcf_file, mustWork = TRUE)
  lines <- read_vcf_lines(vcf_file)
  header_index <- which(startsWith(lines, "#CHROM"))
  if (length(header_index) != 1L) {
    stop("VCF header line starting with `#CHROM` was not found exactly once.",
         call. = FALSE)
  }

  header <- strsplit(lines[[header_index]], "\t", fixed = TRUE)[[1L]]
  header[[1L]] <- sub("^#", "", header[[1L]])
  data_lines <- if (header_index < length(lines)) {
    lines[seq.int(header_index + 1L, length(lines))]
  } else {
    character()
  }
  data_lines <- data_lines[nzchar(data_lines)]

  if (length(data_lines) == 0L) {
    return(empty_vcf_data_frame(header))
  }

  vcf <- utils::read.delim(
    text = paste(data_lines, collapse = "\n"),
    header = FALSE,
    col.names = header,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  vcf <- normalize_vcf_fixed_columns(vcf)

  fixed_cols <- intersect(
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"),
    names(vcf)
  )
  has_info <- "INFO" %in% names(vcf)
  has_format <- "FORMAT" %in% names(vcf)
  format_index <- if (has_format) match("FORMAT", names(vcf)) else NA_integer_
  sample_cols <- if (has_format && format_index < length(names(vcf))) {
    names(vcf)[seq.int(format_index + 1L, length(names(vcf)))]
  } else {
    character()
  }
  sample_cols <- sample_cols[nzchar(sample_cols)]

  if (!is.null(samples)) {
    missing_samples <- setdiff(samples, sample_cols)
    if (length(missing_samples) > 0L) {
      stop("Sample column(s) not found in VCF: ",
           paste(missing_samples, collapse = ", "), call. = FALSE)
    }
    sample_cols <- samples
  }

  out <- vcf[, fixed_cols, drop = FALSE]
  if (isTRUE(add_variant_type) && all(c("REF", "ALT") %in% names(out))) {
    out$Type <- vcf_variant_type(out$REF, out$ALT)
  }

  if (isTRUE(keep_info) && has_info) {
    out$INFO <- vcf$INFO
  }
  info_for_depth <- NULL
  if (isTRUE(parse_info) && has_info) {
    info <- parse_vcf_info(vcf$INFO, prefix = info_prefix, simplify = simplify)
    if (isTRUE(add_allele_depth)) {
      info_for_depth <- parse_vcf_info(vcf$INFO, simplify = FALSE)
    }
    out <- cbind(out, info)
  } else if (isTRUE(add_allele_depth) && has_info) {
    info_for_depth <- parse_vcf_info(vcf$INFO, simplify = FALSE)
  }
  if (isTRUE(add_allele_depth) && !is.null(info_for_depth)) {
    out <- add_vcf_allele_depth_columns(out, info_for_depth)
  }

  if (!isTRUE(parse_genotypes) || !has_format || length(sample_cols) == 0L) {
    if (isTRUE(keep_format) && has_format) {
      out <- cbind(out, vcf[, c("FORMAT", sample_cols), drop = FALSE])
    }
    return(out)
  }

  if (identical(sample_format, "long")) {
    return(parse_vcf_genotypes_long(
      fixed_info = out,
      format = vcf$FORMAT,
      sample_data = vcf[, sample_cols, drop = FALSE],
      keep_format = keep_format,
      simplify = simplify
    ))
  }

  if (isTRUE(keep_format)) {
    out <- cbind(out, vcf[, c("FORMAT", sample_cols), drop = FALSE])
  }
  if (length(sample_cols) == 1L) {
    out$sampleID <- sample_cols[[1L]]
  }

  genotype <- parse_vcf_genotypes_wide(
    format = vcf$FORMAT,
    sample_data = vcf[, sample_cols, drop = FALSE],
    simplify = simplify,
    genotype_prefix = genotype_prefix
  )
  if (is.null(genotype_prefix) && length(sample_cols) == 1L) {
    names(genotype) <- disambiguate_vcf_genotype_names(names(genotype), names(out))
  }
  cbind(out, genotype)
}

#' Infer VCF variant type from REF and ALT alleles
#'
#' `vcf_variant_type()` classifies variants from the allele strings in the same
#' broad way used by wf-amplicon reports, with an extra normalized-indel guard:
#' single-base substitutions are `SNP`, unequal-length alleles with a shared
#' anchor base are `INDEL`, equal-length multi-base substitutions are `MNP`,
#' unequal-length alleles without a shared anchor are `OTHER`, and reference or
#' missing calls are `REF` or `NA`.
#'
#' @param ref Character vector of VCF `REF` alleles.
#' @param alt Character vector of VCF `ALT` alleles. Multi-allelic values such as
#'   `"A,AT"` are supported.
#' @param collapse_multiallelic Logical. If `TRUE`, a multi-allelic `ALT` value
#'   is summarized into one value per VCF row. Rows with mixed non-reference
#'   variant types are reported as `"MIXED"`. If `FALSE`, a list-column with one
#'   type per ALT allele is returned.
#'
#' @return A character vector when `collapse_multiallelic = TRUE`; otherwise a
#'   list of character vectors.
#'
#' @examples
#' vcf_variant_type(c("A", "AT", "A", "AC"), c("G", "A", "AT", "GT"))
#'
#' @export
vcf_variant_type <- function(ref, alt, collapse_multiallelic = TRUE) {
  if (!is.character(ref)) ref <- as.character(ref)
  if (!is.character(alt)) alt <- as.character(alt)
  check_logical_scalar(collapse_multiallelic, "collapse_multiallelic")

  n <- max(length(ref), length(alt))
  ref <- rep_len(ref, n)
  alt <- rep_len(alt, n)

  types <- Map(function(r, a) {
    infer_vcf_variant_type_one(r, a)
  }, ref, alt)

  if (!isTRUE(collapse_multiallelic)) {
    return(types)
  }

  unname(vapply(types, summarize_vcf_variant_types, character(1)))
}

add_vcf_allele_depth_columns <- function(out, info) {
  if ("SR" %in% names(info)) {
    sr <- parse_vcf_sr(info$SR, out$ALT)
    out$ref_fwd_depth <- sr$ref_fwd_depth
    out$ref_rev_depth <- sr$ref_rev_depth
    out$alt_fwd_depth <- sr$alt_fwd_depth
    out$alt_rev_depth <- sr$alt_rev_depth
    out$ref_depth <- sr$ref_depth
    out$alt_depth <- sr$alt_depth
    out$variant_percent <- sr$variant_percent
  }
  if ("AR" %in% names(info)) {
    ar <- parse_vcf_ar(info$AR)
    out$ambiguous_fwd_depth <- ar$ambiguous_fwd_depth
    out$ambiguous_rev_depth <- ar$ambiguous_rev_depth
    out$ambiguous_depth <- ar$ambiguous_depth
  }
  out
}

parse_vcf_sr <- function(sr, alt) {
  parsed <- Map(function(value, alt_value) {
    values <- parse_integer_list(value)
    if (length(values) < 2L) {
      return(rep(NA_integer_, 7L))
    }

    n_alt <- count_alt_alleles(alt_value)
    if (n_alt == 0L) {
      ref_fwd <- values[[1L]]
      ref_rev <- values[[2L]]
      ref_depth <- sum(c(ref_fwd, ref_rev), na.rm = TRUE)
      if (all(is.na(c(ref_fwd, ref_rev)))) ref_depth <- NA_integer_
      return(c(
        ref_fwd, ref_rev, NA_integer_, NA_integer_,
        ref_depth, NA_integer_, NA_real_
      ))
    }
    expected <- 2L + 2L * n_alt
    if (length(values) < expected) {
      length(values) <- expected
    }

    ref_fwd <- values[[1L]]
    ref_rev <- values[[2L]]
    alt_values <- values[seq.int(3L, expected)]
    alt_fwd <- sum(alt_values[c(TRUE, FALSE)], na.rm = TRUE)
    alt_rev <- sum(alt_values[c(FALSE, TRUE)], na.rm = TRUE)
    if (all(is.na(alt_values[c(TRUE, FALSE)]))) alt_fwd <- NA_integer_
    if (all(is.na(alt_values[c(FALSE, TRUE)]))) alt_rev <- NA_integer_

    ref_depth <- sum(c(ref_fwd, ref_rev), na.rm = TRUE)
    alt_depth <- sum(c(alt_fwd, alt_rev), na.rm = TRUE)
    if (all(is.na(c(ref_fwd, ref_rev)))) ref_depth <- NA_integer_
    if (all(is.na(c(alt_fwd, alt_rev)))) alt_depth <- NA_integer_

    denominator <- ref_depth + alt_depth
    variant_percent <- if (is.na(denominator) || denominator == 0L) {
      NA_real_
    } else {
      alt_depth / denominator * 100
    }

    c(
      ref_fwd, ref_rev, alt_fwd, alt_rev,
      ref_depth, alt_depth, variant_percent
    )
  }, sr, alt)

  mat <- do.call(rbind, parsed)
  colnames(mat) <- c(
    "ref_fwd_depth", "ref_rev_depth",
    "alt_fwd_depth", "alt_rev_depth",
    "ref_depth", "alt_depth", "variant_percent"
  )
  as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
}

parse_vcf_ar <- function(ar) {
  parsed <- lapply(ar, function(value) {
    values <- parse_integer_list(value)
    length(values) <- max(length(values), 2L)
    ambiguous_fwd <- values[[1L]]
    ambiguous_rev <- values[[2L]]
    ambiguous_depth <- sum(c(ambiguous_fwd, ambiguous_rev), na.rm = TRUE)
    if (all(is.na(c(ambiguous_fwd, ambiguous_rev)))) ambiguous_depth <- NA_integer_
    c(ambiguous_fwd, ambiguous_rev, ambiguous_depth)
  })

  mat <- do.call(rbind, parsed)
  colnames(mat) <- c(
    "ambiguous_fwd_depth",
    "ambiguous_rev_depth",
    "ambiguous_depth"
  )
  as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
}

parse_integer_list <- function(x) {
  if (is.na(x) || !nzchar(x) || identical(x, ".")) {
    return(integer())
  }
  values <- strsplit(x, ",", fixed = TRUE)[[1L]]
  values[values == "."] <- NA_character_
  suppressWarnings(as.integer(values))
}

count_alt_alleles <- function(alt) {
  if (is.na(alt) || !nzchar(alt) || identical(alt, ".")) {
    return(0L)
  }
  length(strsplit(alt, ",", fixed = TRUE)[[1L]])
}

#' Parse VCF INFO strings
#'
#' @param x Character vector of VCF `INFO` values.
#' @param prefix Prefix added to output column names.
#' @param simplify Logical. If `TRUE`, convert parsed scalar numeric fields to
#'   numeric columns when possible.
#'
#' @return A data frame with one row per input value.
#'
#' @export
parse_vcf_info <- function(x, prefix = "", simplify = TRUE) {
  if (!is.character(x)) x <- as.character(x)
  check_optional_prefix(prefix, "prefix")
  check_logical_scalar(simplify, "simplify")

  parsed <- lapply(x, parse_vcf_info_one)
  keys <- unique(unlist(lapply(parsed, names), use.names = FALSE))
  if (length(keys) == 0L) {
    return(data.frame(row.names = seq_along(x)))
  }

  values <- lapply(keys, function(key) {
    vapply(
      parsed,
      function(item) {
        value <- item[[key]]
        if (is.null(value)) NA_character_ else value
      },
      character(1)
    )
  })
  names(values) <- paste0(prefix, make.names(keys, unique = TRUE))
  out <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  if (isTRUE(simplify)) out <- simplify_vcf_columns(out)
  out
}

read_vcf_lines <- function(vcf_file) {
  con <- if (grepl("[.]gz$", vcf_file, ignore.case = TRUE)) {
    gzfile(vcf_file, open = "rt")
  } else {
    file(vcf_file, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

empty_vcf_data_frame <- function(header) {
  out <- as.data.frame(
    stats::setNames(rep(list(logical()), length(header)), header),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if ("POS" %in% names(out)) out$POS <- integer()
  if ("QUAL" %in% names(out)) out$QUAL <- numeric()
  out
}

normalize_vcf_fixed_columns <- function(vcf) {
  if ("POS" %in% names(vcf)) {
    vcf$POS <- suppressWarnings(as.integer(vcf$POS))
  }
  if ("QUAL" %in% names(vcf)) {
    vcf$QUAL[vcf$QUAL == "."] <- NA_character_
    vcf$QUAL <- suppressWarnings(as.numeric(vcf$QUAL))
  }
  vcf
}

parse_vcf_info_one <- function(info) {
  if (is.na(info) || !nzchar(info) || identical(info, ".")) {
    return(list())
  }

  fields <- strsplit(info, ";", fixed = TRUE)[[1L]]
  fields <- fields[nzchar(fields)]
  out <- list()
  for (field in fields) {
    if (grepl("=", field, fixed = TRUE)) {
      key <- sub("=.*$", "", field)
      value <- sub("^[^=]*=", "", field)
      out[[key]] <- value
    } else {
      out[[field]] <- "TRUE"
    }
  }
  out
}

parse_vcf_genotypes_wide <- function(format,
                                     sample_data,
                                     simplify,
                                     genotype_prefix) {
  sample_cols <- names(sample_data)
  pieces <- lapply(sample_cols, function(sample) {
    parsed <- parse_vcf_sample_fields(format, sample_data[[sample]])
    if (isTRUE(simplify)) parsed <- simplify_vcf_columns(parsed)

    if (length(sample_cols) > 1L || !is.null(genotype_prefix)) {
      prefix <- if (is.null(genotype_prefix)) paste0(sample, "_") else genotype_prefix
      names(parsed) <- paste0(prefix, names(parsed))
    }
    parsed
  })

  do.call(cbind, pieces)
}

parse_vcf_genotypes_long <- function(fixed_info,
                                     format,
                                     sample_data,
                                     keep_format,
                                     simplify) {
  sample_cols <- names(sample_data)
  rows <- lapply(sample_cols, function(sample) {
    parsed <- parse_vcf_sample_fields(format, sample_data[[sample]])
    if (isTRUE(simplify)) parsed <- simplify_vcf_columns(parsed)
    names(parsed) <- disambiguate_vcf_genotype_names(names(parsed), names(fixed_info))
    out <- cbind(
      fixed_info,
      data.frame(sampleID = sample, stringsAsFactors = FALSE)
    )
    if (isTRUE(keep_format)) {
      out$FORMAT <- format
      out$sample_value <- sample_data[[sample]]
    }
    cbind(out, parsed)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

parse_vcf_sample_fields <- function(format, sample_value) {
  parsed <- Map(function(fmt, value) {
    if (is.na(fmt) || !nzchar(fmt) || is.na(value) || !nzchar(value) ||
        identical(value, ".")) {
      return(list())
    }
    keys <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
    values <- strsplit(value, ":", fixed = TRUE)[[1L]]
    length(values) <- length(keys)
    values[is.na(values) | values == "."] <- NA_character_
    stats::setNames(as.list(values), keys)
  }, format, sample_value)

  keys <- unique(unlist(lapply(parsed, names), use.names = FALSE))
  if (length(keys) == 0L) {
    return(data.frame(row.names = seq_along(format)))
  }

  values <- lapply(keys, function(key) {
    vapply(
      parsed,
      function(item) {
        value <- item[[key]]
        if (is.null(value)) NA_character_ else value
      },
      character(1)
    )
  })
  names(values) <- make.names(keys, unique = TRUE)
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

simplify_vcf_columns <- function(data) {
  for (nm in names(data)) {
    value <- data[[nm]]
    non_missing <- value[!is.na(value)]
    if (length(non_missing) == 0L) next
    if (all(grepl("^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$", non_missing))) {
      data[[nm]] <- suppressWarnings(as.numeric(value))
    } else if (all(non_missing %in% c("TRUE", "FALSE"))) {
      data[[nm]] <- value == "TRUE"
    }
  }
  data
}

check_optional_prefix <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be a single character string.", call. = FALSE)
  }
}

disambiguate_vcf_genotype_names <- function(genotype_names, existing_names) {
  conflicts <- genotype_names %in% existing_names
  genotype_names[conflicts] <- paste0("FORMAT_", genotype_names[conflicts])
  genotype_names
}

infer_vcf_variant_type_one <- function(ref, alt) {
  if (is.na(ref) || is.na(alt) || !nzchar(ref) || !nzchar(alt) ||
      identical(ref, ".") || identical(alt, ".")) {
    return(NA_character_)
  }

  alts <- strsplit(alt, ",", fixed = TRUE)[[1L]]
  alts <- alts[nzchar(alts)]
  unname(vapply(alts, function(one_alt) {
    if (is.na(one_alt) || identical(one_alt, ".")) {
      return(NA_character_)
    }
    if (identical(ref, one_alt)) {
      return("REF")
    }
    if (nchar(ref) == 1L && nchar(one_alt) == 1L) {
      return("SNP")
    }
    if (nchar(ref) != nchar(one_alt)) {
      if (has_shared_indel_anchor(ref, one_alt)) {
        return("INDEL")
      }
      return("OTHER")
    }
    "MNP"
  }, character(1)))
}

has_shared_indel_anchor <- function(ref, alt) {
  substr(ref, 1L, 1L) == substr(alt, 1L, 1L) ||
    substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))
}

summarize_vcf_variant_types <- function(types) {
  types <- types[!is.na(types)]
  if (length(types) == 0L) {
    return(NA_character_)
  }
  non_ref_types <- types[types != "REF"]
  if (length(non_ref_types) == 0L) {
    return("REF")
  }
  unique_types <- unique(non_ref_types)
  if (length(unique_types) == 1L) {
    return(unique_types[[1L]])
  }
  "MIXED"
}
