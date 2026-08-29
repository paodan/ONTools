#' Parse base-pair size strings
#'
#' `parse_size()` converts strings such as `"180k"`, `"180kbp"`, `"2.5Mbp"`,
#' or `"1000 bp"` to numeric base-pair counts.
#'
#' @param x Character vector of size strings.
#' @param units Named numeric vector defining unit multipliers. An empty unit is
#'   treated as `"bp"`.
#'
#' @return A numeric vector with parsed sizes in base pairs.
#'
#' @examples
#' parse_size(c("180k", "2.5Mbp", "1000 bp"))
#'
#' @export
parse_size <- function(x,
                       units = c(
                         bp = 1,
                         b = 1,
                         k = 1000,
                         kb = 1000,
                         kbp = 1000,
                         Kbp = 1000,
                         m = 1e6,
                         mb = 1e6,
                         Mb = 1e6,
                         mbp = 1e6,
                         Mbp = 1e6,
                         g = 1e9,
                         gb = 1e9,
                         Gb = 1e9,
                         gbp = 1e9,
                         Gbp = 1e9
                       )) {
  if (!is.character(x)) {
    stop("`x` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(units) || is.null(names(units)) ||
      anyNA(units) || any(!is.finite(units)) || any(!nzchar(names(units)))) {
    stop("`units` must be a named numeric vector with finite values.",
         call. = FALSE)
  }

  compact <- gsub("[[:space:]]+", "", x)
  pattern <- "^([+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?)([[:alpha:]]*)$"
  matched <- !is.na(compact) & grepl(pattern, compact, perl = TRUE)
  if (any(!matched)) {
    stop("Invalid size string: ", x[which(!matched)[1L]], call. = FALSE)
  }

  value <- as.numeric(sub(pattern, "\\1", compact, perl = TRUE))
  unit <- sub(pattern, "\\2", compact, perl = TRUE)
  unit[unit == ""] <- "bp"

  unknown <- !unit %in% names(units)
  if (any(unknown)) {
    stop(
      "Unknown size unit `",
      unit[which(unknown)[1L]],
      "`. Supported units: ",
      paste(names(units), collapse = ", "),
      call. = FALSE
    )
  }

  unname(value * units[unit])
}

#' @rdname parse_size
#' @export
parseSize <- function(x,
                      units = c(
                        bp = 1,
                        b = 1,
                        k = 1000,
                        kb = 1000,
                        kbp = 1000,
                        Kbp = 1000,
                        m = 1e6,
                        mb = 1e6,
                        Mb = 1e6,
                        mbp = 1e6,
                        Mbp = 1e6,
                        g = 1e9,
                        gb = 1e9,
                        Gb = 1e9,
                        gbp = 1e9,
                        Gbp = 1e9
                      )) {
  parse_size(x, units = units)
}
