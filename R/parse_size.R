#' Parse base-pair size strings
#'
#' `parse_size()` converts strings such as `"180k"`, `"180kbp"`, `"2.5Mbp"`,
#' or `"1000 bp"` to numeric base-pair counts.
#'
#' @param x Character vector of size strings.
#' @param units Named numeric vector defining unit multipliers. An empty unit is
#'   treated as `"bp"`.
#' @param na_on_error Logical. If `TRUE`, invalid size strings or unknown units
#'   are returned as `NA_real_` instead of throwing an error. Input type errors
#'   and invalid `units` still throw errors.
#'
#' @return A numeric vector with parsed sizes in base pairs.
#'
#' @examples
#' parse_size(c("180k", "2.5Mbp", "1000 bp"))
#' parse_size(c("180k", "bad", "10tb"), na_on_error = TRUE)
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
                       ),
                       na_on_error = FALSE) {
  if (!is.character(x)) {
    stop("`x` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(units) || is.null(names(units)) ||
      anyNA(units) || any(!is.finite(units)) || any(!nzchar(names(units)))) {
    stop("`units` must be a named numeric vector with finite values.",
         call. = FALSE)
  }
  if (!is.logical(na_on_error) || length(na_on_error) != 1L ||
      is.na(na_on_error)) {
    stop("`na_on_error` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  compact <- gsub("[[:space:]]+", "", x)
  pattern <- "^([+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?)([[:alpha:]]*)$"
  matched <- !is.na(compact) & grepl(pattern, compact, perl = TRUE)
  if (any(!matched)) {
    if (!isTRUE(na_on_error)) {
      stop("Invalid size string: ", x[which(!matched)[1L]], call. = FALSE)
    }
  }

  value <- rep(NA_real_, length(x))
  unit <- rep(NA_character_, length(x))
  value[matched] <- as.numeric(sub(pattern, "\\1", compact[matched], perl = TRUE))
  unit[matched] <- sub(pattern, "\\2", compact[matched], perl = TRUE)
  unit[unit == ""] <- "bp"

  unknown <- matched & !unit %in% names(units)
  if (any(unknown)) {
    if (!isTRUE(na_on_error)) {
      stop(
        "Unknown size unit `",
        unit[which(unknown)[1L]],
        "`. Supported units: ",
        paste(names(units), collapse = ", "),
        call. = FALSE
      )
    }
  }

  parsed <- rep(NA_real_, length(x))
  valid <- matched & !unknown
  parsed[valid] <- value[valid] * units[unit[valid]]
  unname(parsed)
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
                      ),
                      na_on_error = FALSE) {
  parse_size(x, units = units, na_on_error = na_on_error)
}
