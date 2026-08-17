#' Export a Sanger/reference alignment to HTML and PDF
#'
#' `export_sanger_alignment_pdf()` renders a Sanger/reference multiple sequence
#' alignment as a color-coded HTML file with [DECIPHER::BrowseSeqs()], and can
#' optionally convert that HTML file to PDF with [pagedown::chrome_print()].
#'
#' The `alignment` input can be the `alignment` element returned by
#' [align_sanger_with_full_reference()], or the full
#' `"SangerMultipleAlignment"` result object returned by that function.
#'
#' @param alignment A [Biostrings::XStringSet()] alignment object, or a
#'   `"SangerMultipleAlignment"` list containing an `alignment` element.
#' @param pdf_file Output PDF path. Ignored when `render_pdf = FALSE`, except
#'   for deriving `html_file` when `html_file = NULL`.
#' @param html_file Output HTML path. When `NULL`, it is derived from
#'   `pdf_file` by replacing a trailing `.pdf` extension with `.html`.
#' @param block_width Positive integer column width passed to
#'   [DECIPHER::BrowseSeqs()]. It must be a multiple of 20.
#' @param highlight_reference Logical. If `TRUE`, highlight the first sequence,
#'   which is expected to be the reference sequence.
#' @param render_pdf Logical. If `TRUE`, convert the HTML output to PDF with
#'   [pagedown::chrome_print()]. This requires the suggested `pagedown` package
#'   and an available Chromium/Chrome installation.
#' @param open_pdf Logical. If `TRUE` and `render_pdf = TRUE`, open the PDF after
#'   it is written. Defaults to `FALSE` for non-interactive and batch use.
#' @param overwrite Logical. If `FALSE`, stop when an output file already exists.
#' @param quiet Logical. If `FALSE`, print output paths after rendering.
#'
#' @return Invisibly returns a list with `html`, `pdf`, `render_pdf`, and
#'   `alignment_names`.
#'
#' @examples
#' alignment <- Biostrings::DNAStringSet(c(
#'   reference = "AAACCCGGGTTT",
#'   read1 = "AAACCCGGGTTT",
#'   read2 = "AAA---GGGTTT"
#' ))
#' html_file <- tempfile(fileext = ".html")
#'
#' res <- export_sanger_alignment_pdf(
#'   alignment,
#'   html_file = html_file,
#'   render_pdf = FALSE
#' )
#' res$html
#'
#' @export
export_sanger_alignment_pdf <- function(alignment,
                                        pdf_file = "sanger_alignment.pdf",
                                        html_file = NULL,
                                        block_width = 80,
                                        highlight_reference = FALSE,
                                        render_pdf = TRUE,
                                        open_pdf = FALSE,
                                        overwrite = FALSE,
                                        quiet = FALSE) {
  alignment <- extract_sanger_alignment(alignment)

  check_scalar_character(pdf_file, "pdf_file")
  if (is.null(html_file)) {
    html_file <- sub("\\.pdf$", ".html", pdf_file, ignore.case = TRUE)
    if (identical(html_file, pdf_file)) {
      html_file <- paste0(pdf_file, ".html")
    }
  } else {
    check_scalar_character(html_file, "html_file")
  }

  block_width <- validate_positive_integer(block_width, "block_width")
  if (block_width %% 20L != 0L) {
    stop("`block_width` must be a multiple of 20, for example 60, 80, or 100.",
         call. = FALSE)
  }

  check_logical_scalar(highlight_reference, "highlight_reference")
  check_logical_scalar(render_pdf, "render_pdf")
  check_logical_scalar(open_pdf, "open_pdf")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(quiet, "quiet")

  if (length(alignment) == 0L) {
    stop("`alignment` must contain at least one sequence.", call. = FALSE)
  }

  html_file <- normalize_output_path(html_file)
  pdf_file <- normalize_output_path(pdf_file)

  existing_outputs <- c(html_file, if (isTRUE(render_pdf)) pdf_file else character())
  existing_outputs <- existing_outputs[file.exists(existing_outputs)]
  if (length(existing_outputs) > 0L && !isTRUE(overwrite)) {
    stop(
      "Output file already exists. Set `overwrite = TRUE` to replace: ",
      paste(existing_outputs, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(render_pdf) && !requireNamespace("pagedown", quietly = TRUE)) {
    stop(
      "`pagedown` is required when `render_pdf = TRUE`. Install it with ",
      "`install.packages(\"pagedown\")`, or set `render_pdf = FALSE`.",
      call. = FALSE
    )
  }

  highlight_value <- if (isTRUE(highlight_reference)) 1L else NA_integer_

  DECIPHER::BrowseSeqs(
    alignment,
    htmlFile = html_file,
    openURL = FALSE,
    colWidth = block_width,
    highlight = highlight_value,
    patterns = c("-", "A", "C", "G", "T", "N"),
    colors = c(
      "#D9D9D9",
      "#72D572",
      "#70A1FF",
      "#FFD166",
      "#FF7675",
      "#C39BD3"
    )
  )

  pdf_output <- NA_character_
  if (isTRUE(render_pdf)) {
    pagedown::chrome_print(
      input = html_file,
      output = pdf_file,
      format = "pdf",
      options = list(
        landscape = TRUE,
        printBackground = TRUE,
        preferCSSPageSize = FALSE,
        marginTop = 0.3,
        marginBottom = 0.3,
        marginLeft = 0.3,
        marginRight = 0.3
      )
    )
    pdf_output <- pdf_file
  }

  if (!isTRUE(quiet)) {
    message("HTML file: ", html_file)
    if (isTRUE(render_pdf)) {
      message("PDF file: ", pdf_file)
    }
  }

  if (isTRUE(open_pdf) && isTRUE(render_pdf)) {
    utils::browseURL(pdf_file)
  }

  invisible(list(
    html = html_file,
    pdf = pdf_output,
    render_pdf = render_pdf,
    alignment_names = names(alignment)
  ))
}

extract_sanger_alignment <- function(alignment) {
  if (inherits(alignment, "SangerMultipleAlignment")) {
    alignment <- alignment$alignment
  }

  if (!inherits(alignment, "XStringSet")) {
    stop(
      "`alignment` must be an XStringSet object or a SangerMultipleAlignment result.",
      call. = FALSE
    )
  }

  alignment
}

normalize_output_path <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}
