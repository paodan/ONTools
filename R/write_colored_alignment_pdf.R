#' Write a color-coded DNA alignment to PDF
#'
#' `write_colored_alignment_pdf()` writes an aligned
#' [Biostrings::DNAStringSet()] to a multi-page vector PDF. Each base is drawn
#' separately so that nucleotide colors, character spacing, and line spacing
#' remain configurable. Long alignments are split into blocks, and fully
#' conserved columns are marked with an asterisk.
#'
#' The input is expected to be aligned already: every sequence must have the
#' same width. This function does not perform sequence alignment.
#'
#' @param alignment A [Biostrings::DNAStringSet()] containing aligned DNA
#'   sequences of equal width.
#' @param pdf_file Output PDF path. Parent directories are created when needed.
#' @param bases_per_block Positive integer giving the maximum number of
#'   alignment columns displayed in each block.
#' @param blocks_per_page Positive integer giving the number of alignment
#'   blocks displayed on each PDF page.
#' @param font_size Positive numeric font size in points.
#' @param character_spacing Positive numeric multiplier controlling horizontal
#'   space between bases. Values greater than 1 increase spacing.
#' @param line_spacing Positive numeric multiplier controlling vertical space
#'   between sequences. Values greater than 1 increase spacing.
#' @param coordinate_start Positive integer used as the displayed coordinate of
#'   the first column in `alignment`. This is useful when `alignment` is a
#'   subsequence of a larger alignment.
#' @param base_colors Named character vector mapping bases to valid R colors.
#'   Unlisted IUPAC ambiguity symbols are drawn with `ambiguity_color`.
#' @param ambiguity_color Color used for bases not present in `base_colors`.
#' @param show_conservation Logical. If `TRUE`, mark columns that are identical
#'   in every sequence with `"*"`.
#' @param overwrite Logical. If `FALSE`, stop when `pdf_file` already exists.
#' @param quiet Logical. If `FALSE`, report the output path after writing.
#'
#' @return Invisibly returns the normalized path to the PDF file.
#'
#' @examples
#' alignment <- Biostrings::DNAStringSet(c(
#'   Sample_1 = "ACTGAA--CTGA",
#'   Sample_2 = "ACTGAA--CTGA",
#'   Sample_3 = "ACTGAATTCT-A"
#' ))
#'
#' pdf_file <- tempfile(fileext = ".pdf")
#' write_colored_alignment_pdf(
#'   alignment,
#'   pdf_file,
#'   bases_per_block = 10,
#'   blocks_per_page = 2,
#'   character_spacing = 1.2
#' )
#'
#' # Preserve genomic coordinates when plotting a sliced alignment.
#' region <- Biostrings::subseq(alignment, start = 3, end = 10)
#' region_pdf <- tempfile(fileext = ".pdf")
#' write_colored_alignment_pdf(
#'   region,
#'   region_pdf,
#'   coordinate_start = 3
#' )
#'
#' @export
write_colored_alignment_pdf <- function(
    alignment,
    pdf_file,
    bases_per_block = 60L,
    blocks_per_page = 3L,
    font_size = 8,
    character_spacing = 1,
    line_spacing = 1,
    coordinate_start = 1L,
    base_colors = c(
      A = "#009E73",
      C = "#0072B2",
      G = "#E69F00",
      T = "#D55E00",
      U = "#D55E00",
      `-` = "#999999",
      N = "#000000",
      `?` = "#000000"
    ),
    ambiguity_color = "#984EA3",
    show_conservation = TRUE,
    overwrite = FALSE,
    quiet = FALSE) {
  if (!inherits(alignment, "DNAStringSet")) {
    stop("`alignment` must be a DNAStringSet object.", call. = FALSE)
  }
  if (length(alignment) == 0L) {
    stop("`alignment` must contain at least one sequence.", call. = FALSE)
  }

  alignment_widths <- Biostrings::width(alignment)
  if (length(unique(alignment_widths)) != 1L) {
    stop(
      "All sequences in `alignment` must have the same width.",
      call. = FALSE
    )
  }
  if (alignment_widths[[1L]] == 0L) {
    stop("Sequences in `alignment` must not be empty.", call. = FALSE)
  }

  check_scalar_character(pdf_file, "pdf_file")
  bases_per_block <- validate_positive_integer(
    bases_per_block,
    "bases_per_block"
  )
  blocks_per_page <- validate_positive_integer(
    blocks_per_page,
    "blocks_per_page"
  )
  coordinate_start <- validate_positive_integer(
    coordinate_start,
    "coordinate_start"
  )
  font_size <- validate_alignment_pdf_positive_number(font_size, "font_size")
  character_spacing <- validate_alignment_pdf_positive_number(
    character_spacing,
    "character_spacing"
  )
  line_spacing <- validate_alignment_pdf_positive_number(
    line_spacing,
    "line_spacing"
  )
  check_logical_scalar(show_conservation, "show_conservation")
  check_logical_scalar(overwrite, "overwrite")
  check_logical_scalar(quiet, "quiet")
  validate_alignment_pdf_colors(base_colors, ambiguity_color)

  pdf_file <- normalize_output_path(pdf_file)
  if (file.exists(pdf_file) && !isTRUE(overwrite)) {
    stop(
      "Output file already exists. Set `overwrite = TRUE` to replace: ",
      pdf_file,
      call. = FALSE
    )
  }

  sequence_names <- names(alignment)
  if (is.null(sequence_names)) {
    sequence_names <- paste0("Seq_", seq_along(alignment))
  } else {
    missing_names <- is.na(sequence_names) | !nzchar(sequence_names)
    sequence_names[missing_names] <- paste0("Seq_", which(missing_names))
  }

  sequences <- toupper(as.character(alignment))
  alignment_width <- alignment_widths[[1L]]
  number_of_sequences <- length(alignment)
  block_starts <- seq.int(1L, alignment_width, by = bases_per_block)

  # Physical dimensions are used so character spacing remains meaningful in
  # the PDF instead of being cancelled by automatic plot-coordinate scaling.
  character_step <- 0.085 * character_spacing
  row_step <- 0.19 * line_spacing
  left_margin <- 0.3
  right_margin <- 0.3
  top_margin <- 0.3
  bottom_margin <- 0.3
  name_width <- max(1, max(nchar(sequence_names)) * 0.085 + 0.25)
  title_height <- 0.28
  coordinate_height <- 0.38
  block_height <- title_height +
    number_of_sequences * row_step +
    row_step +
    coordinate_height
  page_width <- max(
    8,
    left_margin + name_width + bases_per_block * character_step + right_margin
  )
  page_height <- max(
    6,
    top_margin + blocks_per_page * block_height + bottom_margin
  )

  grDevices::pdf(
    file = pdf_file,
    width = page_width,
    height = page_height,
    onefile = TRUE,
    useDingbats = FALSE
  )
  device_open <- TRUE
  on.exit({
    if (isTRUE(device_open)) {
      grDevices::dev.off()
    }
  }, add = TRUE)

  block_groups <- split(
    seq_along(block_starts),
    ceiling(seq_along(block_starts) / blocks_per_page)
  )

  for (page_blocks in block_groups) {
    grid::grid.newpage()

    for (slot in seq_along(page_blocks)) {
      block_index <- page_blocks[[slot]]
      start_position <- block_starts[[block_index]]
      end_position <- min(
        start_position + bases_per_block - 1L,
        alignment_width
      )
      fragments <- substr(sequences, start_position, end_position)
      base_matrix <- do.call(rbind, strsplit(fragments, split = ""))
      block_width <- ncol(base_matrix)

      block_top <- page_height -
        top_margin -
        (slot - 1L) * block_height
      sequence_start_x <- left_margin + name_width
      x_positions <- sequence_start_x +
        (seq_len(block_width) - 1L) * character_step
      first_row_y <- block_top - title_height
      y_positions <- first_row_y -
        (seq_len(number_of_sequences) - 1L) * row_step
      displayed_start <- coordinate_start + start_position - 1L
      displayed_end <- coordinate_start + end_position - 1L

      grid::grid.text(
        label = sprintf(
          "Alignment positions %d-%d",
          displayed_start,
          displayed_end
        ),
        x = grid::unit(left_margin, "inches"),
        y = grid::unit(block_top - 0.03, "inches"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = font_size + 1, fontface = "bold")
      )

      for (i in seq_len(number_of_sequences)) {
        grid::grid.text(
          label = sequence_names[[i]],
          x = grid::unit(sequence_start_x - 0.12, "inches"),
          y = grid::unit(y_positions[[i]], "inches"),
          just = "right",
          gp = grid::gpar(
            fontsize = font_size,
            fontfamily = "mono",
            fontface = "bold",
            col = "black"
          )
        )

        bases <- base_matrix[i, ]
        colors <- unname(base_colors[bases])
        colors[is.na(colors)] <- ambiguity_color
        grid::grid.text(
          label = bases,
          x = grid::unit(x_positions, "inches"),
          y = grid::unit(rep(y_positions[[i]], block_width), "inches"),
          just = "centre",
          gp = grid::gpar(
            fontsize = font_size,
            fontfamily = "mono",
            fontface = "bold",
            col = colors
          )
        )
      }

      conserved_y <- first_row_y - number_of_sequences * row_step
      if (isTRUE(show_conservation)) {
        conserved <- apply(base_matrix, 2L, function(bases) {
          if (length(unique(bases)) == 1L) "*" else ""
        })
        grid::grid.text(
          label = conserved,
          x = grid::unit(x_positions, "inches"),
          y = grid::unit(rep(conserved_y, block_width), "inches"),
          gp = grid::gpar(
            fontsize = font_size,
            fontfamily = "mono",
            col = "#555555"
          )
        )
      }

      intermediate_ticks <- if (block_width >= 10L) {
        seq.int(10L, block_width, by = 10L)
      } else {
        integer()
      }
      tick_positions <- unique(c(1L, intermediate_ticks, block_width))
      tick_x <- x_positions[tick_positions]
      tick_labels <- coordinate_start +
        start_position +
        tick_positions -
        2L
      axis_y <- conserved_y - 0.22

      grid::grid.segments(
        x0 = grid::unit(tick_x, "inches"),
        x1 = grid::unit(tick_x, "inches"),
        y0 = grid::unit(axis_y + 0.05, "inches"),
        y1 = grid::unit(axis_y + 0.10, "inches"),
        gp = grid::gpar(col = "#555555")
      )
      grid::grid.text(
        label = tick_labels,
        x = grid::unit(tick_x, "inches"),
        y = grid::unit(rep(axis_y, length(tick_x)), "inches"),
        just = "top",
        gp = grid::gpar(
          fontsize = max(font_size - 2, 5),
          fontfamily = "mono",
          col = "#555555"
        )
      )
    }
  }

  grDevices::dev.off()
  device_open <- FALSE

  if (!isTRUE(quiet)) {
    message("PDF file: ", pdf_file)
  }

  invisible(pdf_file)
}

validate_alignment_pdf_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0) {
    stop("`", name, "` must be a single positive number.", call. = FALSE)
  }
  as.numeric(x)
}

validate_alignment_pdf_colors <- function(base_colors, ambiguity_color) {
  if (!is.character(base_colors) || length(base_colors) == 0L ||
      is.null(names(base_colors)) || anyNA(base_colors) ||
      any(!nzchar(names(base_colors))) || any(!nzchar(base_colors))) {
    stop(
      "`base_colors` must be a non-empty named character vector of colors.",
      call. = FALSE
    )
  }
  check_scalar_character(ambiguity_color, "ambiguity_color")

  tryCatch(
    grDevices::col2rgb(c(base_colors, ambiguity_color)),
    error = function(error) {
      stop(
        "`base_colors` and `ambiguity_color` must contain valid R colors.",
        call. = FALSE
      )
    }
  )
  invisible(TRUE)
}
