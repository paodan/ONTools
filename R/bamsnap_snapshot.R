#' Export a BamSnap alignment snapshot from R
#'
#' `bamsnap_snapshot()` runs the command-line BamSnap viewer for a reference
#' FASTA, BAM file, genomic region, and output image. Unlike IGV Desktop, this
#' does not require opening a graphical application.
#'
#' @param genome_fasta Path to the reference FASTA file passed to BamSnap.
#' @param bam Path to the coordinate-sorted BAM file passed to BamSnap.
#' @param chr Reference sequence, chromosome, or contig name.
#' @param start Optional start coordinate of the region. When both `start` and
#'   `end` are `NULL`, `chr` is passed as the whole target region.
#' @param end Optional end coordinate of the region. When both `start` and
#'   `end` are `NULL`, `chr` is passed as the whole target region.
#' @param out_dir Directory where the snapshot image is written.
#' @param format Snapshot output format. Supported values are `"png"` and
#'   `"jpg"`.
#' @param snapshot_name Output snapshot filename. The extension should match
#'   `format`. When `NULL`, a filename is generated from `chr`, `start`, `end`,
#'   and `format`.
#' @param bamsnap Command or full path used to launch BamSnap.
#' @param title Optional label shown for the BAM track.
#' @param draw Character vector of BamSnap tracks passed after `-draw`.
#' @param bamplot Character vector of BAM subtracks passed after `-bamplot`.
#' @param width,height Optional image width and height in pixels.
#' @param margin Optional genomic margin size passed to BamSnap.
#' @param no_title Logical. If `TRUE`, add `-no_title`.
#' @param no_target_line Logical. If `TRUE`, add `-no_target_line`.
#' @param save_image_only Logical. If `TRUE`, add `-save_image_only`.
#' @param process Number of BamSnap worker processes.
#' @param extra_args Character vector of additional raw BamSnap arguments.
#' @param auto_index_fasta Logical. If `TRUE`, create a missing FASTA `.fai`
#'   index with [Rsamtools::indexFa()] before running BamSnap.
#' @param dry_run Logical. If `TRUE`, return the command without launching
#'   BamSnap.
#' @param quiet Logical. If `TRUE`, suppress console output from FASTA indexing
#'   and the BamSnap system command.
#' @param overwrite Logical. If `FALSE`, stop when the output snapshot already
#'   exists.
#'
#' @details
#' BamSnap must be installed separately, for example with `pip install bamsnap`.
#' BAM files should be coordinate sorted and indexed with a standard index name
#' such as `sample.bam.bai`.
#'
#' @return Invisibly returns a list containing the command, arguments, exit
#'   status, region, output snapshot path, and FASTA indexing metadata.
#'
#' @examples
#' ref <- tempfile(fileext = ".fasta")
#' bam <- tempfile(fileext = ".bam")
#' writeLines(c(">barcode09", "ACGTACGTACGT"), ref)
#' file.create(bam)
#' file.create(paste0(bam, ".bai"))
#'
#' res <- bamsnap_snapshot(
#'   genome_fasta = ref,
#'   bam = bam,
#'   chr = "barcode09",
#'   start = 1,
#'   end = 12,
#'   dry_run = TRUE
#' )
#' res$args
#'
#' @export
bamsnap_snapshot <- function(genome_fasta,
                             bam,
                             chr,
                             start = NULL,
                             end = NULL,
                             out_dir = "alignment",
                             format = "png",
                             snapshot_name = NULL,
                             bamsnap = "bamsnap",
                             title = NULL,
                             draw = c("coordinates", "bamplot", "base"),
                             bamplot = c("coverage", "base", "read"),
                             width = 1000,
                             height = NULL,
                             margin = NULL,
                             no_title = FALSE,
                             no_target_line = FALSE,
                             save_image_only = TRUE,
                             process = 1L,
                             extra_args = NULL,
                             auto_index_fasta = TRUE,
                             dry_run = FALSE,
                             quiet = TRUE,
                             overwrite = TRUE) {
  check_file_arg(genome_fasta, "genome_fasta")
  check_file_arg(bam, "bam")
  check_scalar_character(chr, "chr")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(bamsnap, "bamsnap")
  format <- match.arg(format, c("png", "jpg"))

  if (xor(is.null(start), is.null(end))) {
    stop("`start` and `end` must both be supplied or both be `NULL`.",
         call. = FALSE)
  }

  if (!is.null(start)) {
    check_region_coordinate(start, "start")
    check_region_coordinate(end, "end")

    if (start > end) {
      stop("`start` must be less than or equal to `end`.", call. = FALSE)
    }
  }

  if (!is.null(title)) {
    check_scalar_character(title, "title")
  }

  check_character_vector_arg(draw, "draw", allow_null = TRUE)
  check_character_vector_arg(bamplot, "bamplot", allow_null = TRUE)
  check_optional_numeric_arg(width, "width")
  check_optional_numeric_arg(height, "height")
  check_optional_numeric_arg(margin, "margin")
  check_positive_integer_arg(process, "process")
  check_logical_arg(no_title, "no_title")
  check_logical_arg(no_target_line, "no_target_line")
  check_logical_arg(save_image_only, "save_image_only")
  check_logical_arg(auto_index_fasta, "auto_index_fasta")
  check_logical_arg(dry_run, "dry_run")
  check_logical_arg(quiet, "quiet")
  check_logical_arg(overwrite, "overwrite")
  check_character_vector_arg(extra_args, "extra_args", allow_null = TRUE)

  locus <- if (is.null(start)) {
    chr
  } else {
    paste0(chr, ":", start, "-", end)
  }

  genome_fasta <- normalizePath(genome_fasta, mustWork = TRUE)
  bam <- normalizePath(bam, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)

  fasta_index <- paste0(genome_fasta, ".fai")
  fasta_index_created <- FALSE
  if (!file.exists(fasta_index) && isTRUE(auto_index_fasta)) {
    if (isTRUE(quiet)) {
      suppressMessages(suppressWarnings(Rsamtools::indexFa(genome_fasta)))
    } else {
      Rsamtools::indexFa(genome_fasta)
    }
    fasta_index_created <- file.exists(fasta_index)
    if (!isTRUE(fasta_index_created)) {
      stop("Failed to create FASTA index with Rsamtools::indexFa(): ",
           fasta_index, call. = FALSE)
    }
  }

  if (!file.exists(fasta_index)) {
    warning(
      "No FASTA index found at `", fasta_index,
      "`. BamSnap may fail to load the reference; set ",
      "`auto_index_fasta = TRUE` or create one with `samtools faidx`.",
      call. = FALSE
    )
  }

  if (!has_standard_bam_index(bam)) {
    warning(
      "No standard BAM index found next to `bam`. BamSnap may fail to load the ",
      "BAM; create one with `samtools index`.",
      call. = FALSE
    )
  }

  if (is.null(snapshot_name)) {
    snapshot_name <- paste0(sanitize_igv_snapshot_name(locus), ".", format)
  } else {
    check_scalar_character(snapshot_name, "snapshot_name")
    snapshot_ext <- tolower(tools::file_ext(snapshot_name))
    if (!nzchar(snapshot_ext)) {
      snapshot_name <- paste0(snapshot_name, ".", format)
    } else if (!identical(snapshot_ext, format)) {
      stop(
        "`snapshot_name` extension must match `format`. Got `.",
        snapshot_ext,
        "` and format `",
        format,
        "`.",
        call. = FALSE
      )
    }
  }

  snapshot <- file.path(out_dir, snapshot_name)
  if (file.exists(snapshot) && !isTRUE(overwrite)) {
    stop("`snapshot` already exists and `overwrite` is FALSE: ", snapshot,
         call. = FALSE)
  }

  args <- c(
    "-bam", bam,
    "-pos", locus,
    "-ref", genome_fasta,
    "-out", snapshot,
    "-imagetype", format,
    "-process", as.character(process)
  )
  if (!is.null(title)) {
    args <- c(args, "-title", title)
  }
  if (!is.null(draw) && length(draw) > 0L) {
    args <- c(args, "-draw", draw)
  }
  if (!is.null(bamplot) && length(bamplot) > 0L) {
    args <- c(args, "-bamplot", bamplot)
  }
  if (!is.null(width)) {
    args <- c(args, "-width", as.character(width))
  }
  if (!is.null(height)) {
    args <- c(args, "-height", as.character(height))
  }
  if (!is.null(margin)) {
    args <- c(args, "-margin", as.character(margin))
  }
  if (isTRUE(no_title)) {
    args <- c(args, "-no_title")
  }
  if (isTRUE(no_target_line)) {
    args <- c(args, "-no_target_line")
  }
  if (isTRUE(save_image_only)) {
    args <- c(args, "-save_image_only")
  }
  if (isTRUE(quiet)) {
    args <- c(args, "-silence")
  }
  args <- c(args, extra_args)

  status <- 0L
  system_stdout <- if (isTRUE(quiet)) FALSE else ""
  system_stderr <- if (isTRUE(quiet)) FALSE else ""
  if (!isTRUE(dry_run)) {
    status <- system2(
      bamsnap,
      args = args,
      stdout = system_stdout,
      stderr = system_stderr
    )
    if (!identical(status, 0L)) {
      stop("BamSnap command failed with exit status ", status, ".",
           call. = FALSE)
    }
  }

  invisible(list(
    command = bamsnap,
    args = args,
    status = status,
    region = locus,
    format = format,
    quiet = quiet,
    stdout = system_stdout,
    stderr = system_stderr,
    fasta_index = fasta_index,
    fasta_index_created = fasta_index_created,
    snapshot = snapshot
  ))
}

check_logical_arg <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
}

check_character_vector_arg <- function(x, name, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(TRUE))
  }

  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    stop("`", name, "` must be a character vector",
         if (isTRUE(allow_null)) " or `NULL`" else "",
         ".", call. = FALSE)
  }

  invisible(TRUE)
}

check_optional_numeric_arg <- function(x, name) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }

  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0) {
    stop("`", name, "` must be a single positive numeric value or `NULL`.",
         call. = FALSE)
  }

  invisible(TRUE)
}

check_positive_integer_arg <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}
