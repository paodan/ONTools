#' Export an IGV alignment snapshot from R
#'
#' `igv_snapshot()` writes an IGV batch script for a reference FASTA, BAM file,
#' genomic region, and output snapshot, then optionally runs IGV from R.
#'
#' @param genome_fasta Path to the reference FASTA file loaded by IGV.
#' @param bam Path to the coordinate-sorted BAM file loaded by IGV.
#' @param chr Reference sequence, chromosome, or contig name passed to IGV.
#' @param start Start coordinate of the region passed to IGV.
#' @param end End coordinate of the region passed to IGV.
#' @param out_dir Directory where IGV writes snapshot images.
#' @param format Snapshot output format used when `snapshot_name` is `NULL`.
#'   Supported values are `"png"`, `"pdf"`, `"svg"`, `"jpg"`, and `"jpeg"`.
#' @param snapshot_name Output snapshot filename. The extension should be
#'   match `format`. When `NULL`, a filename is generated from `chr`, `start`,
#'   `end`, and `format`.
#' @param igv Command or full path used to launch IGV. On Linux this is usually
#'   `igv.sh`; on macOS command-line IGV installations also provide `igv.sh`.
#' @param batch_file Path where the temporary IGV batch script is written. When
#'   `NULL`, a temporary `.igv` file is used.
#' @param sort Alignment sorting mode passed to IGV `sort`. Use `NULL` to skip
#'   sorting.
#' @param display Alignment display mode. One of `"collapse"`, `"expand"`,
#'   `"squish"`, or `NULL`.
#' @param max_panel_height Optional value for IGV `maxPanelHeight`.
#' @param index Optional BAM index path. When supplied, the batch script uses
#'   `load bam index=index`.
#' @param auto_index_fasta Logical. If `TRUE`, create a missing FASTA `.fai`
#'   index with [Rsamtools::indexFa()] before writing the IGV batch script.
#' @param extra_commands Character vector of additional IGV batch commands added
#'   after `load` and before `goto`.
#' @param use_xvfb Logical. If `TRUE`, run IGV through `xvfb-run -a`, useful on
#'   headless Linux servers.
#' @param dry_run Logical. If `TRUE`, write and return the batch script and
#'   command without launching IGV.
#' @param quiet Logical. If `TRUE`, suppress console output from FASTA indexing
#'   and the IGV system command. Set to `FALSE` to show IGV logs while debugging.
#' @param overwrite Logical. If `FALSE`, stop when `batch_file` already exists.
#'
#' @details
#' IGV must be installed separately. Download IGV Desktop from
#' \url{https://igv.org/doc/desktop/DownloadPage/}. The command-line launcher is
#' documented at \url{https://igv.org/doc/desktop/UserGuide/advanced/command_line/}
#' and batch commands are documented at
#' \url{https://igv.org/doc/desktop/UserGuide/tools/batch/}. On Linux, the IGV
#' archive includes `igv.sh`; make sure it is on your `PATH`, or pass the full
#' path to `igv`.
#'
#' On headless Linux servers, install a virtual framebuffer and use
#' `use_xvfb = TRUE`. For Debian/Ubuntu, this is usually:
#' `sudo apt update && sudo apt install xvfb`. Conda users can install IGV and
#' Xvfb in an isolated environment, for example with:
#' `mamba install -c bioconda igv -c conda-forge xorg-x11-server-Xvfb`.
#'
#' BAM files should be coordinate sorted and indexed. If `index` is not supplied,
#' IGV will look for a standard index such as `sample.bam.bai`. The FASTA should
#' also have a `.fai` index. By default, `igv_snapshot()` creates a missing
#' FASTA index with [Rsamtools::indexFa()], the R equivalent of
#' `samtools faidx ref.fasta`.
#'
#' @return Invisibly returns a list containing the batch file path, batch
#'   commands, system command, arguments, exit status, and expected snapshot path.
#'
#' @examples
#' ref <- tempfile(fileext = ".fasta")
#' bam <- tempfile(fileext = ".bam")
#' writeLines(c(">barcode09", "ACGTACGTACGT"), ref)
#' file.create(bam)
#' file.create(paste0(bam, ".bai"))
#'
#' res <- igv_snapshot(
#'   genome_fasta = ref,
#'   bam = bam,
#'   chr = "barcode09",
#'   start = 1,
#'   end = 12,
#'   dry_run = TRUE
#' )
#' res$commands
#'
#' @export
igv_snapshot <- function(genome_fasta,
                         bam,
                         chr,
                         start,
                         end,
                         out_dir = "alignment",
                         format = "png",
                         snapshot_name = NULL,
                         igv = "igv.sh",
                         batch_file = NULL,
                         sort = "base",
                         display = "collapse",
                         max_panel_height = NULL,
                         index = NULL,
                         auto_index_fasta = TRUE,
                         extra_commands = NULL,
                         use_xvfb = FALSE,
                         dry_run = FALSE,
                         quiet = TRUE,
                         overwrite = TRUE) {
  check_file_arg(genome_fasta, "genome_fasta")
  check_file_arg(bam, "bam")
  check_scalar_character(chr, "chr")
  check_region_coordinate(start, "start")
  check_region_coordinate(end, "end")
  check_scalar_character(out_dir, "out_dir")
  check_scalar_character(igv, "igv")
  format <- match.arg(format, c("png", "pdf", "svg", "jpg", "jpeg"))

  if (start > end) {
    stop("`start` must be less than or equal to `end`.", call. = FALSE)
  }

  locus <- paste0(chr, ":", start, "-", end)

  if (!is.null(index)) {
    check_file_arg(index, "index")
  }

  if (!is.logical(auto_index_fasta) || length(auto_index_fasta) != 1L ||
      is.na(auto_index_fasta)) {
    stop("`auto_index_fasta` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(batch_file)) {
    check_scalar_character(batch_file, "batch_file")
    if (file.exists(batch_file) && !isTRUE(overwrite)) {
      stop("`batch_file` already exists and `overwrite` is FALSE: ", batch_file,
           call. = FALSE)
    }
  } else {
    batch_file <- tempfile(fileext = ".igv")
  }

  if (!is.null(sort)) {
    check_scalar_character(sort, "sort")
  }

  if (!is.null(display)) {
    display <- match.arg(display, c("collapse", "expand", "squish"))
  }

  if (!is.null(max_panel_height)) {
    if (!is.numeric(max_panel_height) || length(max_panel_height) != 1L ||
        is.na(max_panel_height)) {
      stop("`max_panel_height` must be a single numeric value or `NULL`.",
           call. = FALSE)
    }
  }

  if (!is.null(extra_commands) && !is.character(extra_commands)) {
    stop("`extra_commands` must be a character vector or `NULL`.",
         call. = FALSE)
  }

  genome_fasta <- normalizePath(genome_fasta, mustWork = TRUE)
  bam <- normalizePath(bam, mustWork = TRUE)
  if (!is.null(index)) {
    index <- normalizePath(index, mustWork = TRUE)
  }
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
      "`. IGV may fail to load the reference; set `auto_index_fasta = TRUE` ",
      "or create one with `samtools faidx`.",
      call. = FALSE
    )
  }

  if (is.null(index) && !has_standard_bam_index(bam)) {
    warning(
      "No standard BAM index found next to `bam`. IGV may fail to load the BAM; ",
      "create one with `samtools index` or pass `index`.",
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

  load_command <- paste("load", bam)
  if (!is.null(index)) {
    load_command <- paste(load_command, paste0("index=", index))
  }

  commands <- c(
    "new",
    paste("genome", genome_fasta),
    paste("snapshotDirectory", out_dir),
    if (!is.null(max_panel_height)) {
      paste("maxPanelHeight", max_panel_height)
    },
    load_command,
    extra_commands,
    paste("goto", locus),
    if (!is.null(sort)) {
      paste("sort", sort)
    },
    display,
    paste("snapshot", snapshot_name),
    "exit"
  )
  commands <- commands[!vapply(commands, is.null, logical(1))]

  writeLines(commands, batch_file)

  if (isTRUE(use_xvfb)) {
    command <- "xvfb-run"
    args <- c("-a", igv, "-b", batch_file)
  } else {
    command <- igv
    args <- c("-b", batch_file)
  }

  status <- 0L
  system_stdout <- if (isTRUE(quiet)) FALSE else ""
  system_stderr <- if (isTRUE(quiet)) FALSE else ""
  if (!isTRUE(dry_run)) {
    status <- system2(
      command,
      args = args,
      stdout = system_stdout,
      stderr = system_stderr
    )
    if (!identical(status, 0L)) {
      stop("IGV command failed with exit status ", status, ".", call. = FALSE)
    }
  }

  invisible(list(
    batch_file = batch_file,
    commands = commands,
    command = command,
    args = args,
    status = status,
    format = format,
    quiet = quiet,
    stdout = system_stdout,
    stderr = system_stderr,
    fasta_index = fasta_index,
    fasta_index_created = fasta_index_created,
    snapshot = file.path(out_dir, snapshot_name)
  ))
}

check_file_arg <- function(x, name) {
  check_scalar_character(x, name)
  if (!file.exists(x)) {
    stop("`", name, "` does not exist: ", x, call. = FALSE)
  }
}

check_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a single non-empty character string.",
         call. = FALSE)
  }
}

check_region_coordinate <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be a single finite numeric coordinate.",
         call. = FALSE)
  }

  if (x < 1 || x != as.integer(x)) {
    stop("`", name, "` must be a positive integer coordinate.", call. = FALSE)
  }
}

has_standard_bam_index <- function(bam) {
  file.exists(paste0(bam, ".bai")) ||
    file.exists(sub("\\.bam$", ".bai", bam, ignore.case = TRUE))
}

sanitize_igv_snapshot_name <- function(locus) {
  gsub("[^A-Za-z0-9_.-]+", "_", locus)
}
