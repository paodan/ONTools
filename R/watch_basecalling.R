#' List sequencing files in a directory
#'
#' `watch_file()` lists files with a requested extension and reports the file
#' stem, size, and absolute path. It is primarily used by
#' [watch_basecalling()] to compare POD5 inputs with completed BAM outputs.
#'
#' @param path Directory to scan.
#' @param extension File extension, including the leading dot.
#'
#' @return A data frame with columns `name`, `size`, and `path`. `name` is the
#'   basename without `extension` and `size` is measured in bytes.
#'
#' @examples
#' path <- tempfile("watch-file-")
#' dir.create(path)
#' file.create(file.path(path, c("batch01.pod5", "batch02.pod5")))
#' watch_file(path)
#'
#' @export
watch_file <- function(path, extension = ".pod5") {
  check_dir_arg(path, "path")
  check_scalar_character(extension, "extension")

  if (!startsWith(extension, ".") || nchar(extension) < 2L) {
    stop("`extension` must start with '.' and include a suffix.", call. = FALSE)
  }

  path <- normalizePath(path, mustWork = TRUE)
  files <- list.files(path, full.names = TRUE, no.. = TRUE)
  files <- files[endsWith(basename(files), extension)]

  if (length(files) == 0L) {
    return(data.frame(
      name = character(),
      size = numeric(),
      path = character(),
      stringsAsFactors = FALSE
    ))
  }

  info <- file.info(files)
  keep <- !is.na(info$isdir) & !info$isdir
  files <- files[keep]
  info <- info[keep, , drop = FALSE]
  file_names <- basename(files)
  stems <- substr(file_names, 1L, nchar(file_names) - nchar(extension))

  data.frame(
    name = stems,
    size = info$size,
    path = normalizePath(files, mustWork = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Basecall one POD5 file with Dorado
#'
#' `basecall_one()` runs `dorado basecaller` for one POD5 file. Dorado writes
#' BAM data to standard output, so the function first captures it in a partial
#' file and only renames it to `bam_file` after a successful, non-empty run.
#'
#' @param pod5_file Input POD5 file.
#' @param bam_file Output BAM file. An existing file is never overwritten.
#' @param model Dorado model name, model alias such as `"fast"`, `"hac"`, or
#'   `"sup"`, or a path to a downloaded model.
#' @param kit_name Optional Dorado barcode kit name. When supplied, barcode
#'   classification is performed during basecalling.
#' @param sample_sheet Optional MinKNOW-compatible sample sheet. Requires
#'   `kit_name` and applies barcode filtering and aliases during basecalling.
#' @param dorado Command name or executable path for Dorado.
#' @param conda_env Optional Conda environment name. When supplied, Dorado is
#'   run using `conda run --no-capture-output`.
#' @param conda Conda command name or executable path.
#' @param dry_run Logical. If `TRUE`, return the planned command without
#'   running Dorado or creating the BAM file.
#' @param echo Logical. If `TRUE`, print the planned command.
#' @param stderr Passed to [system2()] for Dorado's standard error stream. Use
#'   `""` to show progress in the R console or a file path to save a log.
#'
#' @return Invisibly returns a list with `status`, `command`, `args`,
#'   `command_string`, `paths`, and `conda_env`. In dry-run mode, `status` is
#'   `NA_integer_`.
#'
#' @examples
#' pod5 <- tempfile(fileext = ".pod5")
#' file.create(pod5)
#' basic <- basecall_one(
#'   pod5,
#'   tempfile(fileext = ".bam"),
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' basic$args
#'
#' # Plan inline barcode classification and aliasing with a sample sheet.
#' sample_sheet <- tempfile(fileext = ".csv")
#' writeLines(
#'   c(
#'     "experiment_id,flow_cell_id,kit,barcode,alias",
#'     "run1,FLOWCELL01,SQK-NBD114-96,barcode01,sampleA",
#'     "run1,FLOWCELL01,SQK-NBD114-96,barcode02,sampleB"
#'   ),
#'   sample_sheet
#' )
#' classified <- basecall_one(
#'   pod5_file = pod5,
#'   bam_file = tempfile(fileext = ".bam"),
#'   model = "hac",
#'   kit_name = "SQK-NBD114-96",
#'   sample_sheet = sample_sheet,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#' classified$command_string
#'
#' @export
basecall_one <- function(pod5_file,
                         bam_file,
                         model = "hac",
                         dorado = "dorado",
                         conda_env = NULL,
                         conda = "conda",
                         dry_run = FALSE,
                         echo = TRUE,
                         stderr = "",
                         kit_name = NULL,
                         sample_sheet = NULL) {
  check_file_arg(pod5_file, "pod5_file")
  check_scalar_character(bam_file, "bam_file")
  check_scalar_character(model, "model")
  if (!is.null(kit_name)) check_scalar_character(kit_name, "kit_name")
  if (!is.null(sample_sheet)) {
    if (is.null(kit_name)) {
      stop("`sample_sheet` requires `kit_name`.", call. = FALSE)
    }
    check_file_arg(sample_sheet, "sample_sheet")
    sample_sheet <- normalizePath(sample_sheet, mustWork = TRUE)
  }
  check_scalar_character(dorado, "dorado")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  validate_system2_stderr(stderr)

  pod5_file <- normalizePath(pod5_file, mustWork = TRUE)
  bam_dir <- dirname(bam_file)
  dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)
  bam_file <- normalizePath(bam_file, mustWork = FALSE)

  if (file.exists(bam_file)) {
    stop("`bam_file` already exists: ", bam_file, call. = FALSE)
  }

  dorado_args <- c("basecaller", model, pod5_file)
  if (!is.null(kit_name)) {
    dorado_args <- c(dorado_args, "--kit-name", kit_name)
  }
  if (!is.null(sample_sheet)) {
    dorado_args <- c(dorado_args, "--sample-sheet", sample_sheet)
  }
  if (is.null(conda_env)) {
    call <- list(command = dorado, args = dorado_args)
  } else {
    call <- list(
      command = conda,
      args = c(
        "run", "--no-capture-output", "-n", conda_env, "--",
        dorado, dorado_args
      )
    )
  }

  command_string <- paste(
    c(shQuote(call$command), shQuote(call$args), ">", shQuote(bam_file)),
    collapse = " "
  )
  if (is.character(stderr) && length(stderr) == 1L && nzchar(stderr)) {
    command_string <- paste(command_string, "2>", shQuote(stderr))
  }
  if (isTRUE(echo)) message(command_string)

  paths <- list(
    pod5_file = pod5_file,
    bam_file = bam_file,
    sample_sheet = sample_sheet
  )
  result <- list(
    status = NA_integer_,
    command = call$command,
    args = call$args,
    command_string = command_string,
    paths = paths,
    conda_env = conda_env
  )
  if (isTRUE(dry_run)) return(invisible(result))

  if (is.null(conda_env)) {
    require_external_command(dorado)
  } else {
    require_external_command(conda)
  }

  partial_bam <- tempfile(
    pattern = paste0(basename(bam_file), "."),
    tmpdir = dirname(bam_file),
    fileext = ".partial"
  )
  on.exit(unlink(partial_bam), add = TRUE)

  status <- suppressWarnings(system2(
    command = call$command,
    args = shQuote(call$args),
    stdout = partial_bam,
    stderr = stderr
  ))

  if (!identical(status, 0L)) {
    stop("dorado basecaller failed with exit status: ", status,
         call. = FALSE)
  }
  if (!file.exists(partial_bam) || is.na(file.info(partial_bam)$size) ||
      file.info(partial_bam)$size == 0) {
    stop("dorado basecaller did not create a non-empty BAM: ", bam_file,
         call. = FALSE)
  }
  if (!file.rename(partial_bam, bam_file)) {
    stop("Unable to move the completed BAM to: ", bam_file, call. = FALSE)
  }

  result$status <- status
  invisible(result)
}

#' Watch a POD5 directory and basecall new files
#'
#' `watch_basecalling()` repeatedly compares POD5 files with same-stem BAM
#' files. For each new POD5 file, it waits until the file size is unchanged
#' across two checks and then calls [basecall_one()].
#'
#' @param pod5_path Directory containing POD5 files.
#' @param bam_path Directory in which same-stem BAM files are written.
#' @param max_run_time Maximum watcher runtime in seconds. A basecalling job
#'   already in progress is allowed to finish after this limit.
#' @param check_interval Seconds between consecutive POD5 size checks.
#' @param idle_interval Maximum seconds to wait when no new POD5 file exists.
#' @inheritParams basecall_one
#'
#' @return Invisibly returns a list containing elapsed time and the results
#'   returned by [basecall_one()] for files processed during this invocation.
#'
#' @details
#' Existing same-stem BAM files are treated as completed outputs. In dry-run
#' mode, commands are planned immediately for the POD5 files present at the
#' time of the call; the watcher does not sleep or wait for file stability.
#'
#' @examples
#' pod5_dir <- tempfile("pod5-")
#' bam_dir <- tempfile("bam-")
#' dir.create(pod5_dir)
#' file.create(file.path(pod5_dir, "batch01.pod5"))
#' watch_basecalling(
#'   pod5_dir,
#'   bam_dir,
#'   dry_run = TRUE,
#'   echo = FALSE
#' )
#'
#' @export
watch_basecalling <- function(pod5_path,
                              bam_path,
                              max_run_time = 10,
                              check_interval = 3,
                              idle_interval = 30,
                              model = "hac",
                              dorado = "dorado",
                              conda_env = NULL,
                              conda = "conda",
                              dry_run = FALSE,
                              echo = TRUE,
                              stderr = "") {
  check_dir_arg(pod5_path, "pod5_path")
  check_scalar_character(bam_path, "bam_path")
  max_run_time <- validate_watch_duration(max_run_time, "max_run_time")
  check_interval <- validate_watch_duration(check_interval, "check_interval")
  idle_interval <- validate_watch_duration(idle_interval, "idle_interval")
  check_scalar_character(model, "model")
  check_scalar_character(dorado, "dorado")
  check_scalar_character(conda, "conda")
  check_logical_scalar(dry_run, "dry_run")
  check_logical_scalar(echo, "echo")
  if (!is.null(conda_env)) check_scalar_character(conda_env, "conda_env")
  validate_system2_stderr(stderr)

  pod5_path <- normalizePath(pod5_path, mustWork = TRUE)
  dir.create(bam_path, recursive = TRUE, showWarnings = FALSE)
  bam_path <- normalizePath(bam_path, mustWork = TRUE)

  start_time <- Sys.time()
  elapsed_time <- function() {
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  }
  results <- list()

  repeat {
    elapsed <- elapsed_time()
    if (elapsed >= max_run_time && !isTRUE(dry_run)) break
    if (isTRUE(echo)) message("Watcher elapsed time: ", round(elapsed, 1), " seconds")

    pod5_files <- watch_file(pod5_path, ".pod5")
    bam_files <- watch_file(bam_path, ".bam")
    new_names <- setdiff(pod5_files$name, bam_files$name)

    if (isTRUE(dry_run)) {
      for (file_name in new_names) {
        input <- pod5_files$path[pod5_files$name == file_name][[1L]]
        output <- file.path(bam_path, paste0(file_name, ".bam"))
        results[[file_name]] <- basecall_one(
          input, output, model, dorado, conda_env, conda,
          dry_run = TRUE, echo = echo, stderr = stderr
        )
      }
      break
    }

    if (length(new_names) == 0L) {
      remaining <- max_run_time - elapsed_time()
      wait_time <- min(idle_interval, max(0, remaining))
      if (wait_time > 0) {
        if (isTRUE(echo)) message("No new POD5 files; waiting ", round(wait_time, 1), " seconds")
        Sys.sleep(wait_time)
      }
      next
    }

    for (file_name in new_names) {
      if (elapsed_time() >= max_run_time) break

      input <- pod5_files$path[pod5_files$name == file_name][[1L]]
      last_size <- file.info(input)$size
      stable <- FALSE

      while (elapsed_time() < max_run_time) {
        remaining <- max_run_time - elapsed_time()
        Sys.sleep(min(check_interval, max(0, remaining)))
        current_size <- file.info(input)$size

        if (length(current_size) == 0L || is.na(current_size)) {
          if (isTRUE(echo)) message("POD5 file disappeared: ", input)
          break
        }
        if (isTRUE(current_size == last_size)) {
          stable <- TRUE
          break
        }

        if (isTRUE(echo)) message("POD5 file is still growing: ", input)
        last_size <- current_size
      }

      if (!stable || elapsed_time() >= max_run_time) next

      output <- file.path(bam_path, paste0(file_name, ".bam"))
      results[[file_name]] <- basecall_one(
        input, output, model, dorado, conda_env, conda,
        dry_run = FALSE, echo = echo, stderr = stderr
      )
    }
  }

  invisible(list(
    elapsed_time = elapsed_time(),
    basecalls = results,
    paths = list(pod5_path = pod5_path, bam_path = bam_path)
  ))
}

validate_watch_duration <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0) {
    stop("`", name, "` must be a single positive number.", call. = FALSE)
  }
  x
}

validate_system2_stderr <- function(x) {
  valid <- (is.logical(x) && length(x) == 1L && !is.na(x)) ||
    (is.character(x) && length(x) == 1L && !is.na(x))
  if (!valid) {
    stop("`stderr` must be a single character path or logical value.",
         call. = FALSE)
  }
  invisible(x)
}
