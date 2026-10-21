test_that("create_standard_directories creates work and delivery layouts", {
  root <- tempfile("standard-directories-")
  paths <- standard_paths(
    input_path_root = file.path(root, "input"),
    delivery_path_root = file.path(root, "delivery-root"),
    delivery_folder = "project01"
  )

  result <- create_standard_directories(paths)

  expect_true(all(result$created))
  expect_true(all(result$exists_after))
  expect_true(all(result$status == "created"))
  expect_true(all(dir.exists(result$path)))
  expect_false(any(result$name == "input_path_root"))

  second <- create_standard_directories(paths)
  expect_true(all(second$existed_before))
  expect_false(any(second$created))
  expect_true(all(second$status == "exists"))
})

test_that("create_standard_directories supports dry-run planning", {
  root <- tempfile("standard-directories-dry-")
  paths <- standard_paths(
    input_path_root = file.path(root, "input"),
    delivery_path_root = file.path(root, "delivery-root"),
    delivery_folder = "project01"
  )

  result <- create_standard_directories(paths, dry_run = TRUE)

  expect_false(any(result$existed_before))
  expect_false(any(result$created))
  expect_false(any(result$exists_after))
  expect_true(all(result$status == "planned"))
  expect_false(dir.exists(root))
})

test_that("create_standard_directories can create only one path group", {
  root <- tempfile("standard-directories-work-")
  paths <- standard_paths_gene_editing(
    input_path_root = file.path(root, "input"),
    delivery_path_root = file.path(root, "delivery-root")
  )

  result <- create_standard_directories(
    paths,
    create_delivery = FALSE
  )

  expect_true(all(startsWith(result$name, "work")))
  expect_true(all(dir.exists(result$path)))
  expect_false(dir.exists(paths$delivery_path_root))
  expect_false(any(result$name %in% c(
    "ref_sanitized", "xlsx_file_full", "xlsx_file_filtered"
  )))
})

test_that("create_standard_directories checks all file collisions first", {
  root <- tempfile("standard-directories-collision-")
  dir.create(root)
  paths <- standard_paths(
    input_path_root = file.path(root, "input"),
    delivery_path_root = file.path(root, "delivery-root"),
    delivery_folder = "project01"
  )
  dir.create(dirname(paths$work_qc_path), recursive = TRUE)
  file.create(paths$work_qc_path)

  expect_error(
    create_standard_directories(paths),
    "exists as a file.*QC"
  )
  expect_false(dir.exists(paths$work_bam_path))
  expect_false(dir.exists(paths$delivery_path))
})

test_that("create_standard_directories validates its inputs", {
  expect_error(
    create_standard_directories(list()),
    "named list"
  )
  expect_error(
    create_standard_directories(
      list(work_path = "work"),
      create_delivery = FALSE
    ),
    "missing required directory"
  )

  paths <- standard_paths("input", "delivery", "project01")
  expect_error(
    create_standard_directories(
      paths,
      create_work = FALSE,
      create_delivery = FALSE
    ),
    "At least one"
  )
  expect_error(
    create_standard_directories(paths, dry_run = NA),
    "dry_run"
  )
})
