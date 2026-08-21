get_project_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

  if (length(script_arg) > 0L) {
    script_path <- normalizePath(
      sub("^--file=", "", script_arg[[1L]]),
      mustWork = TRUE
    )
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

assert_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) > 0L) {
    stop(
      "Missing R packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running the workflow.",
      call. = FALSE
    )
  }
}

assert_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Missing input file(s):\n- ",
      paste(missing, collapse = "\n- "),
      "\nSee data/README.md for the expected inputs.",
      call. = FALSE
    )
  }
}

ensure_parent_dirs <- function(paths) {
  invisible(lapply(unique(dirname(paths)), dir.create, recursive = TRUE, showWarnings = FALSE))
}

PROJECT_ROOT <- get_project_root()
source(file.path(PROJECT_ROOT, "config", "analysis_config.R"))
