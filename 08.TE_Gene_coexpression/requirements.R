required_packages <- c(
  "data.table",
  "doParallel",
  "dplyr",
  "foreach",
  "ggplot2",
  "patchwork",
  "stringr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}
