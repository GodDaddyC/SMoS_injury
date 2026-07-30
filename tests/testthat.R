library(testthat)

proj_root <- normalizePath(getwd())

suppressMessages({
  lapply(list.files(
    file.path(proj_root, "functions"), pattern = "\\.R$", full.names = TRUE
  ), source)
})

test_dir(file.path(proj_root, "tests", "testthat"))
