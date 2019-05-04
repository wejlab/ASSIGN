if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free(linters = lintr::with_defaults(
      camel_case_linter = NULL,
      snake_case_linter = NULL,
      commas_linter = NULL,
      line_length_linter = NULL,
      spaces_left_parentheses_linter = NULL,
      multiple_dots_linter = NULL))
  })
}
