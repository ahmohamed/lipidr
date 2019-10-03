library(lipidr)
d <- read_skyline(c("A1.csv", "F2.csv"))

expect_assay_equal <- function(e1, e2, measure) {
  mat <- assay(e1, measure)
  rownames(mat) <- rowData(e1)$Molecule

  sum_mat <- assay(e2, measure)
  rownames(sum_mat) <- rowData(e2)$Molecule
  expect_equal(mat, sum_mat[rownames(mat), ] )
}

context("test-mw")
test_that("can list studies for a specific keyword", {
  skip_if_offline("www.metabolomicsworkbench.org")
  res <- list_mw_studies(keyword = "lipidomics")
  expect_gt(nrow(res), 30)
  expect_equal(ncol(res), 9)
  na_cols <- sapply(res, function(x) all(is.na(x)))
  expect_false(any(na_cols))
})
test_that("can fetch study by id", {
  skip_if_offline("www.metabolomicsworkbench.org")
  expect_warning(
    res <- fetch_mw_study("ST001111"),
    "Some lipid names couldn't be parsed"
  )
  expect_valid_lipidex(c(777, 118))
  row_data <- rowData(res)
  expect_equal(unique(row_data$filename), c("AN001805", "AN001806"))
  expect_lt(sum(row_data$not_matched), 10)
})

