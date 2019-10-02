library(lipidr)
d <- read_skyline(c("A1.csv", "F2.csv"))

expect_assay_equal <- function(e1, e2, measure) {
  mat <- assay(e1, measure)
  rownames(mat) <- rowData(e1)$Molecule

  sum_mat <- assay(e2, measure)
  rownames(sum_mat) <- rowData(e2)$Molecule
  expect_equal(mat, sum_mat[rownames(mat), ] )
}

context("test-update_molnames")
test_that("updates molecules names and reparses annotation", {
})

test_that("Preserves annotation for nonupdated molecules", {
})
