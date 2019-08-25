library(lipidr)
d <- read_skyline(c("A1.csv", "F2.csv"))

expect_assay_equal <- function(e1, e2, measure) {
  mat <- assay(e1, measure)
  rownames(mat) <- rowData(e1)$Molecule

  sum_mat <- assay(e2, measure)
  rownames(sum_mat) <- rowData(e2)$Molecule
  expect_equal(mat, sum_mat[rownames(mat), ] )
}

context("test-summarize_transitiions")
test_that("does not have an effect when no multi-transition are present", {
  expect_s4_class(d, "LipidomicsExperiment")
  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")
  expect_false(any(duplicated(rowData(d)$Molecule)))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)

  d_summarized <- summarize_transitions(d)
  expect_equal(dim(d), dim(d_summarized))
  expect_equal(assayNames(d_summarized), c("Retention.Time", "Area", "Background"))
  expect_equal(metadata(d_summarized)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d_summarized)$summarized)
  expect_false(any(unlist(mcols(assays(d_summarized)))))

  # Check assay matrices
  expect_assay_equal(d, d_summarized, "Area")
  expect_assay_equal(d, d_summarized, "Retention.Time")
})

test_that("Can summarize multi-transition molecules", {
  dup <- d[1,]
  rowData(dup)$Product.Mz <- rowData(dup)$Product.Mz + 1

  area <- assay(dup, "Area")
  assay(dup, "Area") <- area * 3
  rownames(dup) <- rownames(d[2,])
  d[2,] <- dup

  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")

  d_summarized <- summarize_transitions(d, method = "max")
  expect_equal(dim(d)[[1]] - 1, dim(d_summarized)[[1]])
  expect_equal(assayNames(d_summarized), c("Retention.Time", "Area", "Background"))
  expect_equal(metadata(d_summarized)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d_summarized)$summarized)
  expect_false(any(unlist(mcols(assays(d_summarized)))))

  # Check assay matrices
  # Non-duplicated
  dup_mols <- rowData(dup)$Molecule
  d_sum_non_dup <- d_summarized[!rowData(d_summarized)$Molecule %in% dup_mols, ]
  d_non_dup <- d[!rowData(d)$Molecule %in% dup_mols, ]
  expect_assay_equal(d_non_dup, d_sum_non_dup, "Area")
  expect_assay_equal(d_non_dup, d_sum_non_dup, "Retention.Time")

  d_sum_dup <- d_summarized[rowData(d_summarized)$Molecule %in% dup_mols, ]
  area_dup <- assay(d_sum_dup, "Area")
  expect_equal(as.numeric(area_dup), as.numeric(area) * 3)

  d_summarized <- summarize_transitions(d, method = "average")
  d_sum_dup <- d_summarized[rowData(d_summarized)$Molecule %in% dup_mols, ]
  area_dup <- assay(d_sum_dup, "Area")
  expect_equal(as.numeric(area_dup), as.numeric(area) * 2)
})

