library(lipidr)
d <- read_skyline(c("A1.csv", "F2.csv"))
expect_valid_lipidex <- function(d, dim) {
  act <- quasi_label(rlang::enquo(d), arg = "object")
  act_val_lab <- paste(methods::is(d), collapse = "/")
  if (!is(act$val, "LipidomicsExperiment")) {
    fail(sprintf("%s inherits from `%s` not `%s`.",
      act$lab, act_val_lab, "LipidomicsExperiment"))
  }

  if (!is(act$val, "SummarizedExperiment")) {
    fail(sprintf("%s inherits from `%s` not `%s`.",
      act$lab, act_val_lab, "SummarizedExperiment"))
  }

  if (!validObject(d)) {
    fail(sprintf("%s is not a valid LipidomicsExperiment object.", act$lab))
  }
  expect_equal(dim(d), dim)
}
expect_assay_equal <- function(e1, e2, measure) {
  mat <- assay(e1, measure)
  rownames(mat) <- rowData(e1)$Molecule

  sum_mat <- assay(e2, measure)
  rownames(sum_mat) <- rowData(e2)$Molecule
  expect_equal(mat, sum_mat[rownames(mat), ] )
}

context("test-normalize")
test_that("Can be applied to both summarized and unsuammrized datasets", {
  expect_valid_lipidex(d, c(22,11))
  expect_valid_lipidex(
    suppressWarnings(normalize_pqn(d, measure = "Area", exclude=NULL, log = TRUE)),
    c(22,11)
  )
  expect_valid_lipidex(
    suppressWarnings(normalize_istd(d, measure = "Area", exclude=NULL, log = TRUE)),
    c(22,11)
  )

  d_sum <- summarize_transitions(d)
  expect_valid_lipidex(
    suppressWarnings(normalize_pqn(d_sum, measure = "Area", exclude=NULL, log = TRUE)),
    c(22,11)
  )
  expect_valid_lipidex(
    suppressWarnings(normalize_istd(d_sum, measure = "Area", exclude=NULL, log = TRUE)),
    c(22,11)
  )
})

test_that("Can normalize the selected measure only", {
  expect_valid_lipidex(d, c(22,11))
  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", exclude=NULL, log = TRUE))
  expect_assay_equal(d, d_norm, "Retention Time")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Area"))

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", exclude=NULL, log = TRUE))
  expect_assay_equal(d, d_norm, "Retention Time")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Area"))

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Retention Time", exclude=NULL))
  expect_assay_equal(d, d_norm, "Area")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Retention Time"))

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Retention Time", exclude=NULL))
  expect_assay_equal(d, d_norm, "Area")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Retention Time"))


  ###### ISTD norm
  d_norm <- suppressWarnings(normalize_istd(d, measure = "Area", exclude=NULL, log = TRUE))
  expect_assay_equal(d, d_norm, "Retention Time")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Area"))

  d_norm <- suppressWarnings(normalize_istd(d, measure = "Area", exclude=NULL, log = TRUE))
  expect_assay_equal(d, d_norm, "Retention Time")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Area"))

  d_norm <- suppressWarnings(normalize_istd(d, measure = "Retention Time", exclude=NULL))
  expect_assay_equal(d, d_norm, "Area")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Retention Time"))

  d_norm <- suppressWarnings(normalize_istd(d, measure = "Retention Time", exclude=NULL))
  expect_assay_equal(d, d_norm, "Area")
  expect_assay_equal(d, d_norm, "Background")
  expect_failure(expect_assay_equal(d, d_norm, "Retention Time"))
})

test_that("Excludes blanks by default", {
  expect_valid_lipidex(d, c(22,11))
  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", exclude=NULL, log = TRUE))
  expect_valid_lipidex(d_norm, c(22,11))
  expect_true(all(c("Blank_1", "Blank_2") %in% colnames(d_norm)))

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", log = TRUE))
  expect_valid_lipidex(d_norm, c(22,9))
  expect_false(any(c("Blank_1", "Blank_2") %in% colnames(d_norm)))
})

test_that("Excludes blanks for other measures", {
  d2 <- d
  assay(d2, "newAssay") <- assay(d2, "Area")
  d2 <- d2 %>% set_normalized("newAssay", FALSE) %>% set_logged("newAssay", FALSE)
  d_norm <- suppressWarnings(normalize_pqn(d2, measure = "newAssay", log = TRUE))
  expect_valid_lipidex(d_norm, c(22,9))
  expect_false(any(c("Blank_1", "Blank_2") %in% colnames(d_norm)))
})

test_that("Gives error if measure is not in dataset", {
  expect_error(
    normalize_pqn(d, measure = "newAssay", log = TRUE),
    'newAssay is not in the dataset.'
  )
})
test_that("Can exclude selected samples", {
  expect_valid_lipidex(d, c(22,11))
  excluded <- c("TQC_1", "TQC_2", "TQC_3")

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", exclude=excluded, log = TRUE))
  expect_valid_lipidex(d_norm, c(22,8))
  expect_false(any(excluded %in% colnames(d_norm)))
  expect_true(all(c("Blank_1", "Blank_2") %in% colnames(d_norm)))
})

test_that("Values in excluded samples donot affect normlization", {
  d2 <- d
  expect_valid_lipidex(d2, c(22,11))
  excluded <- c("TQC_1", "TQC_2", "TQC_3")
  assay(d2)[,excluded] <- assay(d2)[,excluded] * 10

  d_norm <- suppressWarnings(normalize_pqn(d, measure = "Area", exclude=excluded, log = TRUE))
  d2_norm <- suppressWarnings(normalize_pqn(d2, measure = "Area", exclude=excluded, log = TRUE))
  expect_assay_equal(d_norm, d2_norm, "Area")
})

test_that("Gives error when all samples are excluded", {
  excluded <- colnames(d)

  expect_error(
    normalize_pqn(d, measure = "Area", exclude=excluded, log = TRUE),
    'You cannot exclude all samples'
  )
})

test_that("PQN can correct for dilution factors", {
  d2 <- cbind(d, d[,1])
  colnames(d2)[[12]] <- 'Snew'

  m <- assay(d2, "Area")
  m[, 12] <- m[, 12] * 10
  assay(d2, "Area") <- m
  expect_valid_lipidex(d2, c(22,12))

  m <- assay(d2, "Area")
  expect_equal(m[, "S1A"], m[, "Snew"] / 10)

  d_norm <- suppressWarnings(normalize_pqn(d2, measure = "Area", log = TRUE))
  expect_valid_lipidex(d_norm, c(22,10))
  m <- assay(d_norm, "Area")
  expect_equal(m[, "S1A"], m[, "Snew"])
})

test_that("ISTD can correct for specific classes", {
  d2 <- cbind(d, d[,1])
  colnames(d2)[[12]] <- 'Snew'

  istd <- rowData(d)$istd
  m <- assay(d2, "Area")
  m[istd, ] <- 1
  assay(d2, "Area") <- m

  d_norm <- suppressWarnings(normalize_istd(d2, measure = "Area", log = FALSE))
  expect_valid_lipidex(d_norm, c(22,10))
  expect_assay_equal(d_norm, d2[, colnames(d_norm)], "Area")

  LPEs <- rowData(d)$Class == "LPE"
  m[, 1] <- m[, 1] * 10
  m[istd & !LPEs, 1] <- 1
  assay(d2, "Area") <- m
  d_norm <- suppressWarnings(normalize_istd(d2, measure = "Area", log = FALSE))
  m <- assay(d_norm, "Area")
  expect_equal(m[LPEs, "S1A"], m[LPEs, "Snew"]) # corrected
  expect_equal(m[!LPEs & !istd, "S1A"] / 10, m[!LPEs & !istd, "Snew"]) # not corrected
})

test_that("ISTD will error if no istd molecules are found", {
  d2 <- d
  rowData(d2)$istd <- FALSE
  expect_error(normalize_istd(d2), 'No internal standards found in your lipid list.')
})

test_that("Normalization will give error if measure is normalized", {
  d2 <- normalize_pqn(d, measure="Area")
  expect_error(normalize_pqn(d2, measure="Area"), 'Area is already normalized.')
  expect_error(normalize_istd(d2, measure="Area"), 'Area is already normalized.')

  assay(d2, "Area2") <- assay(d2, "Area")
  d2 <- normalize_pqn(d2, measure = "Area2", log = FALSE)
})
