context("test-fromdf")

f1 <- read_skyline("A1.csv") %>%
  summarize_transitions()
f1_matrix <- f1 %>%
  assay("Area") %>%
  as.data.frame() %>%
  `rownames<-`(rowData(f1) %>% as.data.frame() %>% .$Molecule)

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


test_that("Can read dataframe with the first column as features, rest are samples", {
  f1df <- f1_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(Features, everything())
  d <- as_lipidomics_experiment(f1df)
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), f1df$Features)

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read dataframe with rownames as features, rest are samples", {
  f1df <- f1_matrix %>% as.data.frame()
  d <- as_lipidomics_experiment(f1df)
  expect_valid_lipidex(d, c(19, 11))

  expect_equal(rownames(d), rownames(f1df))
  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read dataframe with duplicate molecules", {
  f1df <- f1_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(Features, everything())

  f1df <- f1df %>% bind_rows(f1df[1:5, ])

  d <- as_lipidomics_experiment(f1df)
  expect_valid_lipidex(d, c(24, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
  expect_valid_lipidex(d %>% summarize_transitions(), c(19, 11))
  # expect_equal(dim(d %>% summarize_transitions()), c(19, 11))
})

test_that("Gives error when no lipid names are provided", {
  f1df <- f1_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(-Features)

  expect_error(as_lipidomics_experiment(f1df), "Data frame does not contain valid lipid names.")
})

context("test-fromdf-skyline")
test_that("Can read Skyline export as dataframe", {
  a <- data.table::fread("A1.csv") %>% as.data.frame()
  expect_valid_lipidex(as_lipidomics_experiment(a), c(19, 11))
})
