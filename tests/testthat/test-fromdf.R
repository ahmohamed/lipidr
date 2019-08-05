context("test-fromdf")

f1 = read_skyline("A1.csv") %>%
  summarize_transitions()
f1_matrix <-  f1 %>%
  assay("Area") %>%
  as.data.frame() %>%
  `rownames<-`(rowData(f1) %>% as.data.frame() %>% .$Molecule)


test_that("Can read dataframe with the first column as features, rest are samples", {
  f1df <- f1_matrix %>% as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(Features, everything())
  d <- as_skyline_experiment(f1df)
  expect_s4_class(d, "SkylineExperiment")
  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")
  expect_equal(dim(d), c(19, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d)$summarized)
  expect_false( any(unlist(mcols(assays(d)))) )
})

test_that("Can read dataframe with rownames as features, rest are samples", {
  f1df <- f1_matrix %>% as.data.frame()
  d <- as_skyline_experiment(f1df)
  expect_s4_class(d, "SkylineExperiment")
  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")
  expect_equal(dim(d), c(19, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("MoleculeId", "Sample"))
  expect_true(metadata(d)$summarized)
  expect_false( any(unlist(mcols(assays(d)))) )
})

test_that("Can read dataframe with duplicate molecules", {
  f1df <- f1_matrix %>% as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(Features, everything())

  f1df <- f1df %>% bind_rows(f1df[1:5, ])

  d <- as_skyline_experiment(f1df)
  expect_s4_class(d, "SkylineExperiment")
  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")
  expect_equal(dim(d), c(24, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), "dataframe")
  expect_equal(assayNames(d), "Area")
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false( any(unlist(mcols(assays(d)))) )
  expect_equal(dim(d %>% summarize_transitions()), c(19, 11))
})

test_that("Gives error when no lipid names are provided", {
  f1df <- f1_matrix %>% as.data.frame() %>%
    rownames_to_column("Features") %>%
    select(-Features)

  expect_error(as_skyline_experiment(f1df), "Data frame does not contain valid lipid names.")
})

context("test-fromdf-todo")
test_that("Skyline export as dataframe is not yet supported", {
  a <- read.csv("A1.csv")
  expect_error(as_skyline_experiment(a))
})
