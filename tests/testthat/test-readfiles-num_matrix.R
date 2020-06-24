f1 <- "A1.csv"
f2 <- "F2.csv"
f3 <- "A1_pivot.csv"

context("test-readfiles-nummatrix")
test_that("Can read numeric matrix as LipidomicsExp", {
  d <- read_skyline(f1)

  # rownames are not molecules
  mat <- assay(d, "Area")
  expect_equal(.get_mol_dim(mat), "none")
  expect_equal(.get_mol_dim(t(mat)), "none")
  expect_equal(.get_mol_dim(rownames_to_column(as.data.frame(mat))), "none")

  # set rownames to Molecule names
  rownames(mat) <- rowData(d)$Molecule
  expect_equal(.get_mol_dim(mat), "row_names")
  expect_equal(.get_mol_dim(t(mat)), "col_names")
  expect_equal(.get_mol_dim(rownames_to_column(as.data.frame(mat))), "first_col")



  d <- as_lipidomics_experiment(mat)
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')

  d <- as_lipidomics_experiment(t(mat))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')

  d <- as_lipidomics_experiment(rownames_to_column(as.data.frame(mat)))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')


  ## Data.frame
  d <- as_lipidomics_experiment(mat %>% as.data.frame())
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')

  d <- as_lipidomics_experiment(t(mat) %>% as.data.frame())
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')
})

test_that("Can convert character matrix to numeric", {
  d <- read_skyline(f1)

  # rownames are not molecules
  mat <- assay(d, "Area")
  storage.mode(mat) <- storage.mode(mat) <- 'character'
  rownames(mat) <- rowData(d)$Molecule
  expect_equal(.get_mol_dim(mat), "row_names")
  expect_equal(.get_mol_dim(t(mat)), "col_names")
  expect_equal(.get_mol_dim(rownames_to_column(as.data.frame(mat))), "first_col")


  d <- as_lipidomics_experiment(mat)
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')


  d <- as_lipidomics_experiment(t(mat))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')

  d <- as_lipidomics_experiment(rownames_to_column(as.data.frame(mat)))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')
})

test_that("Can convert character matrix to numeric", {
  d <- read_skyline(f1)

  # rownames are not molecules
  mat <- assay(d, "Area")
  storage.mode(mat) <- storage.mode(mat) <- 'character'
  rownames(mat) <- rowData(d)$Molecule
  mat[1:26] <- letters #replace numerical values with characters
  expect_equal(.get_mol_dim(mat), "row_names")
  expect_equal(.get_mol_dim(t(mat)), "col_names")
  expect_equal(.get_mol_dim(rownames_to_column(as.data.frame(mat))), "first_col")


  expect_warning(d <- as_lipidomics_experiment(mat), 'NAs introduced by coercion')
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')
  expect_true(all(is.na(assay(d, "Area")[1:26])))

  expect_warning(d <- as_lipidomics_experiment(t(mat)), 'NAs introduced by coercion')
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')
  expect_true(all(is.na(assay(d, "Area")[1:26])))

  expect_warning(d <- as_lipidomics_experiment(rownames_to_column(as.data.frame(mat))), 'NAs introduced by coercion')
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(class(c(assay(d, "Area"))), 'numeric')
  expect_true(all(is.na(assay(d, "Area")[1:26])))
})

test_that("Gives error if dataset is not numeric", {
  d <- read_skyline(f1)

  # rownames are not molecules
  mat <- assay(d, "Area")
  storage.mode(mat) <- storage.mode(mat) <- 'character'
  rownames(mat) <- rowData(d)$Molecule
  mat[1:length(mat)] <- paste(c(mat), 's') #replace numerical values with characters
  expect_equal(.get_mol_dim(mat), "row_names")
  expect_equal(.get_mol_dim(t(mat)), "col_names")
  expect_equal(.get_mol_dim(rownames_to_column(as.data.frame(mat))), "first_col")


  expect_error(
    suppressWarnings(as_lipidomics_experiment(mat)),
    'Dataset is not numeric'
  )
  expect_error(
    suppressWarnings(as_lipidomics_experiment(t(mat))),
    'Dataset is not numeric'
  )
  expect_error(
    suppressWarnings(as_lipidomics_experiment(rownames_to_column(as.data.frame(mat)))),
    'Dataset is not numeric'
  )
})
test_that(".check_lipids_molecules can handle drop=FALSE", {
  d <- read_skyline(f1)
  mat <- assay(d, "Area")
  rownames(mat) <- rowData(d)$Molecule
  first_col <- rownames_to_column(as.data.frame(mat))
  expect_true(.have_lipids_molecules(first_col[, 1, drop=FALSE]))
})

test_that("Gives warnings when molecules are duplicated", {
  d <- read_skyline(f1)
  mat <- assay(d, "Area")
  rownames(mat) <- rowData(d)$Molecule
  first_col <- rownames_to_column(as.data.frame(mat))
  dup = rbind(first_col, first_col, first_col)
  expect_true(.have_lipids_molecules(dup$rowname))
  expect_equal(.get_mol_dim(dup), "first_col")
  expect_warning(d <- as_lipidomics_experiment(dup),
    '38 Duplicate lipid names detected.'
  )
  expect_valid_lipidex(d, c(19*3, 11))
})
