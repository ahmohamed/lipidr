library(lipidr)
# Generate subset of data
# rsubset = rowData(d)[rowData(d)$Class %in% c("PE","LPE"),]
# mols = sample(rsubset$Molecule, 20)
# mols = unique(c(mols, rsubset$Molecule[rsubset$istd]))
# mol_list = c('PE 40:0', 'LPE 20:5', 'PE(O-40:6)', 'PE(O-36:2)', 'PE 36:0', 'PE 38:5', 'PE 34:1', 'PE(P-36:1)', 'PE(P-34:1)', 'PE 38:2', 'LPE 16:0', 'PE(O-34:0)', 'PE(O-38:5)', 'PE 36:4', 'PE(P-34:2)', 'PE(P-36:2)', 'PE 36:3', 'PE(O-38:3)', 'PE(P-36:4)', 'PE(O-36:5)', '15:0-18:1(d7) PE', '18:1(d7) Lyso PE')
# read.csv("inst/extdata/A1_data.csv") %>%
#   filter(
#     Peptide %in% mols,
#     grepl("S[1-3][AB]|TQC_[1-3]$|Blank", Replicate)
#   ) %>%
#   write.csv("tests/testthat/A1.csv")
#
# read.csv("inst/extdata/F2_data.csv") %>%
#   filter(
#     Peptide %in% mols,
#     grepl("S[1-3][AB]|TQC_[1-3]$|Blank", Replicate)
#   ) %>%
#   write.csv("tests/testthat/F2.csv")
#
# read.csv("tests/testthat/A1_skyline_export_Transition Results_pivot.csv") %>%
#   filter(
#     Peptide %in% mols
#   ) %>%
#   select(c(1:7, grep("A1_S[1-3][AB]|A1_TQC_[1-3]\\D|A1_Blank", colnames(.)))) %>%
#   write.csv("tests/testthat/A1_pivot.csv")

f1 <- "A1.csv"
f2 <- "F2.csv"
f3 <- "A1_pivot.csv"

save_temp_csv <- function(d) {
  file <- tempfile(fileext = ".csv")
  write.csv(d, file = file, row.names = FALSE)
  file
}

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

context("test-readfiles-nopivot")
test_that("Can read a single non-pivoted skyline file", {
  d <- read_skyline(f1)
  expect_valid_lipidex(d, c(19, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), f1)
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read multiple non-pivoted skyline file", {
  d <- read_skyline(c(f1, f2))

  expect_s4_class(d, "LipidomicsExperiment")
  expect_true(validObject(d))
  expect_s4_class(d, "SummarizedExperiment")

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c(f1, f2))
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read data.frame input", {
  d <- read_skyline(data.table::fread(f1) %>% as.data.frame)
  expect_valid_lipidex(d, c(19, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c("dataset"))
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read multiple data.frame inputs", {
  df1 <- data.table::fread(f1) %>% as.data.frame
  df2 <- data.table::fread(f2) %>% as.data.frame
  d <- read_skyline(list(df1, df2))
  expect_valid_lipidex(d, c(22, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c("dataset 1", "dataset 2"))
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read mixed files and data.frame inputs", {
  df1 <- data.table::fread(f1) %>% as.data.frame
  d <- read_skyline(list(df1, f2))
  expect_valid_lipidex(d, c(22, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c("dataset 1", f2))
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Will used provided names when provided", {
  df1 <- data.table::fread(f1) %>% as.data.frame
  d <- read_skyline(list(first=df1, second=f2))
  expect_valid_lipidex(d, c(22, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c("first", "second"))
  expect_equal(assayNames(d), c("Retention Time", "Area", "Background"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can handle files with different columns", {
  f_nobg <- f1 %>%
    data.table::fread() %>% as.data.frame %>%
    select(-Background) %>% save_temp_csv()
  expect_warning(d <- read_skyline(c(f_nobg, f2)), "Some columns were not available in all files")

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c(basename(f_nobg), f2))
  expect_equal(assayNames(d), c("Retention Time", "Area"))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

context("test-readfiles-error-checking")
test_that("Errors for empty files", {
  fempty <- f1 %>% data.table::fread() %>% as.data.frame %>% filter(FALSE) %>% save_temp_csv()
  expect_error(read_skyline(fempty), "does not have any data.$")
})

test_that("Errors when missing an intensity column", {
  ferror <- f1 %>% data.table::fread() %>% as.data.frame %>% select(-Area) %>% save_temp_csv()
  expect_error(read_skyline(ferror), "Data should have one of these measures exported")
})

test_that("Errors when missing a molecule column", {
  ferror <- f1 %>% data.table::fread() %>% as.data.frame %>% select(-Peptide) %>% save_temp_csv()
  expect_error(read_skyline(ferror), "Data should have one of these columns with molecule names")
})

test_that("Errors when missing a replicate column", {
  ferror <- f1 %>% data.table::fread() %>% as.data.frame %>% select(-Replicate) %>% save_temp_csv()
  expect_error(read_skyline(ferror), "either export Replicate column or pivot by replicates")
})

context("test-readfiles-pivoted")
test_that("Can read a single pivoted skyline file", {
  d <- read_skyline(f3)
  expect_valid_lipidex(d, c(19, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), f3)
  expect_true(all(c("Retention Time", "Area", "Background") %in% assayNames(d)))
  expect_false(any(grepl("^Replicate", assayNames(d))))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can read multiple pivoted skyline file", {
  f3_ <- f3 %>% data.table::fread() %>% as.data.frame %>% save_temp_csv()
  d <- read_skyline(c(f3_, f3))
  expect_valid_lipidex(d, c(38, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c(basename(f3_), f3))
  expect_true(all(c("Retention Time", "Area", "Background") %in% assayNames(d)))
  expect_false(any(grepl("^Replicate", assayNames(d))))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})

test_that("Can handle multiple pivoted files with different columns", {
  f_nobg <- f3 %>% data.table::fread() %>% as.data.frame %>% select(-matches("Background")) %>% save_temp_csv()
  expect_warning(d <- read_skyline(c(f_nobg, f3)), "Some columns were not available in all files")

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), c(basename(f_nobg), f3))
  expect_true(all(c("Retention Time", "Area") %in% assayNames(d)))
  expect_false(any(grepl("^Replicate", assayNames(d))))
  expect_false(any(grepl("Background", assayNames(d))))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
})
test_that("Can handle duplicate molecules", {
  f3_double <- f3 %>% data.table::fread() %>% as.data.frame %>% rbind(., .) %>% save_temp_csv()
  d <- read_skyline(f3_double)
  expect_valid_lipidex(d, c(38, 11))

  row_data <- rowData(d)
  expect_equal(unique(row_data$filename), basename(f3_double))
  expect_true(all(c("Retention Time", "Area", "Background") %in% assayNames(d)))
  expect_false(any(grepl("^Replicate", assayNames(d))))
  expect_equal(metadata(d)$dimnames, c("TransitionId", "Sample"))
  expect_false(metadata(d)$summarized)
  expect_false(any(unlist(mcols(assays(d)))))
  d_sum <- summarize_transitions(d)
  expect_valid_lipidex(d_sum, c(19, 11))
})

context("test-readfiles-sample_annot")
gen_sample_annot <- function(d) {
  samples <- colnames(d)
  groups <- c("A", "A", "A", "B", "B", "B", "Blank", "Blank", "QC", "QC", "QC")
  data.frame(Sample=samples, Group=groups, stringsAsFactors = FALSE)
}
test_that("Can add sample annotations using csv file", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  d <- add_sample_annotation(d, save_temp_csv(df))
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})
test_that("Will give error for non-tabular files", {
  file <- tempfile(fileext = ".csv")
  f2 %>% data.table::fread() %>% as.data.frame %>%
    write.table(file = file, row.names = FALSE)

  expect_error(read_skyline(file), 'Data should be tabular in CSV or tab-delimited format.')
})

test_that("Can handle files with missing values in columns", {
  file <- f1 %>% data.table::fread() %>% as.data.frame %>%
    mutate(newcol <- ifelse(1:nrow(.) %% 2 == 0,  NA, 1)) %>%
    save_temp_csv()

  expect_valid_lipidex(read_skyline(file), c(19, 11))
})

test_that("Can add sample annotations using data.frame", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  d <- add_sample_annotation(d, df)
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})

test_that("Can add sample annotations when samples are not in order", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df_shuffle <- df[sample(1:nrow(df), nrow(df), replace = FALSE), ]
  d <- add_sample_annotation(d, df_shuffle)
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})

test_that("Will generate warning when no column is named Sample", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df_rename <- df %>% rename(SS=Sample)
  expect_warning(d <- add_sample_annotation(d, df_rename), "No column named 'Sample'")
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})

test_that("Can add annotation when Sample column is not 1st", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df_reorder <- df[, c(2,1)]
  d <- add_sample_annotation(d, df_reorder)
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})

test_that("Will generate error when not all samples are present", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df_missing <- df[-1,]
  expect_error(d <- add_sample_annotation(d, df_missing), 'All sample names must be in the first column')
})

test_that("Can handle factors", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df_factor <- df %>% mutate_all(as.factor)
  d <- add_sample_annotation(d, df_factor)
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
  expect_true(is.character(col_data$Group))
})

test_that("Can handle duplicate annotation columns", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  d <- add_sample_annotation(d, df)
  expect_warning(d <- add_sample_annotation(d, df), 'These annotation columns already exist in data')
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
})

test_that("Can add non-duplicate annotation columns", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  d <- add_sample_annotation(d, df)
  df2 <- df %>% mutate(Group2=Group)
  expect_warning(d <- add_sample_annotation(d, df2), 'These annotation columns already exist in data')
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group", "Group2"))
  expect_equal(col_data$Group2, df$Group)
})

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
  expect_true(.have_lipids_molecules(first_col[, 1, drop=FALSE]), "first_col")
})
context("test-readfiles-todo")
test_that("Cannot mix pivoted and nonpivoted files", {
  # expect_error(d <- read_skyline(c(f2, f3)), "can't be converted")
})
