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

test_that("Will give error for non-tabular files", {
  file <- f2 %>% data.table::fread() %>% as.data.frame %>%
    save_temp_csv(sep = " ")

  expect_error(read_skyline(file), 'Data should be tabular in CSV or tab-delimited format.')
})

context("test-readfiles-missingvals")
test_that("Can handle files with missing values in columns", {
  file <- f1 %>% .read_tabular() %>%
    mutate(newcol <- ifelse(1:nrow(.) %% 2 == 0,  NA, 1)) %>%
    save_temp_csv(quote = FALSE, na = "")

  expect_valid_lipidex(read_skyline(file), c(19, 11))
})

test_that("Reads empty columns as missing values", {
  file <- f1 %>% .read_tabular() %>%
    mutate(Peptide=ifelse(1:nrow(.) %% 2 == 0,  NA, Peptide)) %>%
    save_temp_csv(quote=FALSE, na = "")

  d <- .read_tabular(file)
  expect_equal(sum(is.na(d$Peptide)), floor(nrow(d)/2))
  expect_warning(d <- read_skyline(file))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(sum(is.na(rowData(d)$Molecule)), 9)
})

context("test-readfiles-todo")
test_that("Cannot mix pivoted and nonpivoted files", {
  # expect_error(d <- read_skyline(c(f2, f3)), "can't be converted")
})
