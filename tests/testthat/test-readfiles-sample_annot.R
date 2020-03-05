f1 <- "A1.csv"
f2 <- "F2.csv"
f3 <- "A1_pivot.csv"

context("test-readfiles-sample_annot")
test_that("Can add sample annotations using csv file", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  d <- add_sample_annotation(d, save_temp_csv(df))
  expect_valid_lipidex(d, c(19, 11))

  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group"))
  expect_equal(col_data$Group, df$Group)
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

test_that("Can add sample annotations using tibble", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d) %>% as_tibble()
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

test_that("Reads empty columns as missing values", {
  d <- read_skyline(f1)
  df <- gen_sample_annot(d)
  df$newcol <- NA
  d <- add_sample_annotation(d, df)
  col_data <- colData(d)
  expect_equal(colnames(col_data), c("Group", "newcol"))
  expect_valid_lipidex(d, c(19, 11))
  expect_equal(sum(is.na(col_data$newcol)), 11)
})
