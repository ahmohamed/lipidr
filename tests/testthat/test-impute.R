file = .read_tabular("A1.csv") %>%
  mutate(Area=ifelse(1:nrow(.) %% 10 == 0,  NA, Area)) %>%
  save_temp_csv(quote=FALSE, na = "")

d <- read_skyline(file) %>% summarize_transitions()

test_that("Can impute with different methods", {
  expect_valid_lipidex(d, c(19, 11))
  mat <- assay(d, "Area")
  expect_true(any(is.na(mat)))

  methods <- list(
    impute_na(d, "Area", "zero"),
    impute_na(d, "Area", "minProb"),
    impute_na(d, "Area", "minDet"),
    impute_na(d, "Area", "QRILC"),
    impute_na(d, "Area", "mle"),
    impute_na(d, "Area", "svd", 10)
  )
  lapply(methods, function(d_imputed) {
    expect_valid_lipidex(d_imputed, c(19, 11))
    mat_imputed <- assay(d_imputed, "Area")
    expect_false(any(is.na(mat_imputed)))
  })
})
