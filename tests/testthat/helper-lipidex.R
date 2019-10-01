save_temp_csv <- function(d, sep=",",...) {
  file <- tempfile(fileext = ".csv")
  write.table(d, file = file, row.names = FALSE, sep=sep,...)
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

gen_sample_annot <- function(d) {
  samples <- colnames(d)
  groups <- c("A", "A", "A", "B", "B", "B", "Blank", "Blank", "QC", "QC", "QC")
  data.frame(Sample=samples, Group=groups, stringsAsFactors = FALSE)
}
