.check_log <- function(d, measure) {
  if (measure == "Retention.Time") {
    warning("Retention time should not be logged")
  }
  if (is_logged(d, measure)) {
    return(measure)
  }
  paste0("log2(", measure, ")")
}
.copy_attr <- function(d, original) {
  attr(d, "skyline") <- attr(original, "skyline")
  d
}


#' @importFrom rlang syms
.uniform_attrs <- function(datalist) {
  attrlist <- lapply(datalist, attr, "skyline")
  cols <- lapply(datalist, colnames)
  all_cols <- Reduce(union, cols)
  common_cols <- Reduce(intersect, cols)
  omitted_cols <- all_cols[!all_cols %in% common_cols]
  if (length(omitted_cols) > 0) {
    warning(
      "Some columns were not available in all files. ",
      paste(omitted_cols, collapse = ", ")
    )
  }

  ret <- bind_rows(datalist, .id = "filename") %>%
    select(!!!syms(c("filename", common_cols)))

  measure_cols <- Reduce(union, lapply(attrlist, "[[", "measures"))
  intensity_cols <- Reduce(union, lapply(attrlist, "[[", "intensity_cols"))
  attr(ret, "skyline") <- list(
    skyline = TRUE,
    measures = measure_cols[measure_cols %in% common_cols],
    intensity_cols = intensity_cols[intensity_cols %in% common_cols]
  )
  ret
}

.is_blank <- function(data, measure = "Area") {
  tic <- colSums(assay(data, measure), na.rm = TRUE)

  (median(tic) / tic) > 50
}

# colname created in .uniform.attrs
utils::globalVariables(c("filename"))
