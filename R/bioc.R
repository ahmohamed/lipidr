#' SkylineExperiment object
#'
#' @slot attrs Extra slot to hold the workflow stage for the data, whether normalized, summarized or log transformed.
#'
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SkylineExperiment <- setClass(
  "SkylineExperiment",
  slots = list(attrs = "list"),
  contains = "SummarizedExperiment"
)

#' Constructor for Skyline experiment from list of assays
#'
#' @param assay_list A list or SimpleList of matrix-like elements, or a matrix-like object. Passed to [SummarizedExperiment()].
#' @param attrs A list of extra attributes to be saved to SkylineExperiment object.
#' @param colData A DataFrame object describing the rows (contains generated lipid annotations). Row names, if present, become the row names of the SummarizedExperiment object. The number of rows of the DataFrame must be equal to the number of rows of the matrices in assays.
#' @param rowData An optional data.frame describing the samples (contains clinical information). Row names, if present, become the column names of the SkylineExperiment.
#'
#' @return SkylineExperiment object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SkylineExperiment <- function(assay_list, attrs,
                              colData = NULL, rowData = NULL) {
  se <- SummarizedExperiment(assay_list, colData = colData, rowData = rowData)
  ret <- .SkylineExperiment(se)
  ret@attrs <- attrs
  return(ret)
}

#' @importFrom S4Vectors mcols mcols<-
.to_summarized_experiment <- function(d) {
  if (is.null(attr(d, "skyline"))) {
    stop("Data.frame does not have skyline attribute")
  }
  assay_list <- lapply(
    attr(d, "skyline")$measures,
    function(m) to_num_matrix(d, "Sample", "TransitionId", m)
  )
  names(assay_list) <- attr(d, "skyline")$measures
  assay_list <- as(assay_list, "SimpleList")
  mcols(assay_list) <- list(logged = FALSE, normalized = FALSE)

  col_data <- d %>%
    distinct(!!!(syms(c("Sample", attr(d, "skyline")$annot_cols))))
  col_data <- toDataFrame(col_data, row.names.col = "Sample")

  row_data <- d %>%
    select(
      TransitionId,
      matches("filename|^Class|^Molecule|^Precursor|^Product")
    ) %>%
    distinct()
  row_data <- toDataFrame(row_data, row.names.col = "TransitionId")
  row_data <- row_data[ row.names(assay_list[[1]]), ]
  attrs <- list(summarized = FALSE, dimnames = c("TransitionId", "Sample"))

  SkylineExperiment(
    assay_list,
    colData = col_data,
    rowData = row_data,
    attrs = attrs
  )
}

#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr gather
to_long_format <- function(ds, measure = "Area") {
  dims <- ds@attrs$dimnames
  assay(ds, measure) %>%
    as.data.frame() %>%
    rownames_to_column(dims[[1]]) %>%
    gather(key = !!dims[[2]], value = !!measure, -!!dims[[1]]) %>%
    left_join(to_df(ds, "row")) %>%
    left_join(to_df(ds, "col"))
}

#' @importFrom S4Vectors DataFrame
toDataFrame <- function(df, row.names.col = "rowname") {
  nm <- data.frame(df)[, row.names.col]
  df <- df[, !(colnames(df) %in% row.names.col)]
  DataFrame(df, row.names = nm)
}

#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
to_df <- function(d, dim = "row") {
  if (dim == "row") {
    row_data <- rowData(d)
    rownames(row_data) <- rownames(d)
    row_data %>%
      as.data.frame() %>%
      rownames_to_column(d@attrs$dimnames[[1]])
  } else {
    colData(d) %>%
      as.data.frame() %>%
      rownames_to_column(d@attrs$dimnames[[2]])
  }
}
.join_wrapper <- function(f) {
  return(function(r, l, by = NULL) {
    r %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      f(l, by) %>%
      toDataFrame()
  })
}

#' @export
setMethod(
  "colData",
  c(x = "SkylineExperiment"),
  selectMethod("colData", list("SummarizedExperiment"))
)
#' @export
setMethod(
  "rowData",
  c(x = "SkylineExperiment"),
  selectMethod("rowData", list("SummarizedExperiment"))
)

#' @export
left_join.DataFrame <- .join_wrapper(dplyr::left_join)
#' @export
full_join.DataFrame <- .join_wrapper(dplyr::full_join)
#' @export
inner_join.DataFrame <- .join_wrapper(dplyr::inner_join)
