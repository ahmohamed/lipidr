# Defined as dimnames / dplyr .
utils::globalVariables(c(".", "TransitionId", "Sample"))

#' SkylineExperiment object
#'
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SkylineExperiment <- setClass(
  "SkylineExperiment",
  contains = "SummarizedExperiment"
)

setValidity("SkylineExperiment", function(object) {
  errors <- character()
  metadata <- metadata(object)

  dn <- metadata$dimnames
  if (!is.character(dn) || length(dn != 2)) {
    msg <- "metadata should have a 'dimnames' character vector of length 2"
    errors <- c(errors, msg)
  }

  summarized <- metadata$summarized
  if (!is.is.logical(summarized) || length(summarized != 1)) {
    msg <- "metadata should have a 'summarized' logical value"
    errors <- c(errors, msg)
  }

  row_data <- rowData(object)
  annotations <- c("filename", "Molecule", "Class", "itsd")
  if (!all(annotations %in% row_data)) {
    msg <- paste("rowData should have these annotation columns:", annotations)
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
})

#' Constructor for Skyline experiment from list of assays
#'
#' @param assay_list A list or SimpleList of matrix-like elements,
#'   or a matrix-like object. Passed to [SummarizedExperiment()].
#' @param rowData A DataFrame object describing the rows (contains generated
#'   lipid annotations). Row names, if present, become the row names of the
#'   SummarizedExperiment object. The number of rows of the DataFrame
#'   must be equal to the number of rows of the matrices in assays.
#' @param colData An optional DataFrame describing the samples (contains
#'   clinical information). Row names, if present, become the column names of
#'   the SkylineExperiment.
#'
#' @return SkylineExperiment object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SkylineExperiment <- function(assay_list, metadata,
                              colData = NULL, rowData = NULL) {
  se <- SummarizedExperiment(assay_list,
    colData = colData, rowData = rowData, metadata = metadata)
  ret <- .SkylineExperiment(se)
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
  metadata <- list(summarized = FALSE, dimnames = c("TransitionId", "Sample"))

  SkylineExperiment(
    assay_list,
    metadata = metadata,
    colData = col_data,
    rowData = row_data
  )
}

#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr gather
to_long_format <- function(ds, measure = "Area") {
  dims <- metadata(ds)$dimnames
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
      rownames_to_column(metadata(d)$dimnames[[1]])
  } else {
    colData(d) %>%
      as.data.frame() %>%
      rownames_to_column(metadata(d)$dimnames[[2]])
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

#' @importFrom SummarizedExperiment colData rowData
#' @export colData
#' @export rowData
NULL

#' @export
left_join.DataFrame <- .join_wrapper(dplyr::left_join)
#' @export
full_join.DataFrame <- .join_wrapper(dplyr::full_join)
#' @export
inner_join.DataFrame <- .join_wrapper(dplyr::inner_join)
