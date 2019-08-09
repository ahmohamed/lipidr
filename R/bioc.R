# Defined as dimnames / dplyr .
utils::globalVariables(c(".", "TransitionId", "Sample"))

#' SkylineExperiment object
#'
#' @export
#' @import methods
.SkylineExperiment <- setClass(
  "SkylineExperiment",
  contains = "SummarizedExperiment"
)

setValidity("SkylineExperiment", function(object) {
  errors <- character()
  metadata <- metadata(object)

  dn <- metadata$dimnames
  if (!is.character(dn) || length(dn) != 2) {
    print(!is.character(dn))
    print(length(dn != 2))
    msg <- "metadata should have a 'dimnames' character vector of length 2"
    errors <- c(errors, msg)
  }

  summarized <- metadata$summarized
  if (!is.logical(summarized) || length(summarized) != 1) {
    msg <- "metadata should have a 'summarized' logical value"
    errors <- c(errors, msg)
  }

  row_data <- rowData(object)
  annotations <- c("filename", "Molecule", "Class", "istd")
  if (!all(annotations %in% colnames(row_data))) {
    annotations <- paste(annotations, collapse = ", ")
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
#' @param metadata A list containing arbitrary information about the experiment.
#'   It should at least contain 2 elements: \itemize{
#'     \item dimnames    2-element character vector with dimension names
#'     \item summarized   Has transitions been summarized?
#'   }
#'
#' @return SkylineExperiment object
#' @export
SkylineExperiment <- function(assay_list, metadata,
                              colData = NULL, rowData = NULL) {
  stopifnot(length(assay_list) > 0)
  if (is.null(colData)) {
    colData <- DataFrame(row.names = colnames(assay_list[[1]]))
  }
  se <- SummarizedExperiment(assay_list,
    colData = colData, rowData = rowData, metadata = metadata)
  ret <- .SkylineExperiment(se)
  return(ret)
}

#' Convert data.frame/matrix to SkylineExperiment
#'
#' @param df A data.frame or matrix where rows are lipids and columns
#'   are samples. Lipid names should be provided in the first column
#'   or in rownames of `df`. Measurements should be numeric.
#'   The data is considered `summarized` unless at least one lipid
#'   is duplicated (has > 1 row).
#' @param logged Whether the data is log-transformed
#' @param normalized Whether the data is normalized
#'
#' @return SkylineExperiment
#' @export
as_skyline_experiment <- function(df, logged=FALSE, normalized=FALSE) {
  # if (!.is_skyline_export(df)) {
  #   #return(.to_summarized_experiment(d))
  # }
  first_col <- 0
  if (is.factor(df[, 1])) {
    df[, 1] <- as.character(df[, 1])
  }
  if (is.character(df[, 1])) {
    first_col <- sum(!annotate_lipids(df[, 1], no_match = "ignore")$not_matched)
  }
  rows <- 0
  if (is.character(rownames(df))) {
    rows <- sum(!annotate_lipids(rownames(df), no_match = "ignore")$not_matched)
  }
  if (first_col == 0 && rows == 0) {
    stop("Data frame does not contain valid lipid names. ",
      "Lipids features should be in rownames or the first column."
    )
  }

  if (first_col > rows) {
    # Molecules are in the first column
    molecules <- df[, 1]
    if (sum(duplicated(molecules))) {
      row_dimname <- "TransitionId"
      summarized <- FALSE
      df <- df[, -1]
    } else {
      row_dimname <- "MoleculeId"
      summarized <- TRUE
      rownames(df) <- df[, 1]
      df <- df[, -1]
    }
  } else {
    # Molecules are in the rownames
    molecules <- rownames(df)
    row_dimname <- "MoleculeId"
    summarized <- TRUE
  }
  df <- df %>% mutate_if(is.factor, as.character)
  assay_list <- list(Area = data.matrix(df, rownames.force = TRUE))
  assay_list <- as(assay_list, "SimpleList")
  mcols(assay_list) <- list(logged = logged, normalized = normalized)
  row_data <- DataFrame(
    filename="dataframe",
    Molecule = molecules,
    row.names=rownames(df)
  ) %>%
    left_join(annotate_lipids(molecules))
  metadata <- list(
    summarized = summarized,
    dimnames = c(row_dimname, "Sample")
  )
  SkylineExperiment(
    assay_list,
    metadata = metadata,
    rowData = row_data
  )
}

#' @importFrom S4Vectors mcols mcols<- metadata
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
  row_data <- row_data %>% left_join(annotate_lipids(row_data$Molecule))
  metadata <- list(summarized = FALSE, dimnames = c("TransitionId", "Sample"))

  SkylineExperiment(
    assay_list,
    metadata = metadata,
    colData = col_data,
    rowData = row_data
  )
}

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

#' @export
left_join.DataFrame <- .join_wrapper(dplyr::left_join)
#' @export
full_join.DataFrame <- .join_wrapper(dplyr::full_join)
#' @export
inner_join.DataFrame <- .join_wrapper(dplyr::inner_join)
