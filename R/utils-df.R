to_named_list <- function(df, name, value) {
  stopifnot(is.data.frame(df))
  ret <- data.frame(df)[, value]
  names(ret) <- as.character(data.frame(df)[, name])
  ret
}

to_rownames_ <- function(df, col = "rowname") {
  stopifnot(is.data.frame(df))
  nm <- data.frame(df)[, col]
  df <- as.matrix(df[, !(colnames(df) %in% col)])
  rownames(df) <- as.character(nm)
  df
}

#' @importFrom tidyr spread
#' @importFrom rlang sym
to_num_matrix <- function(data, sample, feature, measure) {
  sample <- sym(sample)
  feature <- sym(feature)
  measure <- sym(measure)

  data %>%
    select(!!sample, !!feature, !!measure) %>%
    spread(!!sample, !!measure) %>%
    to_rownames_(as.character(feature))
}

.replace_na_rowmean <- function(m) {
  k <- which(is.na(m), arr.ind = TRUE)
  if (length(k) > 0) {
    m[k] <- rowMeans(m, na.rm = TRUE)[k[, 1]]
  }
  return(m)
}

rownames_to_column <- function(df, var = "rowname") {
  stopifnot(is.data.frame(df))
  df <- cbind(data.frame(rownames(df), stringsAsFactors = FALSE), df)
  colnames(df)[[1]] <- var
  rownames(df) <- NULL
  df
}


laply <- function(l, fun) {
  ret <- lapply(l, fun)
  ret.mat <- ret %>% unlist() %>% matrix(nrow = length(l), byrow = TRUE)
  colnames(ret.mat) <- names(ret[[1]])
  ret.mat
}

.silent <- function(f) {
  return(function(...) suppressWarnings(suppressMessages(f(...))))
}
.left_join_silent <- .silent(dplyr::left_join)
.full_join_silent <- .silent(dplyr::full_join)
.inner_join_silent <- .silent(dplyr::inner_join)
