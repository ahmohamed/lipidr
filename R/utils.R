class_pattern <- "Protein\\.Name|Protein"
class_names <- c("Protein.Name", "Protein")

molecule_pattern <- "Peptide\\.Name|Peptide|Molecule\\.Name|Precursor.Ion.Name"
molecule_names <- c("Peptide.Name", "Peptide", "Molecule.Name", "Precursor.Ion.Name")

replicate_pattern <- "Replicate\\.Name|Replicate"
replicate_names <- c("Replicate.Name", "Replicate")

intensity_pattern <- "Area|Height|Area\\.Normalized"
intensity_names <- c("Area", "Height", "Area.Normalized")

measure_pattern <- "Area|Retention\\.Time|Height|Area\\.Normalized|Background"
measure_names <- c("Area", "Height", "Area.Normalized", "Retention.Time", "Background")

#' @importFrom SummarizedExperiment mcols assays
.check_log <- function(d, measure) {
  if (measure == "Retention.Time") {
    warning("Retention time should not be logged")
  }
  if (mcols(assays(d), use.names = TRUE)[measure, "logged"]) {
    return(measure)
  }
  return(paste0("log2(", measure, ")"))
}
.copy.attr <- function(d, original) {
  attr(d, "skyline") <- attr(original, "skyline")
  return(d)
}


#' @importFrom dplyr select bind_rows
#' @importFrom rlang syms
.uniform.attrs <- function(datalist) {
  attrlist <- lapply(datalist, attr, "skyline")
  cols <- lapply(datalist, colnames)
  all_cols <- Reduce(union, cols)
  common_cols <- Reduce(intersect, cols)
  omitted_cols <- all_cols[!all_cols %in% common_cols]
  if (length(omitted_cols) > 0) {
    warning("Some columns were not available in all files. ", omitted_cols)
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
  # print(common_cols)
  # print(names(datalist))
  return(ret)
}

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

#' @importFrom dplyr select
#' @importFrom tidyr spread
#' @importFrom rlang sym
to_num_matrix <- function(data, sample, feature, measure) {
  sample <- sym(sample)
  feature <- sym(feature)
  measure <- sym(measure)

  data %>%
    select(!!sample, !!feature, !!measure) %>%
    spread(!!sample, !!measure) %>%
    to_rownames_(as.character(feature)) %>%
    return()
}

.replace_na_rowmean <- function(m) {
  k <- which(is.na(m), arr.ind = TRUE)
  if (length(k) > 0) {
    m[k] <- rowMeans(m, na.rm = TRUE)[k[, 1]]
  }
  return(m)
}

.is_blank <- function(data, measure = "Area") {
  tic <- colSums(assay(data, measure), na.rm = TRUE)
  return(median(tic) / tic > 50)
}

rownames_to_column <- function(df, var = "rowname") {
  stopifnot(is.data.frame(df))
  df <- cbind(data.frame(rownames(df)), df)
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
.quos_syms <- function(x) {
  ret <- list()
  if (rlang::is_syntactic_literal(x)) {
    return()
  } else if (is.symbol(x)) {
    return(x)
  } else if (is.call(x)) {
    return(unlist(lapply(x[TRUE][-1], .quos_syms)))
  } else if (is.pairlist(x)) {
    return(unlist(lapply(x[TRUE], .quos_syms)))
  } else if (is.list(x)) {
    return(unlist(lapply(x, .quos_syms)))
  }
}

.silent <- function(f) {
  return(function(...) suppressWarnings(suppressMessages(f(...))))
}
.left_join_silent <- .silent(dplyr::left_join)
.full_join_silent <- .silent(dplyr::full_join)
.inner_join_silent <- .silent(dplyr::inner_join)

#' @importFrom stats sd
.cv <- function(a) {
  (sd(a, na.rm = TRUE) / mean(a, na.rm = TRUE)) * 100
}

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}
