#' Differential analysis of lipids between sample groups
#'
#' @param ... Expressions, or character strings which can be parsed to
#' expressions, specifying contrasts. These are passed to
#' `limma::makeContrasts`.
#'
#' @param data Skyline data.frame created by [read_skyline()],
#'   should be normalized and log2 transformed.
#' @param measure Name of the column containing sample names.
#' @param group_col Name of the column containing sample groups.
#'
#' @importFrom forcats fct_drop
#' @importFrom dplyr one_of
#' @importFrom rlang quos
#' @importFrom stats model.matrix setNames
#' @return TopTable as returned by limma package
#' @export
#'
#' @examples
#' # type ?normalize_pqn to see how to normalize and log2-transform your data
#' data(data_normalized)
#' de_results <- de_analysis(
#'   HighFat_water - NormalDiet_water,
#'   data = data_normalized, measure = "Area"
#' )
de_analysis <- function(..., data, measure = "Area", group_col = NULL) {
  if (is.null(group_col)) {
    if (ncol(colData(data)) > 0) {
      group_col <- names(colData(data))[[1]]
    } else {
      stop("Please add clinical data or specify a group column")
    }
  }

  symbols <- as.character(.quos_syms(quos(...)))

  group <- colData(data)[[group_col]]
  if (!all(symbols %in% as.character(group))) {
    stop("Some of the constrasts variables are not present in group_col")
  }
  data <- data[, group %in% symbols]

  group <- fct_drop(colData(data)[[group_col]])
  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  return(de_design(..., design = design, data = data, measure = measure))
}

#' Differential analysis of lipids between sample groups
#'
#' @param ... Expressions, or character strings which can be parsed to
#'   expressions, specifying contrasts. These are passed to
#'   `limma::makeContrasts`. Ignored if `coef` is provided.
#' @param design Design matrix generated from [model.matrix()],
#'   or a design formula.
#' @param coef Column number or column name specifying which coefficient of
#'   the linear model is of interest.
#' @param data Skyline data.frame created by [read_skyline()],
#'   should be normalized and log2 transformed.
#' @param measure Name of the column containing sample names.
#'
#' @return TopTable as returned by limma package
#' @importFrom rlang is_formula
#' @importFrom limma topTable lmFit makeContrasts contrasts.fit eBayes
#' @export
#' @examples
#' # type ?normalize_pqn to see how to normalize and log-transfome your data
#' data(data_normalized)
#'
#' # Using formula
#' de_results <- de_design(
#'   coef = "groupHighFat_water",
#'   design = ~group,
#'   data = data_normalized,
#'   measure = "Area"
#' )
#'
#' # Using design matrix
#' design <- model.matrix(~group, data = colData(data_normalized))
#' de_results <- de_design(
#'   coef = "groupHighFat_water",
#'   design = design,
#'   data = data_normalized,
#'   measure = "Area"
#' )
de_design <- function(..., coef = NULL, design, data, measure = "Area") {
  if (is_formula(design)) {
    design <- model.matrix(design, data = colData(data))
  } else if (!is.matrix(design)) {
    stop("design should be a matrix or formula")
  }
  vfit <- lmFit(assay(data, measure), design)

  if (is.null(coef)) {
    contr.matrix <- limma::makeContrasts(..., levels = colnames(design))
    vfit <- limma::contrasts.fit(vfit, contrasts = contr.matrix)
    coef <- setNames(seq_len(ncol(contr.matrix)), colnames(contr.matrix))
  } else {
    if (!coef %in% colnames(design)) {
      stop(
        "One or more coefficients is not in the design matrix.",
        " Allowed values are ", colnames(design)
      )
    }
    names(coef) <- coef
  }

  efit <- eBayes(vfit)
  dimname_x <- data@attrs$dimnames[[1]]

  top <- lapply(
    coef, function(x)
      topTable(efit, number = Inf, coef = x) %>% rownames_to_column(dimname_x)
  ) %>%
    bind_rows(.id = "contrast")

  top <- to_df(data, dim = "row") %>%
    select(
      one_of("Molecule", "Class", "total_cl", "total_cs", "itsd", dimname_x)
    ) %>%
    .left_join_silent(top)
  return(top)
}

#' Get a list of significantly changed molecules
#'
#' @param de.results Output of [de_analysis()].
#' @param p.cutoff Significance threshold.
#' @param logFC.cutoff Cutoff limit for log2 fold change.
#'
#' @return A character vector with names of significantly differentially
#'   changed lipids.
#' @importFrom dplyr %>% filter
#' @export
#'
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(
#'   HighFat_water - NormalDiet_water,
#'   data = data_normalized, measure = "Area"
#' )
#' significant_molecules(de_results)
significant_molecules <- function(de.results, p.cutoff = 0.05,
                                  logFC.cutoff = 1) {
  de.results %>%
    filter(adj.P.Val < p.cutoff, abs(logFC) > logFC.cutoff) %>%
    (function(x) split(x$Molecule, x$contrast))
}

