#' Differential analysis of lipids between sample groups
#'
#' @param ... Expressions, or character strings which can be parsed to expressions, specifying contrasts.
#' These are passed to \code{limma::makeContrasts}.
#' @param data Skyline data.frame created by \code{\link{read_skyline}}, should be normalized and log2 transformed.
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
#' de_results <- de_analysis(HighFat_water - NormalDiet_water, data = data_normalized, measure = "Area")
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
#' @param ... Expressions, or character strings which can be parsed to expressions, specifying contrasts.
#' These are passed to \code{limma::makeContrasts}. Ignored if \code{coef} is provided.
#' @param design Design matrix generated from \code{\link{model.matrix}}, or a design formula .
#' @param coef Column number or column name specifying which coefficient of the linear model is of interest.
#' @param data Skyline data.frame created by \code{\link{read_skyline}}, should be normalized and log2 transformed.
#' @param measure Name of the column containing sample names.
#'
#' @return TopTable as returned by limma package
#' @importFrom rlang is_formula
#' @export
#' @examples
#' # type ?normalize_pqn to see how to normalize and log-transfome your data
#' data(data_normalized)
#' 
#' # Using formula
#' de_results <- de_design(coef = "groupHighFat_water", design = ~group, data = data_normalized, measure = "Area")
#' 
#' # Using design matrix
#' design <- model.matrix(~group, data = colData(data_normalized))
#' de_results <- de_design(coef = "groupHighFat_water", design = design, data = data_normalized, measure = "Area")
de_design <- function(..., coef = NULL, design, data, measure = "Area") {
  if (is_formula(design)) {
    design <- model.matrix(design, data = colData(data))
  } else if (!is.matrix(design)) {
    stop("design should be a matrix or formula")
  }
  vfit <- limma::lmFit(assay(data, measure), design)

  if (is.null(coef)) {
    contr.matrix <- limma::makeContrasts(..., levels = colnames(design))
    vfit <- limma::contrasts.fit(vfit, contrasts = contr.matrix)
    coef <- setNames(seq_len(ncol(contr.matrix)), colnames(contr.matrix))
  } else {
    if (!coef %in% colnames(design)) {
      stop("One or more coefficients is not in the design matrix. Allowed values are ", colnames(design))
    }
    names(coef) <- coef
  }

  efit <- limma::eBayes(vfit)
  dimname_x <- data@attrs$dimnames[[1]]

  top <- lapply(coef, function(x)
    limma::topTable(efit, number = Inf, coef = x) %>% rownames_to_column(dimname_x)) %>%
    bind_rows(.id = "contrast")

  top <- to_df(data, dim = "row") %>%
    select(one_of("Molecule", "Class", "total_cl", "total_cs", "itsd", dimname_x)) %>%
    .left_join.silent(top)
  return(top)
}

#' Get a list of significantly changed molecules
#'
#' @param de.results Output of \code{\link{de_analysis}}.
#' @param p.cutoff Significance threshold.
#' @param logFC.cutoff Cutoff limit for log2 fold change.
#'
#' @return A character vector with names of significantly differentially changed lipids.
#' @importFrom dplyr %>% filter
#' @export
#'
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(HighFat_water - NormalDiet_water, data = data_normalized, measure = "Area")
#' significant_molecules(de_results)
significant_molecules <- function(de.results, p.cutoff = 0.05, logFC.cutoff = 1) {
  de.results %>%
    filter(adj.P.Val < p.cutoff, abs(logFC) > logFC.cutoff) %>%
    (function(x) split(x$Molecule, x$contrast))
}


#' Lipid set enrichment analysis
#'
#' @param de.results Output of \code{\link{de_analysis}}.
#' @param rank.by Statistic used to rank the lipid list.
#'
#' @return a data.frame with enrichment results as obtained from \code{\link[fgsea]{fgsea}}.
#'   The results also contain the following attributes:#' \itemize{
#'     \item de.results    Original de.results input.
#'     \item rank.by   Measure used to rank lipid molecules.
#'     \item sets   Lipid sets tested, with their member molecules.
#'   }
#'
#' @importFrom dplyr %>% bind_rows arrange rename
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(HighFat_water - NormalDiet_water, data = data_normalized, measure = "Area")
#' enrich_results <- enrich_lipidsets(de_results, rank.by = "logFC")
enrich_lipidsets <- function(de.results, rank.by = c("logFC", "P.Value", "Adj.P.Val")) {
  rank.by <- match.arg(rank.by)
  rank.by.sym <- rlang::sym(rank.by)

  sets <- gen_lipidsets(de.results$Molecule)

  de.results <- de.results %>%
    arrange(-!!rank.by.sym) %>%
    group_by(contrast, Molecule) %>%
    summarise(!!rank.by.sym := first(!!rank.by.sym))

  contrast.list <- split.data.frame(de.results, de.results$contrast)

  res <- lapply(contrast.list, function(x)
    fgsea::fgsea(
      pathways = sets,
      stats = to_named_list(x, "Molecule", rank.by),
      minSize = 2, nperm = 10000
    )) %>%
    bind_rows(.id = "contrast")

  res <- res %>% arrange(padj) %>% rename(set = pathway)
  attr(res, "de.results") <- de.results
  attr(res, "rank.by") <- rank.by
  attr(res, "sets") <- sets

  return(res)
}

#' Get a list of significantly changed lipid sets
#'
#' @param enrich.results Output of \code{\link{enrich_lipidsets}}.
#' @param p.cutoff Significance threshold.
#' @param size.cutoff Minimum number of lipids in a set tested for enrichment.
#'
#' @return A list of character vectors of significantly enriched sets for each contrast.
#' @importFrom dplyr %>% filter
#' @export
#'
#' @examples
#' data(data_normalized)
#' de_results <- de_analysis(HighFat_water - NormalDiet_water, data = data_normalized, measure = "Area")
#' enrich_results <- enrich_lipidsets(de_results, rank.by = "logFC")
#' significant_lipidsets(enrich_results)
significant_lipidsets <- function(enrich.results, p.cutoff = 0.05, size.cutoff = 2) {
  enrich.results %>%
    filter(padj < p.cutoff, size > size.cutoff) %>%
    (function(x) split(x$set, x$contrast))
}


#' Generate lipid sets from lipid molecule names
#'
#' @param molecules A character vector containing lipid molecule names.
#'
#' @return List of lipid sets
#' @export
#' @importFrom dplyr %>% filter distinct
#' @importFrom tidyr gather unite
#' @examples
#' data(data_normalized)
#' molecules <- rowData(data_normalized)$Molecule
#' gen_lipidsets(molecules)
gen_lipidsets <- function(molecules) {
  if (!all(c("itsd", "Class", "total_cl", "total_cs") %in% colnames(data))) {
    # data = annotate_lipids(data)
  }

  data_ <- annotate_lipids(molecules) %>%
    filter(!itsd) %>%
    distinct(Molecule, Class, total_cl, total_cs) %>%
    gather("collection", "value", Class, total_cl, total_cs) %>%
    unite("set", collection, value, sep = "_")

  split(as.character(data_$Molecule), data_$set)
}
