#' Differential analysis of lipids between sample groups
#'
#' @param ... expressions, or character strings which can be parsed to expressions, specifying contrasts. 
#' These are passed to \code{limma::makeContrasts}.
#' @param data Skyline data.frame created by \code{\link{read.skyline}}, should be normalized and log2 transformed
#' @param measure name of the column indicating sample names
#' @param group_col name of the column indicating sample groups
#'
#' @importFrom forcats fct_drop
#' @importFrom dplyr one_of
#' @importFrom rlang quos
#' @return TopTable as returned by limma package
#' @export
#'
#' @examples
de_analysis = function(..., data, measure="Area", group_col=NULL){
  if(is.null(group_col)) {
    if(ncol(colData(data)) > 0) {
      group_col = names(colData(data))[[1]]
    } else {
      stop('Please add clinical data or specify a group column')
    }
  }
  
  symbols = as.character(.quos_syms(quos(...)))
  
  group = colData(data)[[group_col]]
  if (!all(symbols %in% as.character(group))) {
    stop("Some of the constrasts variables are not present in group_col")
  }
  data = data[, group %in% symbols]
  
  group = fct_drop(colData(data)[[group_col]])
  design = model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  return(de_design(..., design=design, data=data, measure=measure))
}

#' Differential analysis of lipids between sample groups
#' 
#' @param ... expressions, or character strings which can be parsed to expressions, specifying contrasts.
#' These are passed to \code{limma::makeContrasts}. Ignored if \code{coef} is provided.
#' @param coef column number or column name specifying which coefficient of the linear model is of interest.
#' @param data Skyline data.frame created by \code{\link{read_skyline}}, should be normalized and log2 transformed
#' @param measure name of the column indicating sample names
#' 
#' @importFrom rlang is_formula
#' @export
de_design <- function(..., coef=NULL, design, data, measure="Area") {
  if (is_formula(design)) {
    design = model.matrix(design, data=colData(data))
  } else if (! is.matrix(design)) {
    stop("design should be a matrix or formula")
  }
  vfit <- limma::lmFit(assay(data, measure), design)
  
  if (is.null(coef)) {
    contr.matrix = limma::makeContrasts(..., levels=colnames(design))
    vfit <- limma::contrasts.fit(vfit, contrasts=contr.matrix)
    coef = setNames(c(1:ncol(contr.matrix)), colnames(contr.matrix))
  } else {
    if (! coef %in% colnames(design)) {
      stop("One or more coefficients is not in the design matrix. Allowed values are ", colnames(design))
    }
    names(coef) = coef
  }
  
  efit <- limma::eBayes(vfit)
  dimname_x = data@attrs$dimnames[[1]]
  
  top = lapply(coef, function(x) 
    limma::topTable(efit, number = Inf, coef = x) %>% tibble::rownames_to_column(dimname_x)) %>%
    bind_rows(.id="contrast")
  
  top = to_df(data, dim="row") %>% select(one_of("Molecule", "Class", "total_cl", "total_cs", "itsd", dimname_x)) %>% 
    .left_join.silent(top)
  return(top)
}

#' Get a list of significantly changed molecules
#'
#' @param de.results output of \code{\link{de.analysis}}
#' @param p.cutoff 
#' @param logFC.cutoff 
#'
#' @return
#' @importFrom dplyr %>% filter
#' @export
#'
#' @examples
significant_molecules = function(de.results, p.cutoff=0.05, logFC.cutoff=1) {
  de.results %>% filter(adj.P.Val < p.cutoff, abs(logFC) > logFC.cutoff) %>% 
    (function(x) split(x$Molecule, x$contrast))
}


#' Lipid set enrichment analysis
#'
#' @param de.results output of \code{\link{de.analysis}}
#' @param rank.by statistic used to rank the lipid list
#'
#' @return
#' @importFrom dplyr %>% bind_rows arrange rename
#' @export
#'
#' @examples
enrich_lipidsets <- function(de.results, rank.by=c("logFC", "P.Value", "Adj.P.Val")) {
  rank.by = match.arg(rank.by)
  rank.by.sym = rlang::sym(rank.by)
  
  sets = gen_lipidsets(de.results)
  
  de.results = de.results %>% 
    arrange(-!!rank.by.sym) %>% 
    group_by(contrast, Molecule) %>%
    summarise(!!rank.by.sym := first(!!rank.by.sym))
  
  contrast.list = split.data.frame(de.results, de.results$contrast)
  
  res = lapply(contrast.list, function(x) 
    fgsea::fgsea(
      pathways = sets, 
      stats = to_named_list(x, "Molecule", rank.by),
      minSize = 2, nperm = 10000)
    ) %>%
    bind_rows(.id="contrast")
  
  res = res %>% arrange(padj) %>% rename(set=pathway)
  attr(res, "de.results") = de.results
  attr(res, "rank.by") = rank.by
  attr(res, "sets") = sets
  
  return(res)
}

#' Get a list of significantly changed lipid sets
#'
#' @param enrich.results output of \code{\link{enrich.lipid.set}}
#' @param p.cutoff 
#' @param size.cutoff 
#'
#' @return
#' @importFrom dplyr %>% filter
#' @export
#'
#' @examples
significant_lipidsets = function(enrich.results, p.cutoff=0.05, size.cutoff=2) {
  enrich.results %>% filter(padj < p.cutoff, size > size.cutoff) %>% 
    (function(x) split(x$set, x$contrast))
}


#' Generate lipid sets
#'
#' @param data data.frame with a column named Molecule
#'
#' @return list of lipid sets
#' @export
#' @importFrom dplyr %>% filter distinct
#' @importFrom tidyr gather unite
#' @examples
gen_lipidsets <- function(data) {
  if (!all( c("itsd", "Class", "total_cl", "total_cs") %in% colnames(data))) {
    #data = annotate_lipids(data)
  }
  
  data_ = annotate_lipids(data) %>%
    filter(!itsd) %>%
    distinct(Molecule, Class, total_cl, total_cs) %>%
    gather("collection", "value", Class, total_cl, total_cs) %>%
    unite("set", collection, value, sep = "_")
  
  split(as.character(data_$Molecule), data_$set)
}
