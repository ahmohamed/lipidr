#' Perform Probabilistic Quotient Normalization for intensities.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to use as intensity, usually Area, Area.Normalized or Height
#' @param exclude Samples to exclude, can be either: \cr 
#' "blank" - automatically detected blank samples and exclude them
#' logical vector with the same length as samples
#' 
#' @param log whether the normalized values should be log2 transformed
#'
#' @return
#' @importFrom SummarizedExperiment assay<- assays<-
#' @importFrom dplyr %>% vars select group_by mutate
#' @importFrom rlang sym UQ
#' @export
#' @references Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006). 
#' Probabilistic quotient normalization as robust method to account for dilution 
#' of complex biological mixtures. Application in 1H NMR metabonomics. 
#' Analytical chemistry, 78(13), 4281-4290.
#' 
#' @examples
normalize_pqn <- function(data, measure="Area", exclude="blank", log=TRUE) {
  if(mcols(assays(data), use.names = T)[measure, "normalized"]) {
    stop(measure, " is already normalized")
  }
  if(!is.null(exclude)) {
    if(exclude == "blank") {
      data = data[, !.is_blank(data)]
    } else {
      data = data[, exclude]
    }  
  }
  m = assay(data, measure)
  
  # factor_n = median ( lipid_i_n/ avg(lipid_i) )
  assay(data, measure) = m / apply(m / rowMeans(m, na.rm=TRUE), 2, median, na.rm=TRUE)
  mcols(assays(data), use.names = T)[measure, "normalized"] = TRUE
  
  if (log && !mcols(assays(data), use.names = T)[measure, "logged"]) {
    assay(data, measure) = log2(assay(data, measure))
    mcols(assays(data), use.names = T)[measure, "logged"] = TRUE
  }
  
  return(data)
}

#' Normalize each class by its corresponding internal standard(s).
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @param measure which measure to use as intensity, usually Area, Area.Normalized or Height
#' @param exclude Samples to exclude, can be either: \cr 
#' "blank" - automatically detected blank samples and exclude them
#' logical vector with the same length as samples
#' @param log whether the normalized values should be log2 transformed
#'
#' @return
#' 
#' @importFrom dplyr %>% select group_by mutate filter ungroup left_join inner_join
#' @importFrom rlang sym UQ
#' @export
#'
#' @examples
normalize_itsd <- function(data, measure="Area", exclude="blank", log=TRUE) {
  if(mcols(assays(data), use.names = T)[measure, "normalized"]) {
    stop(measure, " is already normalized")
  }
  if(!data@attrs$summarized){
    stop("Data should be summarized using summarize_transitions")
  }
  
  if(!is.null(exclude)) {
    if(exclude == "blank") {
      data = data[, !.is_blank(data)]
    } else {
      data = data[, exclude]
    }  
  }
  itsd = rowData(data)$itsd
  if(sum(itsd) == 0) {
    stop("No internal standards found in your lipid list.")
  }
  m = assay(data, measure)
  mitsd = m[itsd, ]
  
  # itsd_n = itsd_ni / mean(itsd_i)
  mitsd = mitsd / rowMeans(mitsd, na.rm=TRUE)
  
  # per class: 
  itsd_list = to_df(data) %>% group_by(filename, Class) %>% 
    mutate(itsd_list=list(as.character(MoleculeId[itsd]))) %>% .$itsd_list
  
  assay(data, measure) = laply(seq_along(itsd_list), function(i) {
    if(length(itsd_list[[i]]) == 0) {
      f = 1
    } else if (length(itsd_list[[i]]) == 1){
      f = mitsd[itsd_list[[i]], ]
    } else {
      f = colMeans( mitsd[itsd_list[[i]], ], na.rm=TRUE) 
    }
    
    return (m[i, ] / f)
  })
  
  if (log && !mcols(assays(data), use.names = T)[measure, "logged"]) {
    assay(data, measure) = log2(assay(data, measure))
    mcols(assays(data), use.names = T)[measure, "logged"] = TRUE
  }
  
  return(data)
}