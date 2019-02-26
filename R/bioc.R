library(rlang)
library(SummarizedExperiment)

#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SkylineExperiment <- setClass("SkylineExperiment", slots=list(attrs="list"),contains = "SummarizedExperiment")

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
SkylineExperiment <- function(assay_list, attrs, colData=NULL, rowData=NULL, ...) {
  se = SummarizedExperiment(assay_list, colData=colData, rowData=rowData)
  ret = .SkylineExperiment(se)
  ret@attrs = attrs
  return(ret)
}

#' @importFrom S4Vectors mcols mcols<-
.to_summarized_experiment <- function(d) {
  if(is.null(attr(d, "skyline"))) {
    stop("Data.frame does not have skyline attribute")
  }
  assay_list = lapply(attr(d, "skyline")$measures, function(m) to_num_matrix(d, "Sample", "TransitionId", m))
  names(assay_list) = attr(d, "skyline")$measures
  assay_list = as(assay_list, "SimpleList")
  mcols(assay_list) = list(logged=FALSE, normalized=FALSE)
  
  col_data = d %>% distinct(!!!(syms(c("Sample", attr(d, "skyline")$annot_cols))))
  col_data = toDataFrame(col_data, row.names = "Sample")
  
  row_data = d %>% 
    select(TransitionId, matches("filename|^Class|^Molecule|^Precursor|^Product")) %>%
    distinct()
  row_data = toDataFrame(row_data, row.names = "TransitionId")
  row_data = row_data[ row.names(assay_list[[1]]), ]
  attrs = list(summarized = FALSE, dimnames=c("TransitionId", "Sample"))
  SkylineExperiment(assay_list, colData=col_data, rowData=row_data, attrs=attrs)
}

de_join <- function(df, variable) {
  .summarise_fun = function(x) {
    u = unique(x)
    if(length(u) == 1) {
      return(u)
    }
    return(list(u))
  }
  df %>% group_by(!!sym(variable)) %>% 
    summarise_all(function(x) list(unique(x))) %>% 
    select_if(function(x) all(sapply(x, length) == 1)) %>%
    group_by(!!sym(variable)) %>% summarise_all(unlist)
}
group_by.SummarizedExperiment <- function(.data, ..., add=FALSE, rowgroups=FALSE) {
  if(rowgroups) {
    metadata(.data)$rowgroupvars <- group_by_prepare(.data, ..., add=add)$group_names
  }else {
    metadata(.data)$colgroupvars <- group_by_prepare(.data, ..., add=add)$group_names
  }
  return(.data)
}

#' @export
filter.SummarizedExperiment <- function(.data, ...) {
  dots <- rlang::quos(...)
  if (any(rlang::have_name(dots))) {
    bad <- dots[rlang::have_name(dots)]
    bad_eq_ops(bad, "must not be named, do you need `==`?")
  } else if (rlang:::is_empty(dots)) {
    return(.data)
  }
  
  quo <- dplyr:::all_exprs(!!!dots, .vectorised = TRUE)
  d = DataFrame(colData(.data), rowid=c(1:ncol(.data)))
  if(!is.null(metadata(.data)$colgroupvars)) {
    d = group_by_impl(d, metadata(.data)$colgroupvars)
  }
  out <- dplyr:::filter_impl(d, quo)
  .data[, out$rowid]
}

#' @export
select.SummarizedExperiment <- function(.data, ..., .dots = list()) {
  .names = c(as.list(colData(.data)), as.list(rowData(.data)), as.list(assays(.data)))
  select2(.names, ..., .dots)
}


#' @export
select.DataFrame <- function(.data, ..., .dots = list()) {
  vars <- dplyr:::select_vars(names(.data), !!!quos(...))
  dplyr:::select_impl(.data, vars)
}

select2 <- function(.data, ..., .dots = list()) {
  vars <- dplyr:::select_vars(names(.data), !!!quos(...))
  dplyr:::select_impl(.data, vars)
}

summarise.SummarizedExperiment <- function(.data, .assays=vars(), .cols=vars(), .rows=vars()) {
  col_data = DataFrame(colData(.data), rowid=c(1:ncol(.data)))
  if(!is_empty(metadata(.data)$colgroupvars)) {
    col_data = group_by_impl(col_data, metadata(.data)$colgroupvars)
    colgroups = attr(col_data, "indices")
    col_data_sum = dplyr:::summarise_impl(col_data, .cols)
  }else {
    colgroups = TRUE
  }
  
  row_data = DataFrame(rowData(.data), rowid=c(1:nrow(.data)))
  if(!is_empty(metadata(.data)$rowgroupvars)) {
    row_data = group_by_impl(row_data, metadata(.data)$rowgroupvars)
    rowgroups = attr(row_data, "indices")
    row_data_sum = dplyr:::summarise_impl(row_data, .rows)
  } else {
    rowgroups = TRUE
  }
  
  
  
}
group_by_impl <- function(.data, symbols) {
  fun = dplyr:::grouped_df_impl
  if ( length(formals(fun)) == 3) {
    return (fun(.data, symbols, TRUE))
  }
  return (fun(.data, symbols))
}
get_data_mask <- function(ds) {
  assay_env = env(assays(ds) %>% as.list())
  row_env = child_env(assay_env, rowData(ds) %>% as.data.frame())
  col_env = child_env(row_env, colData(ds) %>% as.data.frame())
  
  return(col_env)
}
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr gather
to_long_format <- function(ds, measure="Area") {
  dims = ds@attrs$dimnames
  assay(ds, measure) %>% as.data.frame() %>% 
    tibble::rownames_to_column(dims[[1]]) %>% 
    gather(key=!!dims[[2]], value=!!measure, -!!dims[[1]]) %>%
    left_join(to_df(ds, "row")) %>%
    left_join(to_df(ds, "col"))
}

#' @importFrom S4Vectors DataFrame
toDataFrame <- function(df, row.names.col="rowname") {
  nm <- data.frame(df)[, row.names.col]
  df = df[, !(colnames(df) %in% row.names.col)]
  DataFrame(df, row.names = nm)
}

#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
to_df <- function(d, dim="row") {
  if (dim == "row") {
    row_data = rowData(d)
    rownames(row_data) = rownames(d)
    row_data %>% as.data.frame() %>%
      tibble::rownames_to_column(d@attrs$dimnames[[1]])
  } else {
    colData(d) %>% as.data.frame() %>%
      tibble::rownames_to_column(d@attrs$dimnames[[2]])
  }
}
.join_wrapper <- function(f) {
  return (function(r, l, by=NULL) {
    r %>% as.data.frame %>%
      tibble::rownames_to_column() %>%
      f(l, by) %>%
      toDataFrame()
  })
}

#' @export
colData <- SummarizedExperiment::colData
#' @export
rowData <- SummarizedExperiment::rowData

#' @export
left_join.DataFrame <- .join_wrapper(dplyr::left_join)
#' @export
full_join.DataFrame <- .join_wrapper(dplyr::full_join)
#' @export
inner_join.DataFrame <- .join_wrapper(dplyr::inner_join)

# df %>% group_by(letter) %>% summarise_all(summarise_fun) %>% select_if(function(x) !is.list(x))
# 
# df %>% group_by(letter) %>% summarise_all(function(x) length(unique(x)) == 1) %>% select_if()
# 
# ddply(df, .(letter), )
