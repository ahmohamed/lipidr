#' Read Skyline exported files
#'
#' @param files list of Skyline exported files in CSV format
#' @importFrom dplyr %>% arrange mutate
#' @importFrom forcats fct_inorder
#' @return SummarizedExperiment object.
#' @export
#'
#' @examples
read_skyline <- function(files) {
  names(files) = basename(files)
  datalist = lapply(files, .read_skyline_file) %>% .uniform.attrs()
  
  original_data = datalist %>%
    mutate(Sample = fct_inorder(Sample)) %>%
    group_by(Sample) %>%
    mutate(TransitionId=c(1:n())) %>%
    ungroup() %>%
    .copy.attr(datalist) %>%
    .to_summarized_experiment()
  
  
  rowData(original_data) <- rowData(original_data) %>% 
    left_join(annotate_lipids(rowData(original_data)))
  
    # 
    # 
    # 
    # arrange(Class, filename) %>% 
    # mutate(Molecule=fct_inorder(Molecule)) %>%
    # .copy.attr(datalist)
    # 
  message("Successfully read ", length(files), " methods.\n",
          "Your data contain ",  ncol(original_data)," samples, ",
          length(unique(rowData(original_data)$Class))," lipid classes, ",
          length(unique(rowData(original_data)$Molecule))," lipid molecules."
  )
  
  return(original_data)
}

#' Add sample annotation to Skyline data frame
#'
#' @param data Skyline data.frame created by \code{\link{read.skyline}}.
#' @param annot_file CSV file with at least 2 columns, sample names and group(s).
#'
#' @importFrom dplyr %>% left_join
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#' @return Skyline data.frame with sample group information.
#' @export
#'
#' @examples
add_sample_annotation = function(ds, annot_file) {
  annot = read.csv(annot_file)
  stopifnot(ncol(annot) > 1)
  
  # check if any column is named "sample", otherwise take the first column
  sample_col = grep("Sample", colnames(annot), ignore.case = T)
  if(length(sample_col) > 0) {
    sample_col = colnames(annot)[sample_col][[1]]
  } else {
    warning("No column named 'Sample', taking the first column as sample names")
    sample_col = colnames(annot)[[1]]
  }
  annot_cols = colnames(annot)[colnames(annot) != sample_col]
  col_data = colData(ds)
  
  if (any (annot_cols %in% colnames(col_data)) ) {
    annot_cols_exist = annot_cols[annot_cols %in% colnames(col_data)]
    warning("These annotation columns already exist in data: ", paste(annot_cols_exist, collapse = ", "))
    annot_cols = annot_cols[!annot_cols %in% annot_cols_exist]
    
    if(length(annot_cols) == 0) return(ds)
    
    annot = annot[, !colnames(annot) %in% annot_cols_exist]
  }
  
  colData(ds) <- col_data %>% left_join(annot, by=c(rowname=sample_col))
  
  return(ds)
}

####################################################################################################
.is_pivoted <- function(d) {
  return (!any(grepl("^Area|^Height|^Area\\.Normalized", colnames(d))))
}
#' Internal method to read skyline file
#' @param file skyline exported file in CSV format
#' @importFrom dplyr %>% vars matches mutate_at
#' @importFrom tidyr gather spread separate
#' @return std data.frame
.read_skyline_file <- function(file) {
  original_data = read.csv(file, stringsAsFactors = FALSE)
  
  class_cols = grep(class_pattern, colnames(original_data))
  if (length(class_cols) == 0)
    stop("At least one of these columns should be exported from Skyline report.", class_names)
  
  molecule_cols = grep(molecule_pattern, colnames(original_data))  
  if (length(molecule_cols) == 0)
    stop("At least one of these columns should be exported from Skyline report.", molecule_names)
  
  colnames(original_data)[[ molecule_cols[[1]] ]] = "Molecule"
  
  if ( !any(grepl("Area|Height|Area\\.Normalized", colnames(original_data))) )
    stop("At least one of these columns should be exported from Skyline report.", intensity_names)
  
  if (.is_pivoted(original_data)) {
    return (.read_pivoted(original_data))
  } else {
    return (.read_not_pivoted(original_data))
  }
}

.read_not_pivoted <- function(original_data) {
  intensity_cols = grep("^Area|^Height|^Area\\.Normalized", colnames(original_data))
  replicate_cols = grep(replicate_pattern, colnames(original_data))
  if(length(replicate_cols) == 0) 
    stop("In Skyline report, you should either export Replicate column or pivot by replicates")
  
  colnames(original_data)[[ replicate_cols[[1]] ]] = "Sample"
  ret = original_data %>% mutate_at(vars(matches(measure_pattern)), as.numeric)
  measure_cols = colnames(original_data)[grep(measure_pattern, colnames(original_data))]
  
  attr(ret, "skyline") <- list(skyline=TRUE, measures=measure_cols, intensity_cols=colnames(original_data)[intensity_cols])
  return (ret)
}

.read_pivoted <- function(original_data) {
  intensity_cols = grep("Area|Height|Area\\.Normalized", colnames(original_data))
  intensity_cols = colnames(original_data)[intensity_cols]
  samples = unique( sub("(.*)(Area|Height|Area\\.Normalized)", "\\1", intensity_cols) )
  samples_pattern = paste("^", samples, sep="", collapse = "|")
  sample_cols = grep(samples_pattern, colnames(original_data))
  measures = unique( sub(samples_pattern, "", colnames(original_data)[sample_cols]) )
  measures_pattern = paste(measures, collapse = "|")
  
  colnames(original_data) = sub(paste0("\\.(", measures_pattern, ")"), "###\\1", colnames(original_data))
  
  ret = original_data %>% 
    mutate_at(vars(matches(measures_pattern)), as.numeric) %>% 
    gather(sample.measure, value, sample_cols) %>% 
    separate(sample.measure, into=c("Sample", "measure"), sep="###", extra = "merge") %>%
    spread(measure, value)
  
  attr(ret, "skyline") <- list(skyline=TRUE, measures=measures, intensity_cols=intensity_cols)
  return(ret)
}