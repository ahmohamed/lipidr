# Metabolomics Workbench integraation.

#' Retrieve list of lipidomics studies from Metabolomics Workbench.
#'
#' This functions uses Metabolomics Workbench REST API to retrieve summaries
#' of lipidomics studies.
#'
#' @param keyword A keyword to search for in Metabolomics Workbench studies.
#'
#' @return A data frame with studies matching the keyword. Study ID, title,
#'   author and details are retrieved.
#'
#' @importFrom tidyr unite spread unnest
#' @importFrom rlang syms !!!
#' @importFrom utils read.delim
#' @export
#' @family Metabolomics Workbench
#'
#' @examples
#' # list_mw_studies()
list_mw_studies <- function(keyword = "lipid") {
  cols <- c(
    "study_id", "study_title", "study_type",
    "first_name", "last_name", "institute", "department", "submit_date",
    "subject_species", "study_summary"
  )
  url <- paste0(
    "https://www.metabolomicsworkbench.org/rest/study/study_title/",
    keyword, "/summary/txt"
  )
  read.delim(url, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    filter(V1 %in% cols) %>%
    mutate(cs = cumsum(V1 == "study_id")) %>%
    spread(key = V1, value = V2) %>%
    select(!!!rlang::syms(cols)) %>%
    unite("author", first_name, last_name, sep = " ")
}

# colname created in list_mw_studies
utils::globalVariables(c("V1", "V2", "first_name", "last_name"))


#' Download and parse full data for a study from Metabolomics Workbench.
#'
#' This functions uses Metabolomics Workbench REST API to retrieve study
#' data using a study ID. The function returns a LipidomicsExperiment where users
#' can directly apply `lipidr` analysis workflow.
#'
#' @param study_id The Metabolomics Workbench study ID.
#'
#' @return A LipidomicsExperiment object containing clinical and lipid intensity
#'   data.
#'
#' @export
#' @family Metabolomics Workbench
#'
#' @examples
#' # fetch_mw_study("ST001111")
fetch_mw_study <- function(study_id) {
  url <- paste0(
    "https://www.metabolomicsworkbench.org/rest/study/study_id/",
    study_id, "/mwtab/txt"
  )
  read_mwTab(url)
}

#' Parse mwTab file into a LipidomicsExperiment.
#'
#' @param mwTab File path or url for a mwTab file.
#'
#' @return A LipidomicsExperiment object containing clinical and lipid intensity
#'   data.
#'
#' @export
#' @family Metabolomics Workbench
read_mwTab <- function(mwTab) {
  txt <- readLines(mwTab)
  analyses <- .mw_analysis_lists(txt)
  prased_analyses <- lapply(analyses, .parse_mw_analysis)
  combined <- do.call("rbind", prased_analyses)
  mcols(assays(combined)) <- list(logged = FALSE, normalized = FALSE)
  combined
}

#' Parse a Metabolomics Workbench data matrix into a LipidomicsExperiment.
#'
#' Data matrix downloaded from Metabolomics Workbench are parsed into
#' a LipidomicsExperiment object to enable `lipidr` workflow analysis.
#'
#' @param file File path or url for the file containing the data matrix.
#'
#' @return A LipidomicsExperiment object containing clinical and lipid intensity
#'   data.
#'
#' @export
#' @family Metabolomics Workbench
read_mw_datamatrix <- function(file) {
  .data <- read.delim(
    file,
    sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE
  )

  d <- .data[-1, ]
  colnames(d)[[1]] <- "Molecule"
  original_names <- d$Molecule
  stereo <- "\\((\\d+[ZE]\\.*)+\\)"
  adduct <- " \\[.*$"
  d$Molecule <- sub(adduct, "", d$Molecule)
  d$Molecule <- gsub(stereo, "", d$Molecule)
  d$Molecule <- sub("^(\\w+)-", "\\1", d$Molecule)
  d$Molecule <- sub("\\d+\\:\\d+ \\(", "\\(", d$Molecule)
  d <- as_lipidomics_experiment(d)
  rownames(d) <- original_names

  col_data <- t(.data[1, ]) %>% as.data.frame()
  colnames(col_data)[[1]] <- "Factors"
  col_data <- col_data %>%
    slice(-1) %>%
    rownames_to_column("Sample") %>%
    .unnest_key_value(Factors, kv_sep = "\\:", list_sep = "\\|") %>%
    `rownames<-`(.$Sample) %>%
    DataFrame()
  colData(d) <- col_data[colnames(d), -1, drop = FALSE]
  d
}

.mw_analysis_lists <- function(mwTab) {
  split(mwTab, cumsum(grepl("#METABOLOMICS WORKBENCH", mwTab)))
}

.parse_mw_analysis <- function(mw_analysis) {
  analysis_id <- sub(
    ".*ANALYSIS_ID\\:([[:alnum:]]+).*", "\\1", mw_analysis[[1]]
  )
  sample_data <- .sampledata_from_mw_analysis(mw_analysis)
  ret <- .msdata_from_mw_analysis(mw_analysis)
  colData(ret) <- to_df(ret, "col") %>%
    .left_join_silent(sample_data) %>%
    toDataFrame(row.names.col = "Sample")
  rowData(ret)$filename <- analysis_id
  ret
}

.msdata_from_mw_analysis <- function(mw_analysis) {
  from <- grep("^MS_METABOLITE_DATA_START", mw_analysis) + 1
  to <- grep("^MS_METABOLITE_DATA_END", mw_analysis) - 1
  read_mw_datamatrix(textConnection(mw_analysis[from:to]))
}

.sampledata_from_mw_analysis <- function(mw_analysis) {
  sampledata <- mw_analysis[grep("^SUBJECT_SAMPLE_FACTORS", mw_analysis)]
  sampledata <- read.delim(
    textConnection(sampledata),
    sep = "\t", check.names = FALSE, header = FALSE
  )[, -c(1, 2)]
  colnames(sampledata) <- c("Sample", "Factors", "AdditionalData")
  sampledata %>%
    mutate(Sample = as.factor(Sample)) %>%
    .unnest_key_value(Factors, kv_sep = "\\:", list_sep = "\\|") %>%
    .unnest_key_value(AdditionalData, kv_sep = "=", list_sep = "; ")
}

.unnest_key_value <- function(.data, col, kv_sep, list_sep) {
  col <- enquo(col)
  .data %>%
    mutate(col = strsplit(as.character(!!col), list_sep)) %>%
    unnest(col) %>%
    separate(col, c("annot", "val"), sep = kv_sep) %>%
    mutate(annot = ifelse(annot == "Sample", "SampleType", annot)) %>%
    spread(annot, val) %>%
    select(-!!col)
}

# colname created in .unnest_key_value
utils::globalVariables(c("annot", "val", "Factors", "AdditionalData"))
