#' Perform Probabilistic Quotient Normalization for intensities.
#'
#' Perform Probabilistic Quotient Normalization (PQN) for sample intensities.
#' The PQN method determines a dilution factor for each sample by comparing
#' the distribution of quotients between samples and a reference spectrum, followed
#' by sample normalization using this dilution factor.
#' The reference spectrum in this method is the average lipid abundance of all
#' samples (excluding blanks).
#'
#' @param data SkylineExperiment object created by [read_skyline()].
#' @param measure Which measure to use as intensity, usually Area, Area.Normalized or Height.
#' @param exclude Samples to exclude, can be either: \cr
#' "blank" - automatically detected blank samples and exclude them
#' logical vector with the same length as samples.
#'
#' @param log Whether the normalized values should be log2 transformed.
#'
#' @return A SkylineExperiment object with normalized values
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
#' datadir <- system.file("extdata", package = "lipidr")
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#' clinical_file <- system.file("extdata", "clin.csv", package = "lipidr")
#' d <- add_sample_annotation(d, clinical_file)
#' d_summarized <- summarize_transitions(d, method = "average")
#' 
#' # Normalize data that have been summarized (single value per molecule).
#' data_normalized <- normalize_pqn(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
normalize_pqn <- function(data, measure = "Area", exclude = "blank", log = TRUE) {
  if (mcols(assays(data), use.names = TRUE)[measure, "normalized"]) {
    stop(measure, " is already normalized")
  }
  if (!is.null(exclude)) {
    if (exclude == "blank") {
      data <- data[, !.is_blank(data)]
    } else {
      data <- data[, exclude]
    }
  }
  m <- assay(data, measure)

  # factor_n = median ( lipid_i_n/ avg(lipid_i) )
  assay(data, measure) <- m / apply(m / rowMeans(m, na.rm = TRUE), 2, median, na.rm = TRUE)
  mcols(assays(data), use.names = TRUE)[measure, "normalized"] <- TRUE

  if (log && !mcols(assays(data), use.names = TRUE)[measure, "logged"]) {
    assay(data, measure) <- log2(assay(data, measure))
    mcols(assays(data), use.names = TRUE)[measure, "logged"] <- TRUE
  }

  return(data)
}

#' Normalize each class by its corresponding internal standard(s).
#'
#' Normalize each class by its corresponding internal standard(s).
#' Lipid classes are normalized using corresponding internal standard(s)
#' of the same lipid class. If no corresponding internal standard is found
#' the average of all measured internal standards is used instead.
#'
#' @param data Skyline data.frame created by [read_skyline()].
#' @param measure which measure to use as intensity, usually Area, Area.Normalized or Height.
#' @param exclude Samples to exclude, can be either: \cr
#' "blank" - automatically detected blank samples and exclude them
#' logical vector with the same length as samples.
#' @param log whether the normalized values should be log2 transformed.
#'
#' @return A SkylineExperiment object with normalized values. Each molecule
#'     is normalized against the internal standard from the same class.
#'
#' @importFrom dplyr %>% select group_by mutate filter ungroup left_join inner_join
#' @importFrom rlang sym UQ
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", package = "lipidr")
#' filelist <- list.files(datadir, "data.csv", full.names = TRUE)
#' d <- read_skyline(filelist)
#' clinical_file <- system.file("extdata", "clin.csv", package = "lipidr")
#' d <- add_sample_annotation(d, clinical_file)
#' d_summarized <- summarize_transitions(d, method = "average")
#' 
#' # Normalize data that have been summarized (single value per molecule).
#' data_norm_itsd <- normalize_itsd(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
normalize_itsd <- function(data, measure = "Area", exclude = "blank", log = TRUE) {
  if (mcols(assays(data), use.names = TRUE)[measure, "normalized"]) {
    stop(measure, " is already normalized")
  }
  if (!data@attrs$summarized) {
    stop("Data should be summarized using summarize_transitions")
  }

  if (!is.null(exclude)) {
    if (exclude == "blank") {
      data <- data[, !.is_blank(data)]
    } else {
      data <- data[, exclude]
    }
  }
  itsd <- rowData(data)$itsd
  if (sum(itsd) == 0) {
    stop("No internal standards found in your lipid list.")
  }
  m <- assay(data, measure)
  mitsd <- m[itsd, ]

  # itsd_n = itsd_ni / mean(itsd_i)
  mitsd <- mitsd / rowMeans(mitsd, na.rm = TRUE)

  # per class:
  itsd_list <- to_df(data) %>%
    group_by(filename, Class) %>%
    mutate(itsd_list = list(as.character(MoleculeId[itsd]))) %>%
    .$itsd_list

  assay(data, measure) <- laply(seq_along(itsd_list), function(i) {
    if (length(itsd_list[[i]]) == 0) {
      f <- 1
    } else if (length(itsd_list[[i]]) == 1) {
      f <- mitsd[itsd_list[[i]], ]
    } else {
      f <- colMeans(mitsd[itsd_list[[i]], ], na.rm = TRUE)
    }

    return(m[i, ] / f)
  })

  if (log && !mcols(assays(data), use.names = TRUE)[measure, "logged"]) {
    assay(data, measure) <- log2(assay(data, measure))
    mcols(assays(data), use.names = TRUE)[measure, "logged"] <- TRUE
  }

  return(data)
}
