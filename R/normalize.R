#' Perform Probabilistic Quotient Normalization for intensities.
#'
#' Perform Probabilistic Quotient Normalization (PQN) for sample intensities.
#' The PQN method determines a dilution factor for each sample by comparing
#' the distribution of quotients between samples and a reference spectrum,
#' followed by sample normalization using this dilution factor.
#' The reference spectrum in this method is the average lipid abundance of all
#' samples (excluding blanks).
#'
#' @param data LipidomicsExperiment object.
#' @param measure Which measure to use as intensity, usually Area,
#'   Area Normalized or Height. Default is `Area`.
#' @param exclude Samples to exclude, can be either: \cr
#'   "blank" - automatically detected blank samples and exclude them
#'   logical vector with the same length as samples. Default.
#'
#' @param log Whether the normalized values should be log2 transformed. Default
#'   is `TRUE`.
#'
#' @return A LipidomicsExperiment object with normalized values
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
#' data_normalized <- normalize_pqn(
#'   d_summarized,
#'   measure = "Area", exclude = "blank", log = TRUE
#' )
normalize_pqn <- function(data, measure = "Area",
                          exclude = "blank", log = TRUE) {
  data <- .prenormalize_check(data, measure, exclude)
  m <- assay(data, measure)

  # factor_n = median ( lipid_i_n/ avg(lipid_i) )
  factor_n <- apply(
    m / rowMeans(m, na.rm = TRUE), 2, median,
    na.rm = TRUE
  )
  normalized_m <- apply(m, 1, function(x) x/factor_n) %>% t()
  rownames(normalized_m) <- rownames(m)
  assay(data, measure) <- normalized_m
  # assay(data, measure) <- m / apply(
  #   m / rowMeans(m, na.rm = TRUE), 2, median,
  #   na.rm = TRUE
  # )
  data <- set_normalized(data, measure, TRUE)
  return(.log_data(data, measure, log))
}

#' Normalize each class by its corresponding internal standard(s).
#'
#' Normalize each class by its corresponding internal standard(s).
#' Lipid classes are normalized using corresponding internal standard(s)
#' of the same lipid class. If no corresponding internal standard is found
#' the average of all measured internal standards is used instead.
#'
#' @param data LipidomicsExperiment object.
#' @param measure Which measure to use as intensity, usually Area,
#'   Area Normalized or Height. Default is `Area`.
#' @param exclude Samples to exclude, can be either: \cr
#'   "blank" - automatically detected blank samples and exclude them
#'   logical vector with the same length as samples. Default.
#' @param log whether the normalized values should be log2 transformed. Default
#'   is `TRUE`.
#'
#' @return A LipidomicsExperiment object with normalized values. Each molecule
#'     is normalized against the internal standard from the same class.
#'
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
#' data_norm_istd <- normalize_istd(
#'   d_summarized,
#'   measure = "Area", exclude = "blank", log = TRUE
#' )
normalize_istd <- function(data, measure = "Area",
                           exclude = "blank", log = TRUE) {
  data <- .prenormalize_check(data, measure, exclude)
  istd <- rowData(data)$istd
  if (sum(istd) == 0) {
    stop("No internal standards found in your lipid list.")
  }
  m <- assay(data, measure)
  mistd <- m[istd, ]

  # istd_n = istd_ni / mean(istd_i)
  mistd <- mistd / rowMeans(mistd, na.rm = TRUE)

  # per class:
  x_dimname <- metadata(data)$dimnames[[1]]
  istd_list <- to_df(data) %>%
    group_by(filename, Class) %>%
    mutate(istd_list = list(as.character( (!!sym(x_dimname))[istd] ))) %>%
    .$istd_list

  assay(data, measure) <- laply(seq_along(istd_list), function(i) {
    if (length(istd_list[[i]]) == 0) {
      f <- 1
    } else if (length(istd_list[[i]]) == 1) {
      f <- mistd[istd_list[[i]], ]
    } else {
      f <- colMeans(mistd[istd_list[[i]], ], na.rm = TRUE)
    }

    return(m[i, ] / f)
  })
  return(.log_data(data, measure, log))
}

.prenormalize_check <- function(data, measure, exclude) {
  if (!measure %in% assayNames(data)) {
    stop(measure, " is not in the dataset.")
  }
  if (is_normalized(data, measure)) {
    stop(measure, " is already normalized.")
  }
  if (!is.null(exclude)) {
    if (length(exclude) == 1 && "blank" %in% exclude) {
      data <- data[, !.is_blank(data, measure)]
    } else {
      excluded_cols <- colnames(data[, exclude])
      data <- data[, !colnames(data) %in% excluded_cols]
      if (ncol(data) == 0) {
        stop("You cannot exclude all samples.")
      }
    }
  }
  assay_ <- assay(data, measure)
  if (any(!is.finite(assay_))) {
    warning(measure, " contains missing/non-finite values.",
            " Replacing with mimnum detected value.")
    assay_[!is.finite(assay_)] <- min(assay_, na.rm = TRUE)
  }
  # if (any(assay_ < 1)) {
  #   warning(measure, " contains values < 1. Replacing with 1.")
  #   assay_[assay_ < 1] <- 1
  # }
  assay(data, measure) <- assay_
  data
}


.log_data <- function(data, measure, log) {
  assay_ <- assay(data, measure)
  if (any(assay_ < 1)) {
    warning(measure, " contains values < 1. Replacing with 1.")
    assay_[assay_ < 1] <- 1
  }
  if (log && !is_logged(data, measure)) {
    assay(data, measure) <- log2(assay_)
    data <- set_logged(data, measure, TRUE)
  }
  data
}
