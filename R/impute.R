#' Impute missing values in a LipidomicsExperiment
#'
#' @param data LipidomicsExperiment object.
#' @param measure Which measure to use as intensity, usually Area,
#'   Area Normalized or Height. Default is `Area`.
#' @param method The imputation method to use. All methods are wrappers for
#' [imputeLCMD] package. These include \itemize{
#'   \item knn    Wraps [imputeLCMD:impute.knn]. Default. This requires an
#'   additional argument `K` (Number of neighbors used to infer the missing data).
#'   \item svd    Wraps [imputeLCMD:impute.wrapper.SVD]. This requires an
#'   additional argument `K` (Number of principal components to use).
#'   \item mle    Wraps [imputeLCMD:impute.wrapper.MLE],
#'   \item minDet    Wraps [imputeLCMD:impute.MinDet],
#'   \item minProb    Wraps [imputeLCMD:impute.MinProb],
#'   \item zero    Wraps [imputeLCMD:impute.ZERO],
#'  }

#' @param ... Other arguments passed to the imputation method.
#'
#' @return LipidomicsExperiment object with missing values imputed.
#' @export
#' @importFrom imputeLCMD impute.wrapper.KNN impute.wrapper.SVD impute.wrapper.MLE
#' @importFrom imputeLCMD impute.QRILC impute.MinDet impute.MinProb impute.ZERO
#'
#' @examples
#' data(data_normalized)
#'
#' # Replace with values calculated using K-nearest neighbors
#' impute_na(data_normalized, "Area", "knn", 10)
#'
#' # Replace with values calculated from the first K principal components
#' impute_na(data_normalized, "Area", "svd", 3)
#'
#' # Replace with Maximum likelihood estimates
#' impute_na(data_normalized, "Area", "mle")
#'
#' # Replace with randomly drawn values from a truncated distribution
#' impute_na(data_normalized, "Area", "QRILC")
#'
#' # Replace with a minimal value
#' impute_na(data_normalized, "Area", "minDet")
#'
#' # Replace with randomly drawn values from a Gaussian distribution
#' # cerntered around a minimal value
#' impute_na(data_normalized, "Area", "minProp")
#'
#' # Replace with zero (not recommended)
#' impute_na(data_normalized, "Area", "zero")
impute_na <- function(data, measure = "Area",
  method = c("knn", "svd", "mle", "QRILC", "minDet", "minProb", "zero"), ...) {
  method <- match.arg(method)
  fun <- c(knn=impute.wrapper.KNN, svd=impute.wrapper.SVD, mle=impute.wrapper.MLE,
    QRILC=impute.QRILC, minDet=impute.MinDet, minProb=impute.MinProb,
    zero=impute.ZERO
  )[[method]]

  mat <- assay(data, measure)
  mat <- fun(mat, ...)
  if (method == "knn") {
    mat <- mat$data
  }
  if (method == "QRILC") {
    mat <- mat[[1]]
  }
  assay(data, measure) <- mat
  data
}
