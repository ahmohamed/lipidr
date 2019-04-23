#' Analysis workflow for targeted lipidomics
#'
#' `lipidr` implements a series of functions to facilitate inspection,
#' analysis and visualization of targeted lipidomics datasets. `lipidr`
#' takes exported Skyline CSV as input, allowing for multiple methods to
#' be analyzed together.
#'
#' `lipidr` represents Skyline files as SummarizedExperiment objects, which can
#' easily be integrated with a wide variety of Bioconductor packages.
#' Sample annotations,
#' such as sample group or other clinical information can be loaded. `lipidr`
#' generates various plots, such as PCA score plots and box plots, for quality
#' control of samples and measured lipids. Normalization methods with and
#' without internal standards are also supported.
#'
#' Differential analysis can be performed using any of the loaded clinical
#' variables, which can be readily visualized as volcano plots. A novel lipid
#' set enrichment analysis (LSEA) is implemented to detect preferential
#' enrichment of certain lipid classes, chain lengths or saturation patterns.
#' Plots for the visualization of enrichment results are also implemented.
#'
#' @author Ahmed Mohamed \email{ahmed.mohamed@@qimrberghofer.edu.au}
#' @name lipidr-package
#' @docType package
#' @aliases lipidr lipidr-package
#'
NULL


#' Description of lipidr datasets
#'
#' lipidr-package has 3 datasets: \itemize{
#'   \item data_normalized    Example lipidomics dataset,
#'     normalized & log2-transformed.
#'   \item lipidDefaults   A list of default mappings and annotations for lipids.
#'   \item lipidnames_pattern   A list of patterns used in parsing lipid names.
#'  }
#'  See below for detailed descrition of each dataset.
#'
#' @docType data
#' @name lipidr-data
#' @family lipidr datasets
#' @examples data(data_normalized)
NULL

#' Default values for lipidr internal functions
#' A set of default mappings and annotation used internally to correctly parse
#' lipid molecule names.
#'
#' @docType data
#' @name lipidDefaults
#' @family lipidr datasets
#' @examples data(lipidDefaults)
"lipidDefaults"

#' Example dataset (normalized and log2 transformed)
#'
#' A dataset containing MRM mass spectrometry-based lipidomics data from murine
#' serum samples. Mice were fed a normal or high-fat diet and had access to
#' normal drinking water or drinking water containing the bile acid
#' deoxycholic acid. Lipid peaks were integrated using Skyline and exported
#' results were imported into R using `lipidr`. The dataset has been normalized
#' and log2 transformed. Please see \link{normalize_pqn} for details on how to
#' generate this dataset.
#'
#' @docType data
#' @name data_normalized
#' @family lipidr datasets
#' @examples data(data_normalized)
"data_normalized"

#' Patterns used in parsing lipid names
#'
#' A collection of patterns to extract lipid class and chain information
#' from lipid names. Used internally by the package.
#'
#' @docType data
#' @name lipidnames_pattern
#' @family lipidr datasets
#' @examples data(lipidnames_pattern)
"lipidnames_pattern"

.myDataEnv <- new.env(parent = emptyenv()) # not exported

.onAttach <- function(lib, pkg) {
  .myDataEnv$interactive <- FALSE
}

#' Activate interactive graphics
#'
#' Use this function to turn on/off interactive graphics
#' plotting. Interactive plots require `plotly`
#' to be installed. Interactive graphics are disabled by default.
#'
#' @param interactive Should interactive plots be displayed? Default is TRUE.
#'
#' @return None
#' @export
#'
#' @examples
#' data(data_normalized)
#' use_interactive_graphics()
#'
#' # plot the variation in intensity and retention time of all measured
#' #  lipids in QC samples
#' d_qc <- data_normalized[, data_normalized$group == "QC"]
#' # plot_molecules(d_qc, "cv", "Area")
#'
#' # turn off interactivity
#' use_interactive_graphics(interactive = FALSE)
use_interactive_graphics <- function(interactive = TRUE) {
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' must be installed for interactive graphics")
    }
  }
  .myDataEnv$interactive <- interactive
}
