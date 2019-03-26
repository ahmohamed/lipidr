#' Analysis workflow for targeted lipidomics.
#'
#' \code{lipidr} implements a series of functions to facilitate inspection,
#' analysis and visualization of targeted lipidomics datasets. \code{lipidr}
#' takes exported Skyline CSV as input, allowing for multiple methods to
#' be analyzed togther.
#' 
#' \code{lipidr} represents Skyline files as annotated data.frames, which can
#' be easily manipulated by a wide variety of R packages. Sample annotations, such 
#' as sample group or other clinical information can be easily loaded. \code{lipidr} 
#' generates various plots, such as PCA score plots and box plots, for quality control of samples 
#' and measured lipids. Normalization methods with and without internal standards are also supported.
#' 
#' Differential analysis can be performed using any of the loaded clinical variables,
#' which can be readily visualized as volcano plots. A novel lipid set enrichment analysis (LSEA)
#' is implemented to detect preferential enrichment of certain lipid classes,
#' chain lengths or saturation patterns. Plots for visualization of enrichment results 
#' are also implemented.
#' 
#' @author Ahmed Mohamed \email{ahmed.mohamed@@qimrberghofer.edu.au}
#' @name lipidr-package
#' @docType package
#' 
NULL


#' lipidDefaults.
#' A set of defualt mappings and annoataion used internally to correctly parse
#' lipid molecule names.
#'
#' @docType data
#' @name lipidDefaults
#' @examples
#' data(lipidDefaults)
NULL

#' Example dataset (normalized and log transformed)
#' 
#' Please see \link{normalize_pqn} for details on how to generate this dataset.
#' 
#' @docType data
#' @name data_normalized
#' @examples
#' data(data_nomalized)
NULL


.myDataEnv <- new.env(parent=emptyenv()) # not exported

.onAttach <- function(lib, pkg) {
   utils::data(lipidDefaults, envir=.myDataEnv)
  .myDataEnv$interactive = FALSE
}

#' Activate interactive graphics
#'
#' Use this function to turn on/off interactive graphics
#' plotting. Interactive plots require \code{plotly}
#' to be installed. Interactive graphics are disabled by default.
#' 
#' @param interactive Should interactive plots be displayed? Default is TRUE.
#'
#' @return None
#' @export
#'
#' @examples
#' data(data_normalized)
#' 
#' # plot the variation in intensity and retention time of all measured lipids in QC samples
#' d_qc = data_normalized[, data_normalized$group == "QC"]
#' plot_molecule_cv(d_qc, "Area")
use_interactive_graphics <- function(interactive=TRUE) {
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' must be installed for interactive graphics")
    }
  }
  .myDataEnv$interactive = interactive
}

# Internal function to update the dataset, in case of future changes.
# Should be called once before each package release.
.updateData <- function(pkgpath) {
  datadir = system.file("extdata", package="lipidr")
  filelist = list.files(datadir, "data.csv", full.names = TRUE)
  d = read_skyline(filelist)
  clinical_file = system.file("extdata", "clin.csv", package="lipidr")
  d = add_sample_annotation(d, clinical_file)
  d_summarized = summarize_transitions(d, method = "average")
  data_normalized = normalize_pqn(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
  save(data_normalized, file=file.path(pkgpath, "data/data_normalized.rda"))
}