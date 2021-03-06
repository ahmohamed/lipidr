% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioc.R
\name{LipidomicsExperiment}
\alias{LipidomicsExperiment}
\title{Constructor for Lipidomics experiment from list of assays}
\usage{
LipidomicsExperiment(assay_list, metadata, colData = NULL, rowData = NULL)
}
\arguments{
\item{assay_list}{A list or SimpleList of matrix-like elements,
or a matrix-like object. Passed to \code{\link[=SummarizedExperiment]{SummarizedExperiment()}}.}

\item{metadata}{A list containing arbitrary information about the experiment.
It should at least contain 2 elements: \itemize{
\item dimnames    2-element character vector with dimension names
\item summarized   Has transitions been summarized?
}}

\item{colData}{An optional DataFrame describing the samples (contains
clinical information). Row names, if present, become the column names of
the LipidomicsExperiment.}

\item{rowData}{A DataFrame object describing the rows (contains generated
lipid annotations). Row names, if present, become the row names of the
SummarizedExperiment object. The number of rows of the DataFrame
must be equal to the number of rows of the matrices in assays.}
}
\value{
LipidomicsExperiment object
}
\description{
Constructor for Lipidomics experiment from list of assays
}
