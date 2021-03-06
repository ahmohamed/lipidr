% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manip.R
\name{remove_non_parsed_molecules}
\alias{remove_non_parsed_molecules}
\title{Remove molecules that couldn't be parsed by \code{lipidr} from the dataset}
\usage{
remove_non_parsed_molecules(data)
}
\arguments{
\item{data}{LipidomicsExperiment object.}
}
\value{
A filtered LipidomicsExperiment object.
}
\description{
Remove molecules that couldn't be parsed by \code{lipidr} from the dataset
}
\examples{
data(data_normalized)
remove_non_parsed_molecules(data_normalized)
}
