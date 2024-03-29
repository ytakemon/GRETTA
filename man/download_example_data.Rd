% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_example_data.R
\name{download_example_data}
\alias{download_example_data}
\title{Download minimal example data}
\usage{
download_example_data(path = NULL, test = FALSE)
}
\arguments{
\item{path}{character. Path for example data. Default: NULL.}

\item{test}{logical. Used for examples. Default: FALSE.}
}
\value{
A data frame containing rank at which the threshold should be drawn for positive and negative co-essential genes.
}
\description{
Downloads minimal data for examples shown in the tutorial. This can also be downloaded manually from:
https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2_sample
}
\details{
All files are stored in './GRETTA_example'
}
\examples{
download_example_data(path = ".", test = TRUE)

}
