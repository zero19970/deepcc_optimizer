% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{cross_validation}
\alias{cross_validation}
\title{Cross Validation of DeepCC Model}
\usage{
cross_validation(fs, labels, fold = 5)
}
\arguments{
\item{fs}{a data.frame containing functional spectra of patients (each row presents one sample)}

\item{labels}{a character vector containing training lables}

\item{fold}{a integer indicating the fold number of cross validation (default: 5)}
}
\value{
a numeric indicating error rate in a single run
}
\description{
This function performs cross validation of DeepCC Mode on the training data.
}
\examples{
cross_validation(tcga_fs, tcga_labels)
}
