% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter_zscore}
\alias{filter_zscore}
\title{Filter the metabolites based on their zscore}
\usage{
filter_zscore(zscore, thr = 2, max_n = 10)
}
\arguments{
\item{d_all}{A list containing WT and MUT dataframes}
}
\value{
A means list containing WT and MUT mean dataframes
}
\description{
This function filters metabolites based on their previously calculated zscore and density.
}
