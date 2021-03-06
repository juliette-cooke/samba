% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_multi_distrib}
\alias{plot_multi_distrib}
\title{Plot multiple distributions from raw samples}
\usage{
plot_multi_distrib(sdata_filtered, thr, max_n = 10)
}
\arguments{
\item{sdata}{A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites}

\item{zscore}{A dataframe of zscore calculated using calc_zscore}
}
\description{
Compare WT and MUT flux distributions for the top set of zscore metabolites for multiple diseases.
}
