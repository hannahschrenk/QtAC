% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QtAC_package.R
\name{QtAC_CTE}
\alias{QtAC_CTE}
\title{QtAC_CTE}
\usage{
QtAC_CTE(
  data,
  num_timepoints = 5,
  javapath,
  noise_level = "1e-20",
  num_permcheck = 1000L,
  k = 1L,
  k_tau = 1L,
  l = 1L,
  l_tau = 1L,
  delay = 1L,
  condEmbedDims = 1L,
  cond_taus = 1L,
  cond_delays = 1L,
  mode = c("average", "local"),
  save = FALSE,
  filename = "result_QtAC"
)
}
\arguments{
\item{data}{data array containing time series of the system's components' abundance data}

\item{num_timepoints}{length of the time windows of abundance data serving as basis of the collective transfer entropy estimations}

\item{javapath}{path of the file "MTinfodynamics.jar"}

\item{noise_level}{amount of random Gaussian noise added in the estimation}

\item{num_permcheck}{number of surrogate samples to bootstrap to generate the distribution in the significance test}

\item{k}{embedding length of destination past history to consider}

\item{k_tau}{embedding delay for the destination variable}

\item{l}{embedding length of source past history to consider}

\item{l_tau}{embedding delay for the source variable}

\item{delay}{time lag between last element of source and destination next value}

\item{condEmbedDims}{array of embedding lengths for each conditional variable}

\item{cond_taus}{array of embedding delays for the conditional variables}

\item{cond_delays}{array of time lags between last element of each conditional variable and destination next value}

\item{mode}{Transfer entropy is either averaged over time windows ("average") or local values are estimated ("local"). In the latter case, the value of num_timepoints is ignored.}

\item{save}{If save=TRUE, the output is saved in a file called filename.}
}
\value{
list of two lists (adjacency and corresponding significance matrices)
}
\description{
This function calculates the collective transfer entropy of each species for shifting time windows of fixed length or local values.
The output is a list of adjacency matrices and the corresponding significance matrices.
}
