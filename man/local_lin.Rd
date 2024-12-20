% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_lin.R
\name{local_lin}
\alias{local_lin}
\title{Estimates E(Y|X) at the observed X by local linear regression.
We use a Gaussian kernel, and also output the kernel estimator of the density of X.}
\usage{
local_lin(data, h = 1)
}
\arguments{
\item{data}{is an environment variable containing the relevant data:
- data$Y a matrix of size n x Tmax containing the values of the dependent variable Y.
- data$X an array of size n x Tmax x dimX containing the values of the covariates X.}

\item{h}{(default 1) the bandwidth for the local linear regressions. We increase it at points
where the regression is degenerate otherwise.}
}
\value{
a list containing.
 - condlMeans: a vector of length n containing the
   estimates for P(S = j - 1 | X) at each observation.
 - densityEstimate: a vector of length n containing, at each row (individual),
   the estimated density for having covariates X.
}
\description{
Estimates E(Y|X) at the observed X by local linear regression.
We use a Gaussian kernel, and also output the kernel estimator of the density of X.
}
