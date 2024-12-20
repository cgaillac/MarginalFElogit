% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/choose_bandwidth_and_estim_condl_proba.R
\name{choose_bandwidth_and_estim_condl_proba}
\alias{choose_bandwidth_and_estim_condl_proba}
\title{This function chooses the bandwidths to be used for the estimation through
local linear regressions of the conditional distribution of S, as described
in online appendix B of the DDL paper. It then proceeds to estimate these
conditional probabilities, and saves them in estimators.}
\usage{
choose_bandwidth_and_estim_condl_proba(data, estimators)
}
\arguments{
\item{data}{is an environment variable containing the relevant data,
formatted by format_data (see function documentation):
 - data$Y a matrix of size n x Tmax containing the values of the
dependent variable Y.
 - data$X an array of size n x Tmax x dimX containing the values of
the covariates X.
 - data$clusterIndexes a vector of size n containing the index of the
cluster each observation belongs to. The computed asymptotic variance is
clustered.}

\item{estimators}{is an environment variable containing the results from
logit and FE logit estimation:
   -  estimators$beta_hat$beta_hat is the CMLE estimate of the FE logit
slope parameter.
   -  estimators$alphaFElogit is the estimate of the constant parameter in
a standard logit model using the estimated CMLE slope parameter. If empty,
it is estimated.}
}
\value{
The function does not return anything, but the estimators environment
has the following parameters added:
- estimators$h_local_lin_proba a vector of length (Tmax + 1) containing, in
  j-th position, the bandwidth used to estimate the P(S = j - 1 | X)'s.
- estimators$condlProbas: a matrix of size n x (Tmax + 1) containing,
  in position (i, j), the estimate for P(S = j - 1 | X) at the i-th
  observation.
- estimators$densityEstimate: a matrix of size n x (Tmax + 1) containing,
  at each row (individual), the estimated density for having covariates
  (X_1, ..., X_T). Each column represents the value found using the
  corresponding bandwidths from h.
}
\description{
This function chooses the bandwidths to be used for the estimation through
local linear regressions of the conditional distribution of S, as described
in online appendix B of the DDL paper. It then proceeds to estimate these
conditional probabilities, and saves them in estimators.
}
