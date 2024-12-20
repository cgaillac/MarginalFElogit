#' Computes the bounds on the AME on the dataset, at all periods given in
#' compute_T, using the method specified as entry parameter.
#'
#' @param data is an environment variable containing the data. If not already
#' formatted, it will be formatted by the format_data function (see its
#' documentation for information about the stored variables). It must contain at
#' least:
#' - data$Y is a matrix of size n x Tmax containing the values of the dependent variable Y,
#' - data$X is an array of size n x Tmax x dimX containing the values of the covariates X,
#' - data$clusterIndexes is a vector of size n x 1 that specifies the cluster each
#' observation pertains to. If it does not exist, the function enforces the default
#' setting of i.i.d. observations - the parameter takes value 1:n so that each
#' observation is its own cluster.
#' @param estimators (default empty environment) is an environment in which the
#' results of the CMLE estimation and of the non-parametric estimation of the
#' conditional distribution of S will be stored. That way, the results will be
#' saved during the first call of the compute_AME_t function and will not need
#' to be computed in the subsequent calls. Results are stored as follows:
#' - estimators$beta_hat, a list which contains the results from CMLE
#' estimation:
#' - estimators$beta_hat$beta_hat a vector of length dimX, the
#' estimated value for the slope parameter.
#' - estimators_beta_hat$phi_b a matrix of size n x dimX containing
#' the value of the influence function at each observation (rows) w.r.t.
#' each dimension of the covariates (columns).
#' - estimators$beta_hat$var_b the estimated asymptotic covariance
#' matrix, of size dimX x dimX, for the estimator beta_hat.
#' - stimators$h_local_lin_proba a vector of length (Tmax + 1) containing,
#' in j-th position, the bandwidth used to estimate the P(S = j - 1 | X)'s.
#' - estimators$condlProbas: a matrix of size n x (Tmax + 1) containing,
#' in position (i, j), the estimate for P(S = j - 1 | X) at the i-th
#' observation.
#' - estimators$densityEstimate a matrix of size n x (Tmax + 1) containing,
#' at each row (individual), the estimated density for having covariates
#' (X_1, ..., X_T). Each column represents the value found using the
#' corresponding bandwidths from h.
#' @param selectX (default NULL) a vector containing the indices of the
#' covariates w.r.t. which the AME must be computed. If null, the AME w.r.t. all
#' covariates will be computed. All variables of interest should be given in the
#' same call, as then there will be no additional cost relative to estimating
#' w.r.t. only one covariate.
#' @param compute_T (default "all") is a vector containing all periods at which
#' the AME must be computed. Alternatively, it can be "all", in which case the
#' AME will be computed successively at every period and, on top of that, the
#' average AME across all periods will also be computed using the function
#' compute_average_AMTE. Also note that non-positive values will be counted
#' backwards from the last period at which each individual is observed, as in
#' an event-study.
#' @param Option (default "quick") Estimation method to be used. If "quick" or
#' "outer" (case-insensitive) the outer bounds are computed. Otherwise, the sharp
#' bounds are computed. We recommend using the outer bounds method if the number
#' of covariates is at least three, if the number of periods observed is four or
#' more, or if the sample size is small (less than 500) or large (more than 10^4).
#' @param CIOption (default "CI2") When the outer bounds method is being used,
#' specifies which confidence interval should be used. If "CI2", the CI2
#' confidence interval is being used (see DDL, section 4.2), otherwise the CI3
#' confidence interval will be used (see DDL, appendix C). We recommend using
#' CI3 only if the user suspects the FE logit model may be a severely
#' misspecified model for the data.
#' @param alpha (default 0.05) desired asymptotic level of the estimated
#' confidence intervals
#' @param nbCores (default 4) number of cores to be used for parallel computing,
#' to speed up the estimation of the sharp bounds.
#'
#' @return a list containing for each entry in compute_T the output of
#' compute_AME_T specifying the relevant TEstim. If compute_T = "all", we use
#' compute_T = 1:Tmax and also add an extra entry for the average AME over all
#' periods. We refer to the documentation of the compute_AME_t and
#' compute_average_AME functions for a detailed description of the output of
#' each function.
#' @export
#'
# @examples
compute_AME <- function(data, estimators = env(), selectX = NULL, compute_T = "all", Option  = "quick", CIOption = "CI2", alpha = 0.05, nbCores = 4) {

  # Initialisation ----------------------------------------------------------

  # We create environments to store some variables which are computed in
  # compute_AME_t and can be used several times, to avoid computing them
  # several times
  other_numbers <- env()

  # Format the data
  format_data(data)
  Tobsd <- data$Tobsd

  # Interpret the value of compute_T to know which AMEs to compute
  compute_average <- FALSE
  if ((length(compute_T) == 1) & (compute_T[1] == "all")) {
    compute_average <- TRUE
    compute_T <- 1:(dim(data$X)[2])
  } else if (is.null(compute_T)) {
    compute_T <- 1
  }

  # Compute all requested AMEs ----------------------------------------------

  # We store them in a list with names
  requested_AME_outputs <- list()
  names_requested_AME <- c()

  for (i in 1:length(compute_T)) {
    t <- compute_T[i]
    if (t <= 0) {
      # Negative values (whether t is a scalar or a vector) count periods backwards starting from the last observation for each row
      t <- ifelse(Tobsd + t >= 1, Tobsd + t, NA)
    }
    requested_AME_outputs[[i]] <- compute_AME_t(data, t, estimators, other_numbers, selectX, Option, CIOption, alpha, nbCores)
  }
  names_requested_AME <- paste0("T_", compute_T)

  # Compute the average if we computed all periods
  if (compute_average) {
    requested_AME_outputs[[length(compute_T) + 1]] <- compute_average_AMTE(requested_AME_outputs, data, estimators, selectX, alpha, Option, CIOption)
    names_requested_AME <- c(names_requested_AME, "average")
  }

  # Name the outputs
  names(requested_AME_outputs) <- names_requested_AME

  # Give proper name to the average bounds
  if (compute_average) {
    names(requested_AME_outputs$average)[4] <- "estimatedAMEbounds"
  }

  return(requested_AME_outputs)
}
