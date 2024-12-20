#' Computes the average of the bounds on the AME or ATE over all periods in the whole
#' dataset. The average is computed by first taking the average across all
#' observed period for each individual, and then averaging the result over all
#' individuals.
#'
#' @param output is a list of outputs from the compute_AME_t or compute_ATE_t
#' functions. Each element of the list corresponds to the output for one period
#' at which the effect is estimated (one value of TEstim), and the average is
#' computed over all periods corresponding to elements of outputs.
#' @param data is an environment variable containing the data. If not already
#' formatted, it will be formatted by the format_data function (see its
#' documentation for information about the stored variables). We only use
#' data$clusterIndexes, a vector of size n x 1 that specifies the cluster each
#' observation pertains to.
#' @param estimators is an environment in which the results of the CMLE
#' estimation and of the non-parametric estimation of the conditional
#' distribution of S are stored. The only parameters of interest are in
#' estimators$beta_hat, a list which contains the results from CMLE
#' estimation:
#' - estimators$beta_hat$beta_hat is a vector of length dimX, the
#' estimated value for the slope parameter.
#' - estimators$beta_hat$var_b the estimated asymptotic covariance
#' matrix, of size dimX x dimX, for the estimator beta_hat.
#' @param selectX a vector containing the indices of the covariates w.r.t.
#' which the AME is being computed. If null, the AME w.r.t. all covariates
#' is being computed.
#' @param Option estimation method being used. If "quick" or "outer"
#' (case-insensitive) the outer bounds are computed. Otherwise, the sharp bounds
#' are computed.
#' @param CIOption (default "CI2") When the outer bounds method is being used,
#' specifies which confidence interval should be used. If "CI2", the CI2
#' confidence interval is being used (see DDL, section 4.2), otherwise the CI3
#' confidence interval will be used (see DDL, appendix C).
#' @param alpha (default 0.05) desired asymptotic level of the estimated
#' confidence intervals.
#' @param nbCores (default 4) number of cores to be used for parallel computing,
#' to speed up the estimation of the sharp bounds.
#'
#' @return a list containing, the same variables as for compute_AME_t, where NAs
#' are used when the intended variable is not relevant:
#'  - Option: same as input
#'  - reducedSampleSize: the number of individuals used to estimate the AME.
#'    Here this is the same as the number of individuals in the whole dataset,
#'    after excluding individuals observed at only one or zero period.
#'  - computationTime NA
#'  - estimatedAMEbounds: a matrix of size |selectX| x 2 containing the averaged
#'    out estimated (sharp or outer, depending on Option) bounds on the AME,
#'    for each covariate in selectX.
#'  - CI: a matrix of size |selectX| x 2 containing in each row the estimated
#'    confidence interval for the average AME, for the corresponding covariate
#'    in selectX.
#'  - alt_CI: a matrix of size |selectX| x 2 containing in each row the estimated
#'    confidence interval for the average AME, for the corresponding covariate
#'    in selectX, using the confidence interval option NOT requested by the
#'    user, for the outer bounds. If the sharp bounds were being computed,
#'    this output is NA.
#'  - indl_estimates: NA
#'
#' @export
#'
# @examples
compute_average_AMTE <- function(output, data, estimators, selectX, alpha, Option, CIOption, nbCores) {

  # Initialisation
  if (is.null(selectX)) {
    selectX <- 1:(dim(data$X)[3])
  }
  nbIndls <- dim(data$X)[1]
  clusterIndexes <- data$clusterIndexes
  nb_var <- length(selectX)
  beta_hat <- estimators$beta_hat$beta_hat
  cur_beta_is_not_null <- rep(NA, nb_var)

  # Prepare output
  indl_estimates <- vector("list")
  indl_estimates[["indl_bounds_on_delta"]] <- array(NA, c(nbIndls, 2, nb_var))
  indl_estimates[["indl_infl_func_delta"]] <- array(NA, c(nbIndls, 2, nb_var))
  average_bounds_on_delta <- matrix(NA, nb_var, 2)
  CI <- matrix(NA, nb_var, 2)
  CI2 <- matrix(NA, nb_var, 2)
  CI3 <- matrix(NA, nb_var, 2)

  for (k in 1:nb_var) {

    ### Average the sharp bounds
    indl_inf_bound_on_delta <- lapply(output, function(x) {x[["indl_estimates"]][["indl_bounds_on_delta"]][, 1, k]})
    indl_inf_bound_on_delta <- do.call(cbind, indl_inf_bound_on_delta)
    indl_average_delta_inf_bound_on_delta <- rowMeans(indl_inf_bound_on_delta, na.rm = TRUE)
    average_bounds_on_delta[k, 1] <- mean(indl_average_delta_inf_bound_on_delta, na.rm = TRUE)

    indl_sup_bound_on_delta <- lapply(output, function(x) {x[["indl_estimates"]][["indl_bounds_on_delta"]][, 2, k]})
    indl_sup_bound_on_delta <- do.call(cbind, indl_sup_bound_on_delta)
    indl_average_delta_sup_bound_on_delta <- rowMeans(indl_sup_bound_on_delta, na.rm = TRUE)
    average_bounds_on_delta[k, 2] <- mean(indl_average_delta_sup_bound_on_delta, na.rm = TRUE)

    # We also compute individual averages for the influence functions
    indl_infl_func_inf_bound <- lapply(output, function(x) {x[["indl_estimates"]][["indl_infl_func_delta"]][, 1, k]})
    indl_infl_func_inf_bound <- do.call(cbind, indl_infl_func_inf_bound)
    indl_average_infl_func_inf_bound <- rowMeans(indl_infl_func_inf_bound, na.rm = TRUE)

    indl_infl_func_sup_bound <- lapply(output, function(x) {x[["indl_estimates"]][["indl_infl_func_delta"]][, 1, k]})
    indl_infl_func_sup_bound <- do.call(cbind, indl_infl_func_sup_bound)
    indl_average_infl_func_sup_bound <- rowMeans(indl_infl_func_sup_bound, na.rm = TRUE)

    indl_average_infl_func <- cbind(indl_average_infl_func_inf_bound, indl_average_infl_func_sup_bound)

    # Now we know the influence functions, we can compute the clustered standard deviation
    # For the quick method, std_bounds is the same for both bounds and in fact is the standard
    # deviation on Delta_hat (the mean of the interval)
    clustered_infl_func <- rowsum(indl_average_infl_func, clusterIndexes)
    std_bounds <- apply(clustered_infl_func, 2, sd, na.rm = TRUE) / sqrt(nbIndls)

    if (Option == "quick") { # Average for outer bounds
      # We extract the individual values of the AME and maximal bias for each period
      # and each individual, which we turn into a matrix. We then average over all
      # observed periods for each individual, and then average these values over all
      # individuals. We do this for each variable w.r.t. which we compute AMEs
      average_delta_hat <- mean(average_bounds_on_delta[k,])
      average_bias_sup <- (average_bounds_on_delta[k, 2] - average_bounds_on_delta[k, 1]) / 2

      # We can deduce the relevant confidence interval
      bCI3 <- (abs(beta_hat[selectX[k]]) + qnorm(1 - 4 * alpha / 5) * estimators$beta_hat$var_b[selectX[k], selectX[k]] / nbIndls) * average_bias_sup
      length_CI2 <- 2 * std_bounds[1] * sqrt(qchisq(1 - alpha, df = 1, ncp = (average_bias_sup / std_bounds[1])^2))
      length_CI3 <- 2 * std_bounds[1] * sqrt(qchisq(1 - alpha / 5, df = 1, ncp = (bCI3 / std_bounds[1])^2))

      CI2[k,] <- cbind(average_delta_hat - length_CI2 / 2, average_delta_hat + length_CI2 / 2)
      CI3[k,] <- cbind(average_delta_hat - length_CI3 / 2, average_delta_hat + length_CI3 / 2)

    } else { # Average for sharp bounds
      # Estimate the confidence interval
      IM_quantile <- compute_IM_quantile(alpha, average_bounds_on_delta[k, ], std_bounds)

      ### This is a t-test to see if the beta corresponding to the current AME
      ### is null
      cur_beta_is_not_null[k] <- sqrt(nbIndls * beta_hat[selectX[k]]^2 / estimators$beta_hat$var_b[selectX[k], selectX[k]]) > qnorm(1 - alpha/2)

      CI[k, ] <- average_bounds_on_delta[k, ] + IM_quantile * c(-std_bounds[1], std_bounds[2])
      if (!cur_beta_is_not_null[k]) {
        CI[k, 1] <- min(CI[k, 1], 0)
        CI[k, 2] <- max(CI[k, 2], 0)
      }
    }
  }

  out <- vector("list")
  out[[1]] <- Option
  out[[2]] <- nbIndls
  out[[3]] <- NA
  out[[4]] <- average_bounds_on_delta
  # We return the two confidence intervals for the outer bounds
  if (Option == "quick") {
    if (CIOption == "CI2") {
      out[[5]] <- CI2
      out[[6]] <- CI3
    }
    else {
      out[[5]] <- CI3
      out[[6]] <- CI2
    }
  } else {
    out[[5]] <- CI
    out[[6]] <- NA
  }
  out[[7]] <- NA

  names(out) <- c("Option", "reducedSampleSize", "computationTime", "estimatedbounds", "CI", "alt_CI", "indl_estimates")

  return(out)
}
