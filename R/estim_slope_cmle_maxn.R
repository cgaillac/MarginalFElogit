#' Estimates the slope parameter on the data in a FE logit model, using CMLE
#' maximisation, along with its variance and the values of the influence function
#' at each point.
#'
#' @param data is an environment variable containing the relevant data:
#' - data$Y a matrix of size n x Tmax containing the values of the
#' dependent variable Y.
#' - data$X an array of size n x Tmax x dimX containing the values of
#' the covariates X.
#' - data$clusterIndexes a vector of size n containing the index of the
#' cluster each observation belongs to. The computed asymptotic variance is
#' clustered.
#' @param beta_init (default NULL) starting value for beta in the estimation
#' algorithm. If null, we take it to be the slope in a linear probability
#' model, divided by 4.
#'
#' @return a list containing:
#' - beta_hat: a vector of length dimX, the estimated value for the slope
#'   parameter.
#' - phi_b: a matrix of size n x dimX containing the value of the influence
#'   function at each observation (rows) w.r.t. each dimension of the
#'   covariates (columns).
#' - var_b: the estimated asymptotic covariance matrix, of size dimX x dimX,
#'   for the estimator beta_hat.
#'
#' @export
#'
# @examples
estim_slope_cmle_maxn <- function(data, beta_init = NULL) {

  # we catch and hide errors
  options(warn = -1)

  # If no initial value specified, fix it as a fourth of the linear probability
  # model estimate (linear regression)
  if (is.null(beta_init)) {
    beta_init <- 1/4 * optim(par = rep(0, dim(data$X)[3]), mean_squared_error, Y = data$Y, X = data$X)$par
  }

  ### fnscale = -1 in control ensures maximisation
  beta_hat <- optim(par = beta_init, condl_log_lik_FE_logit, Y = data$Y, X = data$X, control = list("fnscale" = -1))$par

  # Reactivate errors
  options(warn=0)

  # Now that we know beta, derive the corresponding variables (V, C_S(X, beta) ...)
  format_data(data, beta_hat)

  # Compute the influence function of beta_hat. Useful for inference on
  # Delta, at the end.
  phi_b <- infl_func_beta(beta_hat, data)

  # Get the asymptotic variance, accounting for clustering by summing on each cluster
  infl_func_cluster <- rowsum(phi_b, data$clusterIndexes)
  var_b <- t(infl_func_cluster) %*% infl_func_cluster / dim(data$X)[1]

  return(list(
    "beta_hat" = beta_hat,
    "phi_b" = phi_b,
    "var_b" = var_b
  ))

}
