#' Evaluates the conditional log likelihood function of the panel data logit
#' at beta
#'
#' If Y or X contain NA for some individuals at some periods, these
#' individual-period pairs are discarded.
#'
#' @param beta the value at which the loss is computed
#' @param Y a matrix of size n x Tmax containing the values of the dependent
#' variable Y
#' @param X an array of size n x Tmax x dimX containing the values of the
#' covariates Xs
#'
#' @return the value of the conditional log likelihood function at beta
#' @export
#'
# @examples
condl_log_lik_FE_logit <- function(beta, Y, X) {

  # Force X as an array even if dimX = 1
  if (length(dim(X)) == 2) {
    X <- array(X, c(dim(X), 1))
  }
  k <- dim(X)[3]

  # Compute X' * beta for each individual-period pair
  index <- matrix(0, dim(X)[1], dim(X)[2])
  for (i in 1:k) {
    index <- index + X[,, i] * beta[i];
  }

  # Compute C_S(X, beta) for each individual
  C_S_tab <- C_S_fun(exp(index))
  S <- rowSums(Y, na.rm = TRUE)
  C_S_vec <- extractIndlPeriod(C_S_tab, S + 1)

  # Use the formula for the conditional log-likelihood
  condLogLikelihood <- sum(Y * index, na.rm = TRUE) - sum(log(C_S_vec))

  return(condLogLikelihood)
}
