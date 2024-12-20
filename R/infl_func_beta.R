#' Computes the influence function of beta_hat:
#'     E(d^2ln L/dbeta^2)^{-1}  x (dln L/dbeta)(Yi|Xi,beta)).
#'
#' @param beta the parameter at which the influence function is computed.
#' @param data is an environment variable containing the relevant data:
#' - data$Y a matrix of size n x Tmax containing the values of the dependent variable Y.
#' - data$X an array of size n x Tmax x dimX containing the values of the covariates X.
#' - data$Call a vector of size n containing the indexes of the hypothetical clusters for the individuals
#'
#' @return a vector of size n x dimX containing the value of the influence
#' function of beta w.r.t. each dimension of the covariate (columns), for
#' each observation (rows).
#'
# @examples
infl_func_beta <- function(beta, data) {

  # Initialisation ----------------------------------------------------------

  # Extract relevant variables from the data
  C_S_vec <- extractIndlPeriod(data$C_S_tab, data$S + 1)
  C_S_minus_one_one_period_out_array <- data$C_S_minus_one_one_period_out_array

  # Derive number of observations, periods ...
  n <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  dimX <- dim(data$X)[3]


  # Compute the influence function ------------------------------------------

  # Compute the score (dln L/dbeta)(Yi|Xi,beta)
  score <- matrix(0, n, dimX)
  for (t in 1:Tmax) {
    summand <- data$X[, t,] * repmat(matrix(data$Y[, t] - data$V[, t] * ifelse(data$S == 0, 0, extractIndlPeriod(C_S_minus_one_one_period_out_array[,, t], pmax(data$S, 1))) / C_S_vec, ncol = 1), 1, dimX)
    summand[is.na(summand)] <- 0
    score <- score + summand
  }

  # Get the Fisher information matrix
  fisherInfoMat <- t(score) %*% score / n

  # The influence function is the inverse of the Fisher information matrix times
  # the score
  infl_func <- t(solve(fisherInfoMat, t(score)))

  return(infl_func)
}
