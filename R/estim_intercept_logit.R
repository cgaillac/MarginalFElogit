#' Estimates the intercept of a logit model on the data
#'
#' @param data is an environment variable containing the relevant data:
#' - data$Y a matrix of size n x Tmax containing the values of the
#' dependent variable Y.
#' - data$V an array of size n x Tmax containing the values of X'beta
#' for the estimated value of beta.
#'
#' @return the intercept of a logit model of data$Y on data$X with the
#' beta estimated by conditional logit
#' @export
#'
# @examples
estim_intercept_logit <- function(data) {

  # Stopping parameters for the gradient descent
  tol <- 10^-4
  pas <- 1
  maxIter <- 50
  iter <- 1

  # Extract data
  Y <- data$Y
  V <- data$V

  # Initial values
  PY <- mean(Y, na.rm = TRUE)
  curAlpha <- log(PY/(1-PY)) - mean(log(V), na.rm = TRUE)

  while ((pas > tol) & (iter <= 50)) {
    prevAlpha <- curAlpha
    curAlpha <- prevAlpha + (PY - mean(exp(prevAlpha) * V / (1 + exp(prevAlpha) * V), na.rm = TRUE)) / mean(exp(prevAlpha) * V / (1 + exp(prevAlpha) * V)^2, na.rm = TRUE)
    pas <- abs(curAlpha - prevAlpha)
    iter <- iter + 1
  }

  return(curAlpha)
}
