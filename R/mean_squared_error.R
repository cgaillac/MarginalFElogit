#' Evaluates the mean square error loss function of a linear regression.
#'
#' @param beta a vector of size dimX, the value at which the loss is computed
#' (i.e. we predict Y using X'beta)
#' @param X an array of size n x Tmax x dimX containing the values of the covariates X
#' @param Y a matrix of size n x Tmax containing the values of the dependent variable Y
#'
#' @return the value of the mean squared error loss function of the linear regression at beta
#' @export
#'
# @examples
mean_squared_error <- function(beta, X, Y) {

  # Force X as an array even if dimX = 1
  if (length(dim(X)) == 2) {
    X <- array(X, c(dim(X), 1))
  }
  k <- dim(X)[3]

  # Compute X' * beta for each individual-period pair
  index <- matrix(0, dim(X)[1], dim(X)[2])
  for (i in 1:k) {
    index <- X[,, i] * beta[i]
  }

  # Compute MSE on each row
  mse <- sum((Y - index)^2, na.rm = TRUE)

  return(mse)
}
