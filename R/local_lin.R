#' Estimates E(Y|X) at the observed X by local linear regression.
#' We use a Gaussian kernel, and also output the kernel estimator of the density of X.
#'
#' @param data is an environment variable containing the relevant data:
#' - data$Y a matrix of size n x Tmax containing the values of the dependent variable Y.
#' - data$X an array of size n x Tmax x dimX containing the values of the covariates X.
#' @param h (default 1) the bandwidth for the local linear regressions. We increase it at points
#' where the regression is degenerate otherwise.
#'
#' @return a list containing.
#'  - condlMeans: a vector of length n containing the
#'    estimates for P(S = j - 1 | X) at each observation.
#'  - densityEstimate: a vector of length n containing, at each row (individual),
#'    the estimated density for having covariates X.
#'
#' @export
#'
# @examples
local_lin <- function (data, h = 1) {

  # Extract data
  Y <- data$Y
  n <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  dimX <- dim(data$X)[3]

  # Initialisation
  densityEstimate <- rep(NA, n)
  condlMeans <- matrix(NA, n, dim(Y)[2])

  # We regress each element of Y on all of X_1, ..., X_T, so we flatten X
  X <- list()
  for (k in 1:dimX) {
    X[[k + 1]] <- data$X[,, k]
  }
  X <- do.call(cbind, X) # Binds the different dimensions of X on one row for
  # each observation

  # Perform a linear regression for each point
  for (i in 1:n) {

    # Binary variable to see whether we need to increase the bandwidth
    isDegenerate <- TRUE # True to start with so that we enter the loop
    cur_h <- h

    while (isDegenerate) {

      # Compute weight of each observation for the current linear regression
      # If some rows are missing value, we use only the observed values and
      # "rescale" to not penalise observations that are not missing values
      x <- X[i, ]
      dx <- repmat(x, n, 1) - X[, ]
      w <- exp(Tmax * dimX * rowMeans(-(dx/cur_h)^2/2, na.rm = TRUE))
      recentred_X <- cbind(matrix(rep(1, n), ncol = 1), dx)

      # Compute the denominator, X'WX
      denom <- t(recentred_X) %*% (recentred_X * repmat(matrix(w, ncol = 1), 1, dimX * Tmax + 1))

      if (rcond(denom) < 10^-15) { # If matrix is degenerate, increase bandwidth and start over
        cur_h <- cur_h * 1.5
      } else { # Otherwise we can finish the regression and go out the while loop
        isDegenerate <- FALSE

        # Local density estimate
        densityEstimate[i] <- sum(w) / n / (sqrt(2 * pi) * cur_h)^(dimX * Tmax + 1)

        # Local conditional mean estimate
        numer <- t(recentred_X) %*% (Y * repmat(matrix(w, ncol = 1), 1, dim(Y)[2]))
        condlMeans[i, ] <- solve(denom, numer)[1, ]
      }
    }
  }


  out <- vector("list")
  out[["densityEstimate"]] <- densityEstimate
  out[["condlMeans"]] <- condlMeans
  return(out)
}
