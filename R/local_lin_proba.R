#' Estimates probabilities that S = 0, ..., Tobsd using local linear regression.
#' We use a Gaussian kernel, and then constrain probabilities to be positive and
#' to sum up to 1.The kernel estimator of the density is also given in output.
#'
#' @param data is an environment variable containing the relevant data:
#' - data$S a vector of size n counting for each individual the number of periods for
#' which Y = 1.
#' - data$X an array of size n x Tmax x dimX containing the values of the covariates X.
#' @param h (default 1) the bandwidth for the local linear regressions. We increase it at points
#' where the regression is degenerate otherwise. It should have size Tmax + 1, where the j-th value
#' is the bandwidth used for estimating the value of the P(S = j - 1 | X)'s. If h is too short, the
#' last value is repeated for all subsequent j's.
#'
#' @return a list containing:
#'  - condlProbas: a matrix of size n x (Tmax + 1) containing, in position (i, j), the
#'    estimate for P(S = j - 1 | X) at the i-th observation.
#'  - densityEstimate: a matrix of size n x (Tmax + 1) containing, at each row (individual),
#'    the estimated density for having covariates (X_1, ..., X_T). Each column represents the
#'    value found using the corresponding bandwidths from h.
#'
#' @export
#'
# @examples
local_lin_proba <- function (data, h = 1) {

  # Extract data
  S <- data$S
  n <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  dimX <- dim(data$X)[3]

  # Initialisation
  densityEstimate <- matrix(NA, n, Tmax + 1)
  condlProbas <- matrix(NA, n, Tmax + 1)

  # If only one bandwidth is given, use the last bandwidth for all remaining values of S
  if (length(h) < Tmax + 1) {
    h <- c(h, rep(h[length(h)], Tmax + 1 - length(h)))
  }

  # We regress each element of Y on all of X_1, ..., X_T, so we flatten X
  X <- list()
  for (k in 1:dimX) {
    X[[k + 1]] <- data$X[,, k]
  }
  X <- do.call(cbind, X) # Binds the different dimensions of X on one row for
  # each observation

  # Perform a linear regression for each point and each S == t
  print("Non-parametric estimation of the conditional distribution of S:")
  for (i in 1:n) {

    # Uncomment to track progress if you want
    if (i %% max(floor(n/20), 1) == 0) {
      print(paste0(floor(100 * i / n), "% of the local linear regressions performed"))
    }

    for (t in 0:Tmax) {

      # Binary variable to see whether we need to increase the bandwidth
      isDegenerate <- TRUE # TRUE to start with so that we enter the loop
      cur_h <- h[t + 1]
      Y <- as.numeric(S == t)

      while (isDegenerate) {

        # Compute weight of each observation for the current linear regression
        # If some rows are missing value, we use only the observed values and
        # "rescale" to not penalise observations that are not missing values
        x <- X[i, ]
        dx <- repmat(x, n, 1) - X[, ]

        # Ignore observations which have NAs whenever the i-th observation doesn't
        Y <- Y[rowSums(!is.na(dx)) > 0]
        dx <- dx[rowSums(!is.na(dx)) > 0, ]

        w <- exp(Tmax * dimX * rowMeans(-(dx/cur_h)^2/2, na.rm = TRUE))
        recentred_X <- cbind(matrix(rep(1, dim(dx)[1]), ncol = 1), dx)
        recentred_X <- replace(recentred_X, is.na(recentred_X), 0) # For observations other than i, NAs are as good as 0 (ignored in the regression)
        recentred_X <- recentred_X[, c(TRUE, !is.na(x))] # But if the i-th observation has NAs, we should ignore these columns as otherwise the matrix will be degenerate

        # Compute the denominator, X'WX
        denom <- t(recentred_X) %*% (recentred_X * repmat(matrix(w, ncol = 1), 1, dim(recentred_X)[2]))

        if (rcond(denom) < 10^-15) { # If matrix is degenerate, increase bandwidth and start over
          cur_h <- cur_h * 1.5
        } else { # Otherwise we can finish the regression and go out the while loop
          isDegenerate <- FALSE

          # Local density estimate
          densityEstimate[i, t + 1] <- sum(w) / n / (sqrt(2 * pi) * cur_h)^(dimX * Tmax + 1)

          # Local conditional mean estimate
          numer <- t(recentred_X) %*% matrix(Y * w, ncol = 1)
          condlProbas[i, t + 1] <- solve(denom, numer)[1, ]
        }
      }
    }
  }

  # Normalise conditional probabilities to be positive and sum to 1
  condlProbas <- pmax(pmin(condlProbas, 1), 0)
  condlProbas <- condlProbas / repmat(matrix(rowSums(condlProbas, na.rm = TRUE), ncol = 1), 1, Tmax + 1)

  # If a row was all 0s (so is now NA because we divided by 0s), we take the
  # conditional probabilities from the closest X
  undefinedProbas <- rowSums(condlProbas, na.rm = TRUE) == 0
  indexUndefinedProbas <- which(undefinedProbas)

  for (i in indexUndefinedProbas) {
    curObsn <- X[i,]
    distanceFromCurObsn <- rowSums((X[!undefinedProbas,] - repmat(curObsn, sum(!undefinedProbas), 1))^2)
    closestObsn <- which.min(distanceFromCurObsn)
    condlProbas[i,] <- (condlProbas[!undefinedProbas,])[closestObsn,]
  }

  out <- vector("list")
  out[["densityEstimate"]] <- densityEstimate
  out[["condlProbas"]] <- condlProbas
  return(out)
}
