#' Format the data in the appropriate way for the algorithm
#'
#' @param data is an environment containing the observed data:
#' - data$Y is a matrix of size n x Tmax containing the values of the dependent variable Y,
#' - data$X is an array of size n x Tmax x dimX containing the values of the covariates X
#' -data$clusterIndexes is a vector of size n x 1 that specifies the cluster each
#' observation pertains to. If it does not exist, the function enforces the default
#' setting of i.i.d. observations - the parameter takes value 1:n so that each
#' observation is its own cluster.
#' @param beta (default NULL) the slope parameter for the FE logit model. If specified, computes
#' a few extra variables from the data.
#'
#' @return The algorithm does not return anything but data has been modified so that:
#' - data$Y is a matrix of size n x Tmax containing the values of the dependent variable Y,
#' - data$X is an array of size n x Tmax x dimX containing the values of the covariates X
#' - data$S is a vector of size n x 1 that counts, for each row (individual), the number
#'   of columns (periods) for which Y is 1.
#' - data$Tobsd is a vector of size n x 1 that counts, for each row (individual), the
#'   number of columns (periods) that are not NAs (observed).
#' - data$clusterIndexes is a vector of size n x 1 that specifies the cluster each
#'   observation pertains to. If it does not exist, the function enforces the default
#'   setting of i.i.d. observations - the parameter takes value 1:n so that each
#'   observation is its own cluster.
#' - data$V is a matrix of size n x Tmax containing, for each individual-period pair, the
#'   value of X'beta. It is NA for unobserved individual-period pairs. This parameter
#'   is undefined if beta was NULL (or not specified as input).
#' - data$C_S_tab is a matrix of size n x (Tmax + 1). Each row (individual) has, in the j_th
#'   column, the value of C_(j-1)(X, beta). It is NA when j - 1 is larger than Tobsd. This
#'   parameter is undefined if beta was NULL (or not specified as input).
#' - data$C_S_minus_one_one_period_out_array is a matrix of size n x Tmax x Tmax.
#'   Each row (individual) has, in its t-th third-dimensional slice and j-th column (period),
#'   the value of C_(j-1)(Z, beta) where Z is the same as the relevant row of X deprived
#'   from its t-th period.
#'
# @examples
format_data <- function(data, beta = NULL) {


  # Reshape X and Y ---------------------------------------------------------
  if (length(dim(data$X)) == 2) { # Make X an array even if dimX = 1
    data$X <- array(data$X, c(dim(data$X), 1))
  }

  # Compute S, Tobsd create a unique cluster if none are given --------------
  if (!env_has(data, "S")) { # Computes S for each individual
    env_poke(data, "S", rowSums(data$Y, na.rm = TRUE))
  }

  if (!env_has(data, "Tobsd")) {
  # Computes the number of observed periods for each individual
    env_poke(data, "Tobsd", rowSums(!is.na(data$Y)))
  }

  if (!env_has(data, "clusterIndexes")) {
    # If no clusters are give, consider each data point as its own cluster (i.i.d. observations) ?
    env_poke(data, "clusterIndexes", 1:(dim(data$X)[1]))
  }

  # If beta is specified, compute V and C_S ---------------------------------
  if (!is.null(beta)) {
    if (!env_has(data, "V")) { # Computes exp(X't * beta_hat) for each individual-period pair
      logV <- matrix(0, dim(data$X)[1], dim(data$X)[2])
      for (k in 1:dim(data$X)[3]) {
        logV <- logV + data$X[,, k] * beta[k]
      }
      env_poke(data, "V", exp(logV))
    }

    if (!env_has(data, "C_S_tab")) { # Derives C_t(X, beta) for all individuals and all t
      env_poke(data, "C_S_tab", C_S_fun(data$V))
    }

    # We need to know the value taken by C_t(X, beta) when we exclude the
    # observation at period t' = 1, ..., Tmax and taking t = S - 1
    if(!env_has(data, "C_S_minus_one_one_period_out_array")) {
      C_S_minus_one_one_period_out_array <- array(NA, c(dim(data$Y)[1], dim(data$Y)[2], dim(data$Y)[2]))
      for (t in 1:(dim(data$X)[2])) {
        C_S_minus_one_one_period_out_array[,, t] <- C_S_fun(matrix(data$V[, -t], ncol = dim(data$V)[2] - 1))
        # We cannot exclude an observation which was not observed
        C_S_minus_one_one_period_out_array[,, t] <- ifelse(is.na(repmat(matrix(data$V[, t], ncol = 1), 1, dim(data$Y)[2])), NA, C_S_minus_one_one_period_out_array[,, t])
      }
      env_poke(data, "C_S_minus_one_one_period_out_array", C_S_minus_one_one_period_out_array)
    }
  }

  return()
}
