#' This function computes what corresponds to the function under(over)lined
#' delta in Section 3.2.2 in DDL, that is, estimates the sharp bounds on the
#' AME for each individual observation, using parallel computing. It can also
#' be used to compute the sharp bounds on the ATE, accounting for the Y(2X - 1)
#' term
#'
#' @param data is an environment variable containing the relevant data,
#' formatted by format_data (see function documentation). In particular,
#' it contains:
#'  - data$Y a matrix of size n x Tmax containing the values of the
#' dependent variable Y.
#'  - data$X an array of size n x Tmax x dimX containing the values of
#' the covariates X.
#' @param estimators is an environment variable containing some estimates from
#' the data. It contains:
#'  - estimators$beta_hat, a list which contains the results from CMLE
#' estimation:
#'  - estimators$beta_hat$beta_hat a vector of length dimX, the
#' estimated value for the slope parameter.
#'  - estimators_beta_hat$phi_b a matrix of size n x dimX containing
#' the value of the influence function at each observation (rows) w.r.t.
#' each dimension of the covariates (columns).
#'  - estimators$beta_hat$var_b the estimated asymptotic covariance
#' matrix, of size dimX x dimX, for the estimator beta_hat.
#'  - estimators$h_local_lin_proba a vector of length (Tmax + 1) containing,
#' in j-th position, the bandwidth used to estimate the P(S = j - 1 | X)'s.
#'  - estimators$condlProbas: a matrix of size n x (Tmax + 1) containing,
#' in position (i, j), the estimate for P(S = j - 1 | X) at the i-th
#' observation.
#'  - estimators$densityEstimate a matrix of size n x (Tmax + 1) containing,
#' at each row (individual), the estimated density for having covariates
#' (X_1, ..., X_T). Each column represents the value found using the
#' corresponding bandwidths from h.
#' @param other_numbers is an environment variable containing the variable
#' other_numbers$comb_numbers, a matrix of size (Tmax + 1) x (Tmax + 1)
#' containg in position (i, j) the number (i choose j). If j > i the value
#' is NA.
#' @param Vt a vector of length n containing, for each individual, the value of
#' v(X, beta). If computing the AME, v(X, beta) = X_t'beta where t is the period
#' at which the AME is being computed. If computing the ATE, v(X, beta) is as above
#' but the value of X_tk (where k is the variable relative to which we compute the
#' ATE) must be switched (replaced by 1 if 0, and by 1 if 0).
#' @param coeffsOmega a matrix of size n x (Tmax + 2) where each row represents
#' the coefficient of the polynomial Omega(., x, beta). Coefficients start with
#' the constant coefficient and there are (Tobsd[i] + 2) coefficients on each row
#' where Tobsd[i] is the number of periods observed for the relevant
#' individuals. Subsequent entries on the rows are 0s.
#' @param topCoeffOmega a vector of length n containing the top coefficient
#' for the original Omega polynomial for each individual.
#' @param nbCores the number of cores to be used for parallel computing.
#' @param firstTermATE (default NULL) if computing the ATE, a reminder of
#' the vector of length n declaring for each individual the value of
#' Y * (2X - 1) at the period at which the ATE is computed. If NULL,
#' nothing is added, which corresponds to what is needed for the AME.
#'
#' @return returns a list containing:
#' - average_bounds_on_delta a vector of length 2 containing the average
#'   of the estimated sharp lower and upper bounds on the AME.
#' - indl_bounds_on_delta a matrix of size n x 2 containing, on each row,
#'   the estimated sharp lower and upper bounds on the AME for the relevant
#'   individual.
#' - c_mat a matrix of size n x (Tmax + 1) containing, on each row, in the j_th
#'   column the estimated value of c_(j-1)(X) for the relevant individual. When
#'   j - 1 is larger than the number of observed period, the value is NA.
#'
#' @export
#'
# @examples

bound_delta <- function(data, estimators, other_numbers, Vt, coeffsOmega, topCoeffOmega, nbCores, firstTermATE = NULL) {

  # Initialisation
  n <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  dimX <- dim(data$X)[3]
  Tobsd <- data$Tobsd
  C_S_tab <- data$C_S_tab
  comb_numbers <- other_numbers$comb_numbers
  PSt_X <- estimators$condlProbas
  densityEstimate <- estimators$densityEstimate

  # Compute the Z_t(X, s, beta)'s
  # We use an array. First dimension is the number i of the observation, second
  # is t (starting at 0) and last is the value of S we plug in (starting at 0
  # and possibly different from the one observed)
  Zt_mat <- array(NA, c(n, Tmax + 1, Tmax + 1))

  Tobsd_minus_t <- repmat(matrix(Tobsd, ncol = 1), 1, Tmax + 1) - repmat(0:Tmax, n, 1)
  Tobsd_minus_t[Tobsd_minus_t < 0] <- NA

  for (j in 1:(Tmax + 1)) {
    j_minus_t <- j - 1 - repmat(0:Tmax, n, 1)
    j_minus_t[j_minus_t < 0] <- NA

    T_minus_t_choose_j_minus_t <-
      extractIndlPeriod(
        comb_numbers,
        vecIndls = Tobsd_minus_t + 1,
        vecPeriods = j_minus_t + 1
      )

    Zt_mat[,, j] <- T_minus_t_choose_j_minus_t * repmat(matrix(Vt, ncol = 1), 1, Tmax + 1)^(j-1) / C_S_tab[, rep(j, Tmax + 1)]
  }

  Zt_mat[is.na(Zt_mat)] <- 0 # As per the convention in the paper, Zt is null where it is "undefined"

  # Extract the Zt's for the values of S we observe
  Zt_observed <- matrix(NA, n, Tmax + 1)
  for (t in 1:(Tmax + 1)) {
    Zt_observed[, t] <- extractIndlPeriod(Zt_mat[, t,], data$S + 1)
  }

  # Compute the c_t's. Make it NAs for the undefined value (t larger than the
  # number of period observed for the individual)
  c_mat <- matrix(NA, n, Tmax + 1)
  for (t in 0:Tmax) {
    c_mat[, t + 1] <- rowSums(Zt_mat[,t + 1 ,] * PSt_X, na.rm = TRUE)
    c_mat[Tobsd < t, t + 1] <- NA
  }

  # Compute the m_t's. After taking c_t / c_0 we normalise to ensure
  # Ê(m_t Z_0) = Ê(Z_t) which increases the algorithm's robustness
  # The new m_t is the old E(Z_t) m_t / E(m_t Z_0)
  m_hat <- c_mat[, 2:(Tmax + 1)] / repmat(matrix(c_mat[, 1], ncol = 1), 1, Tmax)
  m_hat <-
    m_hat * repmat(colMeans(Zt_observed[, 2:(Tmax + 1)], na.rm = TRUE), n, 1) /
    repmat(colMeans(m_hat * repmat(matrix(Zt_observed[, 1], ncol = 1), 1, Tmax), na.rm = TRUE), n, 1)

  # Compute the (known) first T + 1 terms of the bounds on delta, using Lemma 1 from the paper
  # known_first_terms <- rowSums(c_mat * coeffsOmega[, 1:(Tmax+1)])
  known_first_terms <- rowSums(Zt_observed * coeffsOmega[, 1:(Tmax + 1)])
  if (!is.null(firstTermATE)) {
    known_first_terms <- known_first_terms + firstTermATE
  }

  # Set up parallel computing to speed up computations
  suppressMessages(sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK"))
  suppressMessages(invisible(capture.output(sfLibrary(R.matlab))))
  # suppressMessages(invisible(capture.output(sfLibrary(pracma))))
  sfLibrary(pracma)
  sfExport("PSt_X", "densityEstimate", "Zt_mat", "n", "Tobsd", "dimX", "estimators", "m_hat", "c_mat", "topCoeffOmega")
  sfExport("bound_next_moment", "build_hankel_matrices", "C_S_fun", "choose",
    "compute_asymptotic_variance_moment_vector", "extractIndlPeriod", "truncate_moment_vector_until_valid")

  # Compute individual bounds
  # First, for the given observation, we compute the asymptotic variance of the
  # (conditional) moment vector, then we plug it in the function to bound the
  # next moment. It's used to check determinants of the Hankel matrices are
  # positive in a statistically significant way.
  bounds_on_the_second_term <- sfLapply(1:n, function(i) {
    asymptotic_variance_m <- compute_asymptotic_variance_moment_vector(PSt_X[i, ], densityEstimate[i], Zt_mat[i,,], n, dimX, estimators$h_local_lin_proba, m_hat[i,], c_mat[i, 1])
    indl_bounds_on_next_moment <- bound_next_moment(m_hat[i,], asymptotic_variance_m[1:Tobsd[i], 1:Tobsd[i]], thresholdDeterminant = sqrt(2 * log(log(n))))
    indl_bounds_on_the_second_term <- c_mat[i, 1] * topCoeffOmega[i] * indl_bounds_on_next_moment
    return(sort(indl_bounds_on_the_second_term)) # We need to sort in case the top coefficient is negative
  })
  bounds_on_the_second_term <- do.call(rbind, bounds_on_the_second_term) # Make matrix

  # Close parallel computing
  suppressMessages(sfStop())

  # Add the known first terms to the bounds on the second term
  indl_bounds_on_delta <- repmat(matrix(known_first_terms, ncol = 1), 1, 2) + bounds_on_the_second_term

  # Average out. We'll still need to multiply by beta to get the proper AME
  average_bounds_on_delta <- colMeans(indl_bounds_on_delta)

  return(list(
    "average_bounds_on_delta" = average_bounds_on_delta,
    "indl_bounds_on_delta" = indl_bounds_on_delta,
    "c_mat" = c_mat
  ))
}
