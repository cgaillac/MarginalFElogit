#' This function computes the asymptotic variance on the estimated moment
#' for a given individual, that results from the (non-parametric) estimation
#' of the P(S = j | X)'s.
#'
#' @param PSt_X a vector of length Tmax + 1 containing in j-th position the
#' estimate for P(S = j - 1 | X) for the observation at hand.
#' @param Xdensity the estimated density for the covariates (X_1, ..., X_T) for
#' the observation at hand
#' @param Zt_mat a matrix of size (Tmax + 1) x (Tmax + 1) containg in position
#' (i, j) the value of Z_i(X, j, beta) for the observation at hand. Thus,
#' c = Zt_mat * PSt_X and m = (c_1, ..., c_T) / c_0, from which we deduce the
#' asymptotic variance.
#' @param nbIndls the total number of individuals used to estimate the
#' P(S = j | X)'s.
#' @param dimX the total number of covariates at a given period.
#' @param h_local_lin_proba vector of length Tmax + 1 where the j-th element is
#' the bandwidth used to estimate P(S = j - 1 | X).
#' @param m estimated moment vector (length Tmax) for the observation at hand.
#' Starts with the 1st moment.
#' @param c0 estimated value of c_0.
#'
#' @return a matrix of size Tmax x Tmax, an empirical estimate for the value
#' of the asymptotic variance of the vector of moments m.
#'
#' @export
#'
# @examples
compute_asymptotic_variance_moment_vector <- function(PSt_X, Xdensity, Zt_mat, nbIndls, dimX, h_local_lin_proba, m, c0) {

  # Int K²(u) du from choose_bandwidth_and_estim_condl_proba
  Int_K2 <- 1/sqrt(pi)
  Tmax <- length(PSt_X) - 1

  # Asymptotic (conditional) variance of the gamma_0t = P(S = t | X), see online
  # appendix B from the DDL paper
  asymptotic_variance_PSt_X <-
    (diag(PSt_X) - repmat(PSt_X, Tmax + 1, 1) * repmat(matrix(PSt_X, ncol = 1), 1, Tmax + 1)) *
    (Int_K2)^(dimX * Tmax) / nbIndls / Xdensity /
    repmat(h_local_lin_proba^(dimX * Tmax / 2), Tmax + 1, 1) * repmat(matrix(h_local_lin_proba, ncol = 1)^(dimX * Tmax / 2), 1, Tmax + 1) # We use a different bandwith for each t in PSt_X


  # Deduce the asymptotic (conditional) variance of c = Zt_mat * PSt_X
  # (matrix multiplication)
  asymptotic_variance_c <- Zt_mat %*% asymptotic_variance_PSt_X %*% t(Zt_mat)

  # Deduce the asymptotic (conditional) variance of m = C_1:T / C_0. Because c (and
  # hence m) asymptotically converges a.s. to a point, the asymptotic variance
  # is asymptotically that of g'(m) c
  deriv_m_d_c <- cbind(-m, eye(Tmax)) / c0
  asymptotic_variance_m <- deriv_m_d_c %*% asymptotic_variance_c %*% t(deriv_m_d_c)

  return(asymptotic_variance_m)
}
