#' This function chooses the bandwidths to be used for the estimation through
#' local linear regressions of the conditional distribution of S, as described
#' in online appendix B of the DDL paper. It then proceeds to estimate these
#' conditional probabilities, and saves them in estimators.
#'
#' @param data is an environment variable containing the relevant data,
#' formatted by format_data (see function documentation):
#'  - data$Y a matrix of size n x Tmax containing the values of the
#' dependent variable Y.
#'  - data$X an array of size n x Tmax x dimX containing the values of
#' the covariates X.
#'  - data$clusterIndexes a vector of size n containing the index of the
#' cluster each observation belongs to. The computed asymptotic variance is
#' clustered.
#' @param estimators is an environment variable containing the results from
#' logit and FE logit estimation:
#'    -  estimators$beta_hat$beta_hat is the CMLE estimate of the FE logit
#' slope parameter.
#'    -  estimators$alphaFElogit is the estimate of the constant parameter in
#' a standard logit model using the estimated CMLE slope parameter. If empty,
#' it is estimated.
#'
#' @return The function does not return anything, but the estimators environment
#' has the following parameters added:
#' - estimators$h_local_lin_proba a vector of length (Tmax + 1) containing, in
#'   j-th position, the bandwidth used to estimate the P(S = j - 1 | X)'s.
#' - estimators$condlProbas: a matrix of size n x (Tmax + 1) containing,
#'   in position (i, j), the estimate for P(S = j - 1 | X) at the i-th
#'   observation.
#' - estimators$densityEstimate: a matrix of size n x (Tmax + 1) containing,
#'   at each row (individual), the estimated density for having covariates
#'   (X_1, ..., X_T). Each column represents the value found using the
#'   corresponding bandwidths from h.
#'
#' @export
#'
# @examples
choose_bandwidth_and_estim_condl_proba <- function(data, estimators) {

  # Initialisation ----------------------------------------------------------

  # Data shape
  nbIndls <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  dimX <- dim(data$X)[3]

  # Extract data
  X <- data$X
  V <- data$V
  C_S_tab <- data$C_S_tab
  C_S_minus_one_one_period_out_array <- data$C_S_minus_one_one_period_out_array
  beta_hat <- estimators$beta_hat$beta_hat

  # Parameters used to choose the bandwidth of the local linear regression in the sharp method
  # We use K(u) = exp(-u²/2) / sqrt(2 * pi) (gaussian kernel)
  Rn_choice_bandwidth <- log(nbIndls)
  Int_K2 <- 1/sqrt(pi) # Int K²(u) du, also used in compute_asymptotic_variance_moment_vector

  # Choose the bandwidths for the local regressions -------------------------

  # We will need uniform draws on Supp(X), which we take to be given,
  # along each dimension and for each period, by the maximum and minimum
  # values taken in the sample
  scaling <- matrix(1, Tmax, dimX)
  for (k in 1:dimX) {
    scaling[, k] <- apply(X[,, k], 2, function(lambda) {max(lambda, na.rm = TRUE) - min(lambda, na.rm = TRUE)})
  }
  sizeSuppX <- prod(c(scaling))

  ### Uniform draws on Supp(X)
  nbUniformDraws <- 2000
  uniformDrawsX <- array(runif(nbUniformDraws * Tmax * dimX), c(nbUniformDraws, Tmax, dimX))
  for (k in 1:dimX) {
    uniformDrawsX[,, k] <- apply(X[,, k], 2, function(x) {min(x, na.rm = TRUE)}) + uniformDrawsX[,, k] * repmat(scaling[, k], nbUniformDraws, 1)
  }

  ### Estimate the constant alpha in a logit model using the original data
  ### We will use that alpha to compute the conditional distribution of S
  ### for the uniformly drawn X
  if (!env_has(estimators, "alphaFELogit")) {
    env_poke(estimators, "alphaFELogit", estim_intercept_logit(data))
  }
  alphaFELogit <- estimators$alphaFELogit

  ### Estimate V = X'beta for the uniformly drawn X
  VuniformDrawsX <- matrix(0, nbUniformDraws, Tmax)
  for(k in 1:dimX){
    VuniformDrawsX <- VuniformDrawsX + uniformDrawsX[,, k] * beta_hat[k]
  }
  VuniformDrawsX <- exp(VuniformDrawsX)

  ### Deduce the distribution of S conditional on the uniformly drawn X
  C_S_tab_uniform_draws_X <- C_S_fun(VuniformDrawsX)
  PSt_uniform_draws_X <-
    C_S_tab_uniform_draws_X *
    repmat(exp(alphaFELogit)^(0:Tmax), nbUniformDraws, 1) /
    repmat(matrix(apply(1 + exp(alphaFELogit) * VuniformDrawsX , 1, prod), ncol = 1), 1, Tmax + 1)

  ### Apply the formula from online appendix B of the DDL paper to get sigma_t(h)^2 * h^(Tmax * dimX)
  std_squared_approx_St <- sizeSuppX / nbIndls * Int_K2^(Tmax * dimX) * colMeans(PSt_uniform_draws_X * (1 - PSt_uniform_draws_X))

  ### Compute the distribution of S conditional on the observed X, for constant alpha
  PSt_obsd_X_constant_alpha <-
    C_S_tab *
    repmat(exp(alphaFELogit)^(0:Tmax), nbIndls, 1) /
    repmat(matrix(apply(1 + exp(alphaFELogit) * V, 1, function(x) {prod(x, na.rm = TRUE)}), ncol = 1), 1, Tmax + 1)

  ### Apply the formula from online appendix B of the DDL paper to get B_t(h) / h^2
  ### Note K is Gaussian so Int u² K(u) du = 1
  indl_bias_approx_St <- matrix(0, nbIndls, Tmax + 1)

  indl_bias_approx_St[, 1] <-
    PSt_obsd_X_constant_alpha[, 1] * sum(beta_hat^2) *
    rowSums(
      abs(- V * exp(alphaFELogit) * (1 - exp(alphaFELogit) * V) / (1 + exp(alphaFELogit) * V)^2), # The abs is not in the original formula but it seems to work better that way
      na.rm = TRUE
    )

  for (t in 1:Tmax) {
    indl_bias_approx_St[, t + 1] <-
      PSt_obsd_X_constant_alpha[, t + 1] * sum(beta_hat^2) *
      rowSums(
        C_S_minus_one_one_period_out_array[, t,] / repmat(matrix(C_S_tab[, t + 1], ncol = 1), 1, Tmax) * V * (1 - exp(alphaFELogit) * V) / (1 + exp(alphaFELogit) * V) -
          exp(alphaFELogit) * V * (1 - exp(alphaFELogit) * V) / (1 + exp(alphaFELogit) * V)^2,
        na.rm = TRUE
      )
  }

  indl_bias_approx_St[is.na(indl_bias_approx_St)] <- 0    # NAs represent observations for which Tobs < Tmax so S cannot take some values, i.e. the probability S > Tobsd is 0
  bias_approx_St <- sqrt(colMeans(indl_bias_approx_St^2)) # If we hadn't replaced NAs by 0s, the mean would be on a sample of smaller size (NAs would be ignored, not seen as 0)

  ### Finally, we solve sigma_t(h)² = R_n B_t(h)² as in online appendix B of
  ### the paper, to choose the bandwith h_t for each t
  h_local_lin_proba <- (std_squared_approx_St / bias_approx_St^2 / Rn_choice_bandwidth)^(1 / (4 + Tmax * dimX))

  h_local_lin_proba <- ifelse(h_local_lin_proba == Inf, 20 * (max(X) - min(X)), h_local_lin_proba) # If no individual bias (no difference between large individuals) take a large window

  env_poke(estimators, "h_local_lin_proba", h_local_lin_proba)

  # Perform the estimation of the local regressions and save them
  out <- local_lin_proba(data, h_local_lin_proba)
  env_poke(estimators, "condlProbas", out$condlProbas)
  if (!env_has(estimators, "densityEstimate")) { # Also save the density estimate for the last bandwidth, maybe should take means ?
    env_poke(estimators, "densityEstimate", out$densityEstimate[, 1])
  }
}
