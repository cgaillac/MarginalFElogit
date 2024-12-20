#' Computes the bounds on the AME on the dataset, at the period given in input,
#' and using the method specified as entry parameter.
#'
#' @param data is an environment variable containing the data. If not already
#' formatted, it will be formatted by the format_data function (see its
#' documentation for information about the stored variables). It must contain at
#' least:
#' - data$Y is a matrix of size n x Tmax containing the values of the dependent variable Y,
#' - data$X is an array of size n x Tmax x dimX containing the values of the covariates X
#' - data$clusterIndexes is a vector of size n x 1 that specifies the cluster each
#' observation pertains to. If it does not exist, the function enforces the default
#' setting of i.i.d. observations - the parameter takes value 1:n so that each
#' observation is its own cluster.
#' @param TEstim is a scalar specifying the period at which the effect must be
#' estimated. It can also be a vector of length n if the period at which the
#' AME must be estimated depends on the individual (e.g. in event-studies models).
#' @param estimators (default empty environment) is an environment in which the
#' results of the CMLE estimation and of the non-parametric estimation of the
#' conditional distribution of S will be stored. That way, if an input is
#' specified, these numbers will not need to be computed again if compute_AME_t
#' is called another time, as the same environment can be passed in input again.
#' Results are stored as follows:
#' - estimators$beta_hat, a list which contains the results from CMLE
#' estimation:
#' - estimators$beta_hat$beta_hat a vector of length dimX, the
#' estimated value for the slope parameter.
#' - estimators_beta_hat$phi_b a matrix of size n x dimX containing
#' the value of the influence function at each observation (rows) w.r.t.
#' each dimension of the covariates (columns).
#' - estimators$beta_hat$var_b the estimated asymptotic covariance
#' matrix, of size dimX x dimX, for the estimator beta_hat.
#' - estimators$h_local_lin_proba a vector of length (Tmax + 1) containing,
#' in j-th position, the bandwidth used to estimate the P(S = j - 1 | X)'s.
#' - estimators$condlProbas: a matrix of size n x (Tmax + 1) containing,
#' in position (i, j), the estimate for P(S = j - 1 | X) at the i-th
#' observation.
#' - estimators$densityEstimate a matrix of size n x (Tmax + 1) containing,
#' at each row (individual), the estimated density for having covariates
#' (X_1, ..., X_T). Each column represents the value found using the
#' corresponding bandwidths from h.
#' @param other_numbers (default empty environment) is an environment variable
#' in which the variable other_numbers$comb_numbers is stored. A matrix of size
#' (Tmax + 1) x (Tmax + 1) containg in position (i, j) the number (i choose j).
#' If j > i the value is NA.
#' @param selectX (default NULL) a vector containing the indices of the
#' covariates w.r.t. which the AME must be computed. If null, the AME w.r.t. all
#' covariates will be computed. All variables of interest should be given in the
#' same call, as then there will be no additional cost relative to estimating
#' w.r.t. only one covariate.
#' @param Option (default "quick") Estimation method to be used. If "quick" or
#' "outer" (case-insensitive) the outer bounds are computed. Otherwise, the sharp
#' bounds are computed. We recommend using the outer bounds method if the number
#' of covariates is at least three, if the number of periods observed is four or
#' more, or if the sample size is small (less than 500) or large (more than 10^4).
#' @param CIOption (default "CI2") When the outer bounds method is being used,
#' specifies which confidence interval should be used. If "CI2", the CI2
#' confidence interval is being used (see DDL, section 4.2), otherwise the CI3
#' confidence interval will be used (see DDL, appendix C).
#' @param alpha (default 0.05) desired asymptotic level of the estimated
#' confidence intervals
#' @param nbCores (default 4) number of cores to be used for parallel computing,
#' to speed up the estimation of the sharp bounds.
#'
#' @return a list containing:
#'  - Option: same as input
#'  - reducedSampleSize: the number of individuals used to estimate the AME. The
#'    individuals that were not observed at the period for which the AME should
#'    be computed were excluded.
#'  - computationTime: the time (in seconds) taken to perform the computation
#'  - estimatedAMEbounds: a matrix of size |selectX| x 2 containing the estimated
#'    (sharp or outer, depending on Option) bounds on the AME, for each covariate
#'    in selectX.
#'  - CI: a matrix of size |selectX| x 2 containing in each row the estimated
#'    confidence interval for the AME, for the corresponding covariate in selectX.
#'  - alt_CI: a matrix of size |selectX| x 2 containing in each row the estimated
#'    confidence interval for the AME, for the corresponding covariate in selectX,
#'    using the confidence interval option NOT requested by the user, for the outer
#'    bounds. If the sharp bounds were being computed, this output is NA.
#'  - indl_estimates: a list containing the bounds and influence functions for
#'    individual observations. They will be used to estimate the average AME
#'    over all periods. They are organised as follows:
#'  - indl_estimates$indl_bounds_on_delta is a matrix of size n x 2 x dimX
#'    containing the estimated lower and upper bounds on the AME, for each
#'    individual observation. It is NA for the individuals unobserved at the
#'    requested period.
#'  - indl_estimates$indl_infl_func_delta is a matrix of size n x 2 x dimX
#'    containing the estimated value of the influence function for the lower and
#'    upper bounds on the AME, for each individual observation. It is NA for the
#'    individuals unobserved at the requested period.
#'
#' @export
#'
# @examples
compute_AME_t <- function(data, TEstim, estimators = env(), other_numbers = env(), selectX = NULL, Option = "quick", CIOption = "CI2", alpha = 0.05,  nbCores = 4) {

  # Start monitoring the computation time
  startTime <- Sys.time()

  # Initialisation ----------------------------------------------------------

  # Values derived from the data (shape, number of periods observed ...)
  dimX <- dim(data$X)[3]
  fullSampleSize <- dim(data$X)[1]
  Tmax <- dim(data$X)[2]
  Tobsd <- rowSums(!is.na(data$Y))

  # Number of variables to compute the AME
  if (is.null(selectX)) {
    selectX <- 1:dimX
  }
  nb_var <- length(selectX)

  # Beta estimated from the CMLE
  if (!env_has(estimators, "beta_hat")) {
    beta_hat <- estim_slope_cmle_maxn(data)
    env_poke(estimators, "beta_hat", beta_hat)
  }
  beta_hat <- estimators$beta_hat$beta_hat
  phi_b <- estimators$beta_hat$phi_b

  # If not already computed, we compute the matrix of combinatorial numbers
  # (i, j)-th value is i choose j. We use recursion
  if (!env_has(other_numbers, "comb_numbers")) {
    comb_numbers <- matrix(NA, Tmax + 1, Tmax + 1)
    comb_numbers[, 1] <- 1
    diag(comb_numbers) <- 1
    if (Tmax > 1) {
      for (k in 2:Tmax) {
        comb_numbers[k + 1, 2:k] <- comb_numbers[k, 1:(k-1)] + comb_numbers[k, 2:k]
      }
    }
    env_poke(other_numbers, "comb_numbers", comb_numbers)
  }
  comb_numbers <- other_numbers$comb_numbers

  # Ensure the data has the right format (and format it if not)
  format_data(data, beta_hat)

  # We extract the individuals for which there is an observation at the period for which we want
  # to compute the AME (and for which TEstim is not NA)
  keepIndls <- !is.na(extractIndlPeriod(data$Y, TEstim))
  nbIndls <- sum(keepIndls)
  X <- array(data$X[keepIndls,,], c(nbIndls, Tmax, dimX))
  Y <- data$Y[keepIndls,]
  clusterIndexes <- data$clusterIndexes[keepIndls]
  S <- data$S[keepIndls]
  V <- data$V[keepIndls,]
  C_S_tab <- data$C_S_tab[keepIndls,]
  C_S_minus_one_one_period_out_array <- data$C_S_minus_one_one_period_out_array[keepIndls,,]
  phi_b <- matrix(phi_b[keepIndls,], ncol = dimX)
  if (length(TEstim) > 1) {
    TEstim <- TEstim[keepIndls]
  }
  Tobsd <- Tobsd[keepIndls]

  # We compute the "tilded" values (difference or ratio of variables at each period v.s. at TEstim)
  Xtilde <- array(NA, dim(X))
  for (k in 1:dim(X)[3]) {
    X_tk <- extractIndlPeriod(X[,, k], TEstim)
    Xtilde[,, k] <- X[,, k] - repmat(matrix(X_tk, ncol = 1), 1, dim(X)[2])
  }

  Vt <- extractIndlPeriod(V, TEstim)
  Vtilde <- V / repmat(matrix(Vt, ncol = 1), 1, dim(X)[2])

  # Computation of statistics common to both algorithms  --------

  # Get the approximated polynomial for Omega(·, x)
  out <- best_approx_poly(Vtilde, Xtilde)
  coeffApproxPolyOmega <- out[["bestPoly"]]
  topCoeffOmega <- out[["topCoeffOmega"]]
  derivCoeffsOmegaDbeta <- out[["derivCoeffsOmegaDbeta"]]
  coeffsOmega <- out[["coeffsOmega"]]

  if (Option == "quick") {

    # Perform the quick estimation algorithm (outer bounds) -------------------

    # We compute the terms in the formula for p(X, S, beta), given in section 4.2.
    ### First we compute the matrix (T - t choose S - t) for T the number of period for the i-th observation (rows),
    ### S the number of times Y = 1 for that observation and t = 0, ..., Tmax (columns).
    Tobsd_minus_t <- repmat(matrix(Tobsd, ncol = 1), 1, Tmax + 1) - repmat(0:Tmax, nbIndls, 1)
    Tobsd_minus_t[Tobsd_minus_t < 0] <- NA

    S_minus_t <- repmat(matrix(S, ncol = 1), 1, Tmax + 1) - repmat(0:Tmax, nbIndls, 1)
    S_minus_t[S_minus_t < 0] <- NA

    T_minus_t_choose_S_minus_t_mat <-
      extractIndlPeriod(
        comb_numbers,
        vecIndls = Tobsd_minus_t + 1,
        vecPeriods = S_minus_t + 1
      )

    ### Multiply by the coefficients of the polynomial
    p_terms_mat <- T_minus_t_choose_S_minus_t_mat * coeffApproxPolyOmega

    ### Sum and multiply by exp(s v(x, beta)) / C_s(x, beta)
    p_vec <- rowSums(p_terms_mat, na.rm = TRUE) * Vt^S / extractIndlPeriod(C_S_tab, S + 1)

    ### Multiply by the different components of beta to get the Delta_hat for each characteristic
    Delta_hat <- beta_hat[selectX] * mean(p_vec)

    # We compute the maximal bias estimate and deduce the bounds on the AME
    indl_bias_sup <- 1 / 2 / 4^Tobsd * T_minus_t_choose_S_minus_t_mat[, 1] * Vt^S / extractIndlPeriod(C_S_tab, S + 1) * abs(topCoeffOmega)
    bias_sup <- abs(beta_hat[selectX]) * mean(indl_bias_sup)
    average_bounds_on_delta <- cbind(Delta_hat - bias_sup, Delta_hat + bias_sup)

    # We save the resulting individual bounds from the individual estimates
    # of delta_hat(X) and of the maximal bias at X
    indl_estimates <- vector("list")
    indl_estimates[["indl_bounds_on_delta"]] <- array(NA, c(fullSampleSize, 2, nb_var))
    indl_estimates[["indl_bounds_on_delta"]][keepIndls, 1,] <- repmat(matrix(p_vec - indl_bias_sup, ncol = 1), 1, nb_var) * repmat(beta_hat[selectX], nbIndls, 1)
    indl_estimates[["indl_bounds_on_delta"]][keepIndls, 2,] <- repmat(matrix(p_vec + indl_bias_sup, ncol = 1), 1, nb_var) * repmat(beta_hat[selectX], nbIndls, 1)

    # Compute the confidence intervals for the quick estimation algorithm -----

    # We compute the influence functions for delta, as given page 22 of DDL
    # Each column corresponds to the influence function relative to one of the beta_k
    # We ignore the constant term E(p(X, S, beta)) which has no impact on the variance
    # Note that in DDL, their p is our beta * p here, so that
    # psi_ik = beta_k * p_i + E(p) phi_ik + beta' * phi_i * E(dp/dbeta)
    infl_func_delta <-
      repmat(matrix(p_vec, ncol = 1), 1, nb_var) * repmat(beta_hat[selectX], nbIndls, 1) + # p
      mean(p_vec) * phi_b[, selectX] # E(p) * phi_i

    # We compute and add the last term, using
    # dp/dbeta_k = p * (S X_Testim,k - sum_t X_tk V_tk C_(S-1)(X1, X_t-1, X_t+1, X_T, beta) / C_S(X, beta))
    # + sum_t=0^S Vt^S / C_S(X, beta) d/dbeta(lambda_t + b*_t,T lambda_(T+1))

    ### First we need to extract the C_(S-1)(...) at the observed values of S
    C_S_minus_one_one_period_out <- matrix(NA, nbIndls, Tmax)
    for (t in 1:Tmax) {
      C_S_minus_one_one_period_out[, t] <- extractIndlPeriod(C_S_minus_one_one_period_out_array[,, t], max(S, 1))
      C_S_minus_one_one_period_out[S == 0, t] <- 0
    }

    ### Now apply the formula
    for (k in 1:dimX) {
      infl_func_delta <-
        infl_func_delta +
        repmat(beta_hat[selectX], nbIndls, 1) * repmat(matrix(phi_b[, k], ncol = 1), 1, nb_var) *
          mean(
            p_vec * (
              S * extractIndlPeriod(X[,, k], TEstim) -
              rowSums(X[,, k] * V * C_S_minus_one_one_period_out, na.rm = TRUE) / extractIndlPeriod(C_S_tab, S + 1)
            ) +
            rowSums(T_minus_t_choose_S_minus_t_mat * derivCoeffsOmegaDbeta[,, k], na.rm = TRUE) * Vt^S / extractIndlPeriod(C_S_tab, S + 1)
          )
    }

    ### Save the individual influence functions
    indl_estimates[["indl_infl_func_delta"]] <- array(NA, c(fullSampleSize, 2, nb_var))
    indl_estimates[["indl_infl_func_delta"]][keepIndls, 1,] <- infl_func_delta
    indl_estimates[["indl_infl_func_delta"]][keepIndls, 2,] <- infl_func_delta

    # Now we know the influence functions, we can compute the clustered standard deviation
    infl_func_delta_cluster <- rowsum(infl_func_delta, clusterIndexes)
    std_delta <- apply(infl_func_delta_cluster, 2, sd, na.rm = TRUE) / sqrt(nbIndls)

    # Compute the confidence intervals for each AME
    length_CI2 <- rep(NA, nb_var)
    length_CI3 <- rep(NA, nb_var)
    for (i in 1:nb_var) {
      ### First compute the corrected parameter for CI3.
      ### We take alpha1 = 4 * alpha / 5 and alpha2 = alpha / 5
      ### Note we use the full sample size, and not nbIndls, for the
      ### asymptotic variance of the CMLE as it's estimated on the
      ### whole set.
      bCI3 <- (abs(beta_hat[selectX[i]]) + qnorm(1 - 4 * alpha / 5) * estimators$beta_hat$var_b[selectX[i], selectX[i]] / fullSampleSize) * mean(indl_bias_sup)

      ### Then compute the length of each CI
      length_CI2[i] <- 2 * std_delta[i] * sqrt(qchisq(1 - alpha, df = 1, ncp = (bias_sup[i] / std_delta[i])^2))
      length_CI3[i] <- 2 * std_delta[i] * sqrt(qchisq(1 - alpha / 5, df = 1, ncp = (bCI3 / std_delta[i])^2))
    }
    CI2 <- cbind(matrix(Delta_hat - length_CI2 / 2, ncol = 1), matrix(Delta_hat + length_CI2 / 2, ncol = 1))
    CI3 <- cbind(matrix(Delta_hat - length_CI3 / 2, ncol = 1), matrix(Delta_hat + length_CI3 / 2, ncol = 1))

  } else{

    # Perform the sharp bound estimation --------------------------------------

    # Estimate the conditional probabilities that S = t if not done already
    if (!env_has(estimators, "condlProbas")) {
      choose_bandwidth_and_estim_condl_proba(data, estimators)
    }
    PSt_X <- estimators$condlProbas[keepIndls,]

    # We restrict the datasets to individuals observed at the period of interest
    focused_data <- env("X" = X, "Y" = Y, "clusterIndexes" = clusterIndexes, "S" = S, "V" = V, "C_S_tab" = C_S_tab, "Tobsd" = Tobsd, "C_S_minus_one_one_period_out_array" = C_S_minus_one_one_period_out_array)
    focused_estimators <- env("beta_hat" = estimators$beta_hat, "condlProbas" = estimators$condlProbas[keepIndls,], "densityEstimate" = estimators$densityEstimate, "alphaFElogit" = estimators$alphaFELogit, "h_local_lin_proba" = estimators$h_local_lin_proba)

    # Obtain the bounds on the AME. Note they still need to be multiplied by the
    # relevant beta's
    out <- bound_delta(focused_data, focused_estimators, other_numbers, Vt, coeffsOmega, topCoeffOmega, nbCores)
    average_bounds_on_delta <- out$average_bounds_on_delta
    indl_bounds_on_delta <- out$indl_bounds_on_delta

    # Compute the influence functions -----------------------------------------

    # First we need to (numerically) estimate the derivative of the bounds
    # in beta and gamma.
    h_deriv <- 10^-3
    deriv_beta <- matrix(NA, dimX, 2)
    deriv_gamma <- array(NA, c(nbIndls, Tmax + 1, 2))

    ### Gamma derivative
    for (t in 0:Tmax) {
      estimators_deriv <- env()
      estimators_deriv$condlProbas <- focused_estimators$condlProbas + repmat((0:Tmax == t) * h_deriv, nbIndls, 1)
      estimators_deriv$densityEstimate <- focused_estimators$densityEstimate
      estimators_deriv$h_local_lin_proba <- focused_estimators$h_local_lin_proba

      for_deriv_indl_bounds_on_delta <-
        bound_delta( # The first argument is the same as data but individuals which are not observed at the period of interest are taken out
          env("X" = X, "Y" = Y, "clusterIndexes" = clusterIndexes, "S" = S, "V" = V, "C_S_tab" = C_S_tab, "Tobsd" = Tobsd, "C_S_minus_one_one_period_out_array" = C_S_minus_one_one_period_out_array),
          estimators_deriv, other_numbers,
          Vt,
          coeffsOmega,
          topCoeffOmega,
          nbCores)$indl_bounds_on_delta

      deriv_gamma[, t + 1,] <- (for_deriv_indl_bounds_on_delta - indl_bounds_on_delta) / h_deriv
    }

    ### Beta derivative
    for (k in 1:dimX) {
      # We create a copy of the data, which we format, to recompute the C_S_beta etc
      data_deriv <- env()
      data_deriv$X <- focused_data$X
      data_deriv$Y <- focused_data$Y
      data_deriv$clusterIndexes <- focused_data$clusterIndexes
      format_data(data_deriv, beta_hat + (1:dimX == k) * h_deriv)

      # We compute the new Vt and the coefficients for the new beta
      Vtilde_deriv <- Vtilde * exp(h_deriv * Xtilde[,, k])
      Vt_deriv <- extractIndlPeriod(data_deriv$V, TEstim)

      out <- best_approx_poly(Vtilde_deriv, Xtilde)
      topCoeffOmega_deriv <- out[["topCoeffOmega"]]
      coeffsOmega_deriv <- out[["coeffsOmega"]]

      # We use a small trick with signs to ensure that we respect
      # the sign criterion used in Appendix C.2 of DDL
      changeOfSign <- sign(topCoeffOmega_deriv * topCoeffOmega)


      for_deriv_average_bounds_on_delta <- colMeans(repmat(matrix(changeOfSign, ncol = 1), 1, 2) * bound_delta(data_deriv, focused_estimators, other_numbers, Vt_deriv, repmat(matrix(changeOfSign, ncol = 1), 1, Tmax + 2) * coeffsOmega_deriv, changeOfSign * topCoeffOmega_deriv, nbCores)$indl_bounds_on_delta)

      deriv_beta[k,] <- (for_deriv_average_bounds_on_delta - average_bounds_on_delta) / h_deriv
    }

    # Now we can compute the influence functions.
    # Here we account for the extra multiplication by beta
    # which was not present in the above values of the bounds.
    infl_func_delta <- array(NA, c(fullSampleSize, 2, nb_var))
    indic_S_minus_PSt_X <- (repmat(matrix(S, ncol = 1), 1, Tmax + 1) == repmat(0:Tmax, nbIndls, 1)) - PSt_X
    for (i in 1:nb_var) {
      k <- selectX[i]
      infl_func_delta[keepIndls,, i] <-
        beta_hat[k] * indl_bounds_on_delta +
        repmat(matrix(phi_b[, k], ncol = 1), 1, 2) * repmat(average_bounds_on_delta, nbIndls, 1) +
        beta_hat[k] * phi_b %*% deriv_beta +
        beta_hat[k] * matrix(
          c(
            rowSums(indic_S_minus_PSt_X * deriv_gamma[,, 1], na.rm = TRUE),
            rowSums(indic_S_minus_PSt_X * deriv_gamma[,, 2], na.rm = TRUE)
          ),
          ncol = 2
        )
    }

    # Estimate the asymptotic variance from the influence functions
    # and compute confidence intervals
    infl_func_delta_cluster <- array(NA, c(length(unique(clusterIndexes)), 2, nb_var))
    infl_func_delta_cluster[, 1,] <- rowsum(infl_func_delta[keepIndls, 1,], clusterIndexes, na.rm = TRUE)
    infl_func_delta_cluster[, 2,] <- rowsum(infl_func_delta[keepIndls, 2,], clusterIndexes, na.rm = TRUE)

    CI <- matrix(NA, nb_var, 2)
    cur_beta_is_not_null <- rep(NA, nb_var)
    for (i in 1:nb_var) {
      k <- selectX[i]
      std_bounds <- apply(infl_func_delta_cluster[,, i], 2, sd, na.rm = TRUE) / sqrt(nbIndls)
      # The IM_quantile is symmetric in the two elements of std_bounds so no need to rearrange
      IM_quantile <- compute_IM_quantile(alpha, sort(beta_hat[k] * average_bounds_on_delta), std_bounds)

      ### This is a t-test to see if the beta corresponding to the current AME
      ### is null
      cur_beta_is_not_null[i] <- sqrt(fullSampleSize * beta_hat[k]^2 / estimators$beta_hat$var_b[k, k]) > qnorm(1 - alpha/2)

      if (beta_hat[k] >= 0) {
        CI[i, ] <- sort(beta_hat[k] * average_bounds_on_delta) + IM_quantile * c(-std_bounds[1], std_bounds[2])
      } else {
        CI[i, ] <- sort(beta_hat[k] * average_bounds_on_delta) + IM_quantile * c(-std_bounds[2], std_bounds[1])
      }

      if (!cur_beta_is_not_null[i]) {
        CI[i, 1] <- min(CI[i, 1], 0)
        CI[i, 2] <- max(CI[i, 2], 0)
      }
    }

    # Account for the missing beta factor for the average bounds
    average_bounds_on_delta <- repmat(matrix(beta_hat[selectX], ncol = 1), 1, 2) * repmat(average_bounds_on_delta, nb_var, 1)

    # Save individual estimates
    indl_estimates <- vector("list")
    indl_estimates[["indl_bounds_on_delta"]] <- array(NA, c(fullSampleSize, 2, nb_var))
    for (k in 1:nb_var) {
      indl_estimates[["indl_bounds_on_delta"]][keepIndls,, k] <- beta_hat[selectX[k]] * indl_bounds_on_delta
    }
    indl_estimates[["indl_infl_func_delta"]] <- infl_func_delta

    # Check the sign of each beta, to revert the columns to the correct answer if it is negative
    for (i in 1:length(selectX)) {
      k <- selectX[i]
      if (beta_hat[k] < 0) {
        average_bounds_on_delta[i,] <- average_bounds_on_delta[i, c(2, 1)]
        indl_estimates$indl_bounds_on_delta[,, i] <- indl_estimates$indl_bounds_on_delta[, c(2, 1), i]
        indl_estimates$indl_infl_func_delta[,, i] <- indl_estimates$indl_infl_func_delta[, c(2, 1), i]
      }
    }
  }

  # Stop monitoring the computation time
  endTime <- Sys.time()
  time <- difftime(endTime, startTime, units = "secs")


  out <- vector("list")
  out[[1]] <- Option
  out[[2]] <- nbIndls
  out[[3]] <- as.numeric(time)
  out[[4]] <- average_bounds_on_delta
  # We return the two confidence intervals for the outer bounds
  if (Option == "quick") {
    if (CIOption == "CI2") {
      out[[5]] <- CI2
      out[[6]] <- CI3
    }
    else {
      out[[5]] <- CI3
      out[[6]] <- CI2
    }
  } else {
    out[[5]] <- CI
    out[[6]] <- NA
  }
  out[[7]] <- indl_estimates

  names(out) <- c("Option", "reducedSampleSize", "computationTime", "estimatedAMEbounds", "CI", "alt_CI", "indl_estimates")

  return(out)
}
