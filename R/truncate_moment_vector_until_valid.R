#' This function truncates m_hat until it belongs to the moment space. The
#' truncated moment is put in m_trun.
#' Sigma is an estimator of the "variance" of m_hat, or more accurately, of
#' the asymptotic variance of m_hat, properly normalized.
#'
#' @param m is the estimated vector of moments. NAs at the end are ignored.
#' @param asymptotic_variance_matrix_m (default NULL) is the estimated variance matrix for
#' the vector of moments m. If NULL, the moment vector is assumed to be known
#' with perfect accuracy (no noise).
#' @param thresholdDeterminant (default 0). The threshold used to decide whether the
#' determinant of the Hankel matrices is null. If the estimated determinant is less
#' than that threshold, it is considered null (so the moment vector is on the boundary).
#' Setting the threshold as a small positive value is meant to account for the noise on
#' the moment vector, e.g. to ensure moment vectors on the boundary are rejected exactly
#' 95% of the time.
#'
#'  @return A list containing:
#'      - truncatedMomentVector: m_1, ..., m_k for the largest k such that m_1, ..., m_k
#'      is a valid vector of moments that is not on the boundary, i.e. the determinants
#'      of the upper and lower hankel matrices for all initial subvectors are positive.
#'      - infringedConstraint, a scalar taking value 0 if the input vector of moments was
#'      valid, 1 if the first invalid moment was too large and -1 if the first invalid
#'      moment was too small.
#'
#' @export
#'
# @examples
truncate_moment_vector_until_valid <- function(m, asymptotic_variance_matrix_m = NULL, thresholdDeterminant = 0) {

  # We start by considering only the first moment and check the vector moment
  # is valid. We then consider the first two moments, etc, until we reach a
  # moment where the moment vector is no longer valid, or we reach the end of the
  # moment vector
  maxMoment <- ifelse(all(is.na(m)), 0, max(which(!is.na(m))))
  curMoment <- 1
  momentVectorIsValid <- TRUE
  infringedConstraint <- "None"

  while (momentVectorIsValid && (curMoment <= maxMoment)) {
    # Compute the determinants of the Hankel matrices
    hankelMatrices <- build_hankel_matrices(m[1:curMoment])
    upperHankelDet <- det(hankelMatrices[["upperHankel"]])
    lowerHankelDet <- det(hankelMatrices[["lowerHankel"]])

    # If we are given the variance matrix of m, we use the
    # Jacobi formula to get the resulting variance on the
    # the determinant of the matrix
    # Jacobi formula : d(det(A(t)))/dt = tr(adj(A(t)) * dA(t)/dt)
    # and the asymptotic variance of f(x) is given by
    # t(df/dx) * Sigma * df/dx if Sigma is the asymptotic variance
    # matrix of x
    if (!is.null(asymptotic_variance_matrix_m)) {
      ### First compute, for the upper and lower Hankel matrices,
      ### d(det(H(m)))/dm_t
      derivUpperHankelDet <- rep(NA, maxMoment)
      derivLowerHankelDet <- rep(NA, maxMoment)
      for (t in 1:maxMoment) {
        ### By linearity dA(m(t))/dm_t = A(dm(t)/dm_t)
        derivHankelMatrices <- build_hankel_matrices(as.numeric((1:curMoment) == t), differentiating = TRUE)

        ### Adjoint is product of inverse and determinant (when defined)
        adjointHankelMatrices <- lapply(hankelMatrices, function(A) {if (dim(A)[1] == 1) {return(1)} else{return(matlib::adjoint(as.matrix(A)))}})

        ### Hankel matrices are symmetric so so are their adjoint
        ### Thus tr(adj(A(t)) dA(t)/dt) is the sum of elementwise products
        derivLowerHankelDet[t] <- sum(adjointHankelMatrices[["lowerHankel"]] * derivHankelMatrices[["lowerHankel"]])
        derivUpperHankelDet[t] <- sum(adjointHankelMatrices[["upperHankel"]] * derivHankelMatrices[["upperHankel"]])
      }

      ### Now deduce the asymptotic variance of the determinants
      asymptoticVarianceUpperHankelDet <- derivUpperHankelDet %*% asymptotic_variance_matrix_m %*% derivUpperHankelDet
      asymptoticVarianceLowerHankelDet <- derivLowerHankelDet %*% asymptotic_variance_matrix_m %*% derivLowerHankelDet
    } else {
    # If we're not given the asymptotic variance of m, we just consider the raw
    # values of the determinants, without comparing it to their respective
    # asymptotic variances (or as if these variances were 1)
      asymptoticVarianceUpperHankelDet <- 1
      asymptoticVarianceLowerHankelDet <- 1
    }

    ### If any of the variances is null or negative, we replace it by 10^-20
    if (asymptoticVarianceUpperHankelDet <= 0) {
      asymptoticVarianceUpperHankelDet <- 10^-20
    }
    if (asymptoticVarianceLowerHankelDet <= 0) {
      asymptoticVarianceLowerHankelDet <- 10^-20
    }

    # Now we use these variances to compare whether the determinants are sufficiently
    # positive for the moment vector to be considered valid (see online appendix B
    # of the DDL paper for details)
    if (min(lowerHankelDet / sqrt(asymptoticVarianceLowerHankelDet), upperHankelDet / sqrt(asymptoticVarianceUpperHankelDet)) <= thresholdDeterminant) {
      momentVectorIsValid <- FALSE
      # If one of the constraints is infringed, specify which one
      if (lowerHankelDet / sqrt(asymptoticVarianceLowerHankelDet) <= upperHankelDet / sqrt(asymptoticVarianceUpperHankelDet)) {
        infringedConstraint <- "Lower"
      } else {
        infringedConstraint <- "Upper"
      }
    } else {
      curMoment <- curMoment + 1
    }
  }

  # We return the largest truncated moment vector that is valid
  lastValidMoment <- curMoment - 1
  if (lastValidMoment == 0) {
    truncatedMomentVector <- c()
  } else{
    truncatedMomentVector <- m[1:lastValidMoment]
  }

  return(list(
    "truncatedMomentVector" = truncatedMomentVector,
    "infringedConstraint" = infringedConstraint
  ))
}
