#' Bounds the (Tobsd + 1)-th moment given the first Tobsd moments, m. We
#' first check that the input moment of vectors is valid, and if necessary
#' replace it by the vector on the boundary of the set of admissible vectors
#' of moment (see DDL paper).
#'
#' @param m is the estimated vector of moments. NAs at the end is ignored.
#' @param asymptotic_variance_matrix_m (default NULL) is the estimated variance matrix for
#' the vector of moments m. If NULL, the moment vector is assumed to be known
#' with perfect accuracy (no noise).
#' @param thresholdDeterminant (default 0). The threshold used to decide whether the
#' determinant of the Hankel matrices is null. If the estimated determinant is less
#' than that threshold, it is considered null (so the moment vector is on the boundary).
#' Setting the threshold as a small positive value is meant to account for the noise on
#' the moment vector, e.g. to ensure moment vectors on the boundary are rejected exactly
#' 95\% of the time.
#'
#' @return a list containing:
#' - lowerBound: the lower bound on the (Tobsd + 1)-th moment
#' - upperBound: the upper bound on the (Tobsd + 1)-th moment
#' @export
#'
#@examples
bound_next_moment <- function(m, asymptotic_variance_matrix_m = NULL, thresholdDeterminant = 0) {

  # Desired moment is maxOrderGivenMoment + 1
  maxOrderGivenMoment <-  ifelse(all(is.na(m)), 0, max(which(!is.na(m))))

  # First truncate the moment vector until it is valid
  out <- truncate_moment_vector_until_valid(m, asymptotic_variance_matrix_m, thresholdDeterminant)
  truncated_m <- out[["truncatedMomentVector"]]
  infringedConstraint <- out[["infringedConstraint"]]

  corrected_m <- truncated_m
  maxOrderInnerMoment <- length(truncated_m) # Careful, this is T' - 1, not T', in the notation of Proposition 2 from the DDL paper

  # We compute the bounds on the next moment
  while (length(corrected_m) <= maxOrderGivenMoment) {

    # short_m is the vector of moments from T - T' + 1 to T for corrected_m,
    # as defined in Proposition 2 of the DDL paper
    # At the first step post-truncation, we need to include the zeroeth moment
    if (length(corrected_m) > maxOrderInnerMoment) {
      short_m <- corrected_m[(length(corrected_m) - maxOrderInnerMoment):length(corrected_m)]
    } else {
      short_m <- c(1, corrected_m)
    }

    hankelMatricesWith0 <- build_hankel_matrices(c(short_m, 0), add_zeroeth_moment = FALSE)
    hankelMatricesWith1 <- build_hankel_matrices(c(short_m, 1), add_zeroeth_moment = FALSE)

    detHankelWith0 <- lapply(hankelMatricesWith0, det)
    detHankelWith1 <- lapply(hankelMatricesWith1, det)

    bounds <- list()
    for (index in c("upperHankel", "lowerHankel")) {
      bounds[[index]] <- - detHankelWith0[[index]] / (detHankelWith1[[index]] - detHankelWith0[[index]])
    }

    # If we are still among the truncated moments, because the "upper bound" and
    # the "lower bound" are achieved for a unique moment vectors, we just need
    # to know which constraint was binding to be able to extend the vector
    # N.B. : technically we could avoid computing both matrices and both
    # determinants for these moments in the above, could save some computation
    # time
    if (length(corrected_m) < maxOrderGivenMoment) {
      if (infringedConstraint == "Upper") {
        corrected_m <- c(corrected_m, bounds[["upperHankel"]])
      } else if (infringedConstraint == "Lower") {
        corrected_m <- c(corrected_m, bounds[["lowerHankel"]])
      }
    } else { # Otherwise move out of the loop - if we were on the boundary,
             # ensure the two bounds are the same (note only one of the two
             # equations in point 2 of Proposition 2 of DDL is valid, which
             # is why we need this correction)
      if (infringedConstraint == "Upper") {
        bounds[["lowerHankel"]] <- bounds[["upperHankel"]]
      } else if (infringedConstraint == "Lower") {
        bounds[["upperHankel"]] <- bounds[["lowerHankel"]]
      }
      break
    }

  }

  return(c("lowerBound" = bounds[["lowerHankel"]],"upperBound" = bounds[["upperHankel"]]))
}
