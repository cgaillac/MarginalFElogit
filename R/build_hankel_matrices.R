#' This function returns the upper and lower Hankel matrices given a moment
#' vector. See DDL paper section 3.2.2. The NAs at the tail of the moment vector
#' are ignored.
#'
#' The function may also return the derivatives of the matrices with respect to
#' each moment.
#'
#' @param m a moment vector with tail in NAs, which should be ignored. Starts at
#' the moment of order 1, unless add_zeroeth_moment is FALSE.
#' @param differentiating (default FALSE) if TRUE, the added zeroeth moment is 0.
#' The input moment should then be only 0's, but for the coordinate w.r.t. which
#' we want to differentiate, which should be one (by linearity).
#' @param add_zeroeth_moment (default TRUE) if FALSE, the zeroeth moment (1) is
#' added at the start of the input moment vector.
#'
#' @return A list with two elements, upperHankel and lowerHankel, corresponding to the
#' upper bar and lower bar H matrices in the DDL paper (see section 3.2.2),
#' respectively. If differentiating is TRUE and m is as specified above, the output
#' is also the derivative of the Hankel matrices w.r.t. the relevant coordinate, by
#' linearity.
#' @export
#'
# @examples
build_hankel_matrices <- function(m, differentiating = FALSE, add_zeroeth_moment = TRUE) {

  # Add the initial 1 (0-th moment) and take out the tail of NAs
  # If we are differentiating, the 0-th moment('s derivative) is null,
  # otherwise it's the mass of a probability law, i.e. 1
  if (add_zeroeth_moment) {
    m <- c(ifelse(differentiating, 0, 1), m)
  }
  m <- m[1:max(which(!is.na(m)))]

  # Find the highest order of a moment in the truncated moment vector
  maxMoment <- length(m) - 1

  # Deduce the Hankel matrices, depending on whether maxMoment is even or odd
  # Remember m[i + 1] corresponds to the i-th moment
  if ((maxMoment %% 2) == 0) {
    lowerHankel <- hankel(m[1:(1 + maxMoment/2)], m[(1 + maxMoment/2):(1 + maxMoment)])
    upperHankel <-
      hankel(m[2:(1 + maxMoment/2)], m[(1 + maxMoment/2):maxMoment]) -
      hankel(m[3:(2 + maxMoment/2)], m[(2 + maxMoment/2):(1 + maxMoment)])
  } else {
    lowerHankel <- hankel(m[2:(1 + (maxMoment + 1)/2)], m[(1 + (maxMoment + 1)/2):(1 + maxMoment)])
    upperHankel <-
      hankel(m[1:((maxMoment + 1)/2)], m[((maxMoment + 1)/2):maxMoment]) -
      hankel(m[2:(1 + (maxMoment + 1)/2)], m[(1 + (maxMoment + 1)/2):(1 + maxMoment)])
  }

  return(list(
    "lowerHankel" = lowerHankel,
    "upperHankel" = upperHankel
  ))
}
