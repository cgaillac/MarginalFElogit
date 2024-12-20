#' Computes the coefficients of T_0(2u - 1), ..., T_n(2u-1), where T_k is the
#' usual k-th Chebyshev polynomial of the 1st kind
#'
#' @param n maximal k for which we compute the coefficients of T_k(2u - 1)
#'
#' @return the result as a (n + 1) x (n + 1) matrix where each row is one
#' of the polynomials. The first column corresponds to the constant coefficient,
#' and rows are filled with 0s at the end for T_0, ..., T_(n-1)
#' @export
#'
# @examples
modifiedChebyshevPoly <- function(n) {

  # We use that T_n = 2 X T_(n-1) - T_(n-2), so that
  # T_n(2u-1) = 2 * (2u - 1) T_(n-1)(2u -1) - T_(n-2)(2u-1)
  # with modified starting values

  # Initialisation
  chebyshevPolyMatrix <- matrix(0, nrow = n + 1, ncol = n + 1)
  chebyshevPolyMatrix[1, 1] <- 1
  if (n >= 1) {
    chebyshevPolyMatrix[2, 1] <- -1
    chebyshevPolyMatrix[2, 2] <- 2
  }

  # Computation of the polynomials using the recursive formula of degree 2
  if (n >= 2) {
    for (k in 2:n) {
      chebyshevPolyMatrix[k + 1, 1] <- - 2 * chebyshevPolyMatrix[k, 1] - chebyshevPolyMatrix[k - 1, 1]
      chebyshevPolyMatrix[k + 1, 2:(k+1)] <-
        4 * chebyshevPolyMatrix[k, 1:k] -
        2 * chebyshevPolyMatrix[k, 2:(k+1)] -
        chebyshevPolyMatrix[k - 1, 2:(k+1)]
    }
  }

  return(chebyshevPolyMatrix)
}
