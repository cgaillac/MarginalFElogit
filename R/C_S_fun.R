#' Computes the elementary polynomials for each row of the input, by induction. This corresponds
#' to the what the function C_S does to V in the DDL paper.
#'
#' @param V the matrix of size n x T whose coefficients are
#'          used in the elementary symmetric polynomial. NAs can be used for
#'          unobserved values.
#'
#' @return a matrix of size n x (T + 1) whose value in the (i, j)-th position
#' is C_(j-1)(X, beta) using the values of V = X' * beta in the i-th row. When
#' j - 1 is larger than the number of non-NA values of V in the relevant row, the
#' value NA is returned.
#'
#' @export
#'
# @examples
C_S_fun <- function(V) {

  # Initialisation
  n <- dim(V)[1]
  Tmax <- dim(V)[2]

  # Note that if there are NAs in V, ignoring them is the same as replacing them
  # with 0s and computing the overall symmetric polynomial.
  # So we do this replacement to simplify the code, and remember where the values
  # of S didn't work
  wrongS <- repmat(0:Tmax, n, 1) > rowSums(!is.na(V))
  V[is.na(V)] <- 0

  # We use the recursion formula, where the k-th symmetric polynomial
  # in V_1, ..., V_n is obtained by multiplying V_n by the (k-1)-th
  # symmetric polynomial in V_1, ..., V_(n-1) and adding the k-th symmetric
  # polynomial in V_1, ..., V_(n-1)
  # Thus we start with the matrix which in the i-th row and j-th column has the
  # (j-1)-th elementary symmetric polynomial in the empty set, and iterate by getting
  # (in the i-th row) that for the symmetric polynomials in V[i, 1], ..., V[i, k]
  # for successive k's.
  tab_C_S <- matrix(0, nrow = n, ncol = Tmax + 1)
  tab_C_S[, 1] <- 1

  for (k in 1:Tmax) {
    tab_C_S[, 2:(k+1)] <- tab_C_S[, 1:k] * V[, k] + tab_C_S[, 2:(k+1)]
  }

  # Put NA where the values of S don't make sense
  tab_C_S[wrongS] <- NA

  return(tab_C_S)
}
