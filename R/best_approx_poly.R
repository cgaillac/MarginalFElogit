#' This function returns the coefficients (functions of V = exp(X'beta)) of
#' the best (in the uniform sense) approximation  of Omega(·, x, beta) * h(·, x, beta)
#' by a polynomial of degree at most T. It also returns the coefficients of the original
#' polynomial, its top coefficients, and the derivative of the coefficients of the
#' approximating polynomials with respect to beta.
#' See Assumption 3 in DDL and the discussion under that assumption for the exact formulae
#' of the functions. Here, we ignore the beta factor at the start, and use the h
#' corresponding to the AME, unless forATE is TRUE.
#'
#' @param Vtilde a matrix of size n x Tmax containing the exp((x-x_t)'beta) for
#' each observation, t being the period at which the AME is being computed. For the
#' ATE, we want Vtilde = exp(x'beta - v(x, beta)).
#' @param Xtilde an array of size n x Tmax x dimX containing, for each observation,
#' X - X_t where t is the period at which the AME is being computed. For the ATE we
#' want X - X_t^(1 - X_tk) (i.e. the second term is the same as X_t except the k_th
#' coordinate changes from 0 to 1 or from 1 to 0).
#' @param forATE (default FALSE) If TRUE, we return the formula used for the ATE
#' instead of for the AME.
#' @param signFactorCoeffshATE (default NULL) The sign to multiply the function h by
#' in the ATE case. This parameter is unnecessary for the AME.
#'
#' @return a list containing:
#'  - bestPoly: a matrix of size n x (Tmax + 1) where each row represents the
#'    coefficient of the best approximation polynomial of Omega(., x, beta).
#'    Coefficients start with the constant coefficient, and there are
#'    (Tobsd[i] + 1) coefficients on each row, where Tobsd[i] is the number
#'    of periods observed for the relevant individuals. Subsequent entries on
#'    the rows are 0s
#'  - topCoeffOmega: a vector of length n containing the top coefficient of the
#'    original Omega coefficient for each individual.
#'  - derivCoeffsOmegaDbeta: derivative of the coefficients stored in bestPoly
#'    relative to beta.
#'  - coeffsOmega: a matrix of size n x (Tmax + 2) where each row represents the
#'    coefficient of the polynomial Omega(., x, beta). Coefficients start with
#'    the constant coefficient and there are (Tobsd[i] + 2) coefficients on each row
#'    where Tobsd[i] is the number of periods observed for the relevant
#'    individuals. Subsequent entries on the rows are 0s.
#'
#' @export
#'
# @examples
best_approx_poly <- function(Vtilde, Xtilde, forATE = FALSE, signFactorCoeffshATE = NULL) {

  n <- dim(Vtilde)[1]
  Tmax <- dim(Vtilde)[2]
  dimX <- dim(Xtilde)[3]
  coeffsOmega <- matrix(0, n, Tmax + 2)
  coeffsOmega[, 1] <- 1
  derivCoeffsOmegaDbeta <- array(0, c(n, Tmax + 2, dimX))

  if (forATE && is.null(signFactorCoeffshATE)) {
    print("The sign to multiply the function h for the ATE formula is missing")
  }

  # We multiply by each term of the product of the (1 + u (exp((x't-x'T)*beta) - 1)) (AME)
  # or (1 + u (exp(x_t'beta - v(x, beta)) - 1)) (ATE)
  # sequentially. Note that, for the AME, in each row one of the Vtilde is 1
  # We proceed in a similar way for the derivatives
  for (t in 1:Tmax) {
    # First update the derivatives (because we use the polynomial coefficients from the previous iteration)
    derivCoeffsOmegaDbeta[, 2:(t+1),] <- derivCoeffsOmegaDbeta[, 2:(t+1),] +
      derivCoeffsOmegaDbeta[, 1:t,] * (replace(Vtilde[, t], is.na(Vtilde[, t]), 1) - 1)

    for (k in 1:dimX) {
      derivCoeffsOmegaDbeta[, 2:(t+1), k] <- derivCoeffsOmegaDbeta[, 2:(t+1), k] +
        coeffsOmega[, 1:t] * repmat(matrix(replace(Xtilde[, t, k] * Vtilde[, t], is.na(Vtilde[, t]), 0), ncol = 1), 1, t)
    }

    # Now update the polynomial coefficients
    coeffsOmega[, 2:(t+1)] <- coeffsOmega[, 2:(t+1)] +
      coeffsOmega[, 1:t] * (replace(Vtilde[, t], is.na(Vtilde[, t]), 1) - 1)
  }

  # We multiply by u (1 - u) (AME) or (-1)^Xtk * u (ATE)
  if (!forATE) {
    coeffsOmega[, 3:(Tmax+2)] <- coeffsOmega[, 2:(Tmax+1)] - coeffsOmega[, 1:Tmax]
    coeffsOmega[, 2] <- coeffsOmega[, 1]
    coeffsOmega[, 1] <- 0
    derivCoeffsOmegaDbeta[, 3:(Tmax+2), ] <- derivCoeffsOmegaDbeta[, 2:(Tmax+1), ] - derivCoeffsOmegaDbeta[, 1:Tmax, ]
    derivCoeffsOmegaDbeta[, 2, ] <- derivCoeffsOmegaDbeta[, 1, ]
    derivCoeffsOmegaDbeta[, 1, ] <- 0
  } else {
    coeffsOmega[, 2:(Tmax+2)] <- signFactorCoeffshATE * coeffsOmega[, 1:(Tmax+1)]
    coeffsOmega[, 1] <- 0
    derivCoeffsOmegaDbeta[, 2:(Tmax+2), ] <- signFactorCoeffshATE * derivCoeffsOmegaDbeta[, 1:(Tmax+1), ]
  }

  # Now pluck the top coefficient
  # We check for NAs to get the right number of observed periods. The top coefficient is
  # that for the degree corresponding to that number plus one (remember there's also a constant
  # coefficient)
  Tobsd <- rowSums(!is.na(Vtilde))
  topCoeff <- extractIndlPeriod(coeffsOmega, Tobsd + 2)

  # Get the coefficients of the "modified" Chebyshev polynomial
  coeffModdChebyshevPoly <- modifiedChebyshevPoly(Tmax + 1)

  # Derive the best uniform approximation of u^(T+1)
  bestUnifApprox <- - coeffModdChebyshevPoly / repmat(matrix(diag(coeffModdChebyshevPoly), ncol = 1), 1, Tmax + 2)
  diag(bestUnifApprox) <- 0

  # Replace the top monomial in the polynomial Omega by the corresponding approximation
  bestUnifApproxOmega <-
    coeffsOmega +
    repmat(matrix(topCoeff, ncol = 1), 1, Tmax + 2) * bestUnifApprox[Tobsd + 2, ]
  bestUnifApproxOmega[(Tobsd + 1) * n + 1:n] <- 0 # The dominant coefficient is changed to 0 (as we replaced it by an approximation)
  bestUnifApproxOmega <- bestUnifApproxOmega[, 1:(Tmax + 1)] # We take away the last column (only 0s)

  # Do the same for the derivatives
  derivBestUnifApproxOmega <- array(0, c(n, Tmax + 2, dimX))
  for (k in 1:dimX) {
    topCoeffDeriv <- extractIndlPeriod(derivCoeffsOmegaDbeta[,, k], Tobsd + 2)
    derivBestUnifApproxOmega[,, k] <-
      derivCoeffsOmegaDbeta[,, k] +
      repmat(matrix(topCoeffDeriv, ncol = 1), 1, Tmax + 2) * bestUnifApprox[Tobsd + 2, ]
    derivBestUnifApproxOmega[(Tobsd + 1) * n + 1:n + (k - 1) * n * (Tmax + 2)] <- 0
  }
  derivBestUnifApproxOmega <- array(derivBestUnifApproxOmega[, 1:(Tmax + 1),], c(n, Tmax + 1, dimX))

  out <- vector("list")
  out[["bestPoly"]] <- bestUnifApproxOmega
  out[["topCoeffOmega"]] <- topCoeff
  out[["derivCoeffsOmegaDbeta"]] <- derivBestUnifApproxOmega
  out[["coeffsOmega"]] <- coeffsOmega
  return(out)
}
