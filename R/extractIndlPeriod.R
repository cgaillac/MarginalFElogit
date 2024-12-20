#' This function takes as input a matrix M and two vectors, vecIndls and vecPeriods
#' It returns the vector whose i-th line is given by M[vecIndls[i], vecPeriods[i]].
#'
#' Instead of vectors, vecIndls and vecPeriods can also be matrices, in which case
#' the (i, j)-th element of the output is M[vecIndls[i, j], vecPeriods[i, j]].
#'
#' If one is a vector and the other a matrix, the vector is interpreted as column
#' vector standing for a matrix constant on each row (for vecIndls) or as a row vector
#' standing for a matrix constant for each column (for vecPeriods).
#'
#' @param M n x T matrix from which we want to extract elements.
#' @param vecPeriods vector corresponding to the column numbers of the elements we want to extract.
#' Can also be a constant, it will then be constant for all elements of vecIndls. NAs are passed on
#' to the output.
#' @param vecIndls vector corresponding to the row numbers of the elements we want to extract. Can
#' also be a constant, it will then be constant for all elements of vecPeriods. Default is all rows,
#' once, in order (the vector 1:n).
#'
#' @return If vecPeriods and vecIndls are vectors, returns the vector whose i-th line is given
#' by M[vecIndls[i], vecPeriods[i]]. If either is a matrix, returns the matrix whose (i, j)-th
#' element of the output is M[vecIndls[i, j], vecPeriods[i, j]].
#' @export
#'
# @examples
extractIndlPeriod <- function(M, vecPeriods, vecIndls = NULL) {

  # The default is to take each row in order
  if (is.null(vecIndls)) {
    vecIndls <- 1:(dim(M)[1])
    if ((length(vecPeriods) != dim(M)[1]) & (length(vecPeriods) != 1)) {
      stop("If vecIndls is not specified, vecPeriods must be constant or have the same length as the matrix has rows")
    }
  }

  # We repeat the vectors/matrices to give them compatible sizes
  isVec <- FALSE # Output a vector or matrix ?
  if (is.vector(vecPeriods) & is.vector(vecIndls)) {
    isVec <- TRUE
  } else if (is.vector(vecPeriods)) {
    vecPeriods <-
      repmat(
        matrix(vecPeriods, nrow = 1),
        dim(vecIndls)[1],
        dim(vecIndls)[2] / length(vecPeriods)
      )
  } else if (is.vector(vecIndls)) {
    vecIndls <-
      repmat(
        matrix(vecIndls, ncol = 1),
        dim(vecPeriods)[1] / length(vecIndls),
        dim(vecPeriods)[2]
      )
  }

  # We extract data from the output, and make it a matrix if necessary
  output <- M[c(vecIndls + dim(M)[1] * (vecPeriods - 1))]

  # Replace by NAs where indices are out of bounds
  output[c((vecIndls < 0) | (vecIndls > dim(M)[1]) |(vecPeriods < 0) | (vecPeriods > dim(M)[2]))] <- NA

  if (!isVec) {
    output <- matrix(output, nrow = dim(vecIndls)[1])
  }

  return(output)
}
