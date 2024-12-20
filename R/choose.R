#' Takes input n, k and returns n choose k
#'
#' @param n number of elements to choose from
#' @param k number of elements to choose
#'
#' @return n choose k
#' @export
#'
# @examples
choose <- function(n, k) {

  # We use n choose k = n! / (k! (n-k)!)
  k <- min(k, n - k)
  if ((n < 0) | (k < 0) | (k > n)) {
    stop("Error : n choose k requires k between 0 and n and n positive")
  }
  else if (k == 0) {
    nchoosek <- 1
  }
  else {
    nchoosek <- prod(((n - k + 1):n)/(1:k))
    nchoosek <- round(nchoosek) # In case numerical computation does not give an exact integer
  }

  return(nchoosek)
}

