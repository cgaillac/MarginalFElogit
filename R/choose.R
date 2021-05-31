#' Compute C(n,k) with n and k scalars
#' Allows for negative k or n<k (then returns 0)
#'
#' @param n scalar in the combinatorial numers C(n,k)
#' @param k scalar in the combinatorial numers C(n,k)
#'
#' @return the value of C(n,k)
#' @export
#'
# @examples
choose <- function(n,k){
  num = factorial(n)
  denom1 = factorial(abs(k))
  denom2 = factorial(abs(n-k))

  res = (num/(denom1*denom2))*(k<=n & k>=0)
return(res)
}

