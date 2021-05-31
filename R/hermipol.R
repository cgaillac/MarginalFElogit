#' Function conding the Hermite polynomial of degree n
#'
#' @param n degree of the Hermite polynomial.
#'
#' @return the coefficients of the Hermite polynomial of degree n.
#' @export
#'
# @examples
hermipol <- function(n){

  p= matrix(0,n+1,n+1)
  p[1,1]=1;
  p[2,1:2]=c(2,0);

  for( k in seq(2,n)){
    p[k+1,1:(k+1)]=2*c(p[k,1:k], 0)-2*(k-1)*c(0,0,p[k-1,1:(k-1)])
  }

  return(p)
}
