#' Evaluate the function Omega in DDL (see Section 4.1), with v the vector of v_t for t
#' less than T0
#'
#' @param u the variable u in Omega, which is a polynomial in u
#' @param v v the vector of v_t  (often equals to exp(x_t'beta) or exp((x_t-x_T0)'beta)) for t less than T0
#' @param T0 T0 +1 is the degree in u of the polynomial Omega
#'
#' @return the value of Omega at u,v
#' @export
#'
# @examples
omega <- function(u,v,T0){
  aux = u*(v[1]-1) + 1;
  if(T0>2){
    for (t in 2:(T0-1)){
      aux = aux*(u*(v[t]-1)+1);
    }
  }
  aux = u*(1-u)*aux/(u*(v[T0]-1)+1)
  return(aux)
}
