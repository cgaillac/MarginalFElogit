#' Given nonnegative integer n, compute the Chebyshev polynomial T_n.
#'
#'This code is adapted from ChebyshevPoly.m by David Terr, Raytheon, 5-10-04
#' @param n indexes the Chebyshev polynomial T_n
#'
#' @return  Return the result as a column vector whose mth element is the coefficient of x^(n+1-m).
#' @export
#'
# @examples
ChebyshevPoly <- function(n){

  if( n==0){
    tk = 1
  }else if(n==1){
    tk = matrix(c(0,1)) #[1 0]';
  }else{
      tkm2 = zeros(n+1,1)
      tkm2[n+1] = 1
      tkm1 = zeros(n+1,1)
      tkm1[n] = 1

    for( k in seq(2,n)){

        tk = zeros(n+1,1)

        for(e  in seq(n-k+1,n,by=2)){
              tk[e] = 2*tkm1[e+1] - tkm2[e]
        }

        if( mod(k,2)==0){
              tk[n+1] = (-1)^(k/2)
        }

        if(k<n){
            tkm2 = tkm1
            tkm1 = tk
        }

    }
      tk = rev(tk)
  }
    return(tk)
}

