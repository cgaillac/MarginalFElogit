#' Computes the coefficients of T_k(2u-1), where T_k is the usual Chebyshev
#' polynomial of the 1st kind
#'
#' @param k parameter k in T_k(2u-1), where T_k is the usual Chebyshev
#'
#' @return the coefficients of T_k(2u-1)
#' @export
#'
# @examples
coeff_Cheb01 <- function(k){

  nodes = Re(polyroot(ChebyshevPoly(k)))

  new_nodes = matrix((nodes+1)/2);

  temp=zeros(k+1,1);

  for(i in seq(0,k)){
    if(i==0){
      temp[i+1]= 1
    }else{
      if(k==2 & i==1){
        temp[i+1]=(-1)^i * sum(apply(new_nodes, 1,prod));
      }else{
        temp[i+1]=(-1)^i * sum(apply(t(combnk(new_nodes,i)), 1,prod));
      }
    }
  }

return(2^k * temp)
}
