#' Function to compute the true bounds for the Delta(x) for each individual value of x
#'  in the simulated examples of DDL using parallel computing.
#'
#' @param i indexes the individual
#' @param m_true a vector of true moments m(x) of size T
#' @param lambda_T_plus1 the value of lambda(x) at T+1
#' @param T the number of observed periods for this indvidual
#' @param c_true a vector of moments c(x) of size T
#'
#' @return the bounds on Delta(x)
#' @export
#'
# @examples
boot_MC <- function(i,m_true,lambda_T_plus1,T,c_true){

  true_bounds = zeros(1,2);

  if(T==1){
    out = trunc_m(matrix(m_true[i]))
  }else{
    out = trunc_m(matrix(m_true[i,]))
  }
  mi= out[[1]]
  indic= out[[2]]

  # if(T==1){
  #
  # }else{
  qinf_sup = bound_mom(mi,T,indic)
  # }
  # second term of Lemma 1.
  # term2 = matrix(matrix(beta0,dimX,1) * c_true[i,1] * lambda_T_plus1[i], dimX,1) %*% t(qinf_sup);
  term2 =  c_true[i,1] %*% lambda_T_plus1[i] %*% t(qinf_sup);

  true_bounds = true_bounds + c(min(Re(term2)), max(Re(term2)));

  return(true_bounds )
}
