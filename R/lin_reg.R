#' Evaluate the mean square error loss function of linear regression by pooling, to obtain a starting value for the conditional maximum likelihood estimator.
#'
#' @param beta the value at which the loss is computed
#' @param Y a matrix of size nxT containing the values of the dependent variable Y
#' @param X an array of size nxTxdimX containing the values of the covariates Xs
#'
#' @return the value of the mean square error loss function of linear regression at beta
#' @export
#'
# @examples
lin_reg <- function(beta,Y, X){

  # [n, T, k]=size(X);
  ## change for dim(X) >1
  n <- dim(X)[1]
  T <- dim(X)[2]
  if(length(dim(X))==3){
    k = dim(X)[3]
  }else{
    k=1
  }
  # beta=0.2
  index=zeros(n,T);
  # if(k ==1){
  #   index = index + X*beta[1];
  # }else{
  for( j in 1:k){
    index = index + X[,,j]*beta[j];
  }
  # }
  #
  # Tall =  apply(Y,1,isnot)
  # grid_T = NULL
  # Tmax = dim(Y)[2]
  # # n_dist = NULL
  # ### find the distibution of T in the population and sample size
  # for(t in 1:Tmax){
  #   if(  sum(Tall==t)>0){
  #     grid_T = c(  grid_T , t)
  #     # n_dist = c( n_dist,sum(Tall==t))
  #   }
  # }
  #
  # V = exp(index);
  # S = apply(Y,1,sumNA);
  # denom = matrix(0,n,1)
  # for(ti in grid_T){
  #   denom[Tall==ti,1] = C_S_fun(S[Tall==ti],V[Tall==ti,1:ti]);
  # }



  # % Careful, we actually take - the log, as we minimize the function at the
  # % end.
  f = sum(apply((Y - index)^2, 1,sumNA));

  return(f)
}
