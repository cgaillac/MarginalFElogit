#' Compute the bounds on q_T(x) for each individual value of x for the sharp method in DDL.
#' Allows to use parallel computing.
#'
#' @param i indexes the individual
#' @param m_hat the vector of estimated moments up to T
#' @param mat_C_k_T a matrix containing the combinatorial numbers
#' @param RK the integral of the kernel K^2 intervening in the estimator of the asymptotic variance of the estimator of gamma
#' @param n the sample size
#' @param dimX the number of covariates
#' @param T the number of observed periods for this indvidual
#' @param h the vector of bandwidths for the nonparametric estimator of gamma
#' @param PSt_X the estimated matrix of size n x T of estimated gamma0t(x)
#' @param c0 the vector containing the first component of the moments c_t(x) for all individuals
#' @param C_mat the matrix containing the S-elementary polynomials C_S
#' @param f  an estimator of the density of X at x.
#' @param lambda_T_plus1 the vector containing the coefficients T+1 of lambda for all individuals
#' @param X matrix of containing the values of the covariates
#'
#' @return bounds on q_T(x) for ith individual
#' @export
#'
# @examples
boot_slow <- function(i,m_hat,mat_C_k_T,RK,n,dimX,T, h ,PSt_X,c0,C_mat, f, lambda_T_plus1,X){
  # % Truncation of m_hat. We also set indic=0 if m_hat was truncated,
  # % indic=1 if the last component of m_hat was above the upper bound,
  # % indic=-1 if it was below the lower bound.

  # % First, we estimate the asymptotic variance of mhat(x)=> Sigma_m
  mi = matrix(m_hat[i,], 1,dim(m_hat)[2]);
  i0 = i

  # if(sum(abs(mi)==0)){
  #   m_hat0 = abs(rowSums(m_hat))==0
  #   # sum(m_hat0)
  #   values = matrix(0,dim(X)[1], 1)
  #   for(k in 1:dimX){
  #     values =  values + rowSums((X[,,k]-matrix(1,dim(X)[1],1)%*%X[i,,k])^2)
  #   }
  #
  #   nonconstant  <- cbind(values,!m_hat0, 1:dim(PSt_X)[1])
  #   nonconstant0 <-  nonconstant[ nonconstant[,2]==1,]
  #   ref =  nonconstant0 [which.min( nonconstant0 [,1]),]
  #   i0 = ref[3]
  #   # X[i0,,1]
  #   mi = matrix(m_hat[i0,], 1,dim(m_hat)[2])
  # }

  H_x = cbind(-t(mi), eye(T))/c0[i0];
  M_x = t(mat_C_k_T)/(matrix(rep(1,dim(mat_C_k_T)[2]))%*%C_mat[i0,]);
  bbx = diag(PSt_X[i0,])- matrix(PSt_X[i0,])%*% PSt_X[i0,]
  Sigma_x = (RK^(dimX*T)/f[i0]) * (bbx  / (matrix(n*h)%*%rep(1,dim(h)[2])));
  Sigma_m = H_x%*% M_x %*% Sigma_x %*% t(M_x) %*% t(H_x)

  out = trunc_m(matrix(mi), Sigma_m, n);
  mi = out[[1]]
  indic = out[[2]]


  qinf_sup <- matrix(c(NA,NA),2,1);
  qinf_sup = bound_mom(matrix(mi,length(mi),1),T,indic)

  term2 =  c0[i0] %*% lambda_T_plus1[i0] %*% t(qinf_sup);

  return(c(min(term2),max(term2)))}
