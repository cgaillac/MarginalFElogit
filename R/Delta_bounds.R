#' This function computes what corresponds to the function under(over)lined h in Section 3.2.2 in DDL
#'
#' @param beta1 value of the variable beta of the function h (4th parameter in DDL)
#' @param V matrix of v = exp(x_t'betahat), where betahat is the estimated value of beta0 using the conditional maximum likelihood estimator.
#' @param S matrix of size nx 1 of values of the statistic S
#' @param dimX the number of covariates
#' @param PSt_X matrix of size n x T containing the estimated conditional probabilities of S given X (gamma in DDL) for all individuals
#' @param f vector of size nx 1 containing the estimated values of the density of X evaluated at the observed individual values of X
#' @param h vector containing the bandwidths for the nonparametric estimator of the conditional probabilities of S given X
#' @param mat_C_k_T matrix containing the combinatorial numbers (T-t, S-t) for all values of t.
#' @param RK integral of the kernel squared for the computation of the nonparametric estimator of the conditional probabilities of S given X
#' @param T_minus_t matrix of all values of T-t
#' @param Tall vector of size n x1 containing the specified values of T where we compute the AME.
#' @param X matrix of containing the values of the covariates
#'
#' @return return a vector containing:
#'
#'     - the average of the functions underline h and overline h over all individuals
#'
#'     - a matrix of size n times 2 containing the individual values of underline h and overline h
#'
#' @export
#'
# @examples
Delta_bounds <- function(beta1,  V, S, dimX, PSt_X, f, h, mat_C_k_T, RK, T_minus_t, Tall,X){
# mat_C_k_T= mat_C_k_T[[T]]
# Tall= Tall0
# beta1 = b_hat
# PSt_X = PSt_X+hderiv*ones(n,1)%*%((0:T)==t)
# ! modif dim(X) > 1
n = dim(V)[1]
T = dim(V)[2]

# if(length(dim(X))==3){
#   k = dim(X)[3]
# }else{
#   k=1
# }
grid0 = sort(unique(Tall))
VT = matrix(1, length(Tall),1)
for(ti in grid0){
  VT[Tall==ti,1] <- V[Tall==ti,ti]
}


Vtilde = V/(VT%*%matrix(1,1,dim(V)[2]));
Vtilde_min1 = cbind(Vtilde - matrix(1,dim(Vtilde)[1],dim(Vtilde)[2]), -ones(n,1))

Mat_lambda=zeros(n,T+1);

t=1
for(t in 1:T){
  Mat_lambda[,t+1]=C_S_fun(t-1,Vtilde_min1);
}

if(T==1){
  lambda_T_plus1 =- matrix(1,n,1)
}else{
  lambda_T_plus1 =- apply(matrix(Vtilde_min1[,1:(dim(Vtilde_min1)[2]-2)],dim(Vtilde_min1)[1],(dim(Vtilde_min1)[2]-2)),1,prod);
}


C_S_vec = C_S_fun(S,Vtilde);
S_minus_t = matrix(S)%*%matrix(1,1,T+1) -   matrix(1,dim(matrix(S))[1],1)%*%(0:T);
mat_combin = choose(T_minus_t,S_minus_t)/repmat(C_S_vec,1,T+1);
Mat_fin = Mat_lambda*mat_combin;
moy_fin = (T+1) * mean(Mat_fin);

term1 = beta1 * moy_fin;

C_mat = zeros(n,T+1);
for (t in 0:T){
  C_mat[,t+1] = C_S_fun(t,Vtilde);
}

mat_ratio = PSt_X/ C_mat;

# mat_C_k_T <- mat_C_k_T[[T]]
c_hat = mat_ratio %*% mat_C_k_T;
c0 = c_hat[,1];
m_hat = c_hat[,2:dim(c_hat)[2]]/(matrix(c0)%*%rep(1,dim(c_hat)[2]-1))


bounds_ind = zeros(n,2);
res0 <- sfLapply(1:n,boot_slow,m_hat,mat_C_k_T,RK,n,dimX,T, h ,PSt_X,c0,C_mat, f, lambda_T_plus1,X)

# all0 = cbind(m_hat0 = abs(rowSums(m_hat))<=10^(-4), m_hat)
# all0
# sum(all0[,1])
# i =96
# res0 = zeros(n,1)
# res0 = zeros(n,2)
# for (i in 1:n){
#    res0[i,]<- boot_slow(i,m_hat,mat_C_k_T,RK,n,dimX,T, h ,PSt_X,c0,C_mat, f, lambda_T_plus1)
#   # res0[i]<-test_sigma(i,m_hat,mat_C_k_T,RK,n,dimX,T, h ,PSt_X,c0,C_mat, f, lambda_T_plus1)
# }
# res0[90,]
# sum(res0)
bounds_ind = zeros(n,2);
for(i in 1:n){
  bounds_ind[i,1] <- res0[[i]][1]
  bounds_ind[i,2] <- res0[[i]][2]
  # bounds_ind[i,1] <- res0[i,1]
  # bounds_ind[i,2] <- res0[i,2]
}

bounds = apply(bounds_ind,2,mean);
bounds = matrix(term1,dimX,1)%*%rep(1,2)+ matrix(beta1,dimX,1)%*%bounds


out = vector("list")
out[[1]] <- bounds
out[[2]]  <- bounds_ind
out[[3]]  <- c_hat
return(out)
}
