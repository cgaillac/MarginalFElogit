#' This function computes what corresponds to the equivalent of the function under(over)lined h in Section 3.2.2 in DDL
#' for the ATE (replacing x_kT by x_kT^0 or x_kT^1 depending on the value of x_kT)
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
#' @param XT0 a matrix of size nx dimX containing either X_kT^0 or X_kT^1 depending on the value of X_kT for the selected variable Xk
#' @param XT  a matrix of size nx dimX containing the values of X_T at T defined in Tall
#' @param YT  a vector of size nx 1 containing the values of Y_T at T defined in Tall
#' @param selectX the selected covariate to compute the ATE (one by one)
#' @param Tall vector of size n x1 containing the specified values of T where we compute the AME.
#'
#' @return return a vector containing:
#'
#'     - the average of the functions underline h and overline h over all individuals
#'
#'     - a matrix of size n times 2 containing the individual values of underline h and overline h
#' @export
#'
# @examples
Delta_bounds_ATE <- function(beta1,  V, S, dimX, PSt_X, f, h, mat_C_k_T, RK, T_minus_t, XT0, XT, YT, selectX,Tall){

# ! modif dim(X) > 1
n = dim(V)[1]
T = dim(V)[2]

# if(length(dim(X))==3){
#   k = dim(X)[3]
# }else{
#   k=1
# }

grid0 = sort(unique(Tall))
# VT = matrix(1, length(Tall),1)
# for(ti in grid0){
#   VT[Tall==ti,1] <- V[Tall==ti,ti]
# }

Vtilde =  V/(matrix(exp(XT0%*%beta1),n,1)%*%rep(1,dim(V)[2]))
# Vtilde_T = V/(VT%*%matrix(1,1,dim(V)[2]));

Vtilde_min1 = Vtilde -1

Mat_lambda=zeros(n,T+2);

for(t in 0:T){
  Mat_lambda[,t+2]=C_S_fun(t,Vtilde_min1);
}

lambda_T_plus1 = matrix(Mat_lambda[,T+2], n,1)

# C_S_vec = C_S_fun(S,Vtilde);
# S_minus_t = matrix(S)%*%matrix(1,1,T+1) -   matrix(1,dim(matrix(S))[1],1)%*%(0:T);
# mat_combin = choose(T_minus_t,S_minus_t)/repmat(C_S_vec,1,T+1);
# Mat_fin = Mat_lambda*mat_combin;
# moy_fin = (T+1) * mean(Mat_fin);
# term1 = beta1 * moy_fin;

C_mat = zeros(n,T+1);
for (t in 0:T){
  C_mat[,t+1] = C_S_fun(t,Vtilde);
}

mat_ratio = PSt_X/ C_mat;

c_hat = mat_ratio %*% mat_C_k_T;
c0 = c_hat[,1];
Mat = matrix(Mat_lambda[,2:(T+1)]*c_hat[,2:dim(c_hat)[2]], n , T)

term1= matrix(rowSums(Mat))
term1 =  mean((1-XT[,selectX])*term1 - XT[,selectX]*term1)

m_hat = matrix(c_hat[,2:dim(c_hat)[2]]/(matrix(c0)%*%rep(1,dim(c_hat)[2]-1)),dim(c_hat)[1],dim(c_hat)[2]-1)

# c0bis = choose(T,S)/C_S_vec;
bounds_ind = zeros(n,2);

# mat_C_k_T = mat_C_k_T[[T]]
#i=14

temp = rowSums(matrix(t(apply(m_hat,1,diff)),dim(m_hat)[1],dim(m_hat)[2]-1))
temp0 = rowSums(abs(m_hat))

for (i in 1:n){
    # % Truncation of m_hat. We also set indic=0 if m_hat was truncated,
    # % indic=1 if the last component of m_hat was above the upper bound,
    # % indic=-1 if it was below the lower bound.

    # % First, we estimate the asymptotic variance of mhat(x)=> Sigma_m
  # % First, we estimate the asymptotic variance of mhat(x)=> Sigma_m
  # % First, we estimate the asymptotic variance of mhat(x)=> Sigma_m
  mi = matrix(m_hat[i,], 1,dim(m_hat)[2]);
  i0 = i
  # # abs(temp)[i]
  # m_hat0 = abs(temp)<=10^(-4)
  # if(m_hat0[i]){
  #         m_hat0[i] <- TRUE
  #         values = rowSums((m_hat-matrix(1,dim(m_hat)[1],1)%*%mi)^2)
  #         nonconstant <- cbind(values,!m_hat0, 1:dim(m_hat)[1])
  #         nonconstant0 <- nonconstant[nonconstant[,2]==1,]
  #         ref = nonconstant0 [which.min(nonconstant0 [,1]),]
  #         i0 = ref[3]
  #         mi = matrix(m_hat[i0,], 1,dim(m_hat)[2])
  # }

  H_x = cbind(-t(mi), eye(T))/c0[i0];
  M_x = t(mat_C_k_T)/(matrix(rep(1,dim(mat_C_k_T)[2]))%*%C_mat[i0,]);
  bbx = diag(PSt_X[i0,])- matrix(PSt_X[i0,])%*% PSt_X[i0,]
  Sigma_x = (RK^(dimX*T)/f[i0]) * (bbx  / (matrix(n*h)%*%rep(1,dim(h)[2])));
  Sigma_m = H_x%*% M_x %*% Sigma_x %*% t(M_x) %*% t(H_x)
  # Sigma_m = Sigma_m + eps0*diag(1,c(dim(Sigma_m)));

  out1 <- tryCatch({
    out = trunc_m(matrix(mi), Sigma_m, n);
    if(sum(abs( out[[1]]))!=0){
      mi = out[[1]]
    }
    indic = out[[2]]
  },error=function(cond){
    ### find closest nonconstant neighbour of mi
    mi = matrix(m_hat[i,], 1,dim(m_hat)[2]);

    # temp[i]
    m_hat0 = abs(temp)<=10^(-2) | abs(temp0)<=10^(-1)
    m_hat0[i] <- TRUE
    m_hat0[i0] <- TRUE

    values = rowSums((m_hat-matrix(1,dim(m_hat)[1],1)%*%mi)^2)
    nonconstant <- cbind(values,!m_hat0, 1:dim(m_hat)[1])
    nonconstant0 <- nonconstant[nonconstant[,2]==1,]
    ref = nonconstant0 [which.min(nonconstant0 [,1]),]
    i0 = ref[3]
    mi = matrix(m_hat[i0,], 1,dim(m_hat)[2])
    eps0 = 10^(-9)
    mi[abs(mi) <eps0 ]<- 0

    H_x = cbind(-t(mi), eye(T))/c0[i0];
    M_x = t(mat_C_k_T)/(matrix(rep(1,dim(mat_C_k_T)[2]))%*%C_mat[i0,]);
    bbx = diag(PSt_X[i0,])- matrix(PSt_X[i0,])%*% PSt_X[i0,]
    Sigma_x = (RK^(dimX*T)/f[i0]) * (bbx  / (matrix(n*h)%*%rep(1,dim(h)[2])));
    Sigma_m = H_x%*% M_x %*% Sigma_x %*% t(M_x) %*% t(H_x)

    out = trunc_m(matrix(mi), Sigma_m, n);
    if(sum(abs( out[[1]]))!=0){
      mi = out[[1]]
    }
    indic = out[[2]]
  })

  # qinf_sup = matrix(c(0,0),2,1);
  out3 <- tryCatch({
    qinf_sup = bound_mom(mi,T,indic)
  },error=function(cond){
    mi = matrix(m_hat[i,], 1,dim(m_hat)[2]);

    # temp[i]
    m_hat0 = abs(temp)<=10^(-1) | abs(temp0)<=10^(-1)
    m_hat0[i] <- TRUE
    m_hat0[i0] <- TRUE

    values = rowSums((m_hat-matrix(1,dim(m_hat)[1],1)%*%mi)^2)
    nonconstant <- cbind(values,!m_hat0, 1:dim(m_hat)[1])
    nonconstant0 <- nonconstant[nonconstant[,2]==1,]
    ref = nonconstant0 [which.min(nonconstant0 [,1]),]
    i0 = ref[3]
    mi = matrix(m_hat[i0,], 1,dim(m_hat)[2])
    eps0 = 10^(-9)
    mi[abs(mi) <eps0 ]<- 0

    H_x = cbind(-t(mi), eye(T))/c0[i0];
    M_x = t(mat_C_k_T)/(matrix(rep(1,dim(mat_C_k_T)[2]))%*%C_mat[i0,]);
    bbx = diag(PSt_X[i0,])- matrix(PSt_X[i0,])%*% PSt_X[i0,]
    Sigma_x = (RK^(dimX*T)/f[i0]) * (bbx  / (matrix(n*h)%*%rep(1,dim(h)[2])));
    Sigma_m = H_x%*% M_x %*% Sigma_x %*% t(M_x) %*% t(H_x)

    out = trunc_m(matrix(mi), Sigma_m, n);
    if(sum(abs( out[[1]]))!=0){
      mi = out[[1]]
    }
    indic = out[[2]]

    # qinf_sup = matrix(c(0,0),2,1);
    # qinf_sup = bound_mom(mi,T,indic)
    out2 <- tryCatch({
      qinf_sup = bound_mom(mi,T,indic)
    },error=function(cond){
      # qinf_sup = matrix(c(0,0),2,1);
    })

  })




    term2 = c0[i] %*% lambda_T_plus1[i] %*% t(qinf_sup);

    bounds_ind[i,1] = term1 + min(term2);
    bounds_ind[i,2] = term1 + max(term2);

}

# mean( Y[,T]*(2*XT-1)) + term1 + mean(XT==0)*  colMeans(    true_bounds0[XT==0,]) - mean(XT==1)* rev(colMeans(   true_bounds0[XT==1,])) ;
bounds =  mean( YT*(2*XT[,selectX]-1))+ term1 + mean(XT[,selectX]==0)*  colMeans(     bounds_ind[XT[,selectX]==0,]) - mean(XT[,selectX]==1)* rev(colMeans(   bounds_ind[XT[,selectX]==1,])) ;
# bounds = apply(bounds_ind,2,mean);

out = vector("list")
out[[1]] <-  bounds
out[[2]]  <- bounds_ind
return(out)
}
