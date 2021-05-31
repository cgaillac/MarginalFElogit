#' Estimates E(Y|X) at the points Xb by local linear regression
#' with bandwidth h1 and the Epanechnikov kernel.
#' Second output f: Kernel estimator of the density
#'
#' @param Yb a matrix of size nxT containing the values of the dependent variable Y.
#' @param Xb an array of size nxTxdimX containing the values of the covariates Xs.
#' @param h1 bandwidth for the local linear regression.
#'
#' @return a list containing
#'
#'  - res: the estimates E(Y|X) at the points Xb by local linear regression;
#'
#'  - f: Kernel estimator of the density.
#'
#' @export
#'
# @examples
local_lin <- function (Yb,Xb,h1){

  T <- dim(Xb)[2]
  if(length(dim(Xb))==3){
    dimX = dim(Xb)[3]
  }else{
    dimX=1
  }

  n = dim(Yb)[1];
  res=zeros(n,1);

  if(length(dim(Xb))==3){
    X1 = cbind(ones(n,1),Xb[,,1])
  if(dimX >1){
    for(k in 2:dimX){
      X1 = cbind(X1,Xb[,,2])
    }
  }
  }else{
    X1 =  cbind(ones(n,1),Xb)
  }

  Xbb = matrix(X1[,-c(1)], n, dim(X1)[2]-1)
  p = dim(X1)[2];


  f = zeros(n,1);
  denom_f = n*(sqrt(2*pi)*h1)^p;

  # i=3
  for (i in 1:n){
    x =  Xbb [i,];
    bx = matrix(rep(1,dim(Xbb)[1]))%*%x - Xbb
    w = apply(exp(-0.5*(bx/h1)^2),1,prod);

    bx = cbind(ones(n,1),bx)

    f[i]=sum(w)/denom_f;

    numer = t(X1)%*%(w*Yb);
    denom = t(X1)%*% (X1 *(matrix(w)%*%rep(1,dim(X1)[2])));
    if(det(denom)<=10^(-15)){
     b = ginv(denom)%*%numer;
    }else{
     b = solve(denom,numer);
    }
    res[i] = b[1];
  }


out = vector("list")
out[[1]] <- res
out[[2]] <- f
return(out)
}
