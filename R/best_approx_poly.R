#' This function returns the coefficients (functions of v=exp(x'beta0)) of
#' the best (in the uniform sense) approximation  of Omega(.,x) by a
#' polynomial of degree T.
#'
#' @param Vtilde_min1 is the vector of coefficients exp((x-x_T)'beta0)-1 appearing
#' in Omega
#' @param grid_T1 a vector containing the different values of the maximal number of periods observed
#' @param Tinf0 a vector of size nx1 containing the number of periods observed for each individual
#'
#' @return  the coefficients of the best approximation of Omega(.,x) by polynomial of degree T.
#' @export
#'
# @examples
best_approx_poly <- function (Vtilde_min1, grid_T1, Tinf0 ){

n = dim(Vtilde_min1)[1]
Tmax = dim(Vtilde_min1)[2]
coeffs_Omega = matrix(NA,n,Tmax+3);

# i=1
for (i in 1:n){
    T =Tinf0[i]
    if (T>1){
    v = Vtilde_min1[i,1:(T-1)]
    }
    # % Coefficients of the polynomial corresponding to the product in Omega(.,x)
    temp=zeros(T,1);
    temp[T]=1
    # j=0
    if(T==2){
      temp[1]=v;
    }else if(T>2){
      for (j in 0:(T-2)){
          temp[j+1]=sum(t(apply(combnk(v,j+1),1,prod)));
      }
    }
    # % Coefficients of Omega(.,x) (in decreasing order of exponents)
    coeffs_Omega[i,1:(T+2)]=conv(c(-1,1),c(temp,0));
}


cheb_coeff=vector("list")
for(t in 1:length(grid_T1)){
  T = grid_T1[t]
  # cheb_coeff[[T]]  = fliplr(t(coeff_Cheb01(T+1)));
  cheb_coeff[[T]]  = coeff_Cheb01(T+1);
}

res = matrix(NA,dim( coeffs_Omega)[1],dim(coeffs_Omega)[2]-1)
for(ti in grid_T1){
  res[Tinf0==ti,1:(ti+1)] = coeffs_Omega[Tinf0==ti,2:(ti+2)] - matrix(coeffs_Omega[Tinf0==ti,1])%*%cheb_coeff[[ti]][2:length(cheb_coeff[[ti]])]/(2^(ti+1));
}

# dim(coeffs_Omega[,1]*cheb_coeff[2:length(cheb_coeff)]/(2^(T+1)))

# x11()
# i=5
# x0= seq(0,1,length.out=100)
# plot(x0 , polyval(coeffs_Omega[i,], x0), type="l")
# lines(x0 , polyval(res[i,], x0), col=2,lty="dotted")

# plot(x0 , polyval(rev(res[i,]), x0),type="l")

g = coeffs_Omega[,1];

out = vector("list")
out[[1]] <- res
out[[2]] <- g
return(out)
}
