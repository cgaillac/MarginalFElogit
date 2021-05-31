#' Returns weights and nodes of the Gauss-Hermite method
#' Allows for a weight exp(-(x-mu)^2/(2sigma^2)) instead of the classical exp(-x^2)
#'
#' @param n the maximal degree of the polynomials H_n used in the Gauss-Hermite quadrature method
#' @param mu mean centering parameter for the normal weight
#' @param sigma standard deviation parameter for the normal weight
#'
#' @return a list containing:
#'
#'      -  the weights and
#'
#'      -  the nodes of the Gauss-Hermite method
#' @export
#'
# @examples
gausshermi <- function(n,mu,sigma){

#Hermite polynomial
p=hermipol(n);
#Roots
racine=Re(polyroot(rev(p[n+1,])))
#x=mu+sqrt(2)*sigma*roots(p(n+1,:));

w=zeros(n,1);
#Weights
for( i in seq(1,n)){
  w[i]=(2^(n-1)*(factorial(n)))/(n^2*(polyval(p[n,1:n],racine[i]))^2);
}

x = mu+sqrt(2)*sigma*racine;

out = vector("list")
out[[1]] <- w
out[[2]] <- x
return(out)
}



