#' Computes the quantile suggested by Imbens & Manski 2004 in order to have
#' CI with as. size of 1-alpha. It solves Eq. (2) in Stoye (2009)
#'
#' @param alpha the level for the confidence intervals
#' @param bounds1 the estimated bounds
#' @param std_bounds1  the estimated standard errors for the bounds
#'
#' @return the value of the quantile suggested by Imbens & Manski 2004
#' @export
#'
# @examples
quantile_IM <- function(alpha, bounds1, std_bounds1){

  term = (bounds1[2]-bounds1[1])/max(std_bounds1);
  binf = qnorm(1-alpha);
  bsup = qnorm(1-alpha/2);

  ff <- function(q){return( (pnorm(q+term)-pnorm(-q)-1+alpha)^2)}
  res = optimize(ff,c(binf,bsup))$minimum;

return(res)
}
