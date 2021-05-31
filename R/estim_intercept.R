#' Estimate the intercept of a Logit model of Y on X_T with the
#' beta estimated by conditional logit
#'
#' @param Y a matrix of size nxT to compute the intercept
#' @param Xbeta a matrix of size nxT containing bhat'x_t where bhat is the estimated value of beta0 using the conditional maximum likelihood estimator.
#'
#' @return the intercept of a Logit model of Y on X_T with the
#' beta estimated by conditional logit
#' @export
#'
# @examples
estim_intercept <- function(Y,Xbeta){

tol=1e-4;
pas = 1;
PY = mean(Y);
res = log(PY/(1-PY)) - mean(Xbeta);
iter=1;

while ((pas > tol) && (iter < 50)){
  prev = res ;
  res = prev - (mean(Lambda(prev+Xbeta)) - PY)/mean(Lambda(prev+Xbeta,1));
  pas = abs(res-prev);
  iter=iter+1;
}

return(res)
}
