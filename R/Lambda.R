#' Evaluate the logit function and its derivative
#'
#' @param x value at which the logit is evaluated
#' @param deriv if not null, evaluate the derivative at x instead of the function
#'
#' @return the value of the logit function at x.
#' @export
#'
# @examples
Lambda <- function(x,deriv=NULL){

res=1./(1+exp(-x));

if( !is.null(deriv) && (deriv==1)){
#derivative of Lambda for computation of the AME
res = exp(-x)*(res^2);
}

return(res)
}



