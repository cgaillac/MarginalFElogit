#' Evaluate the derivative of logit function
#'
#' @param x value at which the derivative of logit function is evaluated
#'
#' @return the value of the derivative of logit function at x
#' @export
#'
# @examples
Lambda_prime <- function( x ){

return(exp(x)/(1+exp(x))^2)
}

