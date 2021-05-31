#' Sum only on non NA values
#'
#' @param x a vector, potentially containing NA
#'
#' @return the sum over the non NA coefficients of x
#' @export
#'
# @examples
sumNA <- function(x){

  return(sum(x[!is.na(x)]))

}
