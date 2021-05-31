#' Count the non NA values of a vector
#'
#' @param x a vector, potentially containing NA
#'
#' @return the number of non NA coefficients of x
#' @export
#'
# @examples
isnot <- function(x){

  return(sum(!is.na(x)))

  }
