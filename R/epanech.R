#' The Epanechnikov kernel
#'
#' @param u variable of the kernel
#'
#' @return the value of the kernel at u
#' @export
#'
# @examples
epanech <- function(u){

res = 0.75*pmax(1 - u^2,0);

return(res)
}
