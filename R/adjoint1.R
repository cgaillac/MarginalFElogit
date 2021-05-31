#' Compute the adjoint of a square matrix. The adjoint is the transpose of the cofactor matrix
#'
#' @param A a square matrix
#'
#' @return The value matrix of A
#' @export
#'
# @examples
adjoint1<-function(A){
  return(det(A)*solve(A))
}
