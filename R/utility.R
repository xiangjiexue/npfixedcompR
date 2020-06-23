# other utility functions

#' Extract or returning the lower triangular part of the matrix
#'
#' The function \code{extractlower} is to extract the strict
#' lower triangular part of a squared matrix and the function
#' \code{returnlower} is to return the vector value into a
#' symmetric matrix with diagonal 1.
#'
#' @param A a matrix to be extracted the lower triangular part
#' @param v a vector to be returned to a symmetric matrix with diagonal 1.
#' @examples
#' a = matrix(1:100, 10, 10)
#' b = extractlower(a)
#' d = returnlower(b)
#' @rdname extractlower
#' @export
extractlower = function(A){
  A[lower.tri(A), drop = TRUE]
}

#' @rdname extractlower
#' @export
returnlower = function(v){
  LLL = (1 + sqrt(1 + 8 * length(v))) / 2
  ans = matrix(0, LLL, LLL)
  ans[lower.tri(ans)] = v
  ans + t(ans) + diag(LLL)
}
