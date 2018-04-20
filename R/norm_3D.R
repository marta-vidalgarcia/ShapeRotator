#' norm_3D internal function
#'
#' This function allows you to obtain the length of a vector u
#' @param u vector
#' @export
#' @examples
#' norm_3D()


norm_3D <- function (u)
{
  return( sqrt(u[1]^2 + u[2]^2 + u[3]^2));}
