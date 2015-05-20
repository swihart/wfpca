#' numerical integration based on Simpsons rule 
#' 
#' based on Bolstad:sintegral(); but needed to remove the message.
#' 
#' 
#' @param x a sequence of x values.
#' @param fx the value of the function to be integrated at x or a function
#' @param n.pts  the number of points to be used in the integration. If x contains more than n.pts then n.pts will be set to length(x)
#' @export
#' @return results a data frame with rmse for each approach for that simulated dataset
#' @examples
#' ---
simpson_nint <- function (x, fx, n.pts = max(256, length(n.x))) 
{
  ##message("Note: sintegral's behavior has changed.\nTo get the value of the integral use sintegral(x,fx)$value.\nTo get the cdf use sintegral(x,fx)$cdf")
  if (class(fx) == "function") 
    fx = fx(x)
  n.x = length(x)
  if (n.x != length(fx)) 
    stop("Unequal input vector lengths")
  if (n.pts < 64) 
    n.pts = 64
  ap = approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + 
                    ap$y[2 * (1:n.pts) + 1])/3
  invisible(list(value = sum(integral), cdf = list(x = ap$x[2 * 
                                                              (1:n.pts)], y = cumsum(integral))))
}