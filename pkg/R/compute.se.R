#' @title Standard error of mediating variable
#'
#' @description Compute standard error
#'
#' @param alpha vector of length N
#' @param beta vector of length N
#'
#' @export
#'
compute.se <- function(alpha, beta, sigma_alpha, sigma_beta) {
  se.sobel <- sqrt((alpha*sigma_beta)^2 + (beta*sigma_alpha)^2)
  return(list(se.sobel=se.sobel))
}

