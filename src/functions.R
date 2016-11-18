#' Simulate multivariate data based on given reference data
#'
#' Uses \code{mvrnorm} function to simulate data matrix of given dimensions, 
#' using provided data as reference for the multivariate normal distributions.
#'
#' @param data matrix to be used as reference
#' @param nrow the desired number of rows.
#' @param ncol the desired number of columns.
#' @param ...  additional parameters to be passed to the \code{mvrnorm} function.
#'
#' @return data matrix with the simulated data.
#' @export
#'
#' @examples
#' 
simulate_data <- function(data, nrow, ncol, ...) {
  
  require(MASS)
  
  simulated_data <- matrix(mvrnorm(n=nrow*ncol, mu=colMeans(reference_data), Sigma=cov(data), ...), nrow=nrow, ncol=ncol)

  simulated_data
}