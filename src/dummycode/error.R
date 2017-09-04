# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

###

pc <- pca(miss_data, nPcs=3, method="ppca")
imputed <- completeObs(pc)

data1
#' Title RMSE_simulated
#' Calculates the root mean square error between simulated(reference) data and imputed data
#' @param data1 = simulated data matrix
#' @param data2 = data matrix with missing values
#' @param data3 = imputed data matrix
#'
#' @return
#' @export rmse
#'
#' @examples
RMSE_simulated <- function(data1,data2,data3){
sum((data1[is.na(data2)] - data3[is.na(data2)])^2) / sum(data1[is.na(data2)]^2)

}