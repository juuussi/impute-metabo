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


#' Title compare_results
#' Calculates the correlation between two matrices/ we use it to evaluate the different 
#' impuatation methods
#' @param x data matrix
#' @param y data matrix
#'
#' @return data frame contains results from the correlation and p-values
#' @export
#'
#' @examples compared_results <- compare_results(simulated_data, imputed_data)
compare_results <- function(x, y) {
  c <- cor.test(as.vector(x), as.vector(y))
  data.frame(Cor=c$estimate, CorP=c$p.value)
  
  
}


compare_results2 <- function(x, y, missing) {
  missing_index <- which(is.na(as.vector(missing)))
  x.val <- as.vector(x)[missing_index]
  y.val <- as.vector(y)[missing_index]
  flog.info(paste('missing:', paste(missing_index, collapse=", "), '\norig:', paste(x.val, collapse=", "),'\nimpu:', paste(y.val, collapse=", "),sep=' '))
  
  c <- cor.test(x.val, y.val)
  
  flog.info(paste("Cor:", c$estimate, "P:", c$p.value, sep=" "))
  data.frame(Cor=c$estimate, CorP=c$p.value)
  
  
}




#' Title simulate_missingness
#' Uses the simulated data to create the different types of missingness (MCAR,MNAR,MAR) 
#' using different percentages of missigness
#' @param data matrix with simulated data
#' @param mcar percenatge of missigness in Missing Completly at Random 
#' @param mar percenatge of missigness in Missing at Random 
#' @param mnar percenatge of missignessMisiing not at Random 
#' @param mnar.type type of trancation 'left' or right
#'
#' @return simulated_data
#' @export
#'
#' @examples miss_data <- simulate_missingness(data=simulated_data, mcar=0.01)
simulate_missingness <- function(data, mcar=0, mar=0, mnar=0, mnar.type="left") {
  
  if(class(data) != "matrix") {
    stop("Variable data should be a matrix.")
  }
  simulated_data <- data
  if (mcar > 0){
    mcar_distribution   = runif(nrow(data)*ncol(data), min=0, max=1)
    simulated_data = matrix(ifelse(mcar_distribution<mcar, NA, data), nrow=nrow(data), ncol=ncol(data))
  }
  if (mnar > 0) {
      added_mnar <- 0
      while (added_mnar < mnar) {
        
        # Select random variable
        variable_index <- sample(1:ncol(simulated_data), 1)
        
        # What percentage of variable to set missing
        cut_percentage <- rchisq(1, df=1) / 30
        
        if (cut_percentage > 1) {
          cut_percentage <- 1
        }
        while (cut_percentage == 0) {
          cut_percentage <- rchisq(1, df=1) / 30
          
        }
        
        # How many values to set missing
        cut_index <- floor(cut_percentage * nrow(simulated_data))
        
        sorted_variable <- sort(simulated_data[,variable_index])
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[cut_index]
        
        
        # Set values to missing
        simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
        
        # Counter to check how much MNAR missingness has been added to data
        added_mnar <- added_mnar + (cut_index-1) / (nrow(simulated_data) * ncol(simulated_data))
      }
  }
  
  simulated_data
  
}



#' Title impute
#' This functions contains different imputation methods and imputes the data with all
#' the different imputation methods 
#' @param data data matrix with simulated data
#' @param methods vector containing the method names 
#'
#' @return the matrix containing the imputed data
#' @export
#'
#' @examples imputation_methods <- c("RF", "mean", "min")
#' imputed_data <- impute(data=miss_data, methods=imputation_methods)
impute <- function(data, methods) {
  
  require(missForest)
  
  if (length(methods) != 1 & length(methods) != ncol(data)) {
    stop("Methods needs to be either one value or of the same length as number of columns in data.")
  }
  
  imputed_data <- data
  
  if (length(methods) == 1) {
    methods <- rep(methods, times=ncol(data))
  }
  
  
  if ("RF" %in% methods) {
    rf_imputed_data <- missForest::missForest(xmis = data)$ximp
    index <- which(methods == "RF")
    imputed_data[,index] <- rf_imputed_data[,index]
  }
  
   if ("PPCA" %in% methods) {
     # Do cross validation with ppca for component 2:10
     esti <- kEstimate(data, method = "ppca", evalPcs = 2:10, nruncv=1, em="nrmsep")
     # The best result was obtained for this number of PCs:esti$bestNPcs
     pc <- pcaMethods::pca(data,nPcs=esti$bestNPcs, method="ppca")
     #index <- which(methods == "PPCA")
     imputed_data <- completeObs(pc)
     
    
  }
  
  foreach (data_column=1:ncol(data)) %do% {
    method <- methods[data_column]
    
    if (method == "mean") {
      impu_value <- mean(data[,data_column], na.rm=TRUE)
      imputed_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }
    
    if (method == "min") {
      impu_value <- min(data[,data_column], na.rm=TRUE)
      imputed_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
      
    }
  }
  
  imputed_data
  
}



#' Title MSE_simulated
#' Calculates the  mean square error between simulated(reference) data and imputed data
#' @param data1 = simulated data matrix
#' @param data2 = data matrix with missing values
#' @param data3 = imputed data matrix
#'
#' @return
#' @export error(mse)
#'
#' @examples
RMSE_simulated <- function(data1,data2,data3){
  error <- sqrt(sum((data1[is.na(data2)] - data3[is.na(data2)])^2) /nrow(data1))
  error
}

