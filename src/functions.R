######### Function No1 ##################################
#' Title missingness_proportions
#'equally division of the percentages of missigness between the three types of missigness (MCAR,MNAR,MAR)
#' @param miss_proportions 
#'
#' @return  result.df a dataframe containing the perecenatges per type
#' @export
#'
#' @examples
#' 
#' #'####################################################

missingness_proportions <- function(miss_proportions) {
  result.df <- NULL
  for (prop in miss_proportions) {
    print(prop)
    df <- data.frame(Total=prop, MCAR=round(c(prop,0, 0,prop/2, prop/2, 0, prop/3),digits = 2), MNAR=round(c(0, prop, 0, prop/2, 0, prop/2, prop/3),digits = 2), MAR=round(c(0,0,prop, 0, prop/2, prop/2, prop/3),digits = 2))
    #df <- data.frame(MCAR=c(prop,0, 0,prop/2, prop/2, 0, prop/3), MNAR=c(0, prop, 0, prop/2, 0, prop/2, prop/3), MAR=rep(0, times=7))
    
    if (is.null(result.df)) {
      result.df <- df
    } else {
      result.df <- rbind(result.df, df)
    }
  }
  result.df
  
}

######### Function No2  ##################################

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
#' #'####################################################

simulate_data <- function(data, nrow, ncol, ...) {
  
  require(MASS)
  
  simulated_data <- matrix(mvrnorm(n=nrow*ncol, mu=colMeans(reference_data), Sigma=cov(data), ...), nrow=nrow, ncol=ncol)
  
  simulated_data
}

######### Function No3  ##################################

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
#' #'####################################################

compare_results <- function(x, y) {
  c <- cor.test(as.vector(x), as.vector(y))
  data.frame(Cor=c$estimate, CorP=c$p.value)
  
}

######### Function No4  ##################################

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
#'####################################################


simulate_missingness <- function(data, mcar=0, mar=0, mnar=0, mnar.type="left", mar.type="left") {
  
  if(class(data) != "matrix") {
    stop("Variable data should be a matrix.")
  }
  if ( ((mcar + mar + mnar) > 1) || ((mcar + mar + mnar) < 0)) {
    stop("Sum of mcar, mar and mnar should be between 0 and 1.")
  }
  
  simulated_data <- data
  if (mcar > 0){
    mcar_distribution   = runif(nrow(data)*ncol(data), min=0, max=1)
    simulated_data = matrix(ifelse(mcar_distribution<mcar, NA, data), nrow=nrow(data), ncol=ncol(data))
  }
  
  if (mnar > 0) {
    added_mnar <- 0
    initial_nas <- sum(colSums(is.na(simulated_data)))
    current_nas <- initial_nas
    simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
    # condition of porpotion of missigness
    while (((current_nas - initial_nas) / simulated_data_n) < mnar) {
      
      # Select random variable
      variable_index <- sample(1:ncol(simulated_data), 1)
      
      # What percentage of variable to set missing
      cut_percentage <- 0
      while (cut_percentage <= 0) {
        cut_percentage <- rchisq(1, df=1) / 30
      }
      cut_percentage <- min(cut_percentage, 1)
      
      # How many values to set missing
      cut_index <- floor(cut_percentage * nrow(simulated_data))
      
      sorted_variable <- sort(simulated_data[,variable_index])
      
      # Set values to missing
      if (mnar.type == "right") {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[length(sorted_variable) - cut_index]
        
        simulated_data[simulated_data[,variable_index] > cut_point, variable_index] <- NA
      } else {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[cut_index]
        
        simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
      }
      
      # Counter to check how much MNAR missingness has been added to data
      current_nas <- sum(colSums(is.na(simulated_data)))
    }
  }
  
  if (mar > 0) {
    added_mar <- 0
    initial_nas <- sum(colSums(is.na(simulated_data)))
    current_nas <- initial_nas
    simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
    # condition of porpotion of missigness
    while (((current_nas - initial_nas) / simulated_data_n) < mar) {
      
      # Select random variable
      variable_index <- sample(1:ncol(simulated_data), 1)
      variable_index2 <- sample(setdiff(1:ncol(simulated_data), variable_index), 1)
      
      
      # What percentage of variable to set missing
      cut_percentage <- 0
      while (cut_percentage <= 0) {
        cut_percentage <- rchisq(1, df=1) / 30
      }
      cut_percentage <- min(cut_percentage, 1)
      
      # How many values to set missing
      cut_index <- floor(cut_percentage * nrow(simulated_data))
      
      sorted_variable <- sort(simulated_data[,variable_index])
      
      # Set values to missing
      if (mar.type == "right") {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[length(sorted_variable) - cut_index]
        
        simulated_data[simulated_data[,variable_index] > cut_point, variable_index2] <- NA
      } else {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[cut_index]
        
        simulated_data[simulated_data[,variable_index] < cut_point, variable_index2] <- NA
      }
      
      # Counter to check how much MNAR missingness has been added to data
      current_nas <- sum(colSums(is.na(simulated_data)))
    }
  }
  
  simulated_data
  
}

######### Function No5  ##################################


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
#' #'####################################################

impute <- function(data, methods) {
  
  require(missForest)
  require(pcaMethods)
  require(impute)
  require(PEMM)
  
  if (length(methods) != 1 & length(methods) != ncol(data)) {
    stop("Methods needs to be either one value or of the same length as number of columns in data.")
  }
  
  imputed_data <- data
  
  if (length(methods) == 1) {
    methods <- rep(methods, times=ncol(data))
  }
  if ("RF" %in% methods) {
    rf_imputed_data <- missForest::missForest(xmis = data,maxiter = 10,verbose = TRUE)$ximp
    index <- which(methods == "RF")
    imputed_data[,index] <- rf_imputed_data[,index]
  }
  
  if ("PPCA" %in% methods) {
    # Do cross validation with ppca for component 2:10
    esti <- kEstimate(Matrix= data, method = "ppca", evalPcs = 2:10, nruncv=1, em="nrmsep")
    # The best result was obtained for this number of PCs:esti$bestNPcs
    pc <- pcaMethods::pca(object = data,nPcs=esti$bestNPcs, method="ppca")
    #index <- which(methods == "PPCA")
    imputed_data <- completeObs(pc)
  }
  if("LLS" %in% methods ){
    # local least squares imputation
    lls_esti <- pcaMethods::llsImpute(Matrix = data, k = 50, center = FALSE, completeObs = TRUE, correlation = "kendall", allVariables = TRUE, maxSteps = 100)
    imputed_data <- completeObs(lls_esti)
    
  } 
  
  if("svdImpute" %in% methods ){
    # Singular Value Decomposition
    svd_esti <- pcaMethods::pca(object = data, method="svdImpute", nPcs=10, center = FALSE)
    imputed_data <- completeObs(svd_esti)
  } 
  
  if("KNNImpute" %in% methods ){
    # K-Nearest neighboors
    KNN_esti <- impute::impute.knn(data = data)
    imputed_data <- KNN_esti$data
  } 
  if ("EM" %in% methods){
    # expectation minimazation algorithm
    EM_esti = PEMM::PEMM_fun(X= data, phi=1)
    imputed_data <- EM_esti$Xhat
    
  }
    if ("BPCA" %in% methods){
      # bayesian principal component analysis
      pc <- pca(object = data, method="bpca", nPcs=10)

      ## Get the estimated complete observations
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

######### Function No6  ##################################


#' Title Rsquare_adjusted
#' Calculates the  coefficient of determination  between simulated(reference) data and imputed data
#' @param data1 = simulated data matrix
#' @param data2 = data matrix with missing values
#' @param data3 = imputed data matrix
#'
#' @return
#' @export differences
#'
#' @examples
#' #'####################################################

Rsquare_adjusted <- function(original.data, missing.data, imputed.data){
  rsquare <- ((sum((original.data[is.na(missing.data)] - imputed.data[is.na(missing.data)])^2) / sum(original.data[is.na(missing.data)]^2)))
  adjustment <- (nrow(original.data) -1)/(nrow(original.data)-ncol(original.data)-1)
  rsquare_adjustment <-  1-(rsquare*adjustment)
  rsquare_adjustment
  
}



######### Function No7  ##################################


#' Title Differences_models
#' Calculates the  differences  between simulated(reference) data and imputed data
#' @param data1 = simulated data matrix
#' @param data2 = data matrix with missing values
#' @param data3 = imputed data matrix
#'
#' @return
#' @export differences
#'
#' @examples
#'####################################################
differences_models <- function(original.data, missing.data, imputed.data){
  differences <- sum((original.data[is.na(missing.data)] - imputed.data[is.na(missing.data)])^2) / sum(original.data[is.na(missing.data)]^2)
  differences
}

######### Function No8  ##################################

#' Title replace_NA_to_0
#'A function that replaces missing values to zero in a matrix 
#' when the percentage of missingness in each column is equal or greater than 0.8
#' @param missing.data 
#'
#' @return missing.data
#' @export
#'
#' @examples
#' 
#' ####################################################
replace_NA_to_0 <- function(missing.data) {
  
  # initialize
  col_nan <- NULL
  total_length_col <- NULL
  percentage_nan<- NULL
  # loop every column
  for (i in col(missing.data)){
    # check how many rows in every column have NA
    col_nan[i] <- length(which(is.na(missing.data[,i])))
    total_length_col[i] <- length(missing.data[,i])
    #calculate the percentage of NA of every column
    percentage_nan[i] <- col_nan[i]/total_length_col[i]
    # if the percenatge of NA every column is above 0.8 replace with zero
    if (percentage_nan[i] >= 0.8){
      missing.data[is.na(missing.data[,i]),i] <- 0
    }
  }
  missing.data
}