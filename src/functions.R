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
  
  simulated_data <- matrix(mvrnorm(n=nrow*ncol, mu=colMeans(data), Sigma=cov(data), ...), nrow=nrow, ncol=ncol)
  
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
#' @examples 
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
  imputed_data <- matrix(NA,nrow = nrow(data),ncol = ncol(data))
  results_data <- data
  
  if (length(methods) == 1) {
    methods <- rep(methods, times=ncol(data))
  }
  
  if ("RF" %in% methods) {
    imputed_data <- missForest::missForest(xmis = data,maxiter = 10,verbose = TRUE)$ximp
    index <- which(methods == "RF")
    results_data[,index] <- imputed_data[,index]
    
  }
  
  if ("PPCA" %in% methods) {
    # Do cross validation with ppca for component 2:10
    esti <- kEstimate(Matrix= data, method = "ppca", evalPcs = 2:10, nruncv=1, em="nrmsep")
    # The best result was obtained for this number of PCs:esti$bestNPcs
    pc <- pcaMethods::pca(object = data,nPcs=esti$bestNPcs, method="ppca")
    #index <- which(methods == "PPCA")
    imputed_data <- completeObs(pc)
    
    index <- which(methods == "PPCA")
    results_data[,index] <- imputed_data[,index]
    
  }
  if("LLS" %in% methods ){
    # local least squares imputation
    lls_esti <- pcaMethods::llsImpute(Matrix = data, k = 10, center = TRUE, completeObs = TRUE, correlation = "pearson", allVariables = TRUE, maxSteps = 100)
    imputed_data <- completeObs(lls_esti)
    
    index <- which(methods == "LLS")
    results_data[,index] <- imputed_data[,index]
    
    
  } 
  
  if("svdImpute" %in% methods ){
    # Singular Value Decomposition
    svd_esti <- pcaMethods::pca(object = data, method="svdImpute", nPcs=10, center = FALSE)
    imputed_data <- completeObs(svd_esti)
    
    index <- which(methods == "svdImpute")
    results_data[,index] <- imputed_data[,index]
    
  } 
  
  if("KNNImpute" %in% methods ){
    # K-Nearest neighboors
    KNN_esti <- impute::impute.knn(data = data)
    imputed_data <- KNN_esti$data
    
    index <- which(methods == "KNNImpute")
    results_data[,index] <- imputed_data[,index]
    
  } 
  if ("EM" %in% methods){
    # expectation minimazation algorithm
    EM_esti = PEMM::PEMM_fun(X= data, phi=1)
    imputed_data <- EM_esti$Xhat
    
    index <- which(methods == "EM")
    results_data[,index] <- imputed_data[,index]
    
    
  }
  if ("BPCA" %in% methods){
    # bayesian principal component analysis
    pc <-pcaMethods:: pca(object = data, method="bpca", nPcs=10)
    
    ## Get the estimated complete observations
    imputed_data <- completeObs(pc)
    
    index <- which(methods == "BPCA")
    results_data[,index] <- imputed_data[,index]
    
  }
  if ("min" %in% methods) {
    #foreach (data_column=1:ncol(data)) %do% {
    
    # foreach (data_column=which(methods == "min")) %do% {
    #   method <- methods[data_column]
    #   
    #   impu_value <- min(data[,data_column], na.rm=TRUE)
    #   imputed_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    #   
    # }
    
    
    for (i in 1:ncol(data)){
      data[is.na(data[,i]),i] <- round(min(data[,i],na.rm = TRUE))
    }
    
    index <- which(methods == "min")
    results_data[,index] <- data[,index]
    
  }
  
  if ("mean" %in% methods){
    # 
    # #foreach (data_column=1:ncol(data)) %do% {
    # foreach (data_column=which(methods == "mean")) %do% {
    #   method <- methods[data_column]
    #   impu_value <- mean(data[,data_column], na.rm=TRUE)
    #   imputed_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    # }
    # index <- which(methods == "mean")
    # results_data[,index] <- imputed_data[,index]
    # 
    
    for (i in 1:ncol(data)){
      data[is.na(data[,i]),i] <- round(mean(data[,i],na.rm = TRUE))
    }
    index <- which(methods == "mean")
    results_data[,index] <- data[,index]
    
  }
  
  
  
  results_data
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



######### Function No9  ##################################

#' Title check.miss
#' check the perce
#' @param data  data matrix
#'
#' @return listVar
#'  list of 3 vectors
#'  MissingVar     = vector containing the columns numbers of the data matrix with missing values that are less than 80 % NAs
#'  ExcludedVar    = vector containing the columns numbers of the data matrix with missing values that are more than 80% NAs
#'  CompleteVar    = vector containing the columns numbers of the data matrix without NAs
#' @export
#'
#' @examples

check.miss <- function(data){
  if (is.data.frame(data)) {
    data <- as.matrix(data)
    
  }
  
  AllVar      <- as.numeric(1:ncol(data))
  MissingVar  <- numeric(0)
  ExcludedVar <- numeric(0)
  CompleteVar <- numeric(0)
  listVar <- list()
  
  
  
  perc.col <- round(colMeans(is.na(data)),digits = 2)
  
  
  for (j in 1:length(perc.col)){
    
    if(perc.col [j] > 0.80){
      ExcludedVar <- c(ExcludedVar,j)
      cat('The variable has more than 80% Nas : ',j,"\n" )
    }else if(perc.col [j] == 0){
      
      
      CompleteVar <- c(CompleteVar,j)
    }else{
      MissingVar <- c(MissingVar,j)
    }
  }
  
  listVar <- list( MissingVar = MissingVar , ExcludedVar = ExcludedVar, CompleteVar = CompleteVar )
  return(listVar)
}


######### Function No10  ##################################
#' Title  detect.miss.MNAR.MAR
# It detects if the missigness depends on other varibales
# Correlates an indicator matrix (0 not missing , 1 missing) with  the original data that can help determine
# if variables tend to be missing together (MAR) or not (MCAR).
# Kabacoff, Robert I. R in Action. manning, 2010.

#'
#' @param data ,matrix with missing values
#'
#' @return results 
#' list of vectors :
#' MissingVar  = vector containing the columns numbers of the data matrix with missing values that are less than 80 % NAs,
#' PairsCorVar = data frame containing the pairs of correlated variables in the data matrix,
#' MAR_MNAR    = vector containing the columns numbers of the data matrix with MAR or MNAR missingness
#' ExcludedVar = vector containing the columns numbers of the data matrix with missing values that are more than 80% NAs
#' @export
#'
#' @examples marietta <-detect.miss.MNAR.MAR (miss_data)


detect.miss.MNAR.MAR <- function(data) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # check 80% missigness, exclude variables 
  Miss_Var <- check.miss(data)
  data     <- data[,Miss_Var[[1]]]
  
  colnames(data) <- as.character(Miss_Var[[1]])
  
  #  Elements of x are 1 if a value in the data is missing and 0 if non-missing.
  
  x <- data.frame(abs(is.na(data)))
  colnames(x) <- colnames(data) 
  
  if(any(is.na(data))){
    
    
    #  Extracting variables that have some missing values.
    #  if the sd is zero all elements in column are missing
    cols <- which(sapply(x, sd) > 0)
    y    <- x[,cols]
    y    <- as.matrix(y)
    
    #  the list of the all missing variables
    
    MissingVar <- as.numeric(colnames(y))
    
    #  Now, looking at the relationship between the 
    #  presence of missing values in each variable and the observed values
    #  in other variables:
    
    
    #  correlation matrix
    corr_matrix <- stats:: cor(data, y, use="pairwise.complete.obs" ,method = "pearson")
    
    #   calculate the probability values from the correlations
    CorTest     <- psych::corr.p(r=corr_matrix,n=nrow(data), adjust="fdr",alpha=.05)
    
    #   Correlations between variables in data and y together with confidence interval and pvalues
    table_P_CI  <- round(CorTest$ci,digits = 2)
    
    #   find which variables are significally correlated
    
    SigPvalues  <- which(table_P_CI$p <= 0.05)
    
    #   correlated Varibales pairs
    PairsCorVar <- data.frame(PairVar = rownames(table_P_CI)[SigPvalues])
    
    #   checking which of the columns have missing values from the pairs of correlated variables
    tmp         <- data.frame(do.call('rbind', strsplit(as.character(PairsCorVar$PairVar),'-',fixed=TRUE)))
    CorVar <- sort(union(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2))),decreasing = F)
    #   check which variables are trully missing
    CorVarTF    <-  is.element(CorVar,MissingVar)
    
    # MAR and MAR variables
    MAR.MNAR    <- CorVar[CorVarTF]
    
    #  final results : all missing variables, the pairs of correlated variables and the MNAR_MAR Variables
    results     <- list(MissingVar = MissingVar,PairsCorVar = PairsCorVar,MAR_MNAR = MAR.MNAR ,ExcludedVar = Miss_Var[[2]])
    
    return(results)
    
  }else{
    print("matrix does not contain any missing values")
  }
}

######### Function No11  ##################################

#' Title detect.MCAR.MNAR.MAR
#'
#' @param MissingVar 
#' vector containing the columns numbers of the data matrix with missing values that are less than 80 % (use the detect.miss.MNAR.MAR function)
#' @param MAR_MNAR 
#'  vector containing the columns with missing values that have MAR or MNAR missigness (use the detect.miss.MNAR.MAR function)
#' @param data 
#'
#' @return results
#' list of  vectors :
#' MCAR = vector conatining the column numbers of MCAR variables, 
#' MNAR = vector conatining the column numbers of MNAR variables ,
#' MAR  = vector conatining the column numbers of MAR variables ,
#' Excluded_marmnar  = vector conatining the column numbers of excluded MAR or MNAR variables that have more than 60% NAs  ,
#' ExcludedVar       = vector conatining the column numbers of excluded variables that have more than 80% NAs  ,
#' CompleteVar       = conatining the column numbers of variables without missing values.
#' @export
#'
#' @examples
detect.MCAR.MNAR.MAR <- function(data ,MissingVar, MAR_MNAR ){
  
  MAR  <- numeric(0)
  MNAR <- numeric(0)
  MAR  <- numeric(0)
  NewMAR_MNAR <- numeric(0)
  rm_MAR_MNAR <- numeric(0)
  
  if (is.data.frame(data)) {
    data <- as.matrix(data)
    
  }
  if(length(MissingVar) > 0){
    
    if (length(MAR_MNAR) == 0) {
      
      MCAR <- MissingVar
      
    }else{
      
      MCAR <- setdiff(MissingVar,MAR_MNAR)
      ## left truncation MNAR 
      #  Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
      #  if the distributions are the same then p-values are high and that means they are left trancated if 
      #  they are different then they are MAR
      
      #  Goodness of fit for left truncated data
      models <- list()
      Pval   <- list()
      Padj   <-list()
      perc.col    <- round(colMeans(is.na(data[,MAR_MNAR])),digits = 2)
      marmnar_mat  <- data[,MAR_MNAR]
      colnames(marmnar_mat) <- as.character(MAR_MNAR)
      
      for (i in 1:length(MAR_MNAR)){
        # the KS test doesnt work for the 60 % NAs per column
        if (perc.col[i] < 0.60){
          
          xt     <- na.omit( marmnar_mat[,i])
          threshold <- min(na.omit(marmnar_mat[,i]))
          #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
          models <- c(models, list(truncgof::ks.test(xt, "pnorm",list(mean(data, na.rm = T),  sd(data,na.rm = T)), H = threshold,  alternative ="two.sided")))
          # pvalues from the KS test
          Pval <- c(Pval, tail(models,1)[[1]]$p.value)
          
        }else{
          cat("Too many NAs  variable the MAR_MNAR is excluded: ",MAR_MNAR [i], "\n")
          marmnar_mat[is.na(marmnar_mat[,i]),i] <- 0
          # list of exluded variables above 60% missigness
          rm_MAR_MNAR <- c(rm_MAR_MNAR, MAR_MNAR [i])
        }
      }
      if  (length(rm_MAR_MNAR )== 0){
        # adjust p values obtained from KS test using fdr correction 
        Padj <- p.adjust(Pval, method = "fdr")
        MAR  <- MAR_MNAR[which(Padj <= 0.05)]
        MNAR <- setdiff(MAR_MNAR,MAR)
      }else {
        #  p values adjusted excluding the variables have MNAR or MAR that have been excluded because more than 60% NAs 
        NewMAR_MNAR <- setdiff(MAR_MNAR,rm_MAR_MNAR)
        Padj <- p.adjust(Pval, method = "fdr")
        MAR  <- NewMAR_MNAR[which(Padj <= 0.05)]
        MNAR <- setdiff(NewMAR_MNAR,MAR)
      }
    }
  } 
  Miss_Var <-  check.miss(data)
  results <- list(MCAR = MCAR,MNAR = MNAR ,MAR = MAR, Excluded_marmnar = rm_MAR_MNAR , ExcludedVar = Miss_Var[[2]], CompleteVar = Miss_Var[[3]])
  return(results)
}

######### Function No12  ##################################


#' Title detect_missingness_type
#' Vector conatining the the columns type of missingness
#' @param missigness list of vectors obtained by the detect.MCAR.MNAR.MAR function
#'
#' @return miss_type =  vector conatining the the columns type of missingness
#' @export
#'
#' @examples
detect_missingness_type <- function(missigness){
  
  miss_type <- numeric(0)
  AllVar    <- numeric(0) 
  MCAR <- missigness[[1]]
  MNAR <- missigness[[2]]
  MAR  <- missigness[[3]]
  ExcludedVar <- union(missigness[[4]],missigness[[5]])
  CompleteVar <- missigness [[6]]
  
  
  
  
  
  miss_type <- character(length(MCAR) + length(MNAR) + length(MAR) + length(ExcludedVar) + length(CompleteVar))
  miss_type[MCAR]        <- "MCAR"
  miss_type[MAR]         <- "MAR"
  miss_type[MNAR]        <- "MNAR"
  miss_type[ExcludedVar] <- "EX"
  miss_type[CompleteVar] <- "NONE"
  
  
  return(miss_type)  
}



######### Function No12  ##################################


#' Title select_imputation_method
#'
#' @param types 
#' @param data 
#'
#' @return imputed_data_list
#' @export
#'
#' @examples
select_imputation_method <- function(types,data){
  imputed_data <-matrix(NA,nrow = nrow(data),ncol = ncol(data))
  imputed_data_list <- list()
  
  for (i in 1:length(types)){
    if(types[i]== "MCAR"){
      imputed_data  <- impute(data,"mean")
    }else if(types[i] == "MAR"){
      imputed_data  <- impute(data,"mean") 
      
    }else if (types[i] == "MNAR"){
      imputed_data <- impute(data,"min") 
      
      
    }else if(types[i]== "NONE"){
      imputed_data <- data
      
      
    }else if(types[i]== "EX"){
      imputed_data <- 0
      
    }
    imputed_data_list <- c(imputed_data_list, list(imputed_data))
    
  }
  return(imputed_data_list)
  
}




