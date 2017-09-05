
#' Title check.miss
#' check percenatge of missgness
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
check.miss <- function(data = miss_data){
  listMissVar <- numeric(0)
  perc.col <- round(colMeans(is.na(miss_data)),digits = 2)
  
  for (j in 1:length(perc.col)){
    
    if(perc.col [j]> 0.80){
      
      cat('The variable has more than 80% Nas : ',j,"\n" )
    }else{
      
      listMissVar <-  c(listMissVar,j)
    }
    
    
    
  }
  return(listMissVar)
}






#' Title  Detect MAR and MNAR missigness
# It detects if the missigness depends on other varibales
# Correlates an indicator matrix (0 not missing , 1 missing) with  the original data that can help determine
# if variables tend to be missing together (MAR) or not (MCAR).

#'
#' @param data ,matrix with missing values
#'
#' @return results ,list of all missing variables, the pairs or correlated variables and the MNAR_MAR Variables
#' @export
#'
#' @examples marietta <-detect.miss.MNAR.MAR (miss_data)


detect.miss.MNAR.MAR <- function(data) {
  
  listMissVar <- check.miss(data)
   data <- data[,listMissVar]
  #  Elements of x are 1 if a value in the data is missing and 0 if non-missing.
  
  x <- as.data.frame(abs(is.na(data)))
  
  if(any(is.na(data))){
    
    
    #  Extracting variables that have some missing values.
    #  if the sd is zero all elements in column are missing
    cols <- which(sapply(x, sd) > 0)
    y    <- x[,cols]
    
    #  the list of the all missing variables
    names(cols) <- NULL
    MissingVar  <- data.frame(MissVar = cols)
    
    
    #  Now, looking at the relationship between the 
    #  presence of missing values in each variable and the observed values
    #  in other variables:
    
    
    #  correlation matrix
    
    
    
    corr_matrix <- stats:: cor(data, y, use="pairwise.complete.obs" ,method = "pearson")
    
    # change the row and col names
    colnames(corr_matrix)  <- 1:ncol(corr_matrix)
    row.names(corr_matrix) <- 1:nrow(corr_matrix)
    
    
    #   calculate the probability values from the correlations
    CorTest <- psych::corr.p(r=corr_matrix,n=nrow(data), adjust="fdr",alpha=.05)
    
    #   Correlations between variables in data and y together with confidence interval and pvalues
    table_P_CI   <- CorTest$ci
    table_P_CI$p <- round(table_P_CI$p,digits = 2)
    
    
    #   find which variables are significally correlated
    
    DetectMiss  <- which(table_P_CI$p <= 0.05)
    
    #   correlated Varibales pairs
    PairsCorVar <- data.frame(PairVar = rownames(table_P_CI)[DetectMiss])
    
    #   checking which of the columns have missing values from the pairs of correlated variables
    tmp        <- data.frame(do.call('rbind', strsplit(as.character(PairsCorVar$PairVar),'-',fixed=TRUE)))
    ListCorVar <- data.frame(ListVar=union(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2))))
    #   list of the variables that have missing values (TRUE) and complete variables (FALSE)
    MissVarTF  <- data.frame(listTF = is.element(ListCorVar$ListVar,MissingVar$MissVar))
    MAR.MNAR   <- ListCorVar$ListVar[MissVarTF$listTF]
    
    #  final results : all missing variables, the pairs or correlated variables and the MNAR_MAR Variables
    results <- list(MissingVar = MissingVar,PairsCorVar = PairsCorVar,MAR_MNAR = MAR.MNAR)
    
    return(results)
    
  }else{
    print("matrix does not contain any missing values")
  }
  
  
}


#' Title detect.MCAR.MNAR.MAR
#'
#' @param MissVAr 
#' @param MARMNAR 
#' @param data 
#' @param simulateDATA 
#'
#' @return results
#' @export
#'
#' @examples
detect.MCAR.MNAR.MAR <- function(MissVAr = missvar,MARMNAR = marmnar, data = miss_data,simulateDATA = simulated_data){
  
  if(length(missvar) > 0){
    
    
    
    if (length(marmnar) == 0) {
      
      
      MCAR <- data.frame(ListMCAR = missvar)
      View(MCAR)
      
    }else{
      MAR_MNAR <- sort(marmnar,decreasing = F)
      MCAR <-setdiff(missvar,MAR_MNAR)
      MCAR <- data.frame(ListMCAR = MCAR)
      View(MCAR)
      
      
      # ----------------------------
      
      ## left truncation MNAR 
      # Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
      # if the distributions are the same then p-values are high and that means they are left trancated if they are different then they are MAR
      
      # Goodness of fit for left truncated data
      models <- list()
      Pval <- list()
      Padj <-list()
      perc.col <- round(colMeans(is.na(miss_data[,MAR_MNAR])),digits = 2)
      marmnar_df <- miss_data[,MAR_MNAR]
      list_ex_MAR_MNAR <- numeric(0)
      
      
      for (i in 1:length(MAR_MNAR)){
        
        if (perc.col[i] < 0.60){
          
          
          
          xt <-na.omit( marmnar_df[,i])
          
          
          threshold <- min(na.omit(marmnar_df[,i]))
          #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
          models <- c(models, list(truncgof::ks.test(xt, "pnorm",list(mean(simulated_data),  sd(simulated_data)), H = threshold,  alternative ="two.sided")))
          
          Pval <- c(Pval, tail(models,1)[[1]]$p.value)
          
        }else{
          cat("Too many NAs variable is excluded: ",MAR_MNAR [i], "\n")
          marmnar_df[is.na(marmnar_df[,i]),i] <- 0
          # list of exluded variables above 60% missigness
          list_ex_MAR_MNAR <- c(list_ex_MAR_MNAR, MAR_MNAR [i])
          
          
          
        }
        
      }
      
      
      if  (length(list_ex_MAR_MNAR )== 0){
        Padj <- p.adjust(Pval, method = "fdr")
        Padj <- data.frame(pvalues= Padj,ListVar = MAR_MNAR) 
        MAR <-data.frame(ListMAR = MAR_MNAR[which(Padj$pvalues<=0.05)])
        View(MAR)
        ##
        MNAR <-setdiff(MAR_MNAR,MAR$ListMAR)
        MNAR <- data.frame(ListMNAR = MNAR)
        View(MNAR)
        
      }else {
        NewMAR_MNAR <- setdiff(MAR_MNAR,list_ex_MAR_MNAR)
        Padj <- p.adjust(Pval, method = "fdr")
        
        Padj <- data.frame(pvalues= Padj,ListVar = NewMAR_MNAR) 
        
        MAR <-data.frame(ListMAR = NewMAR_MNAR[which(Padj$pvalues<=0.05)])
        View(MAR)
        ##
        MNAR <-setdiff(NewMAR_MNAR,MAR$ListMAR)
        MNAR <- data.frame(ListMNAR = MNAR)
        View(MNAR)
      }
      
      
    }
    
    
  } else{
    MCAR <- data.frame(ListMCAR = missvar)
    View(MCAR)
  }
  View(MCAR)
  
  
  
  
  
  results <- list(MCAR = MCAR,MNAR = MNAR ,MAR = MAR, Excluded_Var = list_ex_MAR_MNAR )
  return(results)
}











