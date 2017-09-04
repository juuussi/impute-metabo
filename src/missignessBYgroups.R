rm(list = ls())
# Detect missingness
library(psych)
library(dplyr)
library(tidyr)
library(truncgof)


set.seed(1234)
# choose path
path <- "~/projects/impute-metabo/"
source(paste0(path,"src/functions.R"))
# start logging process by creating a logging file

###################################################################################
# use dummy reference data (combination of different metabolomics data)
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
simulated_data <- simulate_data(data=reference_data, nrow=15, ncol=20)
miss_data <- simulate_missingness(data=simulated_data, mcar=0, mnar=0, mar=0.3)
View(miss_data)
X <-dim(miss_data)

labels <- sample(c(0,1), replace=TRUE, size=X[1])


#x <- data.frame(cbind(labels,miss_data))
x <- miss_data
  #check the percentage of missigness for each column
  
  perc.col <- round(colMeans(is.na(x)),digits = 2)
  
  #x %>% summarize_all(funs(sum(is.na(.)) / length(.)))
  
  
  # check percentage of missingness per group
  
  
  #say you want to get missing values from group 0 and 1
  #perc.group <-x %>% group_by(labels) %>% summarise_all(funs(sum(is.na(.)) / length(.)))
  
  #perc.group <- round(perc.group,digits = 2)
  
  missvar <-data.matrix( billy$MissingVar)
  marnar <- billy$MAR_MNAR
  
  
  if(length(missvar) > 0){
    
    
    
   
    if (length(marnar) == 0) {
      MCAR <- data.frame(ListMCAR = missvar)
      View(MCAR)
      
    }else{
      MAR_MNAR <- sort(marnar,decreasing = F)
      
    
      
      
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
        
        if (perc.col[i] < 0.40){
          
          
          
          xt <-na.omit( marmnar_df[,i])
          
          
          threshold <- min(na.omit(marmnar_df[,i]))
          #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
          
          #models[[i]] <- truncgof::ad.test(xt, "pnorm",list(meanlog = mean(simulated_data), sdlog = sd(simulated_data)), H = threshold,alternative ="two.sided")
          models <- c(models, list(truncgof::ks.test(xt, "pnorm",list(mean(simulated_data),  sd(simulated_data)), H = threshold,  alternative ="two.sided")))
          
          Pval <- c(Pval, tail(models,1)[[1]]$p.value)
          

          
        }else{
          cat("Too many NAs variable is excluded: ",MAR_MNAR [i], "\n")
          marmnar_df[is.na(marmnar_df[,i]),i] <- 0
          # list of exluded 
          list_ex_MAR_MNAR <- c(list_ex_MAR_MNAR, MAR_MNAR [i])
          
          
          
        }
        
        
      }
      
     
      
      if  (length(list_ex_MAR_MNAR )== 0){
        Padj <- p.adjust(Pval, method = "fdr")
        Padj <- data.frame(pvalues= Padj,ListVar = MAR_MNAR) 
        #View(Padj)
        #SigpVal <-Padj$pvalues[which(Padj$pvalues <= 0.05)]
        MAR <-data.frame(ListMAR = MAR_MNAR[which(Padj$pvalues<=0.05)])
        #MARvariables <- data.frame(ListMAR =DetectMissMAR2$ListMar)
        View(MAR)
        ##
        MNAR <-setdiff(MAR_MNAR,MAR$ListMAR)
        MNAR <- data.frame(ListMNAR = MNAR)
        View(MNAR)
        
      }else {
        NewMAR_MNAR <- setdiff(MAR_MNAR,list_ex_MAR_MNAR)
        Padj <- p.adjust(Pval, method = "fdr")
        
        Padj <- data.frame(pvalues= Padj,ListVar = NewMAR_MNAR) 
        #View(Padj)
        #SigpVal <-Padj$pvalues[which(Padj$pvalues <= 0.05)]
        MAR <-data.frame(ListMAR = NewMAR_MNAR[which(Padj$pvalues<=0.05)])
        #MARvariables <- data.frame(ListMAR =DetectMissMAR2$ListMar)
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
  
  
  
  
  
  
  