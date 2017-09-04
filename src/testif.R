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



# -----------------------------------------Detect Type of Missigness  ---------------------------------
#
# This method creates binary matrix, in which the dataset consists of indicator variables 
# where a 1 is given if a value isnot present, and 0 if it isinstall.package  .
#Correlating these with each other and the original data can help determine 
#if variables tend to be missing together (MAR) or not (MCAR).

#Elements of x are 1 if a value in the miss_data is missing and 0 if non-missing.

x <- as.data.frame(abs(is.na(miss_data)))

#head(miss_data)
#head(x)

#Extracting variables that have some missing values.
cols <- which(sapply(x, sd) > 0)
y <- x[,cols]

cor(y)


# the list of the all missing variables
MissingVar <-data.frame(MissVar = cols)

#Now, looking at the relationship between the 
#presence of missing values in each variable and the observed values
#in other variables:


corr_matrix <-cor(miss_data, y, use="pairwise.complete.obs")
row.names(corr_matrix) <- 1:nrow(corr_matrix)
colnames(corr_matrix) <- 1:ncol(corr_matrix)
CorTest <- psych::corr.p(r=corr_matrix,n=nrow(miss_data),adjust="fdr",alpha=.05)

correlation_matrix <-CorTest$r
# extract pvalues
pval <- CorTest$p
pval <- round(pval,digits = 2)
CI <- CorTest$ci
CI$p <-round(CI$p,digits = 2)

# find which variables are significally correlated

DetectMiss <-which(CI$p<0.05)
if(length(DetectMiss) > 0){
  
  
  
  CorrelVariables <- data.frame(VarNames =rownames(CI)[DetectMiss])
  
  View(CorrelVariables)
  
  # checking which of the columns have missing values from the pairs of correlated variables
  tmp <- data.frame(do.call('rbind', strsplit(as.character(CorrelVariables$VarNames),'-',fixed=TRUE)))
  find_var <- union(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2)))
  #find_var <-sort(unique(c(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2)))),decreasing = F)
  find_var <-data.frame(ListVar  = find_var)
  find_var2 <- data.frame(listTF = is.element(find_var$ListVar,MissingVar$MissVar))
  
  #View(find_var2)
  # View(MissingVar)
  if (all(find_var2 == FALSE)){
    MCARvariables <- data.frame(ListMCAR = MissingVar$MissVar)
    View(MCARvariables)
    
  }else{
    MAR_MNAR <- find_var$ListVar[find_var2$listTF]
    
    
    MAR_MNAR <- find_var[find_var2$listTF,]
    MAR_MNARvariables <-data.frame(ListMAR_MNAR = MAR_MNAR)
    
    
    MCARvariables <-setdiff(MissingVar$MissVar,MAR_MNARvariables$ListMAR_MNAR)
    MCARvariables <- data.frame(ListMCAR = MCARvariables)
    View(MCARvariables)
    View(MAR_MNARvariables)
    
    
    # ----------------------------
    
    ## left truncation MNAR 
    #Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
    # if the distributions are the same then p-values are high and that means they are left trancated if they are different then they are MAR
    
    # Goodness of fit for left truncated data
    models <- list()
    Pval <- list()
    Padj <-list()
    
    for (i in 1:nrow(MAR_MNARvariables)){
      
      
      xt <-na.omit(miss_data[,i])
      if(length(xt > 3)){
        
        
        
        
        threshold <- min(na.omit((miss_data[,i])))
        #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
        
        #models[[i]] <- truncgof::ad.test(xt, "pnorm",list(meanlog = mean(simulated_data), sdlog = sd(simulated_data)), H = threshold,alternative ="two.sided")
        models[[i]] <- truncgof::ks.test(xt, "pnorm",list(mean(simulated_data),  sd(simulated_data)), H = threshold,  alternative ="two.sided")
        
        Pval[[i]] <-models[[i]]$p.value
        Padj[[i]]<-p.adjust(Pval[[i]], method = "fdr")
        
      
      
      Padj <- as.numeric(as.character(Padj))
      Padj <- data.frame(pvalues= Padj,ListVar =MAR_MNARvariables$ListMAR_MNAR) 
      #View(Padj)
      #SigpVal <-Padj$pvalues[which(Padj$pvalues <= 0.05)]
      MARvariables <-data.frame(ListMAR = MAR_MNARvariables$ListMAR_MNAR[which(Padj$pvalues<=0.05)])
      #MARvariables <- data.frame(ListMAR =DetectMissMAR2$ListMar)
      View(MARvariables)
      ##
      MNARvariables <-setdiff(MAR_MNARvariables$ListMAR_MNAR,MARvariables$ListMAR)
      MNARvariables <- data.frame(ListMNAR = MNARvariables)
      View(MNARvariables)
      
      }
    }
  else{ print("blalalalla")}
    #cat("the vector conatings less than 3 values")}
  }
} else{
  MCARvariables <- data.frame(ListMCAR = MissingVar$MissVar)
  View(MCARvariables)
}
View(MCARvariables)
