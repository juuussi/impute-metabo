
# Detect missingness
library(psych)
library(dplyr)
library(tidyr)



# choose path
path <- "~/projects/impute-metabo/"
source(paste0(path,"src/functions.R"))
# start logging process by creating a logging file

###################################################################################
# use dummy reference data (combination of different metabolomics data)
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
simulated_data <- simulate_data(data=reference_data, nrow=150, ncol=200)
miss_data <- simulate_missingness(data=simulated_data, mcar=0.30, mnar=0, mar=0)



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
y <- x[which(sapply(x, sd) > 0)]

cor(y)


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
 
 CorrelVariables <- data.frame(VarNames =rownames(CI)[DetectMiss])
 View(CorrelVariables)

 tmp <- data.frame(do.call('rbind', strsplit(as.character(CorrelVariables$VarNames),'-',fixed=TRUE)))
 
 View(tmp)
 MARvariables <-sort(unique(c(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2)))),decreasing = F)
 MARvariables <- data.frame(ListMAR =MARvariables)
 View(MARvariables)
 
 colnames(y)<- 1:ncol(y)
 
 MissingVar <-data.frame(MissVar =1:ncol(y))
View(MissingVar)
 MCARvariables <-setdiff(MissingVar$MissVar,MARvariables$ListMAR)

 MCARvariables <- data.frame(ListMCAR = MCARvariables)
View(MCARvariables)
# ----------------------------

## left truncation MNAR 
#Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
# if the distributions are the same then p-values are small and that means they are left trancated if they are different then they are MAR
library(truncgof)

# Goodness of fit for left truncated data
models <- list()
Pval <- list()
for (i in 1:nrow(MARvariables)){
  
xt <-na.omit(as.matrix(miss_data[,i]))
threshold <- min(na.omit(as.matrix(miss_data[,i])))
models[[i]] <- truncgof::ks.test(xt, "plnorm",list(meanlog = 2, sdlog = 1), H = threshold,alternative ="two.sided")
Pval[[i]] <-models[[i]]$p.value

}
Pval <- as.numeric(as.character(Pval))
Pval <- data.frame(pvalues= Pval,ListVar =MARvariables$ListMAR) 
View(Pval)
SigpVal <-Pval$pvalues[which(Pval$pvalues<0.05)]
DetectMissMNAR <-MARvariables$ListMAR[which(Pval$pvalues<0.05)]

MNARvariables <- data.frame(ListMNAR = DetectMissMNAR)
