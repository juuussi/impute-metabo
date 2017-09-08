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
source(paste0(path,"src/testFuncMIss.R"))

# start logging process by creating a logging file

###################################################################################
# use dummy reference data (combination of different metabolomics data)
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
simulated_data <- simulate_data(data=reference_data, nrow=20, ncol=50)
miss_data <- simulate_missingness(data=simulated_data, mcar=0.3, mnar=0, mar=0)
#View(miss_data)


billy <- detect.miss.MNAR.MAR(data.frame(miss_data))
MissingVar <-data.matrix( billy$MissingVar)
MAR_MNAR <- billy$MAR_MNAR
missigness <- detect.MCAR.MNAR.MAR (miss_data,MissingVar,MAR_MNAR)


listvector <- detect_missingness_type(missigness) 

listvector















# 
# 
# listVAR <- check.miss(miss_data)
# df <- miss_data[,listVAR]
# 
# 
# MNAR <- mak$MNAR
# MCAR <- mak$MCAR
# MAR <- mak$MAR
# newdf <-df[,MNAR$ListMNAR]
# imputation_methods <- c("PPCA", "RF","min","mean","LLS","KNNImpute","BPCA","svdImpute")
# 
# nrmse_error <- 0
# nrmse_error1 <- 0
# nrmse_error2 <- 0
# 
# if (length(MNAR) > 0 ||length(MNAR) > 0 ||length(MAR) > 0 ){
#   tmp <- impute(df[,MNAR$ListMNAR],'KNNImpute')
#   nrmse_error <- missForest:: nrmse(tmp, df[,MNAR$ListMNAR], simulated_data[,MNAR$ListMNAR])
#   tmp1 <- impute(df[,MCAR$ListMCAR],'BPCA')
#   nrmse_error1 <- missForest:: nrmse(tmp1, df[,MCAR$ListMCAR], simulated_data[,MCAR$ListMCAR])
#   
#   tmp2 <- impute(df[,MAR$ListMAR],'KNNImpute')
#   nrmse_error2 <- missForest:: nrmse(tmp2, df[,MAR$ListMAR], simulated_data[,MAR$ListMAR])
#   
# }
# 
# 
# hybrid.impute <- function(data,type){
#   
#   switch(type,
#          mnar = impute(data,"min"),
#          mar = impute(data,"mean"),
#          mcar = impute(data,"mean"),
#          marmnar = impute(data,"mean"),
#          marmcar = impute(data,"mean"),
#          mnarmcar = impute(data,"mean"),
#          marmnarmcar = impute(data,"mean")
#   )
#   
#   
# }
# 
# ######################################################################
# 
# library(mvnmle)
# mlest(miss_data)
# 
# 
