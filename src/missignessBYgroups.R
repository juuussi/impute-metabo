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
simulated_data <- simulate_data(data=reference_data, nrow=15, ncol=20)
miss_data <- simulate_missingness(data=simulated_data, mcar=0, mnar=0, mar=0.4)
View(miss_data)


billy <- detect.miss.MNAR.MAR(miss_data)
missvar <-data.matrix( billy$MissingVar)
marmnar <- billy$MAR_MNAR
mak <- detect.MCAR.MNAR.MAR (MissVAr = missvar,MARMNAR = marmnar, data = miss_data,simulateDATA = simulated_data)
