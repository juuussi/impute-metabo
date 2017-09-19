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
simulated_data <- simulate_data(data=reference_data, nrow=20, ncol=50
                                )
miss_data <- simulate_missingness(data=simulated_data, mcar=0.1, mnar=0.3, mar=0)
#View(miss_data)


billy <- detect.miss.MNAR.MAR(data.frame(miss_data))
MissingVar <-data.matrix( billy$MissingVar)
MAR_MNAR <- billy$MAR_MNAR
missigness <- detect.MCAR.MNAR.MAR (miss_data,MissingVar,MAR_MNAR)



listvector <- detect_missingness_type(missigness) 


aaa <- select_imputation_method(listvector)

test_d <- impute(miss_data,aaa)


