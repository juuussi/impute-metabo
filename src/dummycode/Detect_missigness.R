#detect missigness

#1.1 example with MCAR missigness

# list of libraries

library(missForest)
library(doMC)
library(futile.logger)
library(pcaMethods)
library(BaylorEdPsych)
library(mice)
library(MissMech)
#####################
registerDoMC(cores=30)

path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
flog.appender(appender.tee(paste0(path,"results/logFile.log")))

flog.threshold(DEBUG)

reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
simulated_data <- simulate_data(data=reference_data, nrow=10, ncol=40)
miss_data <- simulate_missingness(data=simulated_data, mcar=0.1, mnar=0, mar=0)
miss_details <- md.pattern(miss_data)
miss_details2 <- LittleMCAR(data.frame(miss_data))

## example 2 missMech package
data(agingdata)
miss_struct <- DelLessData(miss_data, ncases = 0)
h <- Hawkins(data=miss_data, miss_struct$spatcnt)
#####