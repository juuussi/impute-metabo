rm(list = ls())

# list of libraries

library(doMC)
library(dplyr)
# Parallel processing settings
registerDoMC(cores=30)
# Path settings
input.path <- "~/projects/metabolomics_kihd_egg_t2d/"
output.path <- "~/projects/impute-metabo/results/SimulateDatasets/"

# Load external functions
source("/home/users/mariekok/projects/visualizations/src/functions.R")
source("/home/users/mariekok/projects/visualizations/src/visualization_functions.R")
source("~/projects/metabolomics_kihd_egg_t2d/src/functions.R")
source("~/projects/impute-metabo/src/functions.R")


# How many traits will be analyzed, NULL for all
n.responses <- NULL

# Should everything be run from the scratch
initialize <- TRUE

# Log settings
logging <- TRUE
log.file <- paste(output.path, "analysis_log.txt", sep="")
init.log(log.file=log.file, logging=logging)



data.files <- c("KIHD_egg_T2D_hilic_neg","KIHD_egg_T2D_hilic_pos","KIHD_egg_T2D_rp_neg","KIHD_egg_T2D_rp_pos")


full.normalized.data <- read.csv(file=paste(input.path, "data/KIHD_egg_T2D_samples.csv", sep="")) %>% filter(GROUP != "QC") %>% dplyr::select(DATAFILE,T2D,EGGS)
combined.data.list <- list()
for (data.file in data.files) {
  log.text(paste("\nPROCESSING: ", data.file, "\n", sep=""), logging=logging, log.file=log.file)
  
  # Load data
  full.compound.data <- load.metabo.data(file=paste(input.path, "data/", data.file, ".csv", sep=""))
  compound.data <- extract.compound.data(full.compound.data,cols = 2:ncol(full.compound.data))
  compound.data$DATAFILE <- as.factor(gsub("^X", "", compound.data$DATAFILE))
  pheno.data <- load.pheno.data(file=paste(input.path, "data/KIHD_egg_T2D_samples.csv", sep=""))
  combined.data <- combine.pheno.and.compound.data(pheno.data=pheno.data, compound.data=compound.data)
  combined.data$T2D <- relevel(combined.data$T2D,ref="no")
  combined.data$EGGS <- relevel(combined.data$EGGS,ref="low")
  combined.data$T2D <- relevel(combined.data$T2D,ref="QC")
  combined.data$EGGS <- relevel(combined.data$EGGS,ref="QC")
  
  log.text(paste("Combined data: ", nrow(combined.data), " rows, ", ncol(combined.data), sep=""), logging=logging, log.file=log.file)
  
  responses <- select.responses(data=combined.data, pattern="^COMPOUND_",  n=n.responses)
  log.text(paste("Responses: ", length(responses), sep=""), logging=logging, log.file=log.file)
  # Remove QC measurements from analysis
  combined.data <- combined.data[combined.data$GROUP != "QC",]
  combined.data$EGGS <- factor(combined.data$EGGS)
  combined.data$T2D <- factor(combined.data$T2D)
  log.text(paste("Combined data (removed QC samples): ", nrow(combined.data), " rows, ", ncol(combined.data), sep=""), logging=logging, log.file=log.file)
  
  combined.data.list[[data.file]] <-combined.data
}
names(combined.data.list) <- data.files

