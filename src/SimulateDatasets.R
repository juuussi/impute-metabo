rm(list = ls())

# list of libraries
library(foreach)

library(dplyr)
library(futile.logger)
library(missForest)
library(pcaMethods)
library(doMC)
library(speedglm)

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
compound.data.list <- list()
for (data.file in data.files) {
  log.text(paste("\nPROCESSING: ", data.file, "\n", sep=""), logging=logging, log.file=log.file)
  
  # Load data
  full.compound.data <- load.metabo.data(file=paste(input.path, "data/", data.file, ".csv", sep=""))
  compound.data <- extract.compound.data(full.compound.data,cols = 2:ncol(full.compound.data))
  compound.data[compound.data  == 1] <-NA
  
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
  compound.data.list [[data.file]]<- compound.data
}
names(combined.data.list) <- data.files
names(compound.data.list) <- data.files



# *********** Comparing Impuation Methods ******************** #



# start logging process by creating a logging file
flog.appender(appender.tee(paste0(output.path,"logFile_simulateDatasets_2.log")))
flog.threshold(DEBUG)
# propotions of missigness
miss_proportions <- c(0.02)

proportions_df <- missingness_proportions(miss_proportions=miss_proportions)

# impuation methods
imputation_methods <- c("PPCA", "RF","min","mean","LLS","KNNImpute","BPCA","svdImpute")
# number of iterations
n_iterations <- 10
iteration_counter <- 1
total_iterations <-  n_iterations * nrow(proportions_df) * length(imputation_methods)

FULL_results <- foreach(iteration=1:n_iterations, .combine="rbind") %dopar% {
  ##  randomly pick one dataset
  
  random.pick <- sample(length(compound.data.list),1,replace = F)
  
  random.dataset <- compound.data.list[[random.pick]]
  random.dataset <- random.dataset[2:length(random.dataset)]
  
  flog.debug(paste('Dataset:', random.pick))
  # remove NAs
  
  remove.NA <- random.dataset[ , apply(random.dataset, 2, function(x) !any(is.na(x)))]
  #### Randomly pick 200 metabolites
  
  pool.metabolites <-remove.NA[,sample(ncol(remove.NA), 200)]
  pool.metabolites <- as.matrix(pool.metabolites)
  
  # choose the percenatge of missigness 
  proportion_results <- foreach(h=1:nrow(proportions_df), .combine="rbind", .inorder=FALSE) %dopar% {
    
    # seperate the percenatge of missigness per type of misssigness 
    mcar_miss  <- proportions_df$MCAR[h]
    mnar_miss  <- proportions_df$MNAR[h]
    mar_miss   <- proportions_df$MAR[h]
    total_miss <- proportions_df$Total[h]
    
    
    # simulate different type of missgness MCAR,MNAR & MAR using the different percenatges of missigness
    miss_data <- simulate_missingness(data=pool.metabolites, mcar=mcar_miss, mnar=mnar_miss, mar=mar_miss)
    # impute missing values using different imputation methods
    method_results <- foreach(j=1:length(imputation_methods), .combine="rbind", .inorder=FALSE) %do% {
      
      method <- imputation_methods[j] 
      
      # log the total iteration ,method , percenatges and types of missgness 
      flog.debug(paste('Total iterations:', iteration_counter, "/", total_iterations, 'Row size:', nrow(pool.metabolites), 'Method:', method, 'Missingness - Total:', total_miss, "MCAR:", mcar_miss, "MAR:", mar_miss, "MNAR:", mnar_miss, sep=' '))

      # calculate the computational time per method
      
      
      # impute the data
      
      tryCatch({
        flog.debug(paste('Starting imputing:','Method:', method, sep=' '))
        start <- Sys.time ()
        imputed_data <- impute(data=miss_data, methods=method)
        
        
        time_diff <- Sys.time () - start
        # log the computational time per method
        flog.debug(paste('Total impute time:',round(time_diff) ,'Method:', method, sep=' '))
        
        #COMPARING RESULTS USING ROOT MEAN SQUARE ERROR and R square adjusted
        
        #results_Rsquare <- Rsquare_adjusted(original.data = simulated_data, missing.data = miss_data, imputed.data = imputed_data)
        nrmse_error <- missForest:: nrmse(pool.metabolites, miss_data, pool.metabolites)
        #nrmse_error <- differences_models (pool.metabolites, miss_data, imputed_data)
        
        miss_model <- ""
        
        if(mnar_miss>0){
          miss_model <- paste0(miss_model,"MNAR")
        }
        if(mar_miss>0){
          miss_model <- paste0(miss_model,"MAR")
        }
        if(mcar_miss>0){
          miss_model <- paste0(miss_model,"MCAR")
        }
        
        time_diff2 <- Sys.time () - start
        flog.debug(paste('Total nrmse time:',round(time_diff2) ,'Method:', method, sep=' '))
        
        # create a data frame with the results
        results_f <- cbind(data.frame(Method=method, Miss=total_miss,MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss, DatasetNumber = random.pick,n_col=nrow(pool.metabolites),Iteration=iteration, Time_total=time_diff), data.frame(NRMSE=nrmse_error))
        
        time_diff3 <- Sys.time () - start
        flog.debug(paste('Total result binding time:',round(time_diff3) ,'Method:', method, sep=' '))
        
        iteration_counter <- iteration_counter + 1
        results_f 
      }, error=function(x) {
        flog.debug(x)
        results_f <- cbind(data.frame(Method=method, Miss=total_miss, MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss,DatasetNumber = random.pick,n_col=nrow(pool.metabolites) ,Iteration=iteration,Time_total=time_diff), data.frame(NRMSE=NA))
        return(results_f)
      })
    }
    method_results

  }
  proportion_results
}

write.csv(x=FULL_results, file=paste0(path, "results/finalres092017.csv"), row.names=FALSE)
head(FULL_results,10)
