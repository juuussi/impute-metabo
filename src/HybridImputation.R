
# clear everything
rm(list = ls())


# list of libraries
library(foreach)
library(futile.logger)
library(missForest)
library(pcaMethods)
library(doMC)
library(speedglm)
# Detect missingness
library(psych)
library(dplyr)
library(tidyr)
library(truncgof)
# Parallel processing settings
registerDoMC(cores=15)
# Path settings
input.path <- "~/projects/impute-metabo/"
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




# ************************** berries ************************** #
n.responses <- NULL
data.file  <-"Berry_intervention"
# Load data
full.compound.data <- read.csv(file=paste(input.path, "data/", data.file, ".csv", sep=""), header=TRUE, stringsAsFactors=FALSE, sep=",")
compound.data <- full.compound.data
# Add COMPOUND to name
compound.data$Compound <- paste("COMPOUND_", compound.data$Compound, sep="")
# Compound column to factor
#compound.data$Compound <- factor(compound.data$Compound)
# Transpose (compound on top, subject on rownames)
compound.data <- t(compound.data)
compound.data <- data.frame(compound.data)

# add first row as colnames
colnames(compound.data) <- NULL
colnames(compound.data) <- as.character(unlist(compound.data[1,]))
# remove first row/name row
compound.data <- compound.data[-1,]
# set columns as numeric
tmp <- apply(compound.data, 2, as.numeric)
rownames(tmp) <- rownames(compound.data)
compound.data <- tmp
# Add rownames as a column named DATAFILE
result.df <- cbind(data.frame(DATAFILE=rownames(compound.data)), compound.data)

# Read phenotype file
pheno.data <- load.pheno.data(file=paste(input.path, "data/", data.file, "_samples.csv", sep=""))
# Combine phenotype and compound data to single data frame
combined.data <- combine.pheno.and.compound.data(pheno.data=pheno.data, compound.data=result.df)
# replace blank with NA
combined.data[combined.data==""] <- NA

## Check NAs
na.test(combined.data)

log.text(paste("Combined data: ", nrow(combined.data), " rows, ", ncol(combined.data), sep=""), logging=logging, log.file=log.file)

responses <- select.responses(data=combined.data, pattern="^COMPOUND_",  n=n.responses)
log.text(paste("Responses: ", length(responses), sep=""), logging=logging, log.file=log.file)

## Fill NAs in "SUBJECT_ID" column with values from qTOF_SAMPLE_ID column
combined.data$qTOF_SAMPLE_ID <- as.character(combined.data$qTOF_SAMPLE_ID)
# convert column to character (col is factor and causes trouble)
combined.data$SUBJECT_ID <- as.character(combined.data$SUBJECT_ID)
combined.data <- replace.na(combined.data, "SUBJECT_ID", "qTOF_SAMPLE_ID", "QC", add.seq.num = TRUE)
## Fill NAs in "SUBJECT_ID" column with values from GROUP column
combined.data$GROUP <- as.character(combined.data$GROUP)
combined.data$TIME <- as.character(combined.data$TIME)
combined.data <- replace.na(combined.data, "GROUP", "qTOF_SAMPLE_ID", "QC", add.seq.num = FALSE)
combined.data <- replace.na(combined.data, "TIME", "qTOF_SAMPLE_ID", "QC", add.seq.num = FALSE)
combined.data$GROUP <- as.factor(combined.data$GROUP)
combined.data$TIME <- as.factor(combined.data$TIME)

# Remove QC measurements from analysis
combined.data <- combined.data[combined.data$GROUP != "QC",]
combined.data$GROUP <- droplevels(combined.data$GROUP)
log.text(paste("Combined data (removed QC samples): ", nrow(combined.data), " rows, ", ncol(combined.data), sep=""), logging=logging, log.file=log.file)

# Cleaned
combined.data <- replace.values(data=combined.data, responses=responses, original.value=1, target.value=NA)


## divide the data per mode and per time beginning and end 
#hilic
hilic.start <- combined.data %>% filter(TIME=="beginning")%>%
  dplyr:: select(DATAFILE,starts_with("COMPOUND_hilic")) 

hilic.pos.start <- hilic.start %>%dplyr:: select(DATAFILE,contains("pos")) 
hilic.neg.start <- hilic.start %>%dplyr:: select(DATAFILE,contains("neg")) 

hilic.end <- combined.data %>% filter(TIME=="end")%>%
  dplyr:: select(DATAFILE,starts_with("COMPOUND_hilic")) 

hilic.pos.end <- hilic.end %>%dplyr:: select(DATAFILE,contains("pos")) 
hilic.neg.end <- hilic.end %>%dplyr:: select(DATAFILE,contains("neg")) 
# RP

rp.start <- combined.data %>% filter(TIME=="beginning")%>%
  dplyr:: select(DATAFILE,starts_with("COMPOUND_RP")) 

rp.pos.start <- rp.start %>%dplyr:: select(DATAFILE,contains("pos")) 
rp.neg.start <- rp.start %>%dplyr:: select(DATAFILE,contains("neg")) 

rp.end <- combined.data %>% filter(TIME=="end")%>%
  dplyr:: select(DATAFILE,starts_with("COMPOUND_RP")) 

rp.pos.end <- rp.end %>%dplyr:: select(DATAFILE,contains("pos")) 
rp.neg.end <- rp.end %>%dplyr:: select(DATAFILE,contains("neg")) 


compound.data.list$hilic.neg.start  <- hilic.neg.start
compound.data.list$hilic.neg.end  <- hilic.neg.end
compound.data.list$hilic.pos.start  <- hilic.pos.start
compound.data.list$hilic.pos.end  <- hilic.pos.end
compound.data.list$rp.neg.start  <- rp.neg.start
compound.data.list$rp.neg.end  <- rp.neg.end
compound.data.list$rp.pos.start  <- rp.pos.start
compound.data.list$rp.pos.end  <- rp.pos.end
# *********** Hybrid Imputation check the error ******************** #



# start logging process by creating a logging file
flog.appender(appender.tee(paste0(output.path,"logFile_hybrid.log")))
flog.threshold(DEBUG)
# propotions of missigness
miss_proportions <- c(0.02)
#miss_proportions <- c(0.1)

proportions_df <- missingness_proportions(miss_proportions=miss_proportions)



# number of iterations
n_iterations <- 1
iteration_counter <- 1
total_iterations <-  n_iterations * nrow(proportions_df)

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
  proportion_results <- foreach(h=1:nrow(proportions_df), .combine="rbind", .inorder=FALSE) %do% {
    
    # seperate the percenatge of missigness per type of misssigness 
    mcar_miss  <- proportions_df$MCAR[h]
    mnar_miss  <- proportions_df$MNAR[h]
    mar_miss   <- proportions_df$MAR[h]
    total_miss <- proportions_df$Total[h]
    
    
    # simulate different type of missgness MCAR,MNAR & MAR using the different percenatges of missigness
    miss_data <- simulate_missingness(data=pool.metabolites, mcar=mcar_miss, mnar=mnar_miss, mar=mar_miss)
    # impute missing values using different imputation methods
    
      
      # impute the data
      
      tryCatch({
        flog.debug(paste('Starting imputing:','DatasetNumber',random.pick, sep=' '))
        start <- Sys.time ()
        imputed_data <-Run.miss.data(miss_data)
        
        
        time_diff <- Sys.time () - start
        # log the computational time per method
        flog.debug(paste('Total impute time:',round(time_diff) ,'DatasetNumber',random.pick,sep=' '))
        
        #COMPARING RESULTS USING ROOT MEAN SQUARE ERROR and R square adjusted
        
        #results_Rsquare <- Rsquare_adjusted(original.data = simulated_data, missing.data = miss_data, imputed.data = imputed_data)
        nrmse_error <- missForest:: nrmse(imputed_data, miss_data, pool.metabolites)
        nrmse_error2 <- differences_models (pool.metabolites, miss_data, imputed_data)
        #nrmse_error2 <- missForest::mixError(imputed_data,miss_data,pool.metabolites)
        
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
        flog.debug(paste('Total nrmse time:',round(time_diff2) ,'DatasetNumber:', random.pick, sep=' '))
        
        # create a data frame with the results
        results_f <- cbind(data.frame(Method="HybridImputation", Miss=total_miss,MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss,DatasetNumber = random.pick, n_row=nrow(pool.metabolites),Iteration=iteration, Time_total=time_diff), data.frame(NRMSE=nrmse_error),data.frame(ModelDiff=nrmse_error2))
        
        time_diff3 <- Sys.time () - start
        flog.debug(paste('Total result binding time:',round(time_diff3) ,'DatasetNumber:', random.pick, sep=' '))
        
        iteration_counter <- iteration_counter + 1
        results_f 
      }, error=function(x) {
        flog.debug(x)
        results_f <- cbind(data.frame(Method="HybridImputation", Miss=total_miss, MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss,DatasetNumber = random.pick,n_row=nrow(pool.metabolites) ,Iteration=iteration,Time_total=time_diff), data.frame(NRMSE=NA),data.frame(ModelDiff=NA))
        return(results_f)
      })
  
  }
  proportion_results
  
}

write.csv(x=FULL_results, file=paste0(output.path, "hubridImputation.csv"), row.names=FALSE)
head(FULL_results,10)
