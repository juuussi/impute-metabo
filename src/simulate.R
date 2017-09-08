# libraries needed
library(missForest)
library(doMC)
library(futile.logger)
library(pcaMethods)
###################################################################################
# start parallel
registerDoMC(cores=35)
# choose path
path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
# start logging process by creating a logging file
flog.appender(appender.tee(paste0(path,"results/logFile_finalres.log")))
flog.threshold(DEBUG)
###################################################################################
# use dummy reference data (combination of different metabolomics data)
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

###################################################################################

#define the sample space
set.seed(1406)
size_iterations <- 1
# # 
seq_rows <- seq(from=20, to=30, by= 2)
seq_cols <-  seq(from=40, to=100, by= 5)

data_rows <- sample(x=seq_rows, size=size_iterations,replace = TRUE)
data_cols <- sample(x=seq_cols, size=size_iterations,replace = TRUE)
n_iterations <- 1

################################################################################

# define the percenatge of missigness
#miss_proportions <- c(0.01, 0.05,0.1,0.3)
#miss_proportions <- c(0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
#miss_proportions <- c(0.1,0.3,0.6,0.8)
miss_proportions <- c(0.1)

proportions_df <- missingness_proportions(miss_proportions=miss_proportions)
################################################################################

# Define the methods

#imputation_methods <- c( "LLS")

imputation_methods <- c("min","mean", "RF")

#imputation_methods <- c( "RF","min","mean","KNNImpute")

#imputation_methods <- c("PPCA", "RF","min","mean","LLS","KNNImpute","BPCA","svdImpute")
################################################################################
# calculate the total time of iterations
#total_iterations <- length(data_rows) * length(data_cols) * n_iterations * nrow(proportions_df) * length(imputation_methods)
total_iterations <- size_iterations * n_iterations * nrow(proportions_df) * length(imputation_methods)

################################################################################

# log the total number of iterations
flog.info(paste('Total number of iterations: ', total_iterations, sep=' '))

iteration_counter <- 1

################################################################################

full_results <- NULL
full_results <- foreach(r=1:length(data_rows), .combine="rbind") %:%
  
  foreach(c=1:length(data_cols), .combine="rbind") %do% {
    n_rows <- data_rows[r]
    n_cols <- data_cols[c]
    
    iteration_results <- foreach(iteration=1:n_iterations, .combine="rbind") %do% {
      
      # simulate the data from the reference data(using multivariate normal distribution)
      simulated_data <- simulate_data(data=reference_data, nrow=n_rows, ncol=n_cols)
      
      
      # choose the percenatge of missigness 
      proportion_results <- foreach(h=1:nrow(proportions_df), .combine="rbind", .inorder=FALSE) %do% {
        
        # seperate the percenatge of missigness per type of misssigness 
        mcar_miss <- proportions_df$MCAR[h]
        mnar_miss <- proportions_df$MNAR[h]
        mar_miss <- proportions_df$MAR[h]
        total_miss <- proportions_df$Total[h]
        
        
        # simulate different type of missgness MCAR,MNAR & MAR using the different percenatges of missigness
        miss_data <- simulate_missingness(data=simulated_data, mcar=mcar_miss, mnar=mnar_miss, mar=mar_miss)
        
        # impute missing values using different imputation methods
        method_results <- foreach(j=1:length(imputation_methods), .combine="rbind", .inorder=FALSE) %do% {
          
          method <- imputation_methods[j] 
          
          # log the total iteration ,method , percenatges and types of missgness 
          flog.debug(paste('Total iterations:', iteration_counter, "/", total_iterations, 'Data size:', n_rows, "x", n_cols, 'Method:', method, 'Missingness - Total:', total_miss, "MCAR:", mcar_miss, "MAR:", mar_miss, "MNAR:", mnar_miss, sep=' '))
          
          # calculate the computational time per method
          
          # replace missing values with zero
          miss_data <- replace_NA_to_0(miss_data)
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
            nrmse_error <- missForest:: nrmse(imputed_data, miss_data, simulated_data)
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
            results_f <- cbind(data.frame(Method=method, Miss=total_miss,MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss, Iteration=iteration, Rows=n_rows, Cols=n_cols, Total_data=n_rows*n_cols,Time_total=time_diff), data.frame(NRMSE=nrmse_error))
            
            time_diff3 <- Sys.time () - start
            flog.debug(paste('Total result binding time:',round(time_diff3) ,'Method:', method, sep=' '))
            
            iteration_counter <- iteration_counter + 1
            results_f 
          }, error=function(x) {
            flog.debug(x)
            results_f <- cbind(data.frame(Method=method, Miss=total_miss, MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss, Iteration=iteration, Rows=n_rows, Cols=n_cols, Total_data=n_rows*n_cols,Time_total=time_diff), data.frame(NRMSE=NA))
            return(results_f)
          })
        }
        method_results
      }  
      
      proportion_results
    }
    
  }




write.csv(x=full_results, file=paste0(path, "results/finalres0817.csv"), row.names=FALSE)
head(full_results,10)
