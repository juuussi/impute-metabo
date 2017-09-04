# libraries
library(missForest)
library(doMC)
library(futile.logger)
library(pcaMethods)
library(impute)
#####################
registerDoMC(cores=30)

path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
flog.appender(appender.tee(paste0(path,"results/logFile3.log")))

flog.threshold(DEBUG)
#data_size <- matrix(c(50,150,200,500,600,1500),nrow=3,ncol = 2)

reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

set.seed(1406)

size_iterations <- 1

data_rows <- sample(x=10:50, size=size_iterations)
data_cols <- sample(x=3000:5000, size=size_iterations)

#data_rows <- c(50, 150, 200)
#data_cols <- c(500, 600, 1500)
n_iterations <- 1
#miss_proportions <- c(0.01, 0.05,0.1,0.3)
miss_proportions <- c(0.1,0.2,0.3)

proportions_df <- missingness_proportions(miss_proportions=miss_proportions)

#imputation_methods <- c( "min","mean","PPCA","RF")
#miss_proportions <- seq(from=0.01, to=0.9, by=0.05)
#imputation_methods <- c( "LLS")
imputation_methods <- c( "min","mean","PPCA","RF","LLS","svdImpute","KNNImpute")

total_iterations <- length(data_rows) * length(data_cols) * n_iterations * nrow(proportions_df) * length(imputation_methods)
flog.info(paste('Total number of iterations: ', total_iterations, sep=' '))

iteration_counter <- 1

full_results <- NULL

#for(n_rows in data_rows) {
#  for (n_cols in data_cols) {

full_results <- foreach(r=1:length(data_rows), .combine="rbind") %:%
  foreach(c=1:length(data_cols), .combine="rbind") %dopar% {
    n_rows <- data_rows[r]
    n_cols <- data_cols[c]
    #flog.info(paste('Data size:', n_rows, "x", n_cols,  sep=' '))
    
    iteration_results <- foreach(iteration=1:n_iterations, .combine="rbind") %dopar% {
      #prop <- miss_proportions[i]
      
      
      simulated_data <- simulate_data(data=reference_data, nrow=n_rows, ncol=n_cols)
      
      # proportion_results <- foreach(i=1:length(miss_proportions), .combine="rbind") %do% {
      proportion_results <- foreach(h=1:nrow(proportions_df), .combine="rbind") %do% {
        
        
        #miss_data <- simulate_missingness(data=simulated_data, mnar=prop)
        mcar_miss <- proportions_df$MCAR[h]
        mnar_miss <- proportions_df$MNAR[h]
        mar_miss <- proportions_df$MAR[h]
        total_miss <- proportions_df$Total[h]
        
        #flog.info(paste('Missingness - Total:', total_miss, "MCAR:", mcar_miss, "MAR:", mar_miss, "MNAR:", mnar_miss,  sep=' '))
        
        miss_data <- simulate_missingness(data=simulated_data, mcar=mcar_miss, mnar=mnar_miss, mar=mar_miss)
        
        method_results <- foreach(j=1:length(imputation_methods), .combine="rbind") %do% {
          method <- imputation_methods[j]
          #flog.info(paste('Method:', method,  sep=' '))
          
          #flog.info(paste('iteration:', iteration, 'miss proportion:', prop, "method:", method, sep=" "))
          
          #proportion_results <- foreach(h=1:nrow(proportions_df), .combine="rbind") %do% {
          flog.debug(paste('Total iterations:', iteration_counter, "/", total_iterations, 'Data size:', n_rows, "x", n_cols, 'Method:', method, 'Missingness - Total:', total_miss, "MCAR:", mcar_miss, "MAR:", mar_miss, "MNAR:", mnar_miss, sep=' '))
          
          start <- Sys.time ()
          
          imputed_data <- impute(data=miss_data, methods=method)
          time_diff <- Sys.time () - start
          flog.debug(paste('Total time:',time_diff ,'Method:', method,  sep=' '))
          
          #COMPARING RESULTS USING ROOT MEAN SQUARE ERROR
          
          results_differences <- differences_models(original.data = simulated_data, missing.data = miss_data, imputed.data = imputed_data)
          nrmse_error <- missForest:: nrmse(imputed_data, miss_data, simulated_data)
          miss_model <- ""
          if(mnar_miss>0){
            miss_model <- paste0(miss_model," MNAR ")
          }
          if(mar_miss>0){
            miss_model <- paste0(miss_model," MAR ")
          }
          if(mcar_miss>0){
            miss_model <- paste0(miss_model," MCAR ")
          }
          
          results_f <- cbind(data.frame(Method=method, Miss=total_miss,MissModel=miss_model, MNAR=mnar_miss, MCAR=mcar_miss, MAR=mar_miss, Iteration=iteration, Rows=n_rows, Cols=n_cols, Total_data=n_rows*n_cols), data.frame(Result_diff=results_differences, NRMSE=nrmse_error))
          
          iteration_counter <- iteration_counter + 1
          results_f
          #}
          #proportion_results
          
        }
        method_results
        
      }  
      
      proportion_results
      
    }
    
    #full_results
    #flog.info("**** full results ****")
    #write.csv(x=full_results, file=paste0(path, "results/result_NeWmnar_",n_rows,"_",n_cols, ".csv"), row.names=FALSE)
    #  if (is.null(full_results)) {
    #    full_results <- iteration_results
    #  } else {
    #    full_results <- rbind(full_results, iteration_results)
    #  }
    #}
  }


write.csv(x=full_results, file=paste0(path, "results/results_all.csv"), row.names=FALSE)
head(full_results,10)
