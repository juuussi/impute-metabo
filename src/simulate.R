
# libraries
library(missForest)
library(doMC)
library(futile.logger)
library(pcaMethods)
#####################
registerDoMC(cores=1)

path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
flog.appender(appender.tee(paste0(path,"results/logFile.log")))

flog.threshold(INFO)


reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

set.seed(1406)

n_iterations <- 2
n_rows <- 150
n_cols <- 10

miss_proportions <- c(0.01, 0.05, 0.1, 0.2, 0.3)
imputation_methods <- c("RF", "BPCA","mean", "min")

full_results <- foreach(iteration=1:n_iterations, .combine="rbind") %do% {
  
  flog.info(paste('ITERATION', iteration, sep=' '))
  
  simulated_data <- simulate_data(data=reference_data, nrow=n_rows, ncol=n_cols)
  
  proportion_results <- foreach(i=1:length(miss_proportions), .combine="rbind") %do% {
    method_results <- foreach(j=1:length(imputation_methods), .combine="rbind") %do% {
      prop <- miss_proportions[i]
      method <- imputation_methods[j]
      flog.info(paste('iteration:', iteration, 'miss proportion:', prop, "method:", method, sep=" "))
      #flog.info(paste('miss proportion', prop, sep=' '))
      #flog.info(paste('method', method, sep=' '))
      
      
      miss_data <- simulate_missingness(data=simulated_data, mcar=prop)
      imputed_data <- impute(data=miss_data, methods=method)
      # #####################################################################
      # # COMPARING RESULTS USING CORRELATION
      # #compared_results <- compare_results(simulated_data, imputed_data)
      # flog.info("**** COMPARING ****")
      # #compared_results <- compare_results2(simulated_data, imputed_data, missing=miss_data)
      # flog.info("**** DONE ****")
      # 
      # #results <- cbind(data.frame(Method=method, Prop=prop, Iteration=iteration), compared_results)
      # flog.info("**** CBIND DONE ****")
      # 
      ######################################################################
      
      #COMPARING RESULTS USING ROOT MEAN SQUARE ERROR
      
      results_rmse <- RMSE_simulated(data1 = simulated_data,data2 = miss_data,data3 = imputed_data)
      results_f <- cbind(data.frame(Method=method, Prop=prop, Iteration=iteration),results_rmse )
      
      flog.info("**** CBIND DONE ****")
      results_f
      
      
    }
    flog.info("combining results")
    method_results
  }  
  flog.info("combining propotions")
  
  proportion_results
  
}


full_results
write.csv(x=full_results, file=paste0(path, "results/result.csv"), row.names=FALSE)
