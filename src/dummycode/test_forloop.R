

# libraries
library(missForest)
library(doMC)
library(futile.logger)
library(pcaMethods)
#####################
registerDoMC(cores=30)

path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
flog.appender(appender.tee(paste0(path,"results/logFile.log")))

flog.threshold(INFO)
data_size <- matrix(c(50,150,200,500,600,1500),nrow=3,ncol = 2)

reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

set.seed(1406)
for (ia in 1:nrow(data_size)){
  
  for (ja in 1:nrow(data_size)) {
    n_iterations <- 3
    
    n_rows <- data_size[ia,1]
    n_cols <- data_size[ja,2]
    #  n_iterations <- 2
    # n_rows <-100
    #  n_cols <- 50
    miss_proportions <- c(0.01,0.1,0.3,0.5,0.6,0.8)
    imputation_methods <- c( "min","mean")
    
    full_results <- foreach(iteration=1:n_iterations, .combine="rbind") %dopar% {
      
      flog.info(paste('ITERATION', iteration, sep=' '))
      
      simulated_data <- simulate_data(data=reference_data, nrow=n_rows, ncol=n_cols)
      
      proportion_results <- foreach(i=1:length(miss_proportions), .combine="rbind") %do% {
        method_results <- foreach(j=1:length(imputation_methods), .combine="rbind") %do% {
          prop <- miss_proportions[i]
          method <- imputation_methods[j]
          flog.info(paste('iteration:', iteration, 'miss proportion:', prop, "method:", method, sep=" "))
          
          miss_data <- simulate_missingness(data=simulated_data, mnar=prop)
          imputed_data <- impute(data=miss_data, methods=method)
          
          
          #COMPARING RESULTS USING ROOT MEAN SQUARE ERROR
          
          results_differences <- Differences_models(data1 = simulated_data,data2 = miss_data,data3 = imputed_data)
          error <- missForest:: nrmse(imputed_data,miss_data,simulated_data)
          results_f <- cbind(data.frame(Method=method, Prop=prop, Iteration=iteration),results_differences, error)
          
          results_f
          
          
        }
        method_results
        
      }  
      
      
      proportion_results
      
    }
    
    
    full_results
    flog.info("**** full results ****")
    write.csv(x=full_results, file=paste0(path, "results/result_NeWmnar_2_",data_size[ia,1],"_",data_size[ja,2], ".csv"), row.names=FALSE)
    
  }
  
  
}