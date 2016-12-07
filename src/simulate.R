library(missForest)
library(doMC)

registerDoMC(cores=1)

path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

set.seed(1406)

n_iterations <- 2
n_rows <- 50
n_cols <- 10

miss_proportions <- c(0.01, 0.05, 0.1, 0.2, 0.3)
imputation_methods <- c("RF", "mean", "min")

full_results <- foreach(iteration=1:n_iterations, .combine="rbind") %dopar% {

  simulated_data <- simulate_data(data=reference_data, nrow=n_rows, ncol=n_cols)

  proportion_results <- foreach(i=1:length(miss_proportions), .combine="rbind") %do% {
    method_results <- foreach(j=1:length(imputation_methods), .combine="rbind") %do% {
      prop <- miss_proportions[i]
      method <- imputation_methods[j]
      
      miss_data <- simulate_missingness(data=simulated_data, mcar=prop)
      imputed_data <- impute(data=miss_data, methods=method)
      
      compared_results <- compare_results(simulated_data, imputed_data)
      results <- cbind(data.frame(Method=method, Prop=prop, Iteration=iteration), compared_results)
    }
    method_results
  }  
}

summary(full_results)
