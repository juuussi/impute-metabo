path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

simulated_data <- simulate_data(data=reference_data, nrow=100, ncol=100)

