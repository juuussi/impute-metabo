# simulate_missingness2 <- function(data, mcar=0, mar=0, mnar=0, mnar.type="left") {
#   
#   if(class(data) != "matrix") {
#     stop("Variable data should be a matrix.")
#   }
#   if ( ((mcar + mar + mnar) > 1) || ((mcar + mar + mnar) < 0)) {
#     stop("Sum of mcar, mar and mnar should be between 0 and 1.")
#   }
#   
#   simulated_data <- data
#   if (mcar > 0){
#     mcar_distribution   = runif(nrow(data)*ncol(data), min=0, max=1)
#     simulated_data = matrix(ifelse(mcar_distribution<mcar, NA, data), nrow=nrow(data), ncol=ncol(data))
#   }
# 
#   if (mnar > 0) {
#     added_mnar <- 0
#     initial_nas <- sum(colSums(is.na(simulated_data)))
#     current_nas <- initial_nas
#     simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
#     # condition of porpotion of missigness
#     while (((current_nas - initial_nas) / simulated_data_n) < mnar) {
#       
#       # Select random variable
#       variable_index <- sample(1:ncol(simulated_data), 1)
#       
#       # What percentage of variable to set missing
#       cut_percentage <- 0
#       while (cut_percentage <= 0) {
#         cut_percentage <- rchisq(1, df=1) / 30
#       }
#       cut_percentage <- min(cut_percentage, 1)
#       
#       # How many values to set missing
#       cut_index <- floor(cut_percentage * nrow(simulated_data))
#       
#       sorted_variable <- sort(simulated_data[,variable_index])
#       
#       # Set values to missing
#       if (mnar.type == "right") {
#         # Corresponding cut-off point for values
#         cut_point <- sorted_variable[length(sorted_variable) - cut_index]
#         
#         simulated_data[simulated_data[,variable_index] > cut_point, variable_index] <- NA
#       } else {
#         # Corresponding cut-off point for values
#         cut_point <- sorted_variable[cut_index]
#         
#         simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
#       }
#       
#       # Counter to check how much MNAR missingness has been added to data
#       current_nas <- sum(colSums(is.na(simulated_data)))
#     }
#   }
#   simulated_data
#   
# }


# simulate_missingness2 <- function(data, mcar=0, mar=0, mnar=0, mnar.type="left") {
#   
#   if(class(data) != "matrix") {
#     stop("Variable data should be a matrix.")
#   }
#   if ( ((mcar + mar + mnar) > 1) || ((mcar + mar + mnar) < 0)) {
#     stop("Sum of mcar, mar and mnar should be between 0 and 1.")
#   }
#   
#   simulated_data <- data
#   if (mcar > 0){
#     mcar_distribution   = runif(nrow(data)*ncol(data), min=0, max=1)
#     simulated_data = matrix(ifelse(mcar_distribution<mcar, NA, data), nrow=nrow(data), ncol=ncol(data))
#   }
#   # MAR , the missigness doent reach 0,5
#   if (mar > 0){
#     initial_nas <- sum(colSums(is.na(simulated_data)))
#     current_nas <- initial_nas
#     simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
#     # current_nan_percentage <- 0
#     columns <- NULL
#     # number of columns to choose
#     support <- 1:(ncol(simulated_data)-1)
#     while(length(support) > 0 && ((current_nas - initial_nas) / simulated_data_n) < mar ) {
#       # try to take unique column each time
#       column <- sample(x=support, size=1)
#       support <- setdiff(1:(ncol(simulated_data)-1), columns)
#       
#       columns <- c(columns, column)
#       #print(columns)
#       # condition: 
#       # calculate the mean of the choosen column
#       mn <- mean(simulated_data[,column])
#       # dependency : remove the values on the next column(from the one it pick) above the mean 
#       simulated_data[simulated_data[,column+1] > mn, column+1] <- NA 
#       # check the current missing porpotion
#       current_nan <- length(which(is.na(simulated_data)))/length(simulated_data)
#       print(current_nan)
#     }
#   }
#   



#   if (mnar > 0) {
#     added_mnar <- 0
#     initial_nas <- sum(colSums(is.na(simulated_data)))
#     current_nas <- initial_nas
#     simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
#     # condition of porpotion of missigness
#     while (((current_nas - initial_nas) / simulated_data_n) < mnar) {
#       
#       # Select random variable
#       variable_index <- sample(1:ncol(simulated_data), 1)
#       
#       # What percentage of variable to set missing
#       cut_percentage <- 0
#       while (cut_percentage <= 0) {
#         cut_percentage <- rchisq(1, df=1) / 30
#       }
#       cut_percentage <- min(cut_percentage, 1)
#       
#       # How many values to set missing
#       cut_index <- floor(cut_percentage * nrow(simulated_data))
#       
#       sorted_variable <- sort(simulated_data[,variable_index])
#       
#       # Set values to missing
#       if (mnar.type == "right") {
#         # Corresponding cut-off point for values
#         cut_point <- sorted_variable[length(sorted_variable) - cut_index]
#         
#         simulated_data[simulated_data[,variable_index] > cut_point, variable_index] <- NA
#       } else {
#         # Corresponding cut-off point for values
#         cut_point <- sorted_variable[cut_index]
#         
#         simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
#       }
#       
#       # Counter to check how much MNAR missingness has been added to data
#       current_nas <- sum(colSums(is.na(simulated_data)))
#     }
#   }
#   
#   simulated_data
#   
# }
# compare_results2 <- function(x, y, missing) {
#   missing_index <- which(is.na(as.vector(missing)))
#   x.val <- as.vector(x)[missing_index]
#   y.val <- as.vector(y)[missing_index]
#   flog.info(paste('missing:', paste(missing_index, collapse=", "), '\norig:', paste(x.val, collapse=", "),'\nimpu:', paste(y.val, collapse=", "),sep=' '))
#   
#   c <- cor.test(x.val, y.val)
#   
#   flog.info(paste("Cor:", c$estimate, "P:", c$p.value, sep=" "))
#   data.frame(Cor=c$estimate, CorP=c$p.value)
#   
#   
# }