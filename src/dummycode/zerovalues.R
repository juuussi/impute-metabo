simulated_data.rep <-replicate(n=4,simulate_data(reference_data,10,100))
simulated_data <- simulate_data(reference_data,10,100)
apply(simulated_data.rep,MARGIN= 2,FUN = mean)
sd(apply(simulated_data.rep,2,mean))

miss_data2 <- simulate_missingness(simulated_data.rep[,,2],mcar = 0.9)

#re <- cbind(data.frame(colNAN = col_nan,perNAN =percentage_nan))

replace_NA_to_0 <- function(missing.data) {
  col_nan <- NULL
  total_length_col <- NULL
  percentage_nan<- NULL
  for (i in col(missing.data)){
    # check how many rows in every column have NA
    col_nan[i] <- length(which(is.na(missing.data[,i])))
    total_length_col[i] <- length(missing.data[,i])
    #calculate the percentage of NA of every column
    percentage_nan[i] <- col_nan[i]/total_length_col[i]
    # if the percenatge of NA every column is above 0.8 replace with zero
    if (percentage_nan[i] >= 0.8){
      missing.data[is.na(missing.data[,i]),i] <- 0
    }
  }
  missing.data
}

missNew_data <- replace_NA_to_0(miss_data)
##################
for (i in col(miss_data)){
  # check how many rows in every column have NA
  col_nan[i] <- length(which(is.na(miss_data[,i])))
  total_length_col[i] <- length(miss_data[,i])
  #calculate the percentage of NA of every column
  percentage_nan[i] <- col_nan[i]/total_length_col[i]
  
  if (percentage_nan[i] >= 0.8){
    miss_data[is.na(miss_data[,i]),i] <- 0
  }
  
  
}




