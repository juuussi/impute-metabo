
# Detect missingness

# This method creates binary matrix, in which the dataset consists of indicator variables 
  # where a 1 is given if a value is present, and 0 if it isn't.
  #Correlating these with each other and the original data can help determine 
  #if variables tend to be missing together (MAR) or not (MCAR).

#Elements of x are 1 if a value in the miss_data is missing and 0 if non-missing.
  

x <- as.data.frame(abs(is.na(miss_data)))

head(miss_data)
head(x)

#Extracting variables that have some missing values.
y <- x[which(sapply(x, sd) > 0)]
cor(y)


#Now, looking at the relationship between the 
#presence of missing values in each variable and the observed values
#in other variables:


corr_matrix <-cor(miss_data, y, use="pairwise.complete.obs")



  