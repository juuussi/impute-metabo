
rm(list = ls())
## list of libraries

library(doMC)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(matrixStats)
registerDoMC(cores=1)
#### path of theproject
path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))



## load the file
fileNames <- Sys.glob(paste0(path,"results/","results_2it_zero.csv"))

# read cvs file
results <- read.csv(fileNames, stringsAsFactors = FALSE)

# create anew data frame with some of the variables from the result file
results_new <- data.frame(Method=results$Method,Type = results$MissModel,Percentage= results$Miss,Size=results$Cols,Error=results$NRMSE,Time_total = results$Time_total,Iteration=results$Iteration)

###### calculate the MEAN error per type an per percentage and per method
new_mean <- results_new %>%
  group_by(Method,Percentage,Type) %>%
  dplyr::summarise(Total_Mean = mean(Error,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Total_Mean)%>% # it spreads the error through the columns
  arrange(Method)
# rename the columns of the MCAR ..etc 
colnames(new_mean)[3:ncol(new_mean)] <- paste("Mean_Error",colnames(new_mean)[3:ncol(new_mean)], sep = "_" )

##### calclulate the mean error for  ALL the percenatges in total per method
new_MEAN_TOTAL <- results_new %>%
  group_by(Method,Percentage) %>%
  dplyr::summarise(Total_Mean = mean(Error,na.rm=TRUE))#%>%

new_mean <-  cbind.data.frame(new_mean,Total_mean=new_MEAN_TOTAL$Total_Mean) 

MEAN_TOTAL <- results_new %>%
  group_by(Method) %>%
  dplyr::summarise(Total_Mean = mean(Error,na.rm=TRUE))
# ## add a column with total percenatge of the mean error
MEAN_TOTAL$Percentage <- "comb_Perc"
MEAN_TOTAL <- MEAN_TOTAL %>% dplyr::select(Method,Percentage,everything())

## create the TOTAL_MEANmatrix with equal number of columns as the new mean
sad_MEAN_TOTAL <- data.frame(matrix(NA, ncol=ncol(new_mean), nrow=nrow(MEAN_TOTAL)))
colnames(sad_MEAN_TOTAL) <- c(colnames(new_mean))
sad_MEAN_TOTAL[,c(1,2,ncol(new_mean))] <- MEAN_TOTAL

ImputationMatrix1 <- rbind(new_mean, sad_MEAN_TOTAL)%>%
  arrange(Method)

###### calculate the SD error per type an per percentage and per method

new_sd <- results_new %>%
  group_by(Method,Percentage,Type) %>%
  dplyr::summarise(Total_sd = sd(Error,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Total_sd)%>% # it spreads the error through the columns
  arrange(Method)
#change the column names
colnames(new_sd)[3:ncol(new_sd)] <- paste("sd_Error",colnames(new_sd)[3:ncol(new_sd)], sep = "_" )
##### calclulate the SD error for  ALL the percenatges in total per method

new_SD_TOTAL <- results_new %>%
  group_by(Method,Percentage) %>%
  dplyr::summarise(Total_sd = sd(Error,na.rm=TRUE))#%>%
  
new_sd <-  cbind.data.frame(new_sd,Total_sd=new_SD_TOTAL$Total_sd) 

# ## add a column with total percenatge of the sd error
SD_TOTAL <- results_new %>%
  group_by(Method) %>%
  dplyr::summarise(Total_sd = sd(Error,na.rm=TRUE))

SD_TOTAL$Percentage <- "comb_Perc"
SD_TOTAL <- SD_TOTAL %>% dplyr::select(Method,Percentage,everything())


# ## create the SD_TOTAL matrix with equal number of columns as the new sd
sad_SD_TOTAL <- data.frame(matrix(NA, ncol=ncol(new_sd), nrow=nrow(SD_TOTAL)))
colnames(sad_SD_TOTAL) <- c(colnames(new_sd))

sad_SD_TOTAL[,c(1,2,ncol(new_sd))] <- SD_TOTAL

ImputationMatrix2 <- rbind(new_sd, sad_SD_TOTAL)%>%
  arrange(Method)



## combine the 2 matrices
tmp <- ImputationMatrix2[,3:ncol(ImputationMatrix2)]

Imputation_Matrix <- cbind(ImputationMatrix1,tmp)











# #### CALCULATION OF TIME####

## Mean time
new_TIMEMean <- results_new %>%
  group_by(Method,Percentage) %>%
  dplyr::summarise(Time_mean = mean(Time_total,na.rm=TRUE))
MEAN_TOTAL_Time <- results_new %>%
  group_by(Method) %>%
  dplyr::summarise(Time_mean = mean(Error,na.rm=TRUE))

MEAN_TOTAL_Time$Percentage <- "comb_Perc"
MEAN_TOTAL_Time <- MEAN_TOTAL_Time %>% dplyr::select(Method,Percentage,everything())

new_TIMEMean$Percentage <- as.character(new_TIMEMean$Percentage)
new_TIMEMean <- results_new %>%
  group_by(Method,Percentage) %>%
  dplyr::summarise(Time_mean = sd(Error,na.rm=TRUE))#%>%

timematrixMEAN<- rbind.data.frame(new_TIMEMean, MEAN_TOTAL_Time)%>%dplyr::arrange(Method)


# SD TIME

new_TimeSD <-results_new %>%
  group_by(Method,Percentage) %>%
  dplyr::summarise(Time_sd = sd(Time_total,na.rm=TRUE))

SD_TOTAL_Time<- results_new %>%
  group_by(Method) %>%
  dplyr::summarise(Time_sd = sd(Error,na.rm=TRUE))

SD_TOTAL_Time$Percentage <- "comb_Perc"
SD_TOTAL_Time <- SD_TOTAL_Time %>% dplyr::select(Method,Percentage,everything())

new_TimeSD$Percentage <- as.character(new_TimeSD$Percentage)
timematrixSD<- rbind.data.frame(new_TimeSD, SD_TOTAL_Time)%>%dplyr::arrange(Method)







### final imputation matrix
types <- as.character(unique(results_new$Type)) 

for (type in types){
  
  Imputation_Matrix <- Imputation_Matrix %>% dplyr:: select(Method,Percentage,Total_mean,Total_sd,ends_with(type),everything())
  
}
Imputation_Matrix <- cbind.data.frame(Imputation_Matrix,timematrixMEAN[,ncol(timematrixMEAN)],timematrixSD[,ncol(timematrixSD)])


