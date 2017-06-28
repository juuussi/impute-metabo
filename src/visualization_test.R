

## list of libraries

library(doMC)
library(ggplot2)
library("reshape2")
library("dplyr")
library(tidyr)
library(matrixStats)
registerDoMC(cores=1)
#### path of theproject
path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))


#fileNames <- Sys.glob(paste0(path,"results/","results_RF_KNN_big_sample.csv"))

## load the file
fileNames <- Sys.glob(paste0(path,"results/","results_ALLMETHODS_ZERO_062017_5_3p.csv"))


 results <- read.csv(fileNames, stringsAsFactors = FALSE)


results_new <- data.frame(Method=results$Method,Type = results$MissModel,Percentage= results$Miss,Size=results$Cols,Error=results$NRMSE,Time_total = results$Time_total)

###### calculate the error per type an per percenatge
new1 <- results_new %>%
  group_by(Method,Percentage,Type) %>%
  dplyr::summarise(Total_mean = mean(Error,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Total_mean)%>% # it spreads the error through the columns
  arrange(Method)

colnames(new1)[3:ncol(new1)] <- paste("Mean_Error",colnames(new1)[3:ncol(new1)], sep = "_" )
# calculate the total sum error of the columns with the mean
#new1$Total_Error <- rowSums(new1 %>% select(starts_with("Mean_")))



##### calclulate the error for  ALL the percenatges in total per method
new2 <- results_new %>%
  group_by(Method,Type) %>%
  dplyr::summarise(meanMethod = mean(Error,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,meanMethod)%>% # it spreads the error through the columns
  arrange(Method)
colnames(new2)[2:ncol(new2)] <- paste("Mean_Error",colnames(new2)[2:ncol(new2)], sep = "_" )

## add a column percenatge in total1
new2$Percentage <- "combined_Perc"
new2 <- new2 %>% select(Method,Percentage,everything())



#### CALCULATION OF TIME
new3 <- results_new %>%
  group_by(Method,Type) %>%
  dplyr::summarise(Time = mean(Time_total,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Time)%>% # it spreads the error through the columns
  arrange(Method)

colnames(new3)[2:ncol(new3)] <- paste("Mean_Error",colnames(new3)[2:ncol(new3)], sep = "_" )
## add a column percenatge in new3
new3$Percentage <- "combined_Time"
new3 <- new3 %>% select(Method,Percentage,everything())



## combine the 3 matrices


Imputation_Matrix <- rbind.data.frame(new1,new2,new3)%>%
  arrange(Method)



View(Imputation_Matrix)





















# # pick the columns with the name mean
# mean_col <- colnames(new2)[startsWith(colnames(new2), "Mean")]
# 
# new2$Total_MEAN_Error <- apply(new2[,mean_col],1, mean , na.rm = TRUE) # same result as before.
# 
# new2$Total_SD_Error <- apply(new2[,mean_col],1, sd , na.rm = TRUE) # same result as before.
# 
# View(new2)
# 
# 
# # calculate the total sum error of the columns with the mean
# new3$Total_TIME <- rowSums(new3 %>% select(starts_with("Time_")),na.rm = TRUE)
# 
# View(new3)
# 
# 
# new2$Total_TIME_min <- round((rowSums(new3 %>% select(starts_with("Time_")),na.rm = TRUE))/60,digits = 3)
# 
# 
# 



##########################
# total2 <- results_new %>%
#   group_by(Method,Type) %>%
#   summarise(meanT = mean(Time_total,na.rm=TRUE))%>%
#   ungroup() %>%
#   spread(Type,meanT)%>% # it spreads the error through the columns
#   arrange(Method)
# 
#  colnames(total2)[2:ncol(total2)] <- paste("Mean_TIME",colnames(total2)[3:ncol(total2)], sep = "_" )
# 
#  View(total1)
#  
#  View(total2)
#   library(gdata)
#  concat_data <- cbindX(total1, total2)
#  library(plyr)
#  test <- rbind.fill(total1,total2)
#  results_merged <- merge(x=total1, y=total2, by.x="Method", by.y="Method", all.x=TRUE, all.y=FALSE)
# results_merged <- cbind(results_merged[,1], data.frame(Percentage=rep("Total", times=nrow(results_merged)), results_merged[,2:ncol(results_merged)]))
# results_merged 
# head(Imputation_Matrix)
#   
ggplot(data=results_new, aes(x=Type,y=Error,fill=Method) ) +
  geom_bar(stat = "identity",position="dodge")+
  facet_grid(.~ Percentage)
# ###
# View(results_new)
# 
# new <- results_new %>%
#   group_by(Method,Type) %>%
#   summarise(Total_mean = mean(Error, na.rm = TRUE) ,Total_sd = sd(Error,na.rm = TRUE))%>%
#   gather(Statistic, Value, c(Total_mean,Total_sd))%>%
#   spread(Statistic, Value)%>% # it spreads the error through the columns
#   arrange(Percentage,Method)
# 
# colnames(new)[3:ncol(new)] <- paste("Mean_Error",colnames(new)[3:ncol(new)])
# 
# results_new <- data.frame(Method=results$Method,Type = results$MissModel,Percentage= results$Miss,Size=results$Cols,Error=results$NRMSE)




# # read data
# 
# results <- read.csv(fileNames, stringsAsFactors = FALSE)
# head(results)
# str(results)
# 
# p <- ggplot(data=results, aes(Method, NRMSE,fill= MissModel))
# p <- p + geom_boxplot()
# p <- p + theme_bw() + facet_grid(.~Miss)
# 
# p
# 

# read data
# # head(results)
# str(results)
# 
# p <- ggplot(data=results, aes(c(MNAR, NRMSE)))
# p <- p + geom_boxplot()
# p <- p + theme_bw() + facet_grid(.~ Method)
# 
# p
#

# library(ggplot2)
# ggplot(data=results_new, aes(x=Method,y=Error ,fill = Type) )+
#   geom_bar(stat = "identity",position="dodge")+
# 
#   facet_grid(Size~ Percentage)

  
   