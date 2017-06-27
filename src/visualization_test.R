
library(doMC)
library(ggplot2)

library("reshape2")
library("dplyr")
library(tidyr)
library(matrixStats)
registerDoMC(cores=1)
####
path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))


#fileNames <- Sys.glob(paste0(path,"results/","results_RF_KNN_big_sample.csv"))


fileNames <- Sys.glob(paste0(path,"results/","results_ALLMETHODS_ZERO_062017_5_3p.csv"))

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

 results <- read.csv(fileNames, stringsAsFactors = FALSE)
# head(results)
# str(results)
# 
# p <- ggplot(data=results, aes(c(MNAR, NRMSE)))
# p <- p + geom_boxplot()
# p <- p + theme_bw() + facet_grid(.~ Method)
# 
# p
#


results_new <- data.frame(Method=results$Method,Type = results$MissModel,Percentage= results$Miss,Size=results$Cols,Error=results$NRMSE,Time_total = results$Time_total)

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


new2 <- results_new %>%
  group_by(Method,Percentage,Type) %>%
  summarise(Total_mean = mean(Error,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Total_mean)%>% # it spreads the error through the columns
  arrange(Method)

colnames(new2)[3:ncol(new2)] <- paste("Mean_Error",colnames(new2)[3:ncol(new2)], sep = "_" )
# calculate the total sum error of the columns with the mean
new2$Total_Error <- rowSums(new2 %>% select(starts_with("Mean_")))

# pick the columns with the name mean
mean_col <- colnames(new2)[startsWith(colnames(new2), "Mean")]

new2$Total_MEAN_Error <- apply(new2[,mean_col],1, mean , na.rm = TRUE) # same result as before.

new2$Total_SD_Error <- apply(new2[,mean_col],1, sd , na.rm = TRUE) # same result as before.

View(new2)

#### CALCULATION OF TIME
new3 <- results_new %>%
  group_by(Method,Percentage,Type) %>%
  summarise(Time = sum(Time_total,na.rm=TRUE))%>%
  ungroup() %>%
  spread(Type,Time)%>% # it spreads the error through the columns
  arrange(Method)

colnames(new3)[3:ncol(new3)] <- paste("Time_",colnames(new3)[3:ncol(new3)], sep = "_" )
# calculate the total sum error of the columns with the mean
new3$Total_TIME <- rowSums(new3 %>% select(starts_with("Time_")),na.rm = TRUE)

View(new3)


new2$Total_TIME_min <- round((rowSums(new3 %>% select(starts_with("Time_")),na.rm = TRUE))/60,digits = 3)


Imputation_Matrix <- new2


#   

ggplot(data=results_new, aes(x=Type,y=Error,fill=Method) ) +
  geom_bar(stat = "identity",position="dodge")+
  facet_grid(.~ Percentage)
# ###
# 
# library(ggplot2)
# ggplot(data=results_new, aes(x=Method,y=Error ,fill = Type) )+
#   geom_bar(stat = "identity",position="dodge")+
# 
#   facet_grid(Size~ Percentage)

  
   