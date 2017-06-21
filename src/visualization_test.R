
library(doMC)
library(ggplot2)
# install.packages("ggplot2")
# install.packages("resha
library("dplyr")
library(tidyr)
# install.packages("devtools")
# library(devtools)
# install_github("wilkox/ggfittext")
# install_github("wilkox/treemapify")

registerDoMC(cores=1)
####
path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))


#fileNames <- Sys.glob(paste0(path,"results/","results_RF_KNN_big_sample.csv"))


fileNames <- Sys.glob(paste0(path,"results/","results_ALLMETHODS_ZERO_062017_4.csv"))

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
results_new <- data.frame(Method=results$Method,Type = results$MissModel,Percentage= results$Miss,Size=results$Cols,Error=results$NRMSE)

View(results_new)

new <- results_new %>%
  group_by(Method,Type,Percentage) %>%
  summarise(Total_mean = mean(Error, na.rm = TRUE ))%>%
  spread(Type,Total_mean)%>%
  arrange(Percentage,Method)

colnames(new)[3:ncol(new)] <- paste("Mean",colnames(new)[3:ncol(new)])






View(new)

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

  
   