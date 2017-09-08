
library(doMC)
library(ggplot2)

registerDoMC(cores=1)

path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))


fileNames <- Sys.glob(paste0(path,"results/","*.csv"))


for (fileName in fileNames){
  
  # read data
  
  results <- read.csv(fileName, stringsAsFactors = FALSE)
  head(results)
  str(results)
  
  p <- ggplot(data=results, aes(Method, NRMSE))
  p <- p + geom_boxplot()
  p <- p + theme_bw() + facet_grid(.~Prop)
  
  
  ggsave(filename = paste0(fileName,'_nrmse',".pdf") ,plot = p,width = 20, height = 20, units = "cm")
  
  
}
for (fileName in fileNames){
  
  # read data
  
  results <- read.csv(fileName, stringsAsFactors = FALSE)
  head(results)
  str(results)
  
  p <- ggplot(data=results, aes(Method, results_differences))
  p <- p + geom_boxplot()
  p <- p + theme_bw() + facet_grid(.~Prop)
  
  
  ggsave(filename = paste0(fileName,'_relative_difference',".pdf") ,plot = p,width = 20, height = 20, units = "cm")
  
  
}


