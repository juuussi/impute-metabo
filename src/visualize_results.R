library(doMC)
library(ggplot2)

registerDoMC(cores=1)

path <- "~/projects/impute-metabo/"

source(paste0(path,"src/functions.R"))

results <- read.csv(paste0(path, "results/result.csv"), stringsAsFactors = FALSE)
head(results)
str(results)

p <- ggplot(data=results, aes(Method, results_rmse))
p <- p + geom_boxplot()
p <- p + theme_bw() + facet_grid(.~Prop)

p

