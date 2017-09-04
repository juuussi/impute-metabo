# some statistics

library(mdatools)

library(plotly)
library(ggplot2)
packageVersion('plotly')

# create boxplots
p <- plot_ly(results_new, x = ~results_new$Method, y = ~results_new$Error, color = ~results_new$Type, type = "box",split =~results_new$Percentage,frame =  ~results_new$Percentage) %>%
  layout(boxmode = "group")
p


p <- ggplot(data=results_new, aes(Method, Error))
p <- p + geom_boxplot()
p <- p + theme_bw() + facet_grid(.~Percentage)
p


results_new$Error[results_new$Error>1000] <- NA 
