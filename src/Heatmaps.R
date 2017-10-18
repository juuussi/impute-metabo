

library(dplyr)
library(ggplot2)
##################################################################################################

## Heat map 1,Total Mean

totalmean <-Imputation_Matrix %>%
  select(Method,Percentage,Total_mean)


totalmean$Method <-factor(totalmean$Method,labels = c("BPCA","KNN","LLS","mean","min","PPCA","RF","SVD"))
totalmean$Percentage <- factor(totalmean$Percentage )

quantiles_vec <- quantile(totalmean$Total_mean,na.rm = T)
Q1 <-round(quantiles_vec[[1]],digits = 2)
Q2 <-round(quantiles_vec[[2]],digits = 2)
Q3 <-round(quantiles_vec[[3]],digits = 2)
Q4 <-round(quantiles_vec[[4]],digits = 2)
Q5 <-round(quantiles_vec[[5]],digits = 2)



p.breaks <- c(Q5,Q4,Q3, Q2, Q1)

totalmean$Breaks <- cut(x=totalmean$Total_mean, breaks=p.breaks, labels=paste0("< ", as.character(p.breaks)[1:(length(p.breaks)-1)]), include.lowest=TRUE)
palette <- c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")
#palette <- c("#ffffff", "#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")

             
p <- ggplot(totalmean, aes(x=totalmean$Method, y=totalmean$Percentage))

p <- p + geom_tile(aes(fill = factor(Breaks)), colour = "white") 
p <- p + scale_fill_manual("Error Range", values=palette)
#p <- p + scale_x_discrete(limits=totalmean$Method)
#p <- p + scale_y_discrete(limits=totalmean$Percentage)
p <- p + ggtitle("Total Mean") + xlab("Methods") + ylab("Percentage of total Missing Values ")
p <- p + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()      , panel.background = element_blank()) 
p <- p + theme(strip.background = element_blank())
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p      



##################################################################################################
## Heat map 2, Mean MCAR



MCARmean <-Imputation_Matrix %>%
  select(Method,Percentage,Mean_Error_MCAR)


MCARmean$Method <-factor(MCARmean$Method,labels = c("BPCA","KNN","LLS","mean","min","PPCA","RF","SVD"))
MCARmean$Percentage <- factor(MCARmean$Percentage )

quantiles_vec <- quantile(MCARmean$Mean_Error_MCAR,na.rm = T)
Q1 <-round(quantiles_vec[[1]],digits = 2)
Q2 <-round(quantiles_vec[[2]],digits = 2)
Q3 <-round(quantiles_vec[[3]],digits = 2)
Q4 <-round(quantiles_vec[[4]],digits = 2)
Q5 <-round(quantiles_vec[[5]],digits = 2)




p.breaks <- c(Q5,Q4,Q3, Q2, Q1)

MCARmean$Breaks <- cut(x=MCARmean$Mean_Error_MCAR, breaks=p.breaks, labels=paste0("< ", as.character(p.breaks)[1:(length(p.breaks)-1)]), include.lowest=TRUE)
palette <- c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")
#palette <- c("#ffffff", "#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")


p1 <- ggplot(MCARmean, aes(x=MCARmean$Method, y=MCARmean$Percentage))

p1 <- p1 + geom_tile(aes(fill = factor(Breaks)), colour = "white") 
p1 <- p1 + scale_fill_manual("Error Range", values=palette)
#p <- p + scale_x_discrete(limits=totalmean$Method)
#p <- p + scale_y_discrete(limits=totalmean$Percentage)
p1 <- p1 + ggtitle("MCAR Mean") + xlab("") + ylab("Percentage of total Missing Values ")
p1 <- p1 + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()      , panel.background = element_blank()) 
p1 <- p1 + theme(strip.background = element_blank())
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p1      

##################################################################################################
## Heat map 3, Mean MAR




MARmean <-Imputation_Matrix %>%
  select(Method,Percentage,Mean_Error_MAR)


MARmean$Method <-factor(MARmean$Method,labels = c("BPCA","KNN","LLS","mean","min","PPCA","RF","SVD"))
MARmean$Percentage <- factor(MARmean$Percentage )

quantiles_vec <- quantile(MARmean$Mean_Error_MAR,na.rm = T)
Q1 <-round(quantiles_vec[[1]],digits = 2)
Q2 <-round(quantiles_vec[[2]],digits = 2)
Q3 <-round(quantiles_vec[[3]],digits = 2)
Q4 <-round(quantiles_vec[[4]],digits = 2)
Q5 <-round(quantiles_vec[[5]],digits = 2)



p.breaks <- c(Q5,Q4,Q3, Q2, Q1)

MARmean$Breaks <- cut(x=MARmean$Mean_Error_MAR, breaks=p.breaks, labels=paste0("< ", as.character(p.breaks)[1:(length(p.breaks)-1)]), include.lowest=TRUE)
palette <- c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")
#palette <- c("#ffffff", "#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")


p2 <- ggplot(MARmean, aes(x=MARmean$Method, y=MARmean$Percentage))

p2 <- p2 + geom_tile(aes(fill = factor(Breaks)), colour = "white") 
p2 <- p2 + scale_fill_manual("Error Range", values=palette)
#p <- p + scale_x_discrete(limits=totalmean$Method)
#p <- p + scale_y_discrete(limits=totalmean$Percentage)
p2 <- p2 + ggtitle("MAR Mean") + xlab("") + ylab("")
p2 <- p2 + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()      , panel.background = element_blank()) 
p2 <- p2 + theme(strip.background = element_blank())
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p2      



##################################################################################################
## Heat map 4, Mean MNAR




MNARmean <-Imputation_Matrix %>%
  select(Method,Percentage,Mean_Error_MNAR)


MNARmean$Method <-factor(MNARmean$Method,labels = c("BPCA","KNN","LLS","mean","min","PPCA","RF","SVD"))
MNARmean$Percentage <- factor(MNARmean$Percentage )

quantiles_vec <- quantile(MNARmean$Mean_Error_MNAR,na.rm = T)
Q1 <-round(quantiles_vec[[1]],digits = 2)
Q2 <-round(quantiles_vec[[2]],digits = 2)
Q3 <-round(quantiles_vec[[3]],digits = 2)
Q4 <-round(quantiles_vec[[4]],digits = 2)
Q5 <-round(quantiles_vec[[5]],digits = 2)



p.breaks <- c(Q5,Q4,Q3, Q2, Q1)

MNARmean$Breaks <- cut(x=MNARmean$Mean_Error_MNAR, breaks=p.breaks, labels=paste0("< ", as.character(p.breaks)[1:(length(p.breaks)-1)]), include.lowest=TRUE)
palette <- c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")
#palette <- c("#ffffff", "#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d")


p3 <- ggplot(MNARmean, aes(x=MNARmean$Method, y=MNARmean$Percentage))

p3 <- p3 + geom_tile(aes(fill = factor(Breaks)), colour = "white") 
p3 <- p3 + scale_fill_manual("Error Range", values=palette)
#p <- p + scale_x_discrete(limits=totalmean$Method)
#p <- p + scale_y_discrete(limits=totalmean$Percentage)
p3 <- p3 + ggtitle("MNAR Mean") + xlab("Methods") + ylab(" ")
p3 <- p3 + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()      , panel.background = element_blank()) 
p3 <- p3 + theme(strip.background = element_blank())
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p3     


####################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




multiplot(p,p1, p2, p3, cols=2)



