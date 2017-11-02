library(dplyr)
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
tmp_data <- Imputation_Matrix %>% select(Method,Percentage,Total_mean)
colnames(tmp_data)[3] <- "Error"
tmp_data$Type <- "Total"
melted_data <- tmp_data
tmp_data <- Imputation_Matrix %>% select(Method,Percentage,Mean_Error_MCAR)
colnames(tmp_data)[3] <- "Error"
tmp_data$Type <- "MCAR"
melted_data <- rbind(melted_data, tmp_data)
tmp_data <- Imputation_Matrix %>% select(Method,Percentage,Mean_Error_MAR)
colnames(tmp_data)[3] <- "Error"
tmp_data$Type <- "MAR"
melted_data <- rbind(melted_data, tmp_data)
tmp_data <- Imputation_Matrix %>% select(Method,Percentage,Mean_Error_MNAR)
colnames(tmp_data)[3] <- "Error"
tmp_data$Type <- "MNAR"
melted_data <- rbind(melted_data, tmp_data)
melted_data$Type <- as.factor(melted_data$Type)
head(melted_data)

melted_data$Method <-factor(melted_data$Method,labels = c("BPCA","1/2min","KNN","mean","min","PPCA","RF","SVD","zero"))

melted_data$Percentage <- factor(melted_data$Percentage )



#p.breaks <-rev( c(Q5,Q4,Q3, Q2, Q1))
#p.breaks <-c(0,1.05,1.25, 1.5, 1.75,2,3,Inf)
p.breaks <-c( 0,0.09,0.1, 0.2 , 0.4,  0.6, 0.9, 1.0, 1.5,2,Inf)

melted_data$Breaks <- cut(x=melted_data$Error, breaks=p.breaks, labels=paste0("< ", as.character(p.breaks)[2:(length(p.breaks))]),include.lowest=T)

#palette <- c("#ffffff", "#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08459")
#palette <- rev(c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"))
#palette <- rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))
palette <- brewer.pal(11,"BrBG")



p <- ggplot(melted_data, aes(x=melted_data$Method, y=melted_data$Percentage))
p <- p + facet_wrap(~ Type, ncol=2)
p <- p + geom_tile(aes(fill = factor(Breaks)), colour = "white") 
p <- p + scale_fill_manual("Error Range", values=palette, drop=FALSE)
p <- p + ggtitle("") + xlab("") + ylab("Percentage of total Missing Values ")
p <- p + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()      , panel.background = element_blank()) 
p <- p + theme(strip.background = element_blank())
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p      

