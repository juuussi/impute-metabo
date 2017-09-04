mnarmar <- marietta$MAR_MNAR

  if (length(mnarmar)==0){
    MCAR <- data.frame(ListMCAR = MissingVar$MissVar)
   
    
  }else{
    
    MCAR <-data.frame(ListMCAR =setdiff(MissingVar$MissVar,mnarmar))
    
    #  left truncation MNAR 
    #  Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
    #  if the distributions are the same then p-values are high and that means they are left trancated if they are different then they are MAR
    
    # Goodness of fit for left truncated data
    models <- list()
    Pval <- list()
    Padj <-list()
    
    for (i in 1:nrow(MAR_MNARvariables)){
      
      
      xt <-na.omit(miss_data[,i])
      
      
      threshold <- min(na.omit((miss_data[,i])))
      #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
      
      #models[[i]] <- truncgof::ad.test(xt, "pnorm",list(meanlog = mean(simulated_data), sdlog = sd(simulated_data)), H = threshold,alternative ="two.sided")
      models[[i]] <- truncgof::ks.test(xt, "pnorm",list(mean(simulated_data),  sd(simulated_data)), H = threshold,  alternative ="two.sided")
      
      Pval[[i]] <-models[[i]]$p.value
      Padj[[i]]<-p.adjust(Pval[[i]], method = "fdr")
      
    }
    
    Padj <- as.numeric(as.character(Padj))
    Padj <- data.frame(pvalues= Padj,ListVar =MAR_MNARvariables$ListMAR_MNAR) 
    #View(Padj)
    #SigpVal <-Padj$pvalues[which(Padj$pvalues <= 0.05)]
    MARvariables <-data.frame(ListMAR = MAR_MNARvariables$ListMAR_MNAR[which(Padj$pvalues<=0.05)])
    #MARvariables <- data.frame(ListMAR =DetectMissMAR2$ListMar)
    View(MARvariables)
    ##
    MNARvariables <-setdiff(MAR_MNARvariables$ListMAR_MNAR,MARvariables$ListMAR)
    MNARvariables <- data.frame(ListMNAR = MNARvariables)
    View(MNARvariables)
    
    
  }
  
  
} else{
  MCARvariables <- data.frame(ListMCAR = MissingVar$MissVar)
  View(MCARvariables)
}
View(MCARvariables)