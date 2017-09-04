# 
# ColNames <-dplyr::select(Imputation_Matrix,starts_with("Mean_Error")) %>%
#   colnames()
# 
# df <-data.frame()
# 
# for(i in ColNames){
#   fin <- Imputation_Matrix[,i][is.finite(Imputation_Matrix[,i])]
#   index <- which(Imputation_Matrix[,i] == min(fin))
#   method <- as.character(Imputation_Matrix$Method[index]) 
#   df.tmp <- data.frame(Mean_Error = i, Winner_method = method)
#   df <- rbind(df, df.tmp)
# }


# fin <- Imputation_Matrix$`Mean_Error_ MNAR `[is.finite(Imputation_Matrix$`Mean_Error_ MNAR `)]
# index <- which(Imputation_Matrix$`Mean_Error_ MNAR ` == min(fin))
# method <- as.character(Imputation_Matrix$Method[index])
# libraries needed
library(missForest)

library(pcaMethods)
###################################################################################

# choose path
path <- "~/projects/impute-metabo/"
#output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
source(paste0(path,"src/functions.R"))
# start logging process by creating a logging file

###################################################################################
# use dummy reference data (combination of different metabolomics data)
reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))

reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
simulated_data <- simulate_data(data=reference_data, nrow=150, ncol=2000)
miss_data <- simulate_missingness(data=simulated_data, mcar=0.30, mnar=0, mar=0)
## model selector
nFeautures <- dim(miss_data)[2] 
nsamples   <- dim(miss_data)[1]

# initialize the model selector variable
model.selector <-rep(0,nFeautures)
# compute the mean values of each variable
mean.vect <- colMeans(miss_data,na.rm = F)    
# estimate the percenatge of Nas
pNAs = length(which(is.na(mean.vect)))/length(mean.vect)
# - estimate the mean and standard deviation of the original
#   distribution of the means using quantile regression
# -----------------------------------------------------------------------------------
  
  upper.q = 0.99
  
  q.normal = qnorm(seq((pNAs+0.001),(upper.q+0.001),(upper.q-pNAs)/(upper.q*100)), 
                   mean = 0, sd = 1)
  
  q.mean.vect = quantile(mean.vect,
                           probs = seq(0.001,(upper.q+0.001),0.01),
                           na.rm = T) 
  
  temp.QR = lm(q.mean.vect ~ q.normal)
  QR.obj = temp.QR
  
  mean.CDD = temp.QR$coefficients[1]
  sd.CDD = as.numeric(temp.QR$coefficients[2])
  
   # _______________________________________________________________________________________
  # estimate the the censoring threshold
  # ---------------------------------------------------------------------------------------
  
  nPoints = 512
  
  # - create the support of the complete data distribution to evaluate the CDFs
  min.support <- mean.CDD - 4*sd.CDD
  max.support <- mean.CDD + 4*sd.CDD
  
  # - empirical estimation of the cdf of the observed data 
  
  mean.vect.sorted <- sort(mean.vect)
  Fn <- ecdf(mean.vect.sorted)
  
  # - discretize the support of the CDFs
  
  support <-c(seq(min.support,
                  min(mean.vect,na.rm = T),
                  length = nPoints),
              mean.vect.sorted,
              seq(max(mean.vect,na.rm = T),
                  max.support,
                  length = nPoints))
  
   # - evaluate the empirical cdf in the discrete points given by support
  
  cdf.OD <- Fn(support)
  
  # - evaluate the theoretical cdf with MLE parameters estimated above 
  #   in the discrete points given by support
  
  cdf.CD <- pnorm(support, 
                 mean = mean.CDD, 
                 sd = sd.CDD)
  
  # - compute the difference between the CDFs: evaluated in a particular point, 
  #   it gives the proportion of true discoveries (an estimation of the proportion
  #   of the true samples smaller than the evaluation point, given the data)
  
  diff.cdf<-(cdf.CD - cdf.OD)
   
  # - the function to be optimize in order to obtain the optimum censoring threshold
  #   at its optimum, the function should maximize diff.cdf and should minimize 
  #   in the same time the cdf.OD
  
  # obj.fnc = diff.cdf-cdf.OD
  # obj.fnc = (diff.cdf)/(cdf.OD-1)
  # obj.fnc[is.infinite(obj.fnc)] = 0 
  # obj.fnc = 10*obj.fnc
  obj.fnc <- (diff.cdf+1)/(cdf.OD+1) - 1
  obj.fnc <- obj.fnc*10
  
  if (max(obj.fnc) > 0){
    censoring.thr <-support[which(obj.fnc == max(obj.fnc))]
  }else{
    censoring.thr <- min.support-1
  }
  
  # - set the model selector variable to "1" for the proteins with a mean value 
  #   higher than the KS statistic
  
  model.selector[which(mean.vect > censoring.thr)] = 1    
  
  result<-list(model.selector,censoring.thr)
  
  View(result)
  ##
  
  install.packages("BaylorEdPsych")
  install.packages("mvnmle")
  
  library(BaylorEdPsych)
  library(mvnmle)
  library(dplyr)
  #<update>
  data(EndersTable1_1) #retrieve the enders dataset
  View(EndersTable1_1) #view the dataset on R's data viewer
  LittleMCAR(EndersTable1_1)
  View(miss_data)

  
  miss_data <- cbind(miss_data, "observation"=1:nrow(miss_data)) 
  
  
  miss_data <- miss_data %>% as.data.frame()%>%dplyr::select(observation,everything())
  
  a <-LittleMCAR(miss_data)
  #</update>
  
  LittleMCAR(miss_data[ !apply(miss_data, 1, function(x) all(is.na(x))), -10])
  ###############################
 
  library(multcompView) 
  ERROR_0_6 <-data.frame(error =results_new$Error[results_new$Percentage == 0.1& results_new$Type =="MAR"])
  
  methdos_o6 <-data.frame(method =results_new$Method[results_new$Percentage == 0.1& results_new$Type =="MAR"])
 dp <- duplicated(ERROR_0_6)
  df6 <-cbind.data.frame(methdos_o6,data.frame(dp),ERROR_0_6)
  View(df6)
  PP1 <-pairwise.wilcox.test(df6$error,df6$method, p.adj = "BH")
  boxplot(error ~ method,data = df6)
  
  TYPE <- unique(results_new$Type)
  PERCENTAGE <- unique(results_new$Percentage)
  
 # df <- list()
  for(i in 1:length(TYPE)){
    for(j in 1:length(PERCENTAGE)){
      ERROR  <- data.frame(error =results_new$Error[results_new$Percentage == PERCENTAGE[j] & results_new$Type==TYPE[i]])
      METHODS  <-data.frame(method = results_new$Method[results_new$Percentage == PERCENTAGE[j] & results_new$Type==TYPE[i]])
      #df <- c(list(cbind.data.frame(ERROR,METHODS)), df)
      #PP <-pairwise.wilcox.test(ERROR$error,METHODS$method, p.adj = "BH")
      boxplot(ERROR$error ~ METHODS$method)
      
      
    }
  }


  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  tri.to.squ<-function(x)
  {
    rn<-row.names(x)
    cn<-colnames(x)
    an<-unique(c(cn,rn))
    myval<-x[!is.na(x)]
    mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
    for(ext in 1:length(cn))
    {
      for(int in 1:length(rn))
      {
        if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
        mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
        mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      }
      
    }
    return(mymat)
  }
  
  mymat<-tri.to.squ(PP1$p.value)
  myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
  boxplot(error ~ method,data = df6)
  text(c(1,2,3),0.2+max(df6$error),c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]))
  
  
  
  
  
  results_new$Type<-trimws(results_new$Type)
  #############################################
  library(VIM)
  
  path <- "~/projects/impute-metabo/"
  #output_path <- '/home/users/mariekok/projects/impute-metabo/results/result.csv'
  source(paste0(path,"src/functions.R"))
  #use dummy reference data (combination of different metabolomics data)
  reference_data <- as.matrix(read.csv(paste0(path, "data/reference_data.csv")))
  
  simulated_data <- simulate_data(data=reference_data, nrow=100, ncol=1000)
  miss_data <- simulate_missingness(data=simulated_data, mcar=0.02, mnar=0, mar=0.3)
  a <- aggr(miss_data)
  #Load dataset
  # data(sleep, package = "VIM")
  # 
  # x <- as.data.frame(abs(is.na(sleep)))
  # 
  # #Elements of x are 1 if a value in the sleep data is missing and 0 if non-missing.
  # head(sleep)
  # head(x)
  # 
  # #Extracting variables that have some missing values.
  # y <- x[which(sapply(x, sd) > 0)]
  # cor(y)
  # 
  # #We see that variables Dream and NonD tend to be missing together. To a lesser extent, this is also true with Sleep and NonD, as well as Sleep and Dream.
  # 
  # #Now, looking at the relationship between the presence of missing values in each variable and the observed values in other variables:
  # cor(sleep, y, use="pairwise.complete.obs")
  # 
  # #NonD is more likely to be missing as Exp, BodyWgt, and Gest increases, suggesting that the missingness for NonD is likely MAR rather than MCAR.
  # 
  

  x <- as.data.frame(abs(is.na(miss_data)))
  
  #Elements of x are 1 if a value in the sleep data is missing and 0 if non-missing.

  head(x)
  
  #Extracting variables that have some missing values.
  y <- x[which(sapply(x, sd) > 0)]
  cor(y)
  
  #We see that variables Dream and NonD tend to be missing together. To a lesser extent, this is also true with Sleep and NonD, as well as Sleep and Dream.
  
  #Now, looking at the relationship between the presence of missing values in each variable and the observed values in other variables:
  b <-cor(miss_data, y, use="pairwise.complete.obs")
  
  #NonD is more likely to be missing as Exp, BodyWgt, and Gest increases, suggesting that the missingness for NonD is likely MAR rather than MCAR.
  
  ##################################
  
  # ++++++++++++++++++++++++++++
  # flattenCorrMatrix
  # ++++++++++++++++++++++++++++
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  library(Hmisc)
  res2<-rcorr(as.matrix(a[,1:20]))
  flattenCorrMatrix(res2$r, res2$P)
  
  symnum(res2$r, abbr.colnames = FALSE)
  
  cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y)[["p.value"]]
    z <- outer(
      colnames(x), 
      colnames(x), 
      Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
  }
  
  cor.test.p(mtcars)
  
  
  
  
  
  
  
  
  
  ########
  
  
  # 
  # PS1 <- matrix(nrow=length(miss_data), ncol=length(y))
  # 
  # for (i in 1:length(miss_data)){
  #   for(k in 1:length(y)){
  #     PS1[i,k] <-biserial.cor(as.matrix(miss_data[,i]),as.matrix( y[,k]),use = "complete.obs")
  #  
  # }
  # }
  
  
  
  # Estimate the corellations between the continuous and dichotomous variable
  
  mixedcorre <-mixed.cor(x= as.matrix(miss_data), d=as.matrix(y), use="pairwise.complete.obs",method="pearson")
  PCorel <- mixedcorre$rx[]
  
  # correlation test reports the probabilities of correlations
  
  # NOTE
  # the corelation matrix need to have row names and column names in order 
  #for the function corr.p to work!
  row.names(PCorel) <- 1:nrow(PCorel)
  colnames(PCorel) <- 1:ncol(PCorel)
  #corr.p(r,n,adjust="holm",alpha=.05)
  
  CorTest <- psych::corr.p(r=PCorel,n=nrow(miss_data),adjust="fdr",alpha=.05)
  correlation_matrix <-CorTest$r
  # extract pvalues
  pval <- CorTest$p
  pval <- round(pval,digits = 2)
  CI <- CorTest$ci
  CI$p <-round(CI$p,digits = 2)
  
  # find which variables are significally correlated
  DetectMiss <-which(CI$p<0.05)
  
  CorrelVariables <- data.frame(VarNames =rownames(CI)[DetectMiss])
  new <-t(data.frame((strsplit(as.character(CorrelVariables$VarNames),'-'))))
  COL1 <-unique(unlist(new[,1]))
  COL2 <-unique(unlist(new[,2]))
  
  